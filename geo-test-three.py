#!/usr/bin/env python3

import os
import operator
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, box
from shapely.validation import make_valid
import json
import math
import pathlib
import numpy as np


area_limits = [0.005, 0.0001]
simp_limits = [0.025, 0.00075]

areas = np.linspace(area_limits[0], area_limits[1], 4)
simps = np.linspace(simp_limits[0], simp_limits[1], 4)

xmin = -7.0
xmax = 42.5
ymin = 30.0
ymax = 48.0
# minx, miny, maxx, maxy
bounds_keys = [0, -1, 0, 0]


def update_goas():
    # slice from 6 to 7 to get the med only!
    #                  name  latitude   longitude  min_Y      min_X      max_Y      max_X      area_km2
    #6 Mediterranean Region  38.13065   19.70067   30.06809   -6.03255   47.37640   42.35496   2988248
    med_marine_regions = gpd.read_file("data/GOaS_v1_20211214/goas_v01.shp", rows=slice(6, 7))
    med_marine_regions.to_file("data/goas_med/med_only.shp")
    print(med_marine_regions.head())
    return


def punch_bosphorus():
    min_Y = 30.06809
    min_X = -6.03255
    max_Y = 47.3764
    max_X = 42.35496

    mask = Polygon([
        (26.4, max_Y),
        (max_X, max_Y),
        (max_X, 40.54),
        (32.5, 40.54),
        (29.4, 41.15),
        (28.8, 41.2),
        (26.4, 41.72),
        (26.4, max_Y)
    ])

    return mask


def alter(the_poly):
    mask = punch_bosphorus()
    med_poly = the_poly.difference(mask).geoms[0]
    med_marine_regions['geometry'][0] = med_poly
    med_marine_regions.to_file("data/goas_med/med_local.shp")
    #plt.plot(*mask.exterior.xy, color="green")


def virgo(v, m):
    ve = math.copysign(1, v)
    return (math.ceil(abs(v)) * ve)+m


def refine_bounds(bounds):
    tb = []
    for i, b in enumerate(bounds):
        tb.append(virgo(b, bounds_keys[i]))
    w = abs(tb[0]) + abs(tb[2])
    h = abs(tb[3]) - abs(tb[1])
    return [tb, w, h]


def deliver(data, prefix, name):
    filename = f'{prefix}data-{name}-test.json'
    with open(filename, 'w') as fp:
        json.dump(data, fp)  #, indent=2)


def prepare_map_points(the_poly, simplification, min_area):
    shapes = []
    interiors = []
    exterior = []

    def qualify(ct, sub_shape, min_a):
        area_test = Polygon(sub_shape.coords)
        shape_simple = sub_shape.simplify(simplification)
        if area_test.area >= min_a:
            return [ct, shape_simple]
        return None

    for c, n in enumerate(the_poly.interiors):
        e = qualify(c, n, min_area)
        if e is not None:
            interiors.append(e[1])
            shapes.append(e)

    e = qualify(0, the_poly.exterior, min_area)
    if e is not None:
        shapes.append(e)
        exterior.append(e[1])

    return shapes, Polygon(exterior[0], interiors)


def coords_to_flat(coords_list):
    coords_points = []
    list(coords_points.extend((round(cx, 4), round(cy, 4))) for cx, cy in coords_list)
    return coords_points


def package_poly(poly):
    poly_lex = {'out': coords_to_flat(poly.exterior.coords)}
    if len(poly.interiors):
        poly_lex['ins'] = [coords_to_flat(p.coords) for p in poly.interiors]
    return poly_lex


def save_zoom_level(master_poly, simplification, min_area, level=None, places=None):
    bounds = refine_bounds(master_poly.bounds)
    width = bounds[1]
    height = bounds[2]
    tiles_count = width * height

    print(tiles_count, bounds)

    all_shapes, simplified_poly = prepare_map_points(master_poly, simplification, min_area)

    if level == 0:
        for shape in all_shapes:
            shape[1] = coords_to_flat(shape[1].coords)

        deliver(all_shapes, f'../gandi-sac-wp/static/datasets/', f'guides')
        return

    simplified_poly = make_valid(simplified_poly)

    tiles_places = [None] * int(tiles_count)
    tiles_lines = [None] * int(tiles_count)
    tiles_fills = [None] * int(tiles_count)
    tiles_primitive_index = [None] * int(tiles_count)

    for n in range(int(tiles_count)):
        tile_line = []
        tile_fill = []

        x = math.floor(n / height)
        y = n % height

        k = box(
            x + bounds[0][0],
            y + bounds[0][1],
            x + bounds[0][0] + 1,
            y + bounds[0][1] + 1)

        test_fill = k.intersects(simplified_poly)
        k = make_valid(k)

        if places is not None:
            test_places = [p for p in places if k.contains(p['loc']) is True]
            if len(test_places):
                save_places = []
                for p in test_places:
                    save_place = p.copy()
                    save_place['loc'] = coords_to_flat(p['loc'].coords)
                    save_places.append(save_place)

                print(save_places)
                tiles_places[n] = save_places

        if test_fill is True:
            shape = k.intersection(simplified_poly)
            if shape.type == 'Polygon':
                tile_fill.append(package_poly(shape))
                # plt.plot(*s['out'].xy, color=(1.0, 0.0, 0.0), linewidth=4)
                # if 'ins' in s:
                #     [plt.plot(*sn.coords.xy, color=(1.0, 0.0, 0.0), linewidth=4) for sn in s['ins']]
            if shape.type == 'MultiPolygon':
                for n_poly in shape.geoms:
                    tile_fill.append(package_poly(n_poly))
                    #s = package_poly(n_poly)
                    # plt.plot(*s['out'].xy, color=(1.0, 0.0, 1.0), linewidth=4)
                    # if 'ins' in s:
                    #     [plt.plot(*sn.coords.xy, color=(1.0, 0.0, 1.0), linewidth=4) for sn in s['ins']]
            tiles_fills[n] = tile_fill

        overlaps = 0
        for i in all_shapes:
            test_line = k.intersects(i[1])
            # test_o = k.intersects(med_poly)
            if test_line is True:
                lk = k.intersection(i[1])
                if lk.type == 'LineString':
                    tile_line.append([i[0], lk.coords])
                if lk.type == 'MultiLineString':
                    for line in lk.geoms:
                        tile_line.append([i[0], line.coords])
                overlaps += 1
                #break

        if overlaps:
            #plt.plot(*k.exterior.xy, color=(0.8, 0.8, 0.8))
            lines_set = []
            for line in tile_line:
                points = coords_to_flat(line[1])
                lines_set.append([line[0], points])
            tiles_lines[n] = lines_set

        if overlaps or test_fill:
            tiles_primitive_index[n] = 1  #.append(n)

    primitive_index = []
    for index, on_off in enumerate(tiles_primitive_index):

        if on_off is not None:
            primitive_index.append(index)
            #print(index)
            blob = {}

            if tiles_lines[index] is not None:
                blob['lines'] = tiles_lines[index]
            if tiles_fills[index] is not None:
                blob['fills'] = tiles_fills[index]
            if tiles_places[index] is not None:
                blob['places'] = tiles_places[index]

            deliver(blob, f'../gandi-sac-wp/static/datasets/zoom-{level}/', f'tile-{index}')

    deliver(primitive_index, '../gandi-sac-wp/static/datasets/', f'med-low-index-{level}')


def test_zoom_level(master_poly, simplification, min_area, level=None, places=None):
    bounds = refine_bounds(master_poly.bounds)
    width = bounds[1]
    height = bounds[2]
    tiles_count = width * height

    print(tiles_count, bounds)

    all_shapes, simplified_poly = prepare_map_points(master_poly, simplification, min_area)

    if level == 0:
        for shape in all_shapes:
            shape[1] = coords_to_flat(shape[1].coords)

        deliver(all_shapes, f'../gandi-sac-wp/static/datasets/', f'guides')
        return

    simplified_poly = make_valid(simplified_poly)

    tiles_lines = [None] * int(tiles_count)
    tiles_primitive_index = [None] * int(tiles_count)

    for n in range(int(tiles_count)):
        tile_line = []
        tile_fill = []

        x = math.floor(n / height)
        y = n % height

        k = box(
            x + bounds[0][0],
            y + bounds[0][1],
            x + bounds[0][0] + 1,
            y + bounds[0][1] + 1)

        test_fill = k.intersects(simplified_poly)
        k = make_valid(k)

        overlaps = 0
        for i in all_shapes:
            test_line = k.intersects(i[1])
            # test_o = k.intersects(med_poly)
            if test_line is True:
                lk = k.intersection(i[1])
                if lk.type == 'LineString':
                    tile_line.append([i[0], lk.coords])
                if lk.type == 'MultiLineString':
                    for line in lk.geoms:
                        tile_line.append([i[0], line.coords])
                overlaps += 1
                #break

        if test_fill:
            plt.plot(*k.exterior.xy, color=(0.5, 0.5, 0.5))
            lines_set = []
            for line in tile_line:
                points = coords_to_flat(line[1])
                lines_set.append([line[0], points])
            tiles_lines[n] = lines_set

        if overlaps or test_fill:
            tiles_primitive_index[n] = 1  #.append(n)

    # primitive_index = []
    # for index, on_off in enumerate(tiles_primitive_index):
    #
    #     if on_off is not None:
    #         primitive_index.append(index)
    #         #print(index)
    #         blob = {}
    #
    #         if tiles_lines[index] is not None:
    #             blob['lines'] = tiles_lines[index]
    #         if tiles_fills[index] is not None:
    #             blob['fills'] = tiles_fills[index]
    #         if tiles_places[index] is not None:
    #             blob['places'] = tiles_places[index]

            #deliver(blob, f'../gandi-sac-wp/static/datasets/zoom-{level}/', f'tile-{index}')

    #deliver(primitive_index, '../gandi-sac-wp/static/datasets/', f'med-low-index-{level}')






if __name__ == '__main__':

    med_marine_regions = gpd.read_file("data/goas_med/med_local.shp")
    print(med_marine_regions.iloc[0])
    med_poly = med_marine_regions['geometry'][0]
    med_poly_name = med_marine_regions['name'][0]

    print(med_poly.bounds)
    zoom = 1

    fig, ax = plt.subplots(figsize=(10, 6))

    test_zoom_level(med_poly, simps[zoom], areas[zoom], zoom + 1)

    #
    dat = np.loadtxt(open("./data/bathy_med.csv", "rb"), delimiter=";", skiprows=1)
    ldy = np.loadtxt(open("./data/lat_vector.csv", "rb"), delimiter=";", skiprows=1)
    ldx = np.loadtxt(open("./data/lon_vector.csv", "rb"), delimiter=";", skiprows=1)

    h, w = dat.shape
    print(w, h)
    print(h, ldy.shape)
    print(w, ldx.shape)

    bounds = refine_bounds(med_poly.bounds)
    width = bounds[1]
    height = bounds[2]
    size = width*height

    print(bounds, width*height)
    minutes = 60

    plt.imshow(dat, extent=(-6, 36.5, 30, 46), interpolation="none")
    plt.plot(*med_poly.exterior.xy, color=(1.0, 0.0, 0.0), linewidth=1)

    #//verts = np.zeros((int(width*3), int(height*3)))
    p_long = np.zeros(int(width * 3))
    p_lati = np.zeros(int(height * 3))

    xrange = np.arange(0, ldx.shape[0], 20)
    xva = np.take(ldx, xrange)
    yrange = np.arange(0, ldy.shape[0], 20)
    yva = np.take(ldy, yrange)

    # print(xrange)
    # print(yrange)
    # # print(xva)
    # exit()





    for n in range(int(size)):
        x = math.floor(n / height)
        y = n % height

        k = box(
            x + bounds[0][0],
            y + bounds[0][1],
            x + bounds[0][0] + 1,
            y + bounds[0][1] + 1)

        subset = 20  #int(minutes/10)
        subset_width = int(minutes/subset)  #this is the effective width

        si = 0
        m = (subset_width, subset_width)
        indices_subset = np.zeros(m, dtype=[('x', 'i4'), ('y', 'i4')])
        depths_subset = np.zeros(m)

        for i in range(3600):
            dx = math.floor(i / minutes)
            dy = i % minutes
            ix = int((x*60)+dx)
            iy = int((y*60)+dy)
            part = (ix, iy)

            if dx % subset == 0 and dy % subset == 0:
                mx = math.floor(si / subset_width)
                my = si % subset_width
                print(ix, iy)

                try:
                    depths_subset[mx, my] = dat[iy, ix]
                    #try:
                    p_long[dx] = ldx[ix]
                    p_lati[dy] = ldy[iy]

                    # verts[ix, iy] = 0.0  #ldx[ix], ldy[iy]
                    # rpoint = Point([ldx[ix], ldy[iy]])
                    #plt.plot(*rpoint.xy, color=(1.0, 0.0, 0.0), linewidth=1)
                    #plt.scatter(gdx, gdy, s=0.01)

                except IndexError:
                    depths_subset[mx, my] = None
                    pass

                si += 1



        # print(n)
        # # print(indices.shape)
        # # print(indices)
        # # print(indices_subset.shape)
        # # print(indices_subset)
        # # print(depths.shape)
        # # print(depths)
        # print(depths_subset.shape)


        # print(depths_subset)

        # try:
        #     lon, lat = indices[0, 0]
        #     print(ldx[lon], ldy[lat], depths[0, 0])
        # except IndexError:
        #     lon, lat = indices[0, 0]
        #     print(lon, lat, 'null', depths[0, 0])
        #     pass

    gdx, gdy = np.meshgrid(xva, yva, indexing='ij')
    plt.scatter(gdx, gdy, s=1)
    plt.show()

    exit()

    print('zero',ldx[0],ldy[0],ldx[-1],ldy[-1])

    for r in range(0, 60):
        print(ldx[r])

    exit()

    print(w, h)
    print(ldy[0], ldy[h - 1])
    print(ldx[0], ldx[w - 2])

    xoffset = [0.11, 0.28]  #0.1

    minx = ldx[0]
    miny = ldy[h - 1]
    maxx = ldx[w - 2]
    maxy = ldy[0]

    w_ext = (abs(minx-maxx))
    h_ext = (maxy-miny)

    print('w', w_ext/w, 'h', h_ext/h)
    print('extents', w_ext, h_ext)
    # let X = 15 Y = 38
    x = math.floor((w * (abs(minx - 15))) / w_ext)
    y = math.ceil((h * (abs(miny - 38))) / h_ext)

    print(x, y, dat.item((y, x)))
    print(minx+(x*(w_ext/w)))
    print(miny+(y*(h_ext/h)))



    exit()
    # px = (ldx > 5) & (ldx < 6)
    # py = (ldy > 35) & (ldy < 36)
    #
    # indicesx = [i for i, x in enumerate(px) if x == True]
    # print(len(indicesx), indicesx)
    #
    # indicesy = [i for i, y in enumerate(py) if y == True]
    # print(len(indicesy), indicesy)
    #


    minx, miny, maxx, maxy = refine_bounds(med_poly.bounds)[0]

    # px = (ldx > 5) & (ldx < 6)
    # py = (ldy > 35) & (ldy < 36)
    print(minx, miny, maxx, maxy)
    #exit()

    # minx = -5
    # px = (ldx >= minx) & (ldx < minx + 1)
    # indicesx = [i for i, x in enumerate(px) if x == True]
    # print(len(indicesx), indicesx)
    # exit()
    f = 0.016666666666667
    print((minx-ldx[0])/f)  # 62
    print((maxy-ldy[0])/f)  # 2

    sm = (60, 60)
    sa = np.zeros(sm)

    #// taking first and last entries is a bad scene:
    #// these are arc minutes — the table has the extents clipped.
    # lat 46 to 30
    # lon -5 to 36.5
    # this data has omission on bounds, which it should but there needs to be an entry for col 0 and row 0 for this to work



    print(sa)
    exit()

    for lon in range(0, 2):
        print(lon, minx+lon, (minx+1)+lon)
        lon_min = minx
        px = (ldx >= minx+lon) & (ldx < (minx+1)+lon)
        indicesx = [i for i, x in enumerate(px) if x == True]
        print('X', len(indicesx), indicesx)
        if len(indicesx):
            print(ldx[indicesx[0]])


        # 0.016666666666667
        for lat in range(0, 2):

            lat_min = miny
            py = (ldy <= maxy - lat) & (ldy > (maxy - 1)-lat)
            indicesy = [i for i, y in enumerate(py) if y == True]
            print('Y', len(indicesy), indicesy)
            print(ldy[indicesy[0]])

    #
    #         # for slat in range(0, len(indicesx)):
    #         #     for slon in range(0, len(indicesy)):
    #         #
    #         #         print(slon, slat)

    # for starting rect:
    #
    # #aminx, miny, amaxx, maxy = refine_bounds(med_poly.bounds)[0]
    #





if __name__ == '__amain__':
    #update_goas()

    populated_places = gpd.read_file("data/ne_10m_populated_places/ne_10m_populated_places.shp")
    print(populated_places.head())  #iloc[0])
    print(len(populated_places))

    med_marine_regions = gpd.read_file("data/goas_med/med_local.shp")
    print(med_marine_regions.iloc[0])
    med_poly = med_marine_regions['geometry'][0]
    med_poly_name = med_marine_regions['name'][0]




    #bnds = [-6.032549115649033, 30.068086780520986, 36.21572947503529, 45.80891370810119]
    minx, miny, maxx, maxy = refine_bounds(med_poly.bounds)[0]
    k = box(minx, miny, maxx, maxy)
    print(k)

    places_smash = []
    for index, loc in enumerate(populated_places.geometry):
        if loc.within(k):
            places_smash.append({
                'scale': int(math.floor((populated_places.SCALERANK[index]*4)/10)),
                'pop': int(populated_places.POP_MAX[index]),
                'class': populated_places.FEATURECLA[index],
                'tz': populated_places.TIMEZONE[index],
                'name': populated_places.NAME[index],
                'country': populated_places.ADM0NAME[index],
                'locale': populated_places.ADM1NAME[index],
                'loc': loc
            })

    # places_smash.sort(key=lambda x: x['scale'])  #, reverse=True)
    # #valid_places = [i for i in places_smash if i['scale'] == 1]
    #
    # for l in places_smash:
    #     print(l['scale'])
    #
    # exit()





    # save_zoom_level(med_poly, simplification=0.0025, min_area=0.001, level=0)
    # exit()

    for zoom in range(0, 4):
        valid_places = [i for i in places_smash if i['scale'] == zoom]
        print(med_poly_name, med_poly.bounds, 'zoom_level', zoom+1, 'valid_places', len(valid_places))

        pathlib.Path(f'../gandi-sac-wp/static/datasets/zoom-{zoom+1}').mkdir(exist_ok=True)
        save_zoom_level(med_poly, simps[zoom], areas[zoom], zoom+1, valid_places)



















    #fig, ax = plt.subplots(figsize=(10, 6))






    # for c, t in enumerate(tiles_lines):
    #     tile_array = []
    #     for ps in t:
    #         points = coords_to_flat(ps[1])
    #         tile_array.append([ps[0], points])
    #
    #
    #
    #         #print(tiles_primitive_index[c], f'shape-id {s[0]}', len(points))
    #     tiles_lines[c] = tile_array
        #print(tiles_primitive_index[c])
        #//deliver(tiles_lines[c], f'../gandi-sac-wp/static/datasets/zoom-{zoom_level}/', f'tile-{tiles_primitive_index[c]}')



    #print(len(tiles), 'tiles')
    #print(tiles_primitive_index)
    #//deliver(tiles_primitive_index, '../gandi-sac-wp/static/datasets/', f'med-low-index-{zoom_level}')



    #print(all_shapes[1])

    # for c, n in enumerate(all_shapes):
    #     points = coords_to_flat(n[1].coords)
    #     plt.plot(*n[1].xy, color="blue")
    #     n[1] = points

    #//print(json.dumps(all_shapes))
    #//deliver(all_shapes, '../gandi-sac-wp/static/datasets/', 'med-low')




    #exterior_coords = med_poly.exterior.coords[:]

    # for c, n in enumerate(med_poly.interiors):
    #     #print(n.geom_type, n.length)
    #     an = Polygon(n.coords)
    #     ns = n.simplify(0.005, True)
    #     # print(c, an.area*100000)
    #     # ns = ns.buffer(0.005)
    #     # nf = ns.exterior.simplify(0.002, True)
    #     if an.area*10000 >= 2:
    #         satisfy += 1
    #         points += len(ns.coords)
    #         plt.plot(*ns.xy, color="blue")
    #
    # print('satisfy', satisfy, points)
    #
    # #mext = med_poly.exterior.simplify(1, True)
    # exterior_coords = med_poly.exterior.coords[:]
    # print(len(exterior_coords))
    #
    # plt.plot(*med_poly.exterior.xy, color="red")

    # med_marine_regions.plot(
    #     column='name',
    #     categorical=True,
    #     legend=True,
    #     figsize=(10, 6),
    #     markersize=25,
    #     cmap="Set2", ax=ax)

    #plt.show()

    #med_isolate =

    #med_marine_regions.to_file('data/goas_med/dataframe.geojson', driver='GeoJSON')

