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
save_path = 'chonk'

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
    return (math.ceil(abs(v)) * ve) + m


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


def normalize_val(val, mi, ma):
    return (val - mi) / (ma - mi)


def get_data_from_sector(sector_tuple, sector_width, data, lons_ax, lats_ax, resolution, el_w):
    mask_lon = (lons_ax >= sector_tuple[0]) & (lons_ax < sector_tuple[0] + sector_width)
    mask_lat = (lats_ax <= sector_tuple[1]) & (lats_ax > sector_tuple[1] - sector_width)

    xe = [np.where(lons_ax >= i)[0][0] for n, i in enumerate(lons_ax[mask_lon]) if n % resolution == 0]
    ye = [np.where(lats_ax <= i)[0][0] for n, i in enumerate(lats_ax[mask_lat]) if n % resolution == 0]

    xm = (np.amin(lons_ax[mask_lon])-sector_tuple[0])
    ym = sector_tuple[1]-(np.amax(lats_ax[mask_lat]))

    xxs = np.sign(xm)
    xxi = int(math.floor(abs(xm) * (el_w/resolution))*xxs)
    yyi = int(math.floor(ym * (el_w/resolution)))

    data_sector = np.zeros([int((el_w*sector_width)/resolution), int((el_w*sector_width)/resolution)])

    data_min = np.nanmin(data)
    data_max = np.nanmax(data)

    for iy, dy in enumerate(ye):
        for ix, dx in enumerate(xe):
            try:
                data_sector[iy+yyi, ix+xxi] = data[dy, dx]  #normalize_val(data[dy, dx], data_min, data_max)
            except IndexError:
                pass

    extents = [sector_tuple[0], sector_tuple[0]+sector_width, sector_tuple[1]-sector_width, sector_tuple[1]]
    return data_sector, extents, [data_min, data_max]


def save_zoom_level(master_poly, simplification, min_area, level=None, places=None, depths=None, datum=None):
    bounds = refine_bounds(master_poly.bounds)
    width = bounds[1]
    height = bounds[2]
    tiles_count = width * height
    print(tiles_count, bounds)

    depths_res = [None, 20, 15, 6, 3]

    all_shapes, simplified_poly = prepare_map_points(master_poly, simplification, min_area)

    if level == 0:
        for shape in all_shapes:
            shape[1] = coords_to_flat(shape[1].coords)
        deliver(all_shapes, f'../{save_path}/static/datasets/', f'guides')
        return

    simplified_poly = make_valid(simplified_poly)

    tiles_places = [None] * int(tiles_count)
    tiles_lines = [None] * int(tiles_count)
    tiles_fills = [None] * int(tiles_count)
    tiles_depths = [None] * int(tiles_count)
    tiles_data = [None] * int(tiles_count)

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
                #print(save_places)
                tiles_places[n] = save_places

        if test_fill is True:
            shape = k.intersection(simplified_poly)
            if shape.type == 'Polygon':
                tile_fill.append(package_poly(shape))
            if shape.type == 'MultiPolygon':
                for n_poly in shape.geoms:
                    tile_fill.append(package_poly(n_poly))
            tiles_fills[n] = tile_fill
            #//now do the depths:
            # minutes = 60
            # subset = depths_res[level]
            # subset_width = int(minutes / subset)
            # sub_ctr = 0
            # depths_subset = np.zeros((subset_width, subset_width))
            #
            # for i in range(3600):
            #     dx = math.floor(i / minutes)
            #     dy = i % minutes
            #     ix = int((x * minutes) + dx)
            #     iy = int((y * minutes) + dy)
            #     if dx % subset == 0 and dy % subset == 0:
            #         mx = math.floor(sub_ctr / subset_width)
            #         my = sub_ctr % subset_width
            #         try:
            #             o = depths.shape[0]
            #             val = depths[o - iy, ix - minutes]
            #             # if np.isnan(val):
            #             #     val = 0.0
            #
            #             depths_subset[mx, my] = val  #depths[o-iy, ix-minutes]
            #         except IndexError:
            #             depths_subset[mx, my] = None  # 0.0  #'no'
            #         sub_ctr += 1
            #

            depth, extents, limits = get_data_from_sector(
                [x, y],
                1.0,
                depths['data'],
                depths['lons'],
                depths['lats'],
                depths_res[level],
                60
            )

            y = np.where(np.isnan(depth), None, depth)
            tiles_depths[n] = y.tolist()
            #print(depths_subset.tolist())

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
            tiles_primitive_index[n] = 1

    primitive_index = []
    for index, on_off in enumerate(tiles_primitive_index):

        if on_off is not None:
            primitive_index.append(index)
            depth_blob = {}
            blob = {}

            if tiles_depths[index] is not None:
                depth_blob['depths'] = tiles_depths[index]
                deliver(depth_blob, f'../{save_path}/static/datasets/depths/zoom-{level}/', f'tile-{index}')

            if tiles_lines[index] is not None:
                blob['lines'] = tiles_lines[index]
            if tiles_fills[index] is not None:
                blob['fills'] = tiles_fills[index]
            if tiles_places[index] is not None:
                blob['places'] = tiles_places[index]

            deliver(blob, f'../{save_path}/static/datasets/tiles/zoom-{level}/', f'tile-{index}')

    deliver(primitive_index, f'../{save_path}/static/datasets/indices/', f'med-low-index-{level}')


#//original étape
if __name__ == '__main__':

    data_map_spec = {
        "levels": 4,
        "rect":
        {
            "min_X": -7,
            "min_Y": 30,
            "max_X": 38,
            "max_Y": 46
        },
        "name": "med_full",
        "map_degrees": 1.0
    }

    depths = {}
    depths['data'] = np.loadtxt(open("./data/bathy_med.csv", "rb"), delimiter=";", skiprows=1)
    depths['lons'] = np.loadtxt(open("./data/lon_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0)
    depths['lats'] = np.loadtxt(open("./data/lat_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0)




    populated_places = gpd.read_file("data/ne_10m_populated_places/ne_10m_populated_places.shp")
    print(populated_places.head())  #iloc[0])
    print(len(populated_places))

    med_marine_regions = gpd.read_file("data/goas_med/med_local.shp")
    print(med_marine_regions.iloc[0])
    med_poly = med_marine_regions['geometry'][0]
    med_poly_name = med_marine_regions['name'][0]

    minx, miny, maxx, maxy = refine_bounds(med_poly.bounds)[0]
    data_map_spec['rect'] = {"min_X": minx, "min_Y": miny, "max_X": maxx, "max_Y": maxy}
    k = box(minx, miny, maxx, maxy)
    print(k)

    places_smash = []
    for index, loc in enumerate(populated_places.geometry):
        if loc.within(k):
            places_smash.append({
                'scale': int(math.floor((populated_places.SCALERANK[index] * 4) / 10)),
                'pop': int(populated_places.POP_MAX[index]),
                'class': populated_places.FEATURECLA[index],
                'tz': populated_places.TIMEZONE[index],
                'name': populated_places.NAME[index],
                'country': populated_places.ADM0NAME[index],
                'locale': populated_places.ADM1NAME[index],
                'loc': loc
            })

    # save_zoom_level(med_poly, simplification=0.0025, min_area=0.001, level=0)
    # exit()

    pathlib.Path(f'../{save_path}/static/datasets/indices/').mkdir(exist_ok=True)
    pathlib.Path(f'../{save_path}/static/datasets/depths/').mkdir(exist_ok=True)
    pathlib.Path(f'../{save_path}/static/datasets/tiles/').mkdir(exist_ok=True)
    pathlib.Path(f'../{save_path}/static/datasets/data/').mkdir(exist_ok=True)

    for zoom in range(0, 4):
        valid_places = [i for i in places_smash if i['scale'] == zoom]
        print(med_poly_name, med_poly.bounds, 'zoom_level', zoom + 1, 'valid_places', len(valid_places))

        pathlib.Path(f'../{save_path}/static/datasets/depths/zoom-{zoom + 1}').mkdir(exist_ok=True)
        pathlib.Path(f'../{save_path}/static/datasets/tiles/zoom-{zoom + 1}').mkdir(exist_ok=True)
        pathlib.Path(f'../{save_path}/static/datasets/data/zoom-{zoom + 1}').mkdir(exist_ok=True)
        save_zoom_level(med_poly, simps[zoom], areas[zoom], zoom + 1, valid_places, depths)



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
