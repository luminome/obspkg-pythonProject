#!/usr/bin/env python3

import os
import operator
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely import ops
from shapely.geometry import Polygon, MultiPolygon, Point, box, LineString
from shapely.validation import make_valid
import json
import math
import pathlib
import numpy as np


area_limits = [0.025, 0.000025]
simp_limits = [0.025, 0.0005]
areas = np.linspace(area_limits[0], area_limits[1], 4)
simps = np.linspace(simp_limits[0], simp_limits[1], 4)

xmin = -7.0
xmax = 42.5
ymin = 30.0
ymax = 48.0
# minx, miny, maxx, maxy
bounds_keys = [-1, -1, 1, 1]
save_path = 'chonk'

indices_aggregate = []


def update_goas(filename):
    # slice from 6 to 7 to get the med only!
    #                  name  latitude   longitude  min_Y      min_X      max_Y      max_X      area_km2
    # 6 Mediterranean Region  38.13065   19.70067   30.06809   -6.03255   47.37640   42.35496   2988248
    region = gpd.read_file("data/GOaS_v1_20211214/goas_v01.shp", rows=slice(6, 7))
    region.to_file(f"data/{filename}.shp")
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


def refine_bounds(bounds):
    tb = []
    for i, b in enumerate(bounds):
        tb.append(round(b)+bounds_keys[i])

    w = abs(tb[0] - tb[2])
    h = abs(tb[1] - tb[3])
    return [tb, w, h]


def deliver(data, prefix, name):
    filename = f'{prefix}data-{name}-test.json'
    with open(filename, 'w') as fp:
        json.dump(data, fp)  #, indent=2)


def prepare_map_points(the_poly, simplification, min_area):
    shapes = []
    interiors = []
    exterior = []
    multipoly = []
    group = None

    def qualify(ct, sub_shape, min_a):
        area_test = Polygon(sub_shape.coords)
        shape_simple = sub_shape.simplify(simplification, True)
        if area_test.area >= min_a:
            return [ct, shape_simple]
        return None

    if the_poly.type == 'Polygon':
        group = [the_poly]
    if the_poly.type == 'MultiPolygon':
        group = the_poly.geoms

    for the_poly in group:
        sub_interiors = []
        for c, n in enumerate(the_poly.interiors):
            e = qualify(c, n, min_area)
            if e is not None:
                sub_interiors.append(e[1])
                interiors.append(e[1])
                e.append('interior')
                shapes.append(e)

        e = qualify(0, the_poly.exterior, min_area)
        if e is not None:
            e.append('exterior')
            shapes.append(e)
            exterior.append(e[1])

            if len(the_poly.interiors):
                multipoly.append(Polygon(e[1], sub_interiors))
            else:
                multipoly.append(Polygon(e[1]))

    if len(group) == 1:
        print('prepare_map_points poly')
        return shapes, Polygon(exterior[0], interiors)
    else:
        print('prepare_map_points multipoly')
        return shapes, MultiPolygon(multipoly)


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
    print(sector_tuple, sector_width)
    # #// TODO: should the data have a "start-position" index instead or alongside of zeros?
    # #// in an effort to save on empty data being transferred?
    mask_lon = (lons_ax >= sector_tuple[0]) & (lons_ax < sector_tuple[0] + sector_width)
    mask_lat = (lats_ax <= sector_tuple[1]) & (lats_ax > sector_tuple[1] - sector_width)

    xe = [np.where(lons_ax >= i)[0][0] for n, i in enumerate(lons_ax[mask_lon]) if n % resolution == 0]
    ye = [np.where(lats_ax <= i)[0][0] for n, i in enumerate(lats_ax[mask_lat]) if n % resolution == 0]

    if len(xe) == 0 or len(ye) == 0:
        return None

    xm = (np.amin(lons_ax[mask_lon])-sector_tuple[0])
    ym = sector_tuple[1]-(np.amax(lats_ax[mask_lat]))

    xxs = np.sign(xm)
    xxi = int(math.floor(abs(xm) * (el_w/resolution))*xxs)
    yyi = int(math.floor(ym * (el_w/resolution)))

    sector_data = np.zeros([int((el_w*sector_width)/resolution), int((el_w*sector_width)/resolution)])
    sector_data[:] = np.nan

    data_min = np.nanmin(data)
    data_max = np.nanmax(data)

    for iy, dy in enumerate(ye):
        for ix, dx in enumerate(xe):
            try:
                sector_data[iy+yyi, ix+xxi] = data[dy, dx]  #normalize_val(data[dy, dx], data_min, data_max)
            except IndexError:
                pass

    extents = [sector_tuple[0], sector_tuple[0]+sector_width, sector_tuple[1]-sector_width, sector_tuple[1]]
    return {'data': sector_data, 'extents': extents, 'limits': [data_min, data_max]}


def save_zoom_level(master_poly, poly_guides, simplification, min_area, data_map_spec, level=None, data_group=None):
    bounds = data_map_spec['bounds']
    deg_sector = data_map_spec['map_degrees']
    width = bounds[1] * (1 / deg_sector)
    height = bounds[2] * (1 / deg_sector)
    sector_count = width * height

    print('level-', level, sector_count, bounds)

    all_shapes, ms_poly = prepare_map_points(master_poly, simplification, min_area)
    simplified_poly = make_valid(ms_poly)

    # #//CORE
    tiles_places = [None] * int(sector_count)
    tiles_lines = [None] * int(sector_count)
    tiles_fills = [None] * int(sector_count)
    # #//secondary
    tiles_data = {}
    for ke in data_group.keys():
        print(ke)
        tiles_data[ke] = [None] * int(sector_count)

    urban_mash = data_group['urban']['data'].copy()

    tiles_primitive_index = [None] * int(sector_count)

    for n in range(int(sector_count)):
        tile_line = []
        tile_fill = []
        y = math.floor(n / width)
        x = n % width
        k = box(
            bounds[0][0] + x * deg_sector,
            bounds[0][3] - y * deg_sector,
            bounds[0][0] + x * deg_sector + deg_sector,
            bounds[0][3] - y * deg_sector - deg_sector)

        test_fill = k.intersects(simplified_poly)
        if test_fill is True:
            aggreg_data = [t for t in indices_aggregate if t['id'] == n]

            if len(aggreg_data):
                index_data = aggreg_data[0]
            else:
                ag_data = {
                    'id': n,
                    'lv': [0, 0, 0, 0],
                    'dt': [],
                    'fl': [],
                    'au': [],
                    'gd': []
                }
                indices_aggregate.append(ag_data)
                index_data = indices_aggregate[-1]

            index_data['lv'][level - 1] = level

            if 'urban' in data_group and level == 1:
                test_urban_places = [p for p in urban_mash if type(p['poly']) != list and k.intersects(p['poly']) is True]
                for u in test_urban_places:
                    u['poly'] = coords_to_flat(u['poly'].exterior.coords)
                tiles_data['urban'][n] = test_urban_places
                index_data['fl'].append("urban") if len(test_urban_places) > 0 else index_data['fl']

            if 'places' in data_group:
                places_mash = data_group['places']['data']
                places = [i for i in places_mash if i['scale'] == level]

                if len(places):
                    test_places = [p for p in places if k.contains(p['loc']) is True]
                    if len(test_places):
                        print('tile', n, len(test_places), 'places')
                        save_places = []
                        for p in test_places:
                            save_place = p.copy()
                            save_place['loc'] = coords_to_flat(p['loc'].coords)
                            save_places.append(save_place)
                        tiles_places[n] = save_places
                        index_data['au'].append("places") if "places" not in index_data['au'] else index_data['au']


            shape = k.intersection(simplified_poly)
            if shape.type == 'Polygon':
                echo_poly(shape, level)
                tile_fill.append(package_poly(shape))
            if shape.type == 'MultiPolygon':
                for n_poly in shape.geoms:
                    echo_poly(n_poly, level)
                    tile_fill.append(package_poly(n_poly))
            tiles_fills[n] = tile_fill

            for datum in (d for d in data_group.items() if d[1]['type'] == 'pointdata'):
                layer = datum[1]['data']
                sector = (bounds[0][0]+(x * deg_sector), bounds[0][3]-(y * deg_sector),)

                get_data = get_data_from_sector(
                    sector,
                    deg_sector,
                    layer['data'],
                    layer['lons'],
                    layer['lats'],
                    layer['reso'][level],
                    layer['degs']
                )

                if get_data is not None:
                    index_data['dt'].append(datum[0]) if datum[0] not in index_data['dt'] else index_data['dt']
                    d_values = get_data['data'].copy()
                    ty = np.where(np.isnan(d_values), None, d_values)
                    tiles_data[datum[0]][n] = ty.tolist()

        overlaps = 0
        for i in all_shapes:
            #  print(i[0])
            test_line = k.intersects(i[1])
            # #// just have to know what original contour the guide is projected from
            if test_line is True:
                lk = k.intersection(i[1])
                guide_id = None
                test_guide = [g for g in poly_guides if k.intersects(g[2]) and g[0] == i[0]]
                if len(test_guide):
                    guide_id = test_guide[0][1]
                    #print(guide_id, end=" ")
                    gstr = f"guide-{guide_id}"
                    index_data['gd'].append(gstr) if gstr not in index_data['gd'] else index_data['gd']
                if lk.type == 'LineString':
                    tf = [i[0], lk.coords, guide_id]
                    tile_line.append(tf)

                if lk.type == 'MultiLineString':
                    for line in lk.geoms:
                        tfm = [i[0], line.coords, guide_id]
                        tile_line.append(tfm)

                overlaps += 1

        if overlaps:
            lines_set = []
            for line in tile_line:
                points = coords_to_flat(line[1])
                lines_set.append([line[0], points, line[2]])
            tiles_lines[n] = lines_set

        if overlaps or test_fill:
            tiles_primitive_index[n] = 1

    if data_map_spec['debug'] == 1:
        return

    primitive_index = []
    for index, on_off in enumerate(tiles_primitive_index):

        if on_off is not None:
            primitive_index.append(index)
            depth_blob = {}
            temps_blob = {}
            urban_blob = {}
            blob = {}

            if tiles_data['depths'][index] is not None:
                depth_blob['depths'] = tiles_data['depths'][index].copy()
                deliver(depth_blob, f'../{save_path}/static/datasets/data/depths/zoom-{level}/', f'tile-{index}')

            if tiles_data['temps_mini'][index] is not None:
                temps_blob['temps'] = tiles_data['temps_mini'][index].copy()
                deliver(temps_blob, f'../{save_path}/static/datasets/data/temps_mini/zoom-{level}/', f'tile-{index}')

            if tiles_data['urban'][index] is not None and level == 1:
                if len(tiles_data['urban'][index]) > 0:
                    urban_blob['areas'] = tiles_data['urban'][index].copy()
                    deliver(urban_blob, f'../{save_path}/static/datasets/data/urban/', f'tile-{index}')

            if tiles_lines[index] is not None:
                blob['lines'] = tiles_lines[index]
            if tiles_fills[index] is not None:
                blob['fills'] = tiles_fills[index]
            if tiles_places[index] is not None:
                blob['places'] = tiles_places[index]

            deliver(blob, f'../{save_path}/static/datasets/tiles/zoom-{level}/', f'tile-{index}')

    #deliver(primitive_index, f'../{save_path}/static/datasets/indices/', f'med-low-index-{level}')


def echo_poly(poly, level):
    if level == 4:
        if len(poly.interiors):
            plt.plot(*poly.exterior.xy, color="pink")
            for sub_poly in poly.interiors:
                plt.plot(*sub_poly.xy, color="pink")
        else:
            plt.plot(*poly.exterior.xy, color="pink")


def plot(shape, r=1):
    print(shape.type, '!')

    def plot_poly(polygon):
        plt.plot(*polygon.exterior.xy, color=(0.5, 1/r, 1/r), linewidth="1")
        if len(polygon.interiors):
            for i in polygon.interiors:
                plt.plot(*i.xy, color="blue", linewidth="1")

    if shape.type == 'MultiPolygon':
        for s in shape.geoms:
            plot_poly(s)
    else:
        plot_poly(shape)


# #//original étape
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
if __name__ == '__main__':

    # #//SHOULD ALL DATA HAVE LEVELS PER SECTOR?
    # #//TODO: make collection of original paths and areas data
    # #//TODO: make the 'reso' variables scale with map_degrees
    # #//TODO: fix indexing!!!

    pathlib.Path(f'../{save_path}/static/datasets/indices/').mkdir(exist_ok=True)
    pathlib.Path(f'../{save_path}/static/datasets/tiles/').mkdir(exist_ok=True)
    pathlib.Path(f'../{save_path}/static/datasets/data/depths/').mkdir(exist_ok=True)
    pathlib.Path(f'../{save_path}/static/datasets/data/temps_mini').mkdir(exist_ok=True)
    pathlib.Path(f'../{save_path}/static/datasets/data/urban').mkdir(exist_ok=True)

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
        "bounds": [],
        "map_degrees": 1,
        "debug": 2,
        "guides": None
    }

    med_marine_regions = gpd.read_file("data/goas_med/med_local.shp")
    med_poly = med_marine_regions['geometry'][0]
    med_poly_name = med_marine_regions['name'][0]
    print(med_marine_regions.iloc[0])
    print(med_poly.bounds)

    populated_places = gpd.read_file("data/ne_10m_populated_places/ne_10m_populated_places.shp")
    print(populated_places.iloc[0])
    print(len(populated_places))

    urban_planning = gpd.read_file("data/ne_10m_urban_areas/ne_10m_urban_areas.shp")
    print(urban_planning.iloc[0])
    print(max(urban_planning['scalerank']), min(urban_planning['scalerank']))

    minx, miny, maxx, maxy = [6, 32, 16, 44]
    k = box(minx, miny, maxx, maxy)

    bounds = refine_bounds(k.bounds)
    #bounds = refine_bounds(med_poly.bounds)

    minx, miny, maxx, maxy = bounds[0]

    k_grow = box(minx - 0.1, miny - 0.1, maxx + 0.1, maxy + 0.1)

    data_map_spec['rect'] = {"min_X": minx, "min_Y": miny, "max_X": maxx, "max_Y": maxy}

    data_map_spec['bounds'] = bounds

    k = box(minx, miny, maxx, maxy)
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    master_shapes, simplified_poly = prepare_map_points(med_poly, simps[3], simps[0])
    #
    buffer_constant = 0.1
    simplification_const = 0.005
    guides = []
    for c, n in enumerate(master_shapes):

        poly_c = k.intersection(Polygon(n[1]))
        poly_shape = k.difference(poly_c)

        if poly_shape.type == 'MultiPolygon':
            kbf = poly_shape.buffer(buffer_constant)
            kbf = kbf.simplify(simplification_const)
            for w, pl in enumerate(kbf.geoms):
                rtl = k.intersection(pl.exterior)
                if rtl.type == 'MultiLineString':
                    dz = ops.linemerge(rtl)
                    print(dz.type)
                    if dz.type == 'LineString':
                        plt.plot(*dz.coords.xy, color="red", linewidth="4")
                        guides.append([n[0], len(guides), dz, None, dz.is_ring])
                    else:
                        for dzi in dz.geoms:
                            plt.plot(*dzi.coords.xy, color="red", linewidth="4")
                            guides.append([n[0], len(guides), dzi, None, dzi.is_ring])
                else:
                    plt.plot(*rtl.coords.xy, color="blue", linewidth="4")
                    guides.append([n[0], len(guides), rtl, None, rtl.is_ring])

        if poly_shape.type == 'Polygon':
            if len(poly_shape.interiors):
                for i in poly_shape.interiors:
                    bf = buffer_constant
                    if n[2] == 'exterior':
                        bf *= -1

                    i_poly = Polygon(i)
                    kbf = i_poly.buffer(bf)
                    kbf = kbf.simplify(simplification_const)

                    if kbf.type == 'MultiPolygon' and n[2] == 'exterior':
                        s_poly = [sp for sp in kbf.geoms if sp.area > 1.0][0]
                        guides.append([n[0], len(guides), s_poly.exterior, n[2],  True])
                    else:
                        guides.append([n[0], len(guides), kbf.exterior, n[2], True])

    guide_hooks = []
    for r in guides:
        print(r)
        gh = r.copy()
        gh[2] = coords_to_flat(gh[2].coords)
        guide_hooks.append(gh)

    deliver(guide_hooks, f'../{save_path}/static/datasets/indices/', 'guides')
    deliver(data_map_spec, f'../{save_path}/static/datasets/indices/', 'map-spec')

    # plt.plot(*k.exterior.xy, color="green")
    # print(len(shapes))
    # plt.show()
    # exit()
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    src_depths = {
        "data": np.loadtxt(open("./data/bathy_med_tt.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
        "lons": np.loadtxt(open("./data/lon_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
        "lats": np.loadtxt(open("./data/lat_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
        "reso": [None, 10, 6, 2, 1],
        "degs": 60
    }
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    src_temps = {
        "data": np.loadtxt(open("./data/thetao_t.csv", "rb"), delimiter=",", skiprows=0),
        "lons": np.loadtxt(open("./lon.csv", "rb"), delimiter=",", skiprows=0),
        "lats": np.loadtxt(open("./lat.csv", "rb"), delimiter=",", skiprows=0),
        "reso": [None, 8, 6, 3, 1],
        "degs": 24
    }
    src_temps['lats'] = np.flip(src_temps['lats'])
    src_temps['lons'] = np.flip(src_temps['lons'])
    src_temps['data'] = np.flip(src_temps['data'], 0)
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    places_smash = []
    for index, loc in enumerate(populated_places.geometry):
        if loc.within(k):
            places_smash.append({
                'scale': int(math.floor((populated_places.SCALERANK[index] * 4) / 10))+1,
                'pop': int(populated_places.POP_MAX[index]),
                'class': populated_places.FEATURECLA[index],
                'tz': populated_places.TIMEZONE[index],
                'name': populated_places.NAME[index],
                'country': populated_places.ADM0NAME[index],
                'locale': populated_places.ADM1NAME[index],
                'loc': loc
            })
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    urban_smash = []
    u_max = max(urban_planning['scalerank'])
    u_min = min(urban_planning['scalerank'])
    for index, loc in enumerate(urban_planning.geometry):
        if loc.within(k):
            sca = int(math.ceil(normalize_val(urban_planning.scalerank[index], u_min, u_max)*4))
            area = urban_planning.area_sqkm[index]
            urban_smash.append({
                'scale': sca,
                'area': area,
                'poly': loc
            })
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # print(urban_smash)
    # exit()




    data_inventory = {
        "depths": {"type": "pointdata", "data": src_depths},
        "temps_mini": {"type": "pointdata", "data": src_temps},
        "places": {"type": "custom-inline", "data": places_smash},
        "urban": {"type": "custom-inline", "data": urban_smash}
    }

    for zoom in range(0, 4):
        print(med_poly_name, med_poly.bounds, 'zoom_level', zoom + 1)

        pathlib.Path(f'../{save_path}/static/datasets/tiles/zoom-{zoom + 1}').mkdir(exist_ok=True)
        pathlib.Path(f'../{save_path}/static/datasets/data/depths/zoom-{zoom + 1}').mkdir(exist_ok=True)
        pathlib.Path(f'../{save_path}/static/datasets/data/temps_mini/zoom-{zoom + 1}').mkdir(exist_ok=True)

        save_zoom_level(med_poly, guides, simps[zoom], areas[zoom], data_map_spec, zoom + 1, data_inventory)

    for index in indices_aggregate:
        print(index)

    deliver(indices_aggregate, f'../{save_path}/static/datasets/indices/', 'index-aggregated')
    print('save completed.')

    #plt.show()
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
