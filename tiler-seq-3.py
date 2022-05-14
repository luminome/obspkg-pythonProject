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


area_limits = [0.05, 0.000025]
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


def virgo(v, m):
    ve = math.copysign(1, v)
    return (math.ceil(abs(v)) * ve) + m


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


def save_zoom_level(master_poly, poly_guides, simplification, min_area, data_map_spec, level=None, places=None, depths=None, datum=None):
    bounds = data_map_spec['bounds']  #refine_bounds(master_poly.bounds)
    deg_sector = data_map_spec['map_degrees']  # 0.5
    #guides = data_map_spec['guides']


    width = bounds[1] * (1 / deg_sector)
    height = bounds[2] * (1 / deg_sector)
    sector_count = width * height

    print(sector_count, bounds)
    depths_res = [None, 15, 10, 4, 2]
    #print(depths.shape)

    all_shapes, simplified_poly = prepare_map_points(master_poly, simplification, min_area)

    # if level == 0:
    #     for shape in all_shapes:
    #         shape[1] = coords_to_flat(shape[1].coords)
    #     deliver(all_shapes, f'../{save_path}/static/datasets/', f'guides')
    #     return

    simplified_poly = make_valid(simplified_poly)

    tiles_places = [None] * int(sector_count)
    tiles_lines = [None] * int(sector_count)
    tiles_fills = [None] * int(sector_count)
    tiles_depths = [None] * int(sector_count)
    tiles_data = [None] * int(sector_count)

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
                echo_poly(shape, level)
                tile_fill.append(package_poly(shape))
            if shape.type == 'MultiPolygon':
                for n_poly in shape.geoms:
                    echo_poly(n_poly, level)
                    tile_fill.append(package_poly(n_poly))
            tiles_fills[n] = tile_fill

            get_data = get_data_from_sector(
                [bounds[0][0]+(x * deg_sector), bounds[0][3]-(y * deg_sector)],
                deg_sector,
                depths['data'],
                depths['lons'],
                depths['lats'],
                depths_res[level],
                60
            )

            if get_data is not None:
                depth = get_data['data']
                y = np.where(np.isnan(depth), None, depth)
                tiles_depths[n] = y.tolist()

        overlaps = 0
        for i in all_shapes:
            #  print(i[0])
            test_line = k.intersects(i[1])
            # #// doest this need a lookup table?
            # #// just have to know what original contour the guide is projected from
            # TODO: test against test_guide from predetrmined group.
            # //test_o = k.intersects(med_poly)
            if test_line is True:

                lk = k.intersection(i[1])
                guide_id = None
                test_guide = [g for g in poly_guides if k.intersects(g[2]) and g[0] == i[0]]
                if len(test_guide):
                    guide_id = test_guide[0][1]
                    print(guide_id)

                if lk.type == 'LineString':
                    #print(n, 'single', i[0])
                    # test_guide = [g for g in guides if lk.within(Polygon(g[2]))]
                    # print(test_guide)
                    tf = [i[0], lk.coords, guide_id]
                    tile_line.append(tf)

                if lk.type == 'MultiLineString':
                    #print(n, 'multi', i[0])
                    for line in lk.geoms:
                        # test_guide = [g for g in guides if line.within(Polygon(g[2]))]
                        # print(test_guide)
                        tfm = [i[0], line.coords, guide_id]
                        tile_line.append(tfm)

                overlaps += 1
                #break

        if overlaps:
            #plt.plot(*k.exterior.xy, color=(0.8, 0.8, 0.8))
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


def echo_poly(poly, level):
    if level == 4:
        if len(poly.interiors):
            plt.plot(*poly.exterior.xy, color="pink")
            for sub_poly in poly.interiors:
                plt.plot(*sub_poly.xy, color="pink")
        else:
            plt.plot(*poly.exterior.xy, color="pink")


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
        "bounds": [],
        "map_degrees": 1,
        "debug": 2,
        "guides": None
    }

    med_marine_regions = gpd.read_file("data/goas_med/med_local.shp")
    print(med_marine_regions.iloc[0])
    med_poly = med_marine_regions['geometry'][0]
    med_poly_name = med_marine_regions['name'][0]

    print(med_poly.bounds)

    # bounds = refine_bounds(med_poly.bounds)
    #minx, miny, maxx, maxy = bounds[0]
    minx, miny, maxx, maxy = [6, 32, 16, 44]
    # data_map_spec['rect'] = {"min_X": 4, "min_Y": 37, "max_X": 8, "max_Y": 43}
    #data_map_spec['rect'] = {"min_X": minx, "min_Y": miny, "max_X": maxx, "max_Y": maxy}
    #data_map_spec['bounds'] = bounds

    k = box(minx, miny, maxx, maxy)


    print(k.bounds)

    #bounds = refine_bounds(k.bounds)
    bounds = refine_bounds(med_poly.bounds)

    minx, miny, maxx, maxy = bounds[0]

    k_grow = box(minx - 0.1, miny - 0.1, maxx + 0.1, maxy + 0.1)

    data_map_spec['rect'] = {"min_X": minx, "min_Y": miny, "max_X": maxx, "max_Y": maxy}

    data_map_spec['bounds'] = bounds

    k = box(minx, miny, maxx, maxy)
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # #//DESTROY the following ways before
    # med_poly = k_grow.intersection(med_poly)

    #
    # def plot(shape, r=1):
    #     print(shape.type, '!')
    #
    #     def plot_poly(polygon):
    #         plt.plot(*polygon.exterior.xy, color=(0.5, 1/r, 1/r), linewidth="1")
    #         if len(polygon.interiors):
    #             for i in polygon.interiors:
    #                 plt.plot(*i.xy, color="blue", linewidth="1")
    #
    #     if shape.type == 'MultiPolygon':
    #         for s in shape.geoms:
    #             plot_poly(s)
    #     else:
    #         plot_poly(shape)
    #
    #
    master_shapes, simplified_poly = prepare_map_points(med_poly, simps[3], simps[0])
    #
    #
    buffer_constant = 0.1
    simplification_const = 0.005
    guides = []
    for c, n in enumerate(master_shapes):

        poly_c = k.intersection(Polygon(n[1]))
        poly_shape = k.difference(poly_c)

        if poly_shape.type == 'MultiPolygon':
            #plot(poly_shape)
            kbf = poly_shape.buffer(buffer_constant)
            kbf = kbf.simplify(simplification_const)
            #plot(kbf, 4)
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
                        #for s_poly in kbf.geoms:
                        s_poly = [sp for sp in kbf.geoms if sp.area > 1.0][0]
                        guides.append([n[0], len(guides), s_poly.exterior, n[2],  True])
                        #plot(s_poly, 4)
                    else:
                        #plt.plot(*kbf.exterior.xy, color="green", linewidth="1")
                        guides.append([n[0], len(guides), kbf.exterior, n[2], True])


    guide_hooks = []
    for r in guides:
        print(r)
        gh = r.copy()
        gh[2] = coords_to_flat(gh[2].coords)
        guide_hooks.append(gh)

    # plt.show()
    # exit()


    deliver(guide_hooks, f'../{save_path}/static/datasets/', f'guides2')

    # if level == 0:
    #     for shape in all_shapes:
    #         shape[1] = coords_to_flat(shape[1].coords)
    #     deliver(all_shapes, f'../{save_path}/static/datasets/', f'guides')
    #     return


    #data_map_spec['guides'] = guides
    #exit()

    # plt.show()
    # exit()
    #
    # poly_c = k.intersection(med_poly)
    #
    # poly_shape = k.difference(poly_c)
    #
    # guide_shapes, simplified_poly = prepare_map_points(poly_shape, simps[0], 0.1)
    #
    #
    #
    # prefilter = []
    # for n in simplified_poly.geoms:
    #     sp = n.buffer(0.2)
    #     sp = sp.simplify(0.05)
    #     prefilter.append(sp)
    #
    # poly_shape = MultiPolygon(prefilter)
    # guide_shapes, simplified_poly = prepare_map_points(poly_shape, simps[0], areas[0])
    #
    # # renumber the shape indices
    # new_guides = []
    # for c, n in enumerate(guide_shapes):
    #     print(c,n)
    #     new_guides.append([c, n[1]])
    #
    # n_count = 0
    # guides = []
    #
    # for n in new_guides:
    #     p = k.intersection(n[1])  ##n[1]?
    #     if p.type == 'MultiLineString':
    #         for i in p.geoms:
    #             guides.append([n_count, i, i.is_ring])
    #             n_count += 1
    #     elif not p.is_empty:
    #         guides.append([n_count, p, p.is_ring])
    #         n_count += 1
    #
    # print([(a[0], a[2]) for a in guides])
    #
    # data_map_spec['guides'] = guides
    #
    # # #lt.show()
    # # exit()
    #
    # def plot(shape, r=1):
    #     print(shape.type, '!')
    #
    #     def plot_poly(polygon):
    #         plt.plot(*polygon.exterior.xy, color=(0.2, 1/r, 0), linewidth="1")
    #         if len(polygon.interiors):
    #             for i in polygon.interiors:
    #                 plt.plot(*i.xy, color="blue", linewidth="1")
    #
    #     if shape.type == 'MultiPolygon':
    #         for s in shape.geoms:
    #             plot_poly(s)
    #     else:
    #         plot_poly(shape)
    #
    # for guide in guides:
    #     plt.plot(*guide[1].xy, color="blue", linewidth="1")
    #
    # # plot(med_poly)
    #
    # #
    # #
    # # poly_c = k.intersection(med_poly)
    # # poly_shape = k.difference(poly_c)
    # # #
    # # guide_shapes, simplified_poly = prepare_map_points(poly_shape, simps[0], 0.1)
    # #
    # # prefilter = []
    # # is_enclosed_group_flag = False
    # #
    # # for n in simplified_poly:
    # #     if len(n.interiors):
    # #         is_enclosed_group_flag = True
    # #
    # #     sp = n.buffer(0.2)
    # #     sp = sp.simplify(0.05)
    # #     prefilter.append(sp)
    # #
    # #     print(n.type, len(n.interiors))
    # #
    # # poly_shape = MultiPolygon(prefilter)
    # # guide_shapes, simplified_poly = prepare_map_points(poly_shape, simps[0], areas[0])
    # #
    # # new_guides = []
    # # for c, n in enumerate(guide_shapes):
    # #     print(c, n[0], n[1].type)
    # #     new_guides.append([c, n[1]])
    # #     plot(Polygon(n[1]), n[0] + 1)
    # # #plot(simplified_poly)
    # #
    # # #if not is_enclosed_group_flag:
    # #     # #ODO: this only happens when the entire map shape is bounded. So must cut and slice the guides.
    # # n_count = 0
    # #
    # # guides = []
    # #
    # # for n in new_guides:
    # #     p = k.intersection(n[1])
    # #
    # #     if p.type == 'MultiLineString':
    # #         for i in p.geoms:
    # #             if i.is_ring:
    # #                 plt.plot(*i.xy, color="blue", linewidth="1")
    # #             else:
    # #                 plt.plot(*i.xy, color="red", linewidth="2")
    # #
    # #             guides.append([n_count, i, i.is_ring])
    # #             print(n_count, i.type, i.is_ring)
    # #             n_count += 1
    # #     else:
    # #         if p.is_empty:
    # #             pass
    # #         else:
    # #             guides.append([n_count, p, p.is_ring])
    # #             print(n_count, p.type, p.is_ring)
    # #             n_count += 1
    # #
    # #         plt.plot(*p.xy, color="black", linewidth="1")
    # #
    # #     #plot(spc)
    # #
    # #     #pass
    # #
    # #
    # #
    # # print([(a[0],a[2]) for a in guides])
    #
    # # for n in guide_shapes:
    # #     plot(Polygon(n[1]), n[0]+1)
    # #
    # #
    #
    # #
    # # guide_shapes, simplified_poly = prepare_map_points(poly_shape, simps[0], areas[0])
    # #
    # # n = guide_shapes[1]
    # # plot(Polygon(n[1]), n[0] + 1)
    # # for n in guide_shapes:
    # #     plot(Polygon(n[1]), n[0]+1)
    #
    #
    # # sp = simplified_poly.buffer(0.075)
    # # sp = sp.simplify(0.05)
    #
    # #guide_shapes, simplified_poly = prepare_map_points(sp, simps[3], areas[0])
    # #
    # # plot(simplified_poly)
    #
    # # plot(med_poly)
    # # print(poly_shape.type)
    # #
    #
    #
    # # depths = np.loadtxt(open("./data/bathy_med.csv", "rb"), delimiter=";", skiprows=1)
    # # exit()
    #
    # # med_poly = k.intersection(med_poly)
    # #med_poly = k.difference(med_poly)
    #
    # # all_shapes, simplified_poly = prepare_map_points(med_poly, simps[3], areas[0])
    # # print(simplified_poly.type)
    # # for s in all_shapes:
    # #     plt.plot(*s[1].xy, color="red", linewidth="1")
    # #     print(s)
    #
    # # plt.show()
    # # exit()
    # #
    # # #poly_one = simplified_poly
    # # for n, s in enumerate(all_shapes):
    # #     plt.plot(*s[1].xy, color="red", linewidth="1")
    # #     print(s)
    # #
    # # plt.show()
    # # exit()
    # # # plt.plot(*poly_one.exterior.xy, color="red", linewidth="1")
    # # #
    # # # if len(poly_one.interiors):
    # # #     for i in poly_one.interiors:
    # # #         plt.plot(*i.xy, color="blue", linewidth="1")
    # #
    # #
    # #
    # #
    # # if poly_one.type == 'MultiPolygon':
    # #     print(len(poly_one.geoms))
    # #     for p in poly_one.geoms:
    # #         print(p.type)
    # #         if len(p.interiors):
    # #             print(len(p.interiors))
    # #             p_buf = p.buffer(1 / 16)
    # #             plt.plot(*p_buf.exterior.xy, color="green", linewidth="1")
    # #
    # #             for i in p_buf.interiors:
    # #                 plt.plot(*i.xy, color="blue", linewidth="1")
    # #
    # #         else:
    # #             plt.plot(*p.exterior.xy, color="red", linewidth="1")
    # #
    # #     pass
    # #
    # # #p_buf = k.difference(simplified_poly)
    # #
    # # # outline = Polygon(simplified_poly.exterior.coords)
    # # #
    # # # poly = k.symmetric_difference(outline)
    # # #
    # # # print(poly.type)
    # # #
    # # # p_buf = poly.buffer(1 / 16)
    # # #
    # # # plt.plot(*p_buf.exterior.xy, color="red", linewidth="1")
    # # #
    # # # if len(p_buf.interiors):
    # # #     for i in p_buf.interiors:
    # # #         plt.plot(*i.xy, color="blue", linewidth="1")
    # # #
    # # #
    # # #
    # # #
    # # #
    # # #
    # # # # # p_buf = poly.buffer(1 / 16)
    # # # #
    # # # # if p_buf.type == 'MultiPolygon':
    # # # #     for i in p_buf.geoms:
    # # # #         print(i.type)
    # # # #         plt.plot(*i.exterior.xy, color="red", linewidth="1")
    # # # #
    # # # #
    # # # #
    # # # # #plt.plot(*k.exterior.xy, color="green")
    # # #
    # # plt.show()
    # # exit()
    # #
    # #
    # # shapes = [poly.exterior]
    # # print(poly.type)
    # #
    # # for sub_poly in poly.interiors:
    # #     shapes.append(sub_poly)
    # #
    # # # these should be simplified before they're cut!
    # # # if this is the outermost exterior shape, the buffer selection should be reversed.
    # #
    # # for n, s in enumerate(shapes):
    # #     if s.intersects(k):
    # #
    # #         p_buf = s.buffer(1/16)
    # #         print(p_buf.type)
    # #         p_buf = p_buf.simplify(0.02)
    # #
    # #         if p_buf.type =='MultiPolygon':
    # #             for i in p_buf.geoms:
    # #                 plt.plot(*i.exterior.xy, color="red", linewidth="1")
    # #         else:
    # #             if n == 0:
    # #                 #p_buf = p_buf.exterior.union(p_buf.interiors)
    # #                 plt.plot(*p_buf.exterior.xy, color="black", linewidth="4")
    # #                 for f in p_buf.interiors:
    # #                     plt.plot(*f.xy, color="green", linewidth="4")
    # #
    # #
    # #             else:
    # #                 plt.plot(*p_buf.exterior.xy, color="pink", linewidth="1")
    # #
    # #
    # #
    # #         # p_buf = s.buffer(0.5)
    # #         # print(p_buf.type)
    # #         #
    # #         sl = LineString(k.exterior.coords)
    # #         p = s.difference(sl)
    # #
    # #         print(n, p.type)
    # #
    # #         if p.type =='MultiLineString':
    # #             for i in p.geoms:
    # #                 plt.plot(*i.xy, color="black", linewidth="1")
    # #         else:
    # #             plt.plot(*p.xy, color="green", linewidth="1")
    # #         print(s.type)
    # #
    #
    # # med_poly = k.intersection(med_poly)
    # #
    plt.plot(*k.exterior.xy, color="green")
    # print(len(shapes))
    # plt.show()
    # exit()

    # if len(poly.interiors):
    #     plt.plot(*poly.exterior.xy, color="red", linewidth="4")
    #     for sub_poly in poly.interiors:
    #         plt.plot(*sub_poly.xy, color="red", linewidth="4")
    # else:
    #     plt.plot(*poly.exterior.xy, color="black", linewidth="8")


    # for shape in all_shapes:
    #     #shape[1] = coords_to_flat(shape[1].coords)
    #     #base = shape[1]
    #
    #     plt.plot(*shape[1].xy, color="red")
    #
    # #//plt.plot(*med_poly.exterior.xy, color="red")
    #


    #exit()

    depths = {}
    depths['data'] = np.loadtxt(open("./data/bathy_med_tt.csv", "rb"), delimiter=";", encoding=None, skiprows=0)
    depths['lons'] = np.loadtxt(open("./data/lon_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0)
    depths['lats'] = np.loadtxt(open("./data/lat_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0)

    populated_places = gpd.read_file("data/ne_10m_populated_places/ne_10m_populated_places.shp")
    print(populated_places.head())  #iloc[0])
    print(len(populated_places))

    deliver(data_map_spec, f'../{save_path}/static/datasets/', 'map-spec')

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
        save_zoom_level(med_poly, guides, simps[zoom], areas[zoom], data_map_spec, zoom + 1, valid_places, depths)

    plt.show()
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
