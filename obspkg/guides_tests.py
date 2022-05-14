#!/usr/bin/env python3

import json
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import pickle
# from obspkg_util import make_guides, refine_bounds, parse_map_shapes, plot_polygon, line_s_to_list, poly_s_to_list

import obspkg_util as utl

from shapely.geometry import box, Polygon, LinearRing
from json_encoder import JsonSafeEncoder

with open('./configuration.json') as jfp:
    cfg = json.load(jfp)

area_limits = cfg['area_limits']
simp_limits = cfg['simp_limits']

areas = np.linspace(area_limits[0], area_limits[1], 4)
simps = np.linspace(simp_limits[0], simp_limits[1], 4)

cfg['criteria'] = {'areas': areas, 'simplify': simps}


if __name__ == '__main__':

    fig, ax = plt.subplots()

    # #//LOAD MAP DATA SHAPEFILE
    # med_marine_regions = gpd.read_file(f"{cfg['data_source_path']}/goas_med/med_local.shp")
    # med_poly = med_marine_regions['geometry'][0]
    # med_poly_name = med_marine_regions['name'][0]

    # my_data = np.genfromtxt('../corsica/contour20m_Corsica.csv', delimiter=',')
    # med_poly = LinearRing(my_data)

    # Load polygon from disc
    with open('../data/map_polygon', "rb") as poly_file:
        med_poly = pickle.load(poly_file)


    if cfg["rect_ignore"] == "default":
        bounds = utl.refine_bounds(med_poly.bounds, cfg, False)
        cfg['rect'] = {"min_X": bounds[0][0],
                               "min_Y": bounds[0][1],
                               "max_X": bounds[0][2],
                               "max_Y": bounds[0][3]}
    else:
        rect = cfg['rect']
        w = rect["max_X"] - rect["min_X"]
        h = rect["max_Y"] - rect["min_Y"]

        bounds = [[rect["min_X"], rect["min_Y"], rect["max_X"], rect["max_Y"]], w, h]

    minx, miny, maxx, maxy = bounds[0]
    bounds_box = box(minx, miny, maxx, maxy)
    print('Essential Bounds (deg):', minx, miny, maxx, maxy)

    # print(med_poly.type)
    # med_poly = Polygon(med_poly)
    # test_master_poly = med_poly.intersection(bounds_box)
    #
    # if cfg['rect_bounded']:
    #     test_master_poly = bounds_box.difference(test_master_poly)
    #
    # med_poly = test_master_poly



    #
    # print(test_master_poly.type)
    #
    # cfe = utl.line_s_to_list(test_master_poly)
    #
    # les = utl.merge_lines(cfe)
    # print(les)
    #
    # ax.plot(*les.coords.xy, color="black", linewidth=2)
    #
    # pl = Polygon(les)
    # utl.plot_polygon(ax, pl)



    #if "new_map" in cfg["includes"]:
    print('simplifications', cfg['criteria']['simplify'])
    print('areas', cfg['criteria']['areas'])

    map_collection = utl.parse_map_shapes(med_poly, bounds_box, cfg, ax)
    print('\nmap_collection')

    def get_map_shapes_guides(collection):
        guides_flat_list = []
        [guides_flat_list.extend(g['guides']) for g in collection]
        return guides_flat_list

    guides = get_map_shapes_guides(map_collection)
    # for n, l_item in enumerate(guides):
    #     out = json.dumps(l_item, indent=2, cls=JsonSafeEncoder)
    #     print(out)
    #     print(n, l_item)

    # "min_X": 6,
    # "min_Y": 36,
    # "max_X": 11,
    # "max_Y": 44
    # sector_box = box(10, 37, 11, 38)
    sector_box = box(6, 36, 11, 44)

    # #// assume level 4, IU think all these parts are in the right place.
    s_lv = 4

    def get_mapped(param_box, collection, guides_collection):
        # get the collection the overlaps this box
        collect = [ref for ref in collection if ref['geom'][s_lv] and ref['geom'][s_lv].intersects(param_box)]
        print('\n')
        for ne, element in enumerate(collect):
            print(ne)

            poly_color = (0.5 * np.random.random_sample(3)).tolist()
            # shape_guides = [l['guides'] for l in g['lines'] if l['lev'] == s_lv]
            #
            #
            # for g_ref in shape_guides:
            #     for g_sub_ref in g_ref:
            #         gd = [ax.plot(*gp['geom'].coords.xy, color=poly_color) for gp in g['guides'] if gp['id'] == g_sub_ref]

            def lplot(line):
                poly_color = (0.5 * np.random.random_sample(3)).tolist()

                # flatten_list = lambda irregular_list: [element for item in irregular_list for element in
                #                                        flatten_list(item)] if type(irregular_list) is list else [
                #     irregular_list]
                #
                #
                for lg in line['guides']:
                    #ince = line['guides'][0]
                    if lg is not None:
                        g = guides_collection[lg]

                        if g['geom'].type == 'LineString':
                            plt.plot(*g['geom'].coords.xy, color=poly_color, linewidth=2)
                        else:
                            ftl = utl.line_s_to_list(g['geom'])
                            for dline in ftl:
                                plt.plot(*dline.coords.xy, color="red", linewidth=4)
                            pass

                # # guides_flat = sum(line['guides'], []) if line['guides'][0] else None
                # # list(line['guides'])  #np.concatenate(line['guides']).flat)
                #
                # out = json.dumps(line, indent=2, cls=JsonSafeEncoder)
                # print(out)



                #print(guides_flat)

                plt.plot(*line['geom'].coords.xy, color=poly_color, linewidth=1)


            shape_lines = [lplot(l) for l in element['lines'] if l['lev'] == s_lv and l['shape'] == element['id']]

            p = param_box.intersection(element['geom'][s_lv])
            #print(p.type)
            #reli = [plot_polygon(ax, pl) for pl in line_s_to_list(p)]
            #reli = [plt.plot(*pl.exterior.coords.xy, color=(1.0, 0.0, 1.0, 0.9), linewidth=1) for pl in poly_s_to_list(p)]

            reli = [utl.plot_polygon(ax, pl) for pl in utl.poly_s_to_list(p)]

        # poly_color = ((0.5 * np.random.random_sample(3)).tolist(), len(ste))
        # #poly_color.append(0.8)

        return collect

    # outcoords = [list(i.coords) for i in inlines]
    #
    # [i for sublist in outcoords for i in sublist]
    res = get_mapped(bounds_box, map_collection, guides)

    # res = [plot_polygon(ax, p) for p in get_mapped(sector_box, map_collection)]

    # out = json.dumps(map_collection, indent=2, cls=JsonSafeEncoder)
    # print(out)

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # if "guides" in cfg["includes"]:
    #     print("Making Guides", bounds_box)
    #     guides, guides_json, shape_data_record = make_guides(med_poly,
    #                                                           bounds_box,
    #                                                           simps,
    #                                                           areas,
    #                                                           ax,
    #                                                           cfg["guide_buffer"],
    #                                                           cfg["guide_simp"])

    xt = np.arange(minx, maxx+1, 1.0)
    yt = np.arange(miny, maxy+1, 1.0)

    plt.xticks(xt)
    plt.yticks(yt)

    plt.grid()
    plt.tight_layout()
    plt.axis('equal')
    plt.show()
