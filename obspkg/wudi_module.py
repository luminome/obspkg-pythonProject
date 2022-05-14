#!/usr/bin/env python3

import numpy as np
from toolkit import path_tools as tool
from shapely.geometry import box, Point, MultiPoint
from shapely.strtree import STRtree

from typing import List, Dict
import matplotlib.pyplot as plt
import config
from obspkg_util import coords_to_flat


def wudi_prepare(bnd: tuple):
    data_path = f'{config.data_path}/wudi/sample_data_summer_days_continent.csv'
    data_path2 = f'{config.data_path}/wudi/sample_data_summer_days_corsica.csv'
    wudi_med_data = np.genfromtxt(data_path, delimiter=',', names=True)
    wudi_data_corsica = np.genfromtxt(data_path2, delimiter=',', names=True)
    data_set = np.concatenate((wudi_med_data, wudi_data_corsica), axis=0)

    # #//pre-pack here: filter cardinal points and return points as dicts
    # #//a little absurd to force-in a point when nans are there but wth
    cloud = []
    for i in zip(data_set['lonmid'], data_set['latmid']):
        if not np.isnan(i[0]) and not np.isnan(i[1]):
            cloud.append(Point(i[0], i[1]))
        else:
            cloud.append(Point(-100, -100))

    wudi_maximums = {'U': np.nanmax(data_set['UPWdays']), 'D': np.nanmax(data_set['DNWdays'])}
    print(wudi_maximums)
    # cloud = [Point(i[0], i[1]) for i in zip(data_set['lonmid'], data_set['latmid']) if not np.isnan(i[0]) and not np.isnan(i[1])]

    points = MultiPoint(cloud)
    tree = STRtree(cloud)
    index_by_id = dict((id(pt), i) for i, pt in enumerate(cloud))

    bounds_box = box(bnd[0], bnd[1], bnd[2], bnd[3])

    collisions = bounds_box.intersection(points)

    these = [int(index_by_id[id(pt)]) for pt in tree.query(collisions)]
    these.sort()

    data_as_dicts = []

    for n, s in enumerate(these):
        A = Point(data_set['lon'][s], data_set['lat'][s])
        B = Point(data_set['lon'][s+1], data_set['lat'][s+1])
        M = Point(data_set['lonmid'][s], data_set['latmid'][s])
        V = {'U': data_set['UPWdays'][s], 'D': data_set['DNWdays'][s]}

        grp = {'a': A, 'b': B, 'm': M, 'v': V, 'i': s, 'lv': n % 3}
        data_as_dicts.append(grp)
        # print(s, M.is_empty)

    # #//THIS IS ALARMING, there is an empty value for Midpoint?
    scoped_cloud = [g['m'] for g in data_as_dicts if not g['m'].is_empty]
    scoped_tree = STRtree(scoped_cloud)
    scoped_points = MultiPoint(scoped_cloud)
    scoped_indices = dict((id(pt), i) for i, pt in enumerate(scoped_cloud))

    return data_as_dicts, scoped_tree, scoped_points, scoped_indices, wudi_maximums


def wudi_get_sector(bnd: tuple, tree: STRtree, points: MultiPoint, indices: Dict) -> List:
    bounds_box = box(bnd[0], bnd[1], bnd[2], bnd[3])
    if bounds_box.intersects(points):
        collisions = bounds_box.intersection(points)
        new_indices = [int(indices[id(pt)]) for pt in tree.query(collisions)]
        new_indices.sort()
    else:
        new_indices = None
    return new_indices


def wudi_filter_dict(element) -> Dict:
    def nano(val):
        return val if not np.isnan(val) else None

    for k, v in element.items():
        try:
            if v.type == 'Point':
                element[k] = coords_to_flat(v.coords)
        except AttributeError:
            if type(v) is list:
                element[k] = [nano(sv) for sv in v]
            pass
    return element


if __name__ == '__main__':

    fig, ax = plt.subplots()
    data_bounds = (-6, 34, 6, 38)
    wudi_data, w_tree, w_points, w_index = wudi_prepare(data_bounds)

    exit()
    #
    # def get_wudi(bnd: tuple, tree: STRtree, points: MultiPoint, indices: Dict) -> List:
    #     bounds_box = box(bnd[0], bnd[1], bnd[2], bnd[3])
    #     if bounds_box.intersects(points):
    #         collisions = bounds_box.intersection(points)
    #         new_indices = [int(indices[id(pt)]) for pt in tree.query(collisions)]
    #         new_indices.sort()
    #     else:
    #         new_indices = None
    #     return new_indices
    #

    test_bnd = (-4, 35, -3, 36)
    selection = wudi_get_sector(test_bnd, w_tree, w_points, w_index)



    for index in selection:

        bar = wudi_filter_dict(wudi_data[index])


        print(bar)
        # lin = LineString([bar['a'], bar['b']])
        # plt.plot(*lin.coords.xy, marker='o')
        # M = bar['m']
        # plt.scatter(*M.coords.xy, s=bar['v'][0], color='green')

    exit()

    #
    # for bar in selection:
    #     print(wudi_data[bar]['i'])
    #
    #
    # print(selection)

    # for n, bar in enumerate(wudi_data):
    #     print(n, bar['i'])

    #
    # for bar in wudi_data:
    #     print(bar['i'])
    #     lin = LineString([bar['a'], bar['b']])
    #     plt.plot(*lin.coords.xy, marker='o')
    #     M = bar['m']
    #     plt.scatter(*M.coords.xy, s=bar['v'][0], color='green')




    #
    # cloud = [Point(i[0], i[1]) for i in zip(wudi_data['lonmid'], wudi_data['latmid'])]
    # points = MultiPoint(cloud)
    # tree = STRtree(cloud)
    # index_by_id = dict((id(pt), i) for i, pt in enumerate(cloud))
    #
    # print(points.bounds)
    #
    # test_box = box(-5, 35, -3, 37)
    #
    # collisions = test_box.intersection(points)
    # print(collisions)
    #
    # these = [int(index_by_id[id(pt)]) for pt in tree.query(collisions)]
    # these.sort()
    #
    # for s in these:
    #     A = Point(wudi_data['lon'][s], wudi_data['lat'][s])
    #     B = Point(wudi_data['lon'][s+1], wudi_data['lat'][s+1])
    #     M = Point(wudi_data['lonmid'][s], wudi_data['latmid'][s])
    #     lin = LineString([A, B])
    #     plt.plot(*lin.coords.xy, marker='o', color='red')
    #     plt.scatter(*M.coords.xy, color='green')
    #     print(A, B, s)
    #
    #
    #
    # print(these)
    #
    test_box = box(main_bnd[0], main_bnd[1], main_bnd[2], main_bnd[3])
    tool.plot_polygon(ax, test_box, color=(0, 0, 0, 0.2))

    # index_by_id = dict((id(pt), i) for i, pt in enumerate(location.geoms))
    #
    # points = [Point]
    #
    #
    #
    # [(index_by_id[id(pt)], pt.wkt) for pt in tree.query(Point(2, 2).buffer(1.0))]



    # print(wudi_data.shape, wudi_data.dtype, wudi_data)
    #
    # print(wudi_data[0]['lon'])

    # plt.plot(wudi_data['lon'], wudi_data['lat'], marker='o', color='blue')
    #
    # plt.scatter(wudi_data['lonmid'], wudi_data['latmid'])
    # out = tool.do_refine(my_data['lon'], my_data['lat'], 1, 3)
    # plt.plot(out[0], out[1], linewidth=1, color=(0.0, 0.0, 0.0))

    plt.grid()
    plt.tight_layout()
    plt.axis('equal')
    plt.show()

"""
sac@luminome wudi % nano ~/.bash_profile
export PYTHONPATH="${PYTHONPATH}:/Users/sac/Sites/pythonProject"
sac@luminome wudi % source ~/.bash_profile
"""