#!/usr/bin/env python3

import numpy as np
from toolkit import path_tools as tool
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon, Point, box, LineString, LinearRing




def get_bathymetry():
    cfg = {'data_source_path': '../data'}
    depth_points = {
        "data": np.loadtxt(open(f"{cfg['data_source_path']}/bathy_med_tt.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
        "lons": np.loadtxt(open(f"{cfg['data_source_path']}/lon_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
        "lats": np.loadtxt(open(f"{cfg['data_source_path']}/lat_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0)
    }
    depth_points['origin'] = (depth_points['lons'][0], depth_points['lats'][0],)
    depth_points['density'] = 60
    return depth_points


if __name__ == '__main__':
    fig, ax = plt.subplots()

    # my_data = np.genfromtxt('./sample_data_summer_days_continent.csv', delimiter=',', names=True)
    my_data = np.genfromtxt('./sample_data_summer_days_corsica.csv', delimiter=',', names=True)
    print(my_data.shape, my_data.dtype, my_data)
    print(my_data[0]['lon'])


    M = np.genfromtxt('./contour20m_Corsica.csv', delimiter=',')

    #
    #
    med_marine_regions = gpd.read_file("../data/goas_med/med_local.shp")
    med_poly = med_marine_regions['geometry'][0]
    tool.plot_polygon(ax, med_poly, color=(0.0, 0.0, 0.0, 0.1), edgecolor='red')
    #
    # # # def get_contour(data, sigma, level, density, pos_origin):
    # bat = get_bathymetry()
    # contours = tool.get_contour(bat['data'], 0.0, 4, 0.0, bat['density'], bat['origin'])
    #
    # # get the longest one!
    # s_contours = sorted(contours, key=lambda d: len(d.coords))
    # s_contours.reverse()

    M = LinearRing(M)
    # M = s_contours[0]
    #
    total_points = my_data.shape[0]
    points = []
    distance_k = M.length
    skipped = []

    for i, cpc in enumerate(my_data):
        lp = Point(cpc['lon'], cpc['lat'])
        distance = M.project(lp)

        if distance > distance_k:
            skipped.append(i)

        distance_k = distance

        cp = M.interpolate(distance)
        points.append((cp.x, cp.y))

    points.reverse()
    print(skipped)


    def get_closest(p1, p2, all, indx):
        for p in range(indx-2, indx+2):
            point = Point(all[p])
            if p1.distance(point) < p1.distance(p2):
                return point, p

        return None, 0

    def pairs(lst):
        for i in range(1, len(lst)):
            yield lst[i - 1], lst[i]

    def progress():
        pointindex = 0
        segment = []
        segments_all = []

        all_pairs = pairs(list(M.coords))
        rlen = list(M.coords)
        buf = 0

        for n, pair in enumerate(all_pairs):
            A = Point(pair[0])
            B = Point(pair[1])
            f_point, buf = get_closest(A, B, points, pointindex)
            if buf != 0 and buf is not None:
                pointindex = buf

            print(f"\r{n}/{len(rlen)} {pointindex} complete", end='')

            if f_point is not None:
                segment.append((f_point.x, f_point.y))
                segments_all.append(LineString(segment))
                segment = []
                segment.append((f_point.x, f_point.y))
            else:
                segment.append(pair[0])

        for p in segments_all:
            plt.plot(*p.coords.xy, linewidth=4)



    aft = np.array(points)
    print(aft.shape)
    plt.scatter(aft[:, 0], aft[:, 1])  # , s=4.0, marker='o')


    # print(M.type)
    #
    #
    plt.plot(*M.coords.xy)

    # plt.plot(M[:, 0], M[:, 1])  # , s=4.0, marker='o')

    plt.scatter(my_data['lon'], my_data['lat'], s=1)

    plt.plot(my_data['lonmid'], my_data['latmid'], marker='+')

    plt.plot(my_data['lon'], my_data['lat'], marker='o')

    out = tool.do_refine(my_data['lon'], my_data['lat'], 1, 3)
    plt.plot(out[0], out[1], linewidth=1, color=(0.0, 0.0, 0.0))


    plt.grid()
    plt.tight_layout()
    plt.axis('equal')
    plt.show()
    pass
