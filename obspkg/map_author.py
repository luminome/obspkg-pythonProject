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
import obspkg_util as utl
import pickle




def preflight(bounds):
    geo_b = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[bounds])
    print(geo_b)
    map_regions = gpd.read_file("../data/GOaS_v1_20211214/goas_v01.shp")
    intersections2 = gpd.overlay(geo_b, map_regions, how='intersection')
    intersections2.to_file("../data/goas_med/med_zone.shp")
    return True



if __name__ == '__main__':
    fig, ax = plt.subplots()
    bounds = [-12, 28, 44, 48]
    minx, miny, maxx, maxy = bounds

    geographic_bounds = box(minx, miny, maxx, maxy)
    # preflight(geographic_bounds)

    # slice from 6 to 7 to get the med only!
    #                  name  latitude   longitude  min_Y      min_X      max_Y      max_X      area_km2
    # 6 Mediterranean Region  38.13065   19.70067   30.06809   -6.03255   47.37640   42.35496   2988248
    map_regions = gpd.read_file("../data/goas_med/med_zone.shp")  # , rows=slice(6, 7))
    print(map_regions.head())

    group = map_regions['geometry']

    map_multi_poly = geographic_bounds

    for g in group:
        for poly in utl.poly_s_to_list(g):
            map_multi_poly = map_multi_poly.difference(poly)
            # utl.plot_polygon(ax, poly, color=(0.5 * np.random.random_sample(3)).tolist())

    with open('../data/map_polygon', "wb") as poly_file:
        pickle.dump(map_multi_poly, poly_file, pickle.HIGHEST_PROTOCOL)
    #
    #
    # d = utl.poly_s_to_list(resr)
    # print(len(d))
    #
    # pg = [utl.plot_polygon(ax, sg, color='black') for sg in d]
    #
    # print(resr.type)



    xt = np.arange(minx, maxx+1, 1.0)
    yt = np.arange(miny, maxy+1, 1.0)

    plt.xticks(xt)
    plt.yticks(yt)

    plt.grid()
    plt.tight_layout()
    plt.axis('equal')
    plt.show()


    #med_poly = med_marine_regions['geometry'][0]
    #print(map_regions.head())
    #return
