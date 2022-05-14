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



if __name__ == '__main__':

    # d = 43.3125 - 43.270832
    # t = 5.0416665 - 5.0
    # c = -5.983333333 - -5.966666667
    #
    #
    # print(1/d, 1/t, 1/c)
    #
    # exit()


    depths = np.loadtxt(open("./data/bathy_med.csv", "rb"), delimiter=";", skiprows=1)
    print(depths.shape)

    depths_lons = np.loadtxt(open("./data/lon_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0)
    print(depths_lons.shape)

    depths_lats = np.loadtxt(open("./data/lat_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0)
    print(depths_lats.shape)

    temps = np.loadtxt(open("./data/thetao_t.csv", "rb"), delimiter=",", skiprows=0)
    print(temps.shape)

    temps_lons = np.loadtxt(open("./lon.csv", "rb"), delimiter=",", skiprows=0)
    print(temps_lons.shape)

    temps_lats = np.loadtxt(open("./lat.csv", "rb"), delimiter=",", skiprows=0)
    print(temps_lats.shape)

    temps_lats = np.flip(temps_lats)
    temps = np.flip(temps, 0)
    #//temps = np.flip(temps, 1)




    #42.895832,42.9375,42.979168,43.020832,43.0625,43.104168,43.145832,43.1875,43.229168,43.270832,43.3125
    #5.0,5.0416665,5.0833335,5.125,5.1666665,5.2083335,5.25,5.2916665,5.3333335,5.375,5.4166665,5.4583335,5.5,5.5416665,5.5833335,5.625,5.6666665,5.7083335,5.75,5.7916665,5.8333335,5.875,5.9166665,5.9583335,6.0,6.0416665,6.0833335,6.125,6.1666665,6.2083335,6.25,6.2916665,6.3333335,6.375,6.4166665,6.4583335,6.5,6.5416665,6.5833335

    kbx = [5, 44]
    # interval = 1/60
    #
    # fer = np.linspace(-6, -5, 61)
    #
    # for p in fer:
    #     result = np.where(lons > p)
    #     print(p, result)
    #
    # exit()





    # xrange = np.arange(0, lons.shape[0], 20)
    # xva = np.take(lons, xrange)
    # yrange = np.arange(0, lats.shape[0], 20)
    # yva = np.take(lats, yrange)
    # print(lons[0])
    #
    # kbx = [-5.92, 36]
    #
    # mask = (lons > kbx[0]) & (lons <= kbx[0]+1)
    # ft = lons[mask]
    # print(lons[mask], ft.shape)

    #fig, ax = plt.subplots(figsize=(10, 6))

    def normalize_val(val, mi, ma):
        return (val - mi) / (ma - mi)

    def get_data_from_sector(sector_tuple, sector_width, data, lons_ax, lats_ax, resolution, el_w):
        mask_lon = (lons_ax >= sector_tuple[0]) & (lons_ax <= sector_tuple[0] + sector_width)
        mask_lat = (lats_ax <= sector_tuple[1]) & (lats_ax >= sector_tuple[1] - sector_width)

        xe = [np.where(lons_ax >= i)[0][0] for n, i in enumerate(lons_ax[mask_lon]) if n % resolution == 0]
        ye = [np.where(lats_ax <= i)[0][0] for n, i in enumerate(lats_ax[mask_lat]) if n % resolution == 0]

        xmin = (np.amin(lons_ax[mask_lon])-sector_tuple[0])
        ymin = sector_tuple[1]-(np.amax(lats_ax[mask_lat]))

        print(np.amin(lons_ax[mask_lon]), np.amax(lats_ax[mask_lat]))
        print(xmin*el_w, ymin*el_w)

        xsi = np.sign(xmin)

        xxi = int(math.floor(abs(xmin) * (el_w/resolution))*xsi)  #int(math.floor(xmin * el_w))
        yyi = int(math.floor(ymin * (el_w/resolution)))

        print(xxi, yyi)

        print(xmin, ymin)

        data_sector = np.zeros([int((el_w*sector_width)/resolution), int((el_w*sector_width)/resolution)])

        print(data_sector.shape)
        #data_sector = np.zeros((len(xe), len(ye)))

        for iy, dy in enumerate(ye):
            for ix, dx in enumerate(xe):
                try:
                    data_sector[iy+yyi, ix+xxi] = data[dy, dx]  #normalize_val(data[dy, dx], 19, 20)
                except IndexError:
                    pass

        extents = [sector_tuple[0], sector_tuple[0]+sector_width, sector_tuple[1]-sector_width, sector_tuple[1]]
        return data_sector, extents

    #-5.983333333, 36.48333333, 30.01666667, 45.98333333

    plt.imshow(depths, extent=(-5.983333333, 36.48333333, 30.01666667, 45.98333333), interpolation="none", alpha=0.4, cmap="terrain")
    #plt.imshow(depths, extent=(-6, 36.5, 30, 46), interpolation="none", alpha=0.9, cmap="twilight")

    #sector = [5, 44]
    #sector = [-6.5, 46]
    sector = [-6.5, 44.5]
    depth, extents = get_data_from_sector(sector, 1.0, depths, depths_lons, depths_lats, 15, 60)
    #stemps, extents = get_data_from_sector(sector, 1.0, temps, temps_lons, temps_lats, 1, 24)

    #print(stemps)

    plt.imshow(depth, extent=extents, interpolation="none", alpha=0.5, cmap='gnuplot')  #, cmap="twilight")

    #
    #
    # #kbx = [36, 39]
    # mask_lon = (lons > kbx[0]) & (lons <= kbx[0] + 1)
    # mask_lat = (lats < kbx[1]) & (lats >= kbx[1] - 1)
    #
    # xe = [np.where(lons >= i)[0][0] for n, i in enumerate(lons[mask_lon]) if n % 2 == 0]
    # ye = [np.where(lats <= i)[0][0] for n, i in enumerate(lats[mask_lat]) if n % 2 == 0]
    #
    # depths_subset = np.zeros((len(xe), len(ye)))
    #
    # for iy, dy in enumerate(ye):
    #     for ix, dx in enumerate(xe):
    #         try:
    #             depths_subset[iy, ix] = depths[dy, dx]
    #         except IndexError:
    #             pass




    # plt.imshow(depths_subset, extent=(5, 6, 43, 44), interpolation="none", alpha=0.5)
    #plt.imshow(temps, extent=(5.0, 6.5833335, 42.895832, 43.3125), interpolation="none", alpha=0.5)

    k = box(-7, 29, 37, 47)
    plt.plot(*k.exterior.xy, color="blue")

    med_marine_regions = gpd.read_file("data/goas_med/med_local.shp")
    print(med_marine_regions.iloc[0])
    med_poly = med_marine_regions['geometry'][0]
    print(med_poly.bounds)

    plt.plot(*med_poly.exterior.xy, color="blue")
    # print(depths_subset)

    #Z1 = np.add.outer(range(8), range(8))



    #extent=(15, 16, 34, 35)
    # print('kbx', kbx, '=>', x, y, len(x), len(y))
    #plt.imshow(depths_subset, interpolation="none")
    plt.axis('scaled')
    plt.colorbar()
    plt.show()