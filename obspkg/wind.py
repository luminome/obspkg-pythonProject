#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap
from scipy import interpolate
from json_encoder import JsonSafeEncoder
import json
import netCDF4 as nc
import numpy as np
import geopandas as gpd
from scipy.ndimage.filters import gaussian_filter
# import skimage
from skimage import transform
import obspkg_util

# #// http://vincepota.com/plotting-satellite-data-with-python.html
# #// https://coastwatch.gitbook.io/satellite-course/tutorials/python-tutorial/1.-how-to-work-with-satellite-data-in-python

def k_to_deg_c(z):
    return z - 273.15


def get_data_scale_extents(data_packet):
    # d = data_packet['data']
    # print(d.shape)
    lo = data_packet['lons'][:]
    # print(lo[0], lo[-1], lo.shape)
    la = data_packet['lats'][:]
    # print(la[0], la[-1], la.shape)
    return lo[0], la[0], lo[-1], la[-1]


def plot_polygon(ax, poly, **kwargs):
    path = Path.make_compound_path(
        Path(np.asarray(poly.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])

    patch = PathPatch(path, **kwargs)
    collection = PatchCollection([patch], **kwargs)

    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection


if __name__ == '__main__':

    fig, ax = plt.subplots()

    med_marine_regions = gpd.read_file("../data/goas_med/med_local.shp")
    med_poly = med_marine_regions['geometry'][0]

    fn = "../data/download.nc"
    ds = nc.Dataset(fn)

    # dat = ds['analysed_sst'][0, :, :]
    # dat = np.flip(dat, axis=0)

    dpk = {
        "Udata": ds['u10'][0, :, :],
        "Vdata": ds['v10'][0, :, :],
        "lons": ds['longitude'],
        "lats": ds['latitude']  #, #np.flip(ds['latitude'], axis=0)
    }

    b = get_data_scale_extents(dpk)
    print('b', b)

    plot_polygon(ax, med_poly, color=(0.0, 0.0, 0.0, 0.1), edgecolor='red')

    M = np.hypot(dpk['Udata'], dpk['Vdata'])

    xva = dpk['lons'] #np.around(dpk['lons'], decimals=3)
    yva = dpk['lats'] #np.around(dpk['lats'], decimals=3)

    xpt = np.arange(xva[0], xva[-1], 0.25)
    ypt = np.arange(yva[0], yva[-1], 0.25)

    print(xpt, ypt)
    print(dpk['Udata'].shape)

    x, y = np.meshgrid(xva, yva)

    plt.quiver(x, y, dpk['Udata'], dpk['Vdata'], M, units='x', pivot='tip', width=0.022, scale=1 / 0.15)

    # xt = np.arange(xpt[0], xpt[-1], 1.0)
    # yt = np.arange(ypt[0], ypt[-1], 1.0)
    #
    # plt.xticks(xt)
    # plt.yticks(yt)

    plt.grid()

    #plt.colorbar()

    plt.tight_layout()
    plt.show()




    #     surface_temps = {
    #         "data": dat,
    #         "lons": ds['longitude'][:],
    #         "lats": ds['latitude'][:],
    #         "indx": [ds['time'][:], [0]],
    #         "reso": levels,
    #         "degs": 4,
    #         "origin": None
    #     }
    #
    #
    #
    #
    # print(json.dumps(elements, indent=1, cls=JsonSafeEncoder))
