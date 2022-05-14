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


# #// http://vincepota.com/plotting-satellite-data-with-python.html
# #// https://coastwatch.gitbook.io/satellite-course/tutorials/python-tutorial/1.-how-to-work-with-satellite-data-in-python

def k_to_deg_c(z):
    return z - 273.15


def get_data_scale_extents(data_packet):
    d = data_packet['data']
    # print(d.shape)
    lo = data_packet['lons'][:]
    # print(lo[0], lo[-1], lo.shape)
    la = data_packet['lats'][:]
    # print(la[0], la[-1], la.shape)
    return lo[0], la[0], lo[-1], la[-1]


def scale(im, nR, nC):
  nR0 = len(im)     # source number of rows
  nC0 = len(im[0])  # source number of columns
  return [[im[int(nR0 * r / nR)][int(nC0 * c / nC)]
             for c in range(nC)] for r in range(nR)]

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
    new_fn = "../data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1645377903582.nc"

    # g_layer = gpd.read_file(shapes_file, layer=layer_name)
    # print('converting crs to 4326')
    # g_layer = g_layer.to_crs(crs=4326)
    # #//news in seconds for time
    # #//have to restructure this for it to work: maybe reverse since dimensions:
    #     time = 31;
    #     lat = 318;
    #     lon = 867;
    ds = nc.Dataset(new_fn)

    #ds = ds.to_crs(crs=4326)


    # print(ds.variables)
    # nds = nc.Dataset(new_fn)
    dat = ds['analysed_sst'][0, :, :]
    #dat = np.nan_to_num(dat)
    #dat = np.where(np.isnan(dat), 32767, dat)
    #
    #dat = np.vectorize(k_to_deg_c)(dat[:, :, :])
    # #dat = np.where(np.isnan(dat), 20.0, dat)
    # dat = np.place(dat, -32767, 0.0)
    #dat = np.place(dat[0], dat[0] == np.nan, 0.0)

    dat = np.flip(dat, axis=0)
    #dat = np.vectorize(k_to_deg_c)(dat)

    # print(dat)
    #
    # # np.nan_to_num(dat[0])
    # epm = dat[0] # np.where(np.isnan(dat[0]), 20.0, dat[0])  #dat[0]  #gaussian_filter(dat[0], sigma=2)
    # # epm = gaussian_filter(dat[0], sigma=0.002)
    # # epm = np.where(np.isnan(epm), 32767, epm)
    # masked = epm.data
    # print(masked)
    # print(hasattr(masked, 'mask'))
    # mean = epm.mean(axis=1)
    # mean[mean.mask] = epm.fill_value
    #
    # print(mean)

    dpk = {
        "data": dat,
        "lons": ds['lon'],  #[:-1],  #[:],
        "lats": ds['lat']  #[:],
    }

    b = get_data_scale_extents(dpk)
    print('b', b)

    #print((((dpk['lons'][-1]-0.00625)-dpk['lons'][0]))*20)
    print((dpk['lons'][-1]-abs(dpk['lons'][0])))

    xva = np.around(dpk['lons'], decimals=3)
    yva = np.around(dpk['lats'], decimals=3)
    # xva, yva = (dpk['lons'], dpk['lats'])

    xpt = np.arange(xva[0], xva[-1], 0.05)  #[:-1]
    print(xpt.shape)
    ypt = np.arange(yva[0], yva[-1], 0.05)  #[:-1]
    print(ypt.shape)

    #

    gdx, gdy = np.meshgrid(xpt, ypt, indexing='ij')
    # plt.scatter(gdx, gdy, s=1)
    # plt.show()
    print("fails", (xva[-1]-(-7))/0.05)

    print(0.05/0.05004692)
    print(0.05004692*0.9995)

    print("float", "{:.4f}".format(dpk['lons'][0]))

    # cd = np.around(dpk['lons'], decimals=3)
    # print(cd)


    print(dpk['lats'].shape)
    print(dpk['lons'].shape)
    print(dpk['data'].shape)


    #exit()

    # print(dpk)
    plot_polygon(ax, med_poly, color=(0.0, 0.0, 0.0, 0.1), edgecolor='red')
    #exit()


    # # print(dat[0].shape)
    #
    # w, h = dat[0].shape
    #
    # #plt.gca().invert_yaxis()
    #
    #
    # masked = np.where(masked == -32767, 284.06039366914854, masked)
    # masked = transform.pyramid_expand(masked, upscale=4, sigma=None, preserve_range=True)
    # masked = gaussian_filter(masked, sigma=2.7)
    #print(epm.shape)

    print("A", dpk['data'].shape)
    # w, h = dpk['data'].shape
    # epm = scale(dpk['data'], w * 2, h * 2)
    # epm = np.array(epm)
    mt = np.mean(dpk['data'])

    epm = dpk['data'].data
    #print(epm)

    print("mean", mt)
    #epm = np.where(np.isnan(epm), 273.15, epm)
    epm = np.where(epm == -32768, mt, epm)
    epm = transform.pyramid_expand(epm, upscale=2, sigma=0, preserve_range=True)
    #epm = gaussian_filter(epm, sigma=0.7)
    print("B", epm.shape)

    #(-7.014798, 30.125, 36.325, 46.025)
    nmi = np.nanmin(dpk['data'])
    nma = np.nanmax(dpk['data'])
    levs = np.arange(nmi, nma, 0.05)

    jet = ["blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"]

    cm = LinearSegmentedColormap.from_list('my_jet', jet, N=len(levs))

    #plt.contourf(xpt, ypt, dpk['data'][:, :], levs, cmap=cm, origin="lower", extend='both')
    #plt.contourf(dpk['lons'][:], dpk['lats'][:], dpk['data'][:, :], levs, cmap=cm, origin="lower", extend='both')
    #(b[0], b[2], b[1], b[3])

    #ext = (xpt[0]-0.025, xpt[-1]+0.025, ypt[0]-0.025, ypt[-1]-0.025)  #(-7.015, b[2], b[1], b[3])  # (-20.0, 40.0, 29.0, 47.0)
    ext = (xpt[0], xpt[-1], ypt[0], ypt[-1])  #(-7.015, b[2], b[1], b[3])  # (-20.0, 40.0, 29.0, 47.0)
    #
    plt.imshow(epm, extent=ext, interpolation="None")

    xt = np.arange(xpt[0], xpt[-1], 1.0)
    yt = np.arange(ypt[0], ypt[-1], 1.0)[:-1]

    #yt = np.arange(np.round(b[1]), np.round(b[3]), 1.0)
    # xt = np.arange((b[0]), (b[2]), 0.5)
    # yt = np.arange((b[1]), (b[3]), 0.5)

    # print(xt)
    # print(yt)

    plt.xticks(xt)
    plt.yticks(yt)

    plt.grid()

    plt.colorbar()

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
