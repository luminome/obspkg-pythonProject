#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter

from scipy import interpolate
from scipy.interpolate import griddata
from json_encoder import JsonSafeEncoder
import json
import netCDF4 as nc
import numpy as np
import geopandas as gpd
from scipy.ndimage.filters import gaussian_filter
# import skimage
from skimage import transform
import obspkg_util

import pyproj
from pyproj import Transformer

import pickle
# https://all-geo.org/volcan01010/2012/11/change-coordinates-with-pyproj/


from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import resample_nearest


wgs84o = pyproj.CRS("EPSG:3857")  # (EPSG:41001:3857) WGS 84 (EPSG 4326)
wgs84 = pyproj.CRS("EPSG:4326")

def interp_b_to_a(a, b):
    def_a = SwathDefinition(lons=a['lons'], lats=a['lats'])
    def_b = SwathDefinition(lons=b['lons'], lats=b['lats'])

    interp_dat = resample_nearest(def_a, a['data'], def_b, radius_of_influence=5500)
    new_b = {'data': interp_dat,
             'lats': b['lats'],
             'lons': b['lons']
            }
    return new_b


# #// http://vincepota.com/plotting-satellite-data-with-python.html
# #// https://coastwatch.gitbook.io/satellite-course/tutorials/
#     python-tutorial/1.-how-to-work-with-satellite-data-in-python
# #// https://stackoverflow.com/questions/58079075/numpy-select-rows-based-on-condition


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

    with open('/Users/sac/Sites/PythonProject/data/map_polygon', "rb") as poly_file:
        med_poly = pickle.load(poly_file)

    new_fn = "../data/download-2.nc"
    ds = nc.Dataset(new_fn)

    """
    :geospatial_lat_min = 30.25; // double
    :geospatial_lat_max = 46.0; // double
    :geospatial_lon_min = -18.125; // double
    :geospatial_lon_max = 36.25; // double
    """

    r = {
        "data": ds['sst'][0, :, :],
        "lons": ds['longitude'][:],
        "lats": ds['latitude'][:],
        "indx": ds['time'][:],
        "U": ds['u10'],
        "V": ds['v10'],
        "degs": 0.25
    }

    # r['data'] = np.flip(r['data'], axis=0)

    csp = 0.125
    csc = 0.0625
    nlim = (np.min(r['lons'])-csc, np.max(r['lons'])+(3*csc), np.min(r['lats'])-csc, np.max(r['lats'])+(3*csc))
    print(nlim)

    mt = np.mean(r['data'])
    print(mt)
    print(ds['sst']._FillValue)

    epm = r['data']
    epm[epm.data == ds['sst']._FillValue] = mt
    #
    # epm = np.where(epm == -32767, mt, epm)

    print(epm.shape)
    print(r['lats'].shape)
    print(r['lons'].shape)
    epm = transform.pyramid_expand(epm, upscale=2, sigma=0, preserve_range=True)
    epm = np.array(epm)

    print(epm.shape)

    # nmi = np.nanmin(epm)
    # nma = np.nanmax(epm)
    # levs = np.arange(nmi, nma, 0.005)
    # jet = ["blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"]
    # cm = LinearSegmentedColormap.from_list('my_jet', jet, N=len(levs))

    plt.imshow(epm, extent=nlim, interpolation="None")
    plt.colorbar()

    # for poly in med_poly.geoms:
    #     plot_polygon(ax, poly, color=(0.0, 0.0, 0.0, 0.1), edgecolor='red')

    gdx, gdy = np.meshgrid(r['lons'], r['lats'], indexing='ij')
    plt.scatter(gdx, gdy, s=1.0, color="black")

    plt.xticks(r['lons'][::4])
    plt.yticks(r['lats'][::4])


    plt.grid()
    plt.tight_layout()
    plt.axis('equal')
    plt.show()
    exit()






















    # new_fn = "../data/download.nc"
    # # //new_fn = "../data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1645013561954.nc"
    # # #//news in seconds for time
    # # #//have to restructure this for it to work: maybe reverse since dimensions:
    # #     time = 31;
    # #     lat = 318;
    # #     lon = 867;
    # # ds = nc.Dataset(new_fn)
    # # new_fn = "../data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1645377903582.nc"
    # ds = nc.Dataset(new_fn)
    # print(ds)
    # # exit()
    # r = {
    #     "data": np.loadtxt(open("../data/bathy_med_tt.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
    #     "lons": np.loadtxt(open("../data/lon_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
    #     "lats": np.loadtxt(open("../data/lat_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
    # }
    # :northernmost_latitude = 46.025
    # f; // float
    # :southernmost_latitude = 30.125
    # f; // float
    # :easternmost_longitude = 36.325
    # f; // float
    # :westernmost_longitude = -18.125
    # f; // float
    # "geographical coordinates, WGS84 projection"; WGS84 / Simple Mercator (EPSG:41001) WGS 84 (EPSG 4326)
    # :geospatial_lat_min = 30.125; // double
    # :geospatial_lat_max = 46.025001525878906; // double
    # :geospatial_lon_min = -7.014798164367676; // double
    # :geospatial_lon_max = 36.32500076293945; // double
    # https://resources.marine.copernicus.eu/product-detail/SST_MED_SST_L4_REP_OBSERVATIONS_010_021/DOCUMENTATION
    new_fn = "../data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1645377903582.nc"
    ds = nc.Dataset(new_fn)
    dat = ds['analysed_sst'][0, :, :]

    #spe = pyproj.transform(wgs84, wgs84o, ds['lat'][0], ds['lon'][0])

    transformer = Transformer.from_crs("epsg:3857", "epsg:4326")
    spe = transformer.transform(ds['lat'][0], ds['lon'][0])

    print(ds['lat'][0], ds['lon'][0])
    print(spe)
    exit()

    r = {
        "data": dat,
        "lons": ds['lon'],
        "lats": ds['lat']
    }

    for e in r:
        print(e, r[e].shape)

    rtx = np.linspace(-7.015, 36.325, r['lons'].shape[0])  #, endpoint=True)  #, endpoint=True)
    rty = np.linspace(30.125, 46.025, r['lats'].shape[0])  #, endpoint=True)

    rfp = {
        "data": dat[:318, :318],
        "lons": ds['lon'][:318],
        "lats": ds['lat'][:318]
    }

    pr = {
        "lons": rtx[:318],
        "lats": rty[:318]
    }

    print(rfp['lons'].shape, rfp['lats'].shape)

    pek = interp_b_to_a(rfp, pr)

    plt.imshow(pek['data'], interpolation="None")

    print(pek['data'])

    plt.grid()
    plt.axis('equal')
    plt.tight_layout()
    #//plt.axis('equal')
    plt.show()
    exit()

    exit()
    # lat_new, lon_new = np.meshgrid(rtx, rty)
    #
    # new_grid = griddata((r['lats'], r['lons']), r['data'][:], (lat_new, lon_new), method='linear')
    #
    # print(new_grid.shape)

    # print(rtx.shape, rtx)
    # print(rty.shape, rty)
    # yva = np.format_float_positional(r['lats'], precision=4)
    # yva = ["%.4f" % i for i in r['lats']]
    # print(yva)
    #
    #
    # exit()

    gdx, gdy = np.meshgrid(pr['lons'], pr['lats'], indexing='ij')
    plt.scatter(gdx, gdy, s=2.0, color="blue")


    yva = r['lats']  #np.around(r['lats'], decimals=2)
    xva = r['lons']  #np.around(r['lons'], decimals=2)

    agdx, agdy = np.meshgrid(xva, yva, indexing='ij')
    plt.scatter(agdx, agdy, s=2.0, color="red")

    # xt = np.arange(xpt[0], xpt[-1], 1.0)
    # yt = np.arange(ypt[0], ypt[-1], 1.0)  #[:-1]
    #ax.set_ylim(start, stop)

    # plt.xticks(xpt)
    # plt.yticks(ypt)

    plt.grid()
    plt.axis('equal')
    plt.tight_layout()
    #//plt.axis('equal')
    plt.show()
    exit()




    print(yva)

    #
    exit()
    # dat = ds['analysed_sst'][0, :, :]
    # dat = np.flip(dat, axis=0)
    # U = ds['u10']
    # V = ds['v10']
    # a = np.array([U], dtype=np.str)
    # b = np.array([V], dtype=np.str)
    # c = np.char.add(a, ' ').astype(str)
    # dat = np.char.add(c, b).astype(str)[0]


    # new_fn = "../data/download-2.nc"
    # ds = nc.Dataset(new_fn)
    #
    # r = {
    #     "data": ds['sst'][0, :, :],
    #     "lons": ds['longitude'][:],
    #     "lats": ds['latitude'][:],
    #     "indx": ds['time'][:],
    #     "U": ds['u10'],
    #     "V": ds['v10'],
    #     "degs": 0.25
    # }

    nmi = np.nanmin(r['data'])
    nma = np.nanmax(r['data'])
    levs = np.arange(nmi, nma, 0.05)
    jet = ["blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"]
    cm = LinearSegmentedColormap.from_list('my_jet', jet, N=len(levs))

    xva = np.around(r['lons'], decimals=3)
    yva = np.around(r['lats'], decimals=3)

    print(yva[-1])

    # c = 0.125
    # ext = (-18-c, 37+c, 29-c, 46+c)
    # plt.imshow(r['data'], cmap=cm, extent=ext, interpolation="None")
    # plt.imshow(r['data'], cmap=cm, interpolation="None")


    xpt = np.linspace(xva[0], xva[-1], num=len(xva), endpoint=True)
    ypt = np.linspace(yva[0], yva[-1], num=len(yva), endpoint=True)


    # xpt = np.arange(xva[0], xva[-1], r['degs'])
    # ypt = np.arange(yva[0], yva[-1], -r['degs'])

    gdx, gdy = np.meshgrid(r['lons'], r['lats'], indexing='ij')
    plt.scatter(gdx, gdy, s=2.0, color="blue")

    agdx, agdy = np.meshgrid(rtx, rty, indexing='ij')
    plt.scatter(agdx, agdy, s=2.0, color="red")
    # xt = np.arange(xpt[0], xpt[-1], 1.0)
    # yt = np.arange(ypt[0], ypt[-1], 1.0)  #[:-1]
    print(xpt)
    print(ypt)
    #ax.set_ylim(start, stop)

    # plt.xticks(xpt)
    # plt.yticks(ypt)

    plt.grid()
    # //plt.colorbar()
    plt.tight_layout()
    #//plt.axis('equal')
    plt.show()
    exit()






    def tick_func(x, pos):
        # 'The two args are the value and tick position'
        if pos is not None:
            tick_locs = ax.yaxis.get_majorticklocs()  # Get the list of all tick locations

            str_tl = str(tick_locs).split()[1:-1]  # convert the numbers to list of strings
            # print(str_tl)

            p = max(len(i) - i.find('.') - 1 for i in str_tl)  # calculate the maximum number of non zero digit after "."
            # p = max(1, p)  # make sure that at least one zero after the "." is displayed
            # print(p)
            # if p == 0:
            try:
                spk = r['lats'][int(x)]
                return "%.2f" % round(spk, 2)  # "{:.2}f".format(spk)
            except IndexError:
                pass
            # else:
            #
            # try:
            #     spk = r['lats'][int(x)]
            #     # print(spk)
            #     return "pos:{0}/x:{1:1.{2}f}".format(pos, spk, p)
            # except IndexError:
            #     return None

                # spk = x

            #p = max(r['data'][i, :] for i in str_tl)



    formatter = FuncFormatter(tick_func)





    for el in r:
        print(r[el].shape)



    # //plt.imshow(r['U'][0, :, :], cmap=cm, interpolation="None")
    plt.imshow(r['data'], cmap=cm, interpolation="None")

    # xt = np.arange(r['lons'][0], r['lons'][-1], 1.0).tolist()
    # yt = np.arange(r['lats'][0], r['lats'][-1], -1.0).tolist() #[:-1]

    # sector_tuple = (-6, 42)
    # sector_width = 16.0
    # ye = np.where((r['lats'] <= sector_tuple[1]) & (r['lats'] > (sector_tuple[1] - sector_width)))[0]  # indices
    # xe = np.where((r['lons'] >= sector_tuple[0]) & (r['lons'] < (sector_tuple[0] + sector_width)))[0]
    #
    # #spec = r['U'][0, :, :]
    # spec = r['data']
    #
    # x_does_work = spec[ye[:, np.newaxis], xe]
    #
    # ext = (sector_tuple[0],
    #        sector_tuple[0] + sector_width,
    #        sector_tuple[1] - sector_width,
    #        sector_tuple[1])
    #
    # plt.imshow(x_does_work, cmap=cm, extent=ext, interpolation="None")
    #
    #
    # print(x_does_work, x_does_work.size, x_does_work.shape)
    # #//exit()
    #
    #


    # xva = np.around(r['lons'], decimals=3)
    # yva = np.around(r['lats'], decimals=3)
    #
    # xpt = np.arange(xva[0], xva[-1], 0.25)
    # ypt = np.arange(yva[0], yva[-1], -0.25)
    #
    # print(yva)
    #
    # gdx, gdy = np.meshgrid(xpt, ypt, indexing='ij')
    # plt.scatter(gdx, gdy, s=2.0, color="blue")
    #
    # xt = np.arange(xpt[0], xpt[-1], 1.0)
    # yt = np.arange(ypt[0], ypt[-1], -1.0)  #[:-1]
    #
    # plt.xticks(xt)
    # plt.yticks(yt)

    #
    # lax.set_xticks(xt)
    # lax.set_yticks(yt)
    # plt.xticks(xt)
    # plt.yticks(yt)

    # plt.yticks(r['lats'])
    # lax = plt.axes()
    # lax.set_yticklabels(r['lats'][:])

    ax.yaxis.set_major_formatter(formatter)

    plt.grid()
    # //plt.colorbar()
    plt.tight_layout()
    #//plt.axis('equal')
    plt.show()
    exit()

    # new_fn = "../data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1645377903582.nc"
    # ds = nc.Dataset(new_fn)
    # dat = ds['analysed_sst'][0, :, :]
    # r = {
    #     "data": dat,
    #     "lons": np.array(ds['lon']),
    #     "lats": np.array(ds['lat'])
    # }

    # r = {
    #     "data": dat,
    #     "lons": np.array(ds['lon']),
    #     "lats": np.array(ds['lat'])
    # }

    # time, lat, lon
    # r = {
    #     "data": np.loadtxt(open("../data/bathy_med_tt.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
    #     "lons": np.loadtxt(open("../data/lon_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
    #     "lats": np.loadtxt(open("../data/lat_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
    # }

    sector_tuple = (5, 34)
    sector_width = 1.0
    ye = np.where((src['lats'] <= sector_tuple[1]) & (src['lats'] > (sector_tuple[1] - sector_width)))[0]  # indices
    xe = np.where((src['lons'] >= sector_tuple[0]) & (src['lons'] < (sector_tuple[0] + sector_width)))[0]
    # #x_indexed = r['data'][(ye, xe)]
    # print(r['data'].shape)
    # print(len(ye))
    # print(xe.shape)
    #
    #
    # print(r['data'][719][119])
    # print(r['lats'][719],r['lats'].shape)
    # print(r['lons'][119],r['lons'].shape)

    # x_does_work = r['data'][ye[0]:ye[-1], xe[0]:xe[-1]]
    # ye[0]:ye[-1], xe[0]:xe[-1]
    # x_does_work = r['data'][np.ix_(ye, xe)]
    # x_indexed = r['data'][xe][ye]
    # x_does_work = r['data'][xe][:, ye]
    # x_does_work = r['data'][:, xe][ye, :]
    # x_does_work = r['data'][ye[:, np.newaxis], xe]
    # x_does_work = r['data'][np.ix_(ye, xe)]

    x_does_work = r['data'][ye[:, np.newaxis], xe]
    print(x_does_work, x_does_work.size, x_does_work.shape)

    rsd = x_does_work[::30, ::30]

    # rsd = np.resize(x_does_work, (4, 4))

    print(rsd.shape, rsd)

    exit()
    #
    # condition = np.mod(x_does_work, 4) == 0
    # print(condition.shape)
    #
    # x_does_mod = x_does_work[condition]
    # print(x_does_mod.size, x_does_mod.shape)
    #
    # print(x_does_work, x_does_work.size, x_does_work.shape)


    exit()

    #
    # #ref = np.column_stack((r['lats'][:, np.newaxis], r['lons']), axis=1)
    #
    # #ref = np.concatenate((r['lats'][:, np.newaxis], r['data']), axis=1)
    #
    # # yx = zip(r['lats'], r['lons'])
    # # for i, coord in enumerate(yx):
    # #
    # #     print(i, coord)
    #
    # print(ref.shape, ref[0][0], r['data'][0][0])
    #
    # # x = [0, 0, 1, 1, 2, 2]
    # # y = [1, 2, 0, 1, 1, 2]
    # # z = [14, 17, 15, 16, 18, 13]
    # #
    # # print(r['data'].shape)
    # #
    # # arr = np.zeros(r['data'].shape)
    # # yx = zip(r['lats'], r['lons'])
    # #
    # # for i, coord in enumerate(yx):
    # #     arr[coord] = r['data'][i]
    #
    #
    #
    #
    exit()



    fig, ax = plt.subplots()

    med_marine_regions = gpd.read_file("../data/goas_med/med_local.shp")
    med_poly = med_marine_regions['geometry'][0]

    fn = "../data/download.nc"
    new_fn = "../data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1645377903582.nc"
    ds = nc.Dataset(new_fn)

    dat = ds['analysed_sst'][:, :, :]
    dat = np.flip(dat, axis=0)

    sector_tuple = (-4, 34)
    sector_width = 1.0

    print(dat.shape)
    print(dat)
    # print(dat)
    # xe = np.where(dat[0] >= sector_tuple[0])  #;//, dat[0] < sector_tuple[0] + sector_width, dat)  #  & dat[1] < sector_tuple[0] + sector_width)]
    #
    # print(xe)
    exit()

    #
    # xe = [np.where(dat[1] >= sector_tuple[0] & dat[1] < sector_tuple[0] + sector_width]
    #
    # mask_lon = (lons_ax >= sector_tuple[0]) & (lons_ax < sector_tuple[0] + sector_width)
    # mask_lat = (lats_ax <= sector_tuple[1]) & (lats_ax > sector_tuple[1] - sector_width)
    #
    #
    # # try:
    # xe = [np.where(lons_ax >= i)[0][0] for n, i in enumerate(lons_ax[mask_lon]) if n % resolution == 0]
    # ye = [np.where(lats_ax <= i)[0][0] for n, i in enumerate(lats_ax[mask_lat]) if n % resolution == 0]
    #





    dpk = {
        "data": dat,
        "lons": ds['lon'],
        "lats": ds['lat']
    }

    b = get_data_scale_extents(dpk)
    print('b', b)
    print(dat.shape)

    xva = np.around(dpk['lons'], decimals=3)
    yva = np.around(dpk['lats'], decimals=3)
    xpt = np.arange(xva[0], xva[-1], 0.05)
    ypt = np.arange(yva[0], yva[-1], 0.05)

    gdx, gdy = np.meshgrid(xpt, ypt, indexing='ij')
    plt.scatter(gdx, gdy, s=0.5)

    plot_polygon(ax, med_poly, color=(0.0, 0.0, 0.0, 0.1), edgecolor='red')

    mt = np.mean(dpk['data'])
    epm = dpk['data'].data
    epm = np.where(epm == -32768, mt, epm)
    epm = transform.pyramid_expand(epm, upscale=2, sigma=0, preserve_range=True)
    print("B", epm.shape)

    nmi = np.nanmin(epm)
    nma = np.nanmax(epm)
    levs = np.arange(nmi, nma, 0.05)
    jet = ["blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"]
    cm = LinearSegmentedColormap.from_list('my_jet', jet, N=len(levs))

    #plt.contourf(xpt, ypt, dpk['data'][:, :], levs, cmap=cm, origin="lower", extend='both')
    #plt.contourf(dpk['lons'][:], dpk['lats'][:], dpk['data'][:, :], levs, cmap=cm, origin="lower", extend='both')
    #(b[0], b[2], b[1], b[3])

    #ext = (xpt[0]-0.025, xpt[-1]+0.025, ypt[0]-0.025, ypt[-1]-0.025)  #(-7.015, b[2], b[1], b[3])  # (-20.0, 40.0, 29.0, 47.0)
    ext = (xpt[0]-0.025, xpt[-1]+0.025, ypt[0], ypt[-1])  #(-7.015, b[2], b[1], b[3])  # (-20.0, 40.0, 29.0, 47.0)
    #ext = (xpt[0]-0.025, xpt[-1]+0.025, ypt[0], ypt[-1])  #(-7.015, b[2], b[1], b[3])  # (-20.0, 40.0, 29.0, 47.0)
    #
    plt.imshow(epm, extent=ext, cmap=cm, interpolation="None")

    # interval = data_packet['degs']
    # data = data_packet['data']
    # resolution = data_packet['reso'][level]
    # lons_ax = data_packet['lons']
    # lats_ax = data_packet['lats']
    # indices = data_packet['indx']

    xptc = np.arange(xva[0], xva[-1], 0.025)
    yptc = np.arange(yva[0], yva[-1], 0.025)

    spk = {
        "data": epm,
        "lons": xptc,
        "lats": np.flip(yptc, axis=0),
        "indx": [[0], [0]],
        "degs": 1,
        "reso": [1]
    }


    f_set = obspkg_util.get_data_from_sector((5,42), 1, spk, 0)
    print(f_set)
    print(f_set['data'].shape)

    #CREATE POLIES HERE

    xt = np.arange(xpt[0], xpt[-1], 1.0)
    yt = np.arange(ypt[0], ypt[-1], 1.0)[:-1]

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
