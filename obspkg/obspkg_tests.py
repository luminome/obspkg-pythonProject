#!/usr/bin/env python3
from obspkg_util import new_data_from_sector
import numpy as np
import netCDF4 as nc
import pickle
from typing import List
from shapely.geometry import Polygon, MultiPolygon


# #//https://nctoolkit.readthedocs.io/en/latest/introduction.html

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

from skimage import transform

# from scipy.ndimage.interpolation import shift
import scipy.ndimage

import pyproj
from pyproj import Transformer


def k_to_deg_c(z):
    return float("%.30f" % round(z - 273.15, 3))  # float('%0.2f' % z - 273.15)


def show(image_data):

    plt.imshow(image_data, interpolation="None", alpha=1)
    plt.grid()
    plt.tight_layout()
    plt.axis('equal')
    plt.show()


def mangle(s):
    x, y = np.meshgrid(s['lons'], s['lats'])
    M = np.hypot(s['U'][0], s['V'][0])
    plt.quiver(x, y, s['U'][0], s['V'][0], M, units='x', width=0.022, scale=1 / 0.15)  #, [C], **kw) pivot='tip',
    plt.grid()
    plt.tight_layout()
    plt.axis('equal')
    plt.show()


def display(image, original=None, range_tuple=(12, 24)):

    fig, axes = plt.subplots(nrows=1, ncols=5)
    cmap = cm.get_cmap('viridis')
    normalizer = Normalize(range_tuple[0], range_tuple[1])
    im = cm.ScalarMappable(norm=normalizer)

    fgroup = [
        image,
        original[1::2, 1::2],
        original[2::4, 2::4],
        original[5::10, 5::10],
        original
    ]


    for i, ax in enumerate(axes.flat):
        ax.imshow(fgroup[i], interpolation="None",cmap=cmap, norm=normalizer)
        ax.set_title(str(fgroup[i].shape))

    fig.colorbar(im, ax=axes.ravel().tolist())
    plt.show()
    #
    # fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5)
    #
    # axlist = [ax1, ax2, ax3, ax4, ax5]
    # fig.suptitle('Horizontally stacked subplots')
    # ax1.grid()
    #
    # first = ax1.imshow(image, interpolation="None", alpha=1)
    #
    # subset = image[1::2, 1::2]
    # ax2.imshow(subset, interpolation="None", alpha=1)  #extent=scope,
    #
    # subset_2 = image[2::4, 2::4]
    # print(subset_2)
    # ax3.imshow(subset_2, interpolation="None", alpha=1)  #extent=scope,
    #
    # subset_3 = image[4::4, 4::4]
    # print(subset_3)
    # ax4.imshow(subset_3, interpolation="None", alpha=1)  #extent=scope,
    #
    # ax5.imshow(original, interpolation="None", alpha=1)
    #
    # #fig.colorbar()
    # fig.colorbar(first, ax=axlist)
    #
    # plt.tight_layout()
    # plt.axis('equal')

    plt.show()


def prepare_mini_temps():
    fn = '../data/1999-6-10_med-cmcc-tem-rean-d_1638647808820.nc'
    ds = nc.Dataset(fn)
    original_dat = ds['thetao']
    dat = np.array(original_dat)
    print(dat.shape)
    # flip the lats here bc inverted.
    dat = np.flip(dat, axis=2)
    dat[dat == original_dat._FillValue] = np.nan
    meantemp = np.nanmean(dat)

    r = {
        "data": dat,
        "time": np.array(ds['time']),
        "lons": np.array(ds['lon']),
        "lats": np.array(ds['lat'])
    }

    print(r['lats'])
    print(dat.shape)
    concrete = new_data_from_sector(r, (5, 44), 2.0, 1)
    print(concrete.size, concrete.shape)

    show(concrete[0, 0, :, :])

    #print(concrete)


def prepare_sst_and_wind(do_resample=None):
    data_path_fn = "../data/download-2.nc"
    ds = nc.Dataset(data_path_fn)
    # sanitize masked data:
    original_dat = ds['sst']
    dat = np.array(original_dat)
    dat[dat == original_dat._FillValue] = np.nan
    mean_temp = np.nanmean(dat)
    dat[np.isnan(dat)] = mean_temp

    r = {
        "data": dat,
        "time": np.array(ds['time']),
        "lons": np.array(ds['longitude']),
        "lats": np.array(ds['latitude']),
        "U": np.array(ds['u10']),
        "V": np.array(ds['v10']),
    }

    if do_resample is not None:
        r['data'] = transform.pyramid_expand(r['data'], upscale=2, sigma=0.1, preserve_range=True)[::2, :-1, :-1]
        r['U'] = transform.pyramid_expand(r['U'], upscale=2, sigma=0.1, preserve_range=True)[::2, :-1, :-1]
        r['V'] = transform.pyramid_expand(r['V'], upscale=2, sigma=0.1, preserve_range=True)[::2, :-1, :-1]
        r['lats'] = np.linspace(r['lats'][0], r['lats'][-1], r['data'].shape[1])  #-1)
        r['lons'] = np.linspace(r['lons'][0], r['lons'][-1], r['data'].shape[2])  #-1)

    r['data'] = np.vectorize(k_to_deg_c)(r['data'])
    return r


"""
prepare sectorized data selon type
    1) type one is easy and not imaginary lon/lat
    2) has extents, and uneven lon/lat

with the (new) scoped sector, run additional transforms on it.

"""


def get_hi_res_sst_data():
    # analysed_sst(time=31, lat=318, lon=1089);
    #//data_path_fn = "../data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1649679937092.nc"  #(one month)
    data_path_fn = "../data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1650708192445.nc"  #(three months)
    print("opening source:", data_path_fn)

    dds = nc.Dataset(data_path_fn)
    original_dat = dds['analysed_sst']
    dat = np.array(original_dat)
    dat[dat == original_dat._FillValue] = np.nan
    mean_temp = np.nanmean(dat)
    dat[np.isnan(dat)] = mean_temp

    dat = np.flip(dat, axis=1)

    r = {
        "data": dat,
        "time": dds['time'][:],
        "lons": np.linspace(-18.1, 36.3, num=dds['lon'].shape[0], endpoint=True),
        "lats": np.linspace(46.0, 30.15, num=dds['lat'].shape[0], endpoint=True),
        "degs": 20,
        "mean": k_to_deg_c(mean_temp)
    }

    return r


def get_lo_res_wind_data():
    data_path_fn = "../data/download-2.nc"
    ds = nc.Dataset(data_path_fn)
    print("opening source:", data_path_fn)

    r = {
        "u10": np.array(ds['u10']),
        "v10": np.array(ds['v10']),
        "time": np.array(ds['time']),
        "lons": np.array(ds['longitude']),
        "lats": np.array(ds['latitude']),
        "degs": 4
    }

    return r





def prepare_hi_res_sst():
    with open('/Users/sac/Sites/PythonProject/data/map_polygon', "rb") as poly_file:
        med_poly = pickle.load(poly_file)

    #data_path_fn = "../data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1645013561954.nc"
    data_path_fn = "../data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1649679937092.nc"
    # #// https://resources.marine.copernicus.eu/product-detail/SST_MED_SST_L4_REP_OBSERVATIONS_010_021/INFORMATION

    ds = nc.Dataset(data_path_fn)
    original_dat = ds['analysed_sst']

    fill_value = original_dat._FillValue  #dat[dat.data == ds['sst']._FillValue] = mt
    dat = np.array(original_dat)
    dat = np.flip(dat, axis=1)

    r = {
        "data": dat,
        "time": np.array(ds['time']),
        "lons": np.array(ds['lon']),
        "lats": np.array(ds['lat']),
    }

    print(original_dat.shape)
    print(dat.shape)
    print(r['lons'].shape)
    #exit()

    interval = 1/20
    """
    :valid_min = -20.0
    f; // float
    :valid_max = 40.0
    f; // float
    
    :valid_min = 29.0f; // float
    :valid_max = 47.0f; // float
    """
    # transformer = Transformer.from_crs("epsg:4326", "epsg:4326")  #3857
    # spe = transformer.transform(r['lons'][0], r['lats'][0])
    # print(r['lats'][0], r['lons'][0])
    # print(spe)
    # exit()
    """
    :northernmost_latitude = 46.025f; // float
    :southernmost_latitude = 30.125f; // float
    :easternmost_longitude = 36.325f; // float
    :westernmost_longitude = -18.125f; // float
    """

    #latref = np.linspace(30.15, 46.0, num=r['lats'].shape[0], endpoint=True)
    # #//THIS IS REVERSED HERE BC OF FLIP
    latref = np.linspace(46.0, 30.15, num=r['lats'].shape[0], endpoint=True)
    print(latref, r['lats'].shape, latref.shape)

    lonref = np.linspace(-18.1, 36.3, num=r['lons'].shape[0], endpoint=True)
    print(lonref, r['lons'].shape, lonref.shape)



    # exit()

    c = 0.0  #0.025
    # 0.05004501 (for lons)
    # 0.05015945 (for lats)
    #0.0505

    # ext = (lonref[0]-c, lonref[-1]+c, latref[0]-c, latref[-1]+c)
    # plt.imshow(r['data'][0], extent=ext, origin='lower', interpolation='None', alpha=0.5)

    # pref = np.array(r['data'])[1, :, :]  #, 1::20, 2::20]
    # pref = r['data'][1, :, :]  #, 1::20, 2::20]

    ext = (-18.125, 36.325, 30.125, 46.025)
    #ext = (r['lons'][0]-c, r['lons'][-1]+c, r['lats'][0]+c, r['lats'][-1]-c )
    plt.imshow(r['data'][0], extent=ext, origin='upper', interpolation='None', alpha=0.5)  #origin='lower',

    # x, y = np.meshgrid(r['lons'], r['lats'])
    # plt.scatter(x, y, s=0.5, color=(0.0, 0.0, 0.5, 0.2))

    stset = [r['lons'], r['lats'], lonref, latref]
    for s in stset:
        print(s.shape)

    r['lons'] = lonref
    r['lats'] = latref

    deg = 1.0
    tpl = (30.0, 37.0)  #(8.5, 43)  #(35, 34)  #(29.0, 37.0)
    off = (-1, 0)
    concrete = new_data_from_sector(r, tpl, deg, off, 1)

    c = 0.025
    pext = (tpl[0]-c, (tpl[0]+deg)-c, (tpl[1]-deg)+c, tpl[1]+c)

    #scipy.ndimage.shift(input_array, (2.4, 0))
    #shifted_digit_image = scipy.ndimage.shift(concrete[0], (0.5, 0.5))  #, mode='nearest', order=2)
    shifted_digit_image = concrete[0]

    shifted_digit_image[shifted_digit_image == fill_value] = np.nan

    shifted_digit_image = np.vectorize(k_to_deg_c)(shifted_digit_image)

    #shifted_digit_image = transform.pyramid_expand(concrete[0], upscale=2, sigma=0, preserve_range=True)
    plt.imshow(shifted_digit_image, extent=pext, origin='upper', interpolation='None', alpha=0.85)  #origin='lower',

    # print(concrete)


    # xt = np.arange( -20.0, 40, interval)
    # yt = np.arange(29, 47, interval)
    # xt = np.arange(r['lons'][0], r['lons'][-1], interval)
    # yt = np.arange(r['lats'][0], r['lats'][-1], interval)

    # x, y = np.meshgrid(lonref, latref)
    # plt.scatter(x, y, s=1.0, color=(0.0, 0.0, 0.0, 0.5))
    #
    for poly in med_poly.geoms:
        plt.plot(*poly.exterior.xy, color=(0.0, 0.0, 0.0, 1.0))

    # plt.xticks(r['lons'])
    # plt.yticks(r['lats'])

    plt.grid()
    plt.tight_layout()
    plt.axis('equal')
    plt.show()


def check_all():
    hi_res_sst = get_lo_res_wind_data()  #get_hi_res_sst_data()
    part = hi_res_sst['data'][0]
    part = np.vectorize(k_to_deg_c)(part)

    plt.imshow(part, origin='upper', interpolation='None', alpha=1.0)
    plt.grid()
    plt.axis('equal')
    plt.show()


def to_np_string(parts_as_tuple) -> np.array:
    a = np.array(parts_as_tuple[0], dtype=np.str)
    b = np.array(parts_as_tuple[1], dtype=np.str)
    c = np.char.add(a, ' ').astype(str)
    return np.char.add(c, b).astype(str)


def get_sector_wind(data_source, sector_deg, sector_tpl) -> List:
    # #// CALLED PER SECTOR
    sector_U = new_data_from_sector(data_source['u10'], data_source, sector_tpl, sector_deg)
    sector_V = new_data_from_sector(data_source['v10'], data_source, sector_tpl, sector_deg)
    time_frames = data_source['time'].shape[0]
    segments = [(1, 2), (2, 4)]  # set here because wind is 4x4 @ 1deg
    levels = [[], [], []]

    for i in range(time_frames):
        native_U = sector_U[i]
        native_V = sector_V[i]
        offset_native_U = transform.pyramid_expand(native_U, upscale=2, sigma=0, preserve_range=True)[1::2, 1::2]
        offset_native_V = transform.pyramid_expand(native_V, upscale=2, sigma=0, preserve_range=True)[1::2, 1::2]
        offset_WIND = to_np_string((offset_native_U, offset_native_V))
        levels[0].append(offset_WIND)

        for j, s in enumerate(segments):
            offset_WIND = to_np_string((native_U[s[0]::s[1], s[0]::s[1]], native_V[s[0]::s[1], s[0]::s[1]]))
            levels[j+1].append(offset_WIND)

    levels_rfx = [lv for lv in levels if len(lv[0])]
    return levels_rfx


def get_sector_sst(data_source, sector_deg, sector_tpl) -> List:
    # #// CALLED PER SECTOR
    sector_SST = new_data_from_sector(data_source['data'], data_source, sector_tpl, sector_deg)
    time_frames = data_source['time'].shape[0]
    segments = [(1, 2), (2, 4), (5, 10), (10, 20)]  # set here because sst is 20x20 @ 1deg
    levels = [[], [], [], [], []]

    for i in range(time_frames):
        native = sector_SST[i]
        offset_native = transform.pyramid_expand(native, upscale=2, sigma=0, preserve_range=True)[1::2, 1::2]
        # vectorize after transform for float length
        offset_native = np.vectorize(k_to_deg_c)(offset_native)
        levels[0].append(offset_native)
        native = np.vectorize(k_to_deg_c)(native)

        for j, s in enumerate(segments):
            offset_sst = native[s[0]::s[1], s[0]::s[1]]
            levels[j+1].append(offset_sst)

    levels_rfx = [lv for lv in levels if len(lv[0])]
    return levels_rfx


def prepare():
    # with open('/Users/sac/Sites/PythonProject/data/map_polygon', "rb") as poly_file:
    #     med_poly = pickle.load(poly_file)
    # hi_res_sst = get_hi_res_sst_data()
    main_data = get_lo_res_wind_data()
    sector_deg = 1.0
    sector_tpl = (29.0, 37.0)  #(30.0, 37.0)  #(8.5, 43)  #(35, 34)  #(29.0, 37.0)
    levels_and_times = get_sector_wind(main_data, sector_deg, sector_tpl)

    for gl in levels_and_times:
        print(len(gl), gl[0])

    # main_data = get_hi_res_sst_data()
    # levels_and_times = get_sector_sst(main_data, sector_deg, sector_tpl)
    #
    # for gl in levels_and_times:
    #     print(len(gl), gl[0].shape, gl[0])
    #




    exit()
    #
    # segments = [(1, 2), (2, 4)]
    # broad_levels = [[], [], []]
    #
    #
    # off = (0, 0)
    #
    # sector_U = new_data_from_sector(main_data['u10'], main_data, tpl, sector_deg, off, 1)
    # sector_V = new_data_from_sector(main_data['v10'], main_data, tpl, sector_deg, off, 1)
    #
    # time_frames = main_data['time'].shape[0]
    # # root -> level -> time
    # for i in range(time_frames):
    #     native_U = sector_U[i]
    #     native_V = sector_V[i]
    #     offset_native_U = transform.pyramid_expand(native_U, upscale=2, sigma=0, preserve_range=True)[1::2, 1::2]
    #     offset_native_V = transform.pyramid_expand(native_V, upscale=2, sigma=0, preserve_range=True)[1::2, 1::2]
    #     offset_WIND = to_np_string((offset_native_U, offset_native_V))
    #     broad_levels[0].append(offset_WIND)
    #
    #     for j, s in enumerate(segments):
    #         offset_WIND = to_np_string((native_U[s[0]::s[1], s[0]::s[1]], native_V[s[0]::s[1], s[0]::s[1]]))
    #         broad_levels[j+1].append(offset_WIND)

    for gl in levels_and_times:
        print(len(gl), gl[0])

    exit()

    concrete = new_data_from_sector(main_data['u10'], main_data, tpl, sector_deg, off, 1)
    print(concrete.shape)
    time_frames = concrete.shape[0]

    # for i in range(time_frames):
    #     print(i, hi_res_sst['time'][i])
    #     native = concrete[i]
    #     native = np.vectorize(k_to_deg_c)(native)
    #     offset_native = transform.pyramid_expand(native, upscale=2, sigma=0, preserve_range=True)[1::2, 1::2]  #[1::2, 1::2]  #[:-1, :-1]

    native = concrete[0]
    offset_native = transform.pyramid_expand(native, upscale=2, sigma=0, preserve_range=True)[1::2, 1::2]
    display(offset_native, native, (np.nanmin(native), np.nanmax(native)))
    exit()

    #shifted_digit_image = transform.pyramid_expand(concrete[0], upscale=2, sigma=0, preserve_range=True)
    plt.imshow(sdi, origin='upper', interpolation='None', alpha=0.5)  #origin='lower',

    # sxt = (0, 40, 40, 0)
    # plt.imshow(sdi[1::2, 1::2], extent=sxt, origin='upper', interpolation='None', alpha=0.85)  #origin='lower',

    print(sdi)

    x_mesh = np.arange(0, 20, 1)  #;np.linspace(0, 39, 20)
    y_mesh = np.arange(0, 20, 1)  #np.linspace(0, 39, 20)
    x, y = np.meshgrid(x_mesh, y_mesh)
    plt.scatter(x, y, s=2.0, color=(0.0, 0.0, 0.0, 0.5))

    # print(concrete)


    # xt = np.arange( -20.0, 40, interval)
    # yt = np.arange(29, 47, interval)
    # xt = np.arange(r['lons'][0], r['lons'][-1], interval)
    # yt = np.arange(r['lats'][0], r['lats'][-1], interval)

    # x, y = np.meshgrid(lonref, latref)
    # plt.scatter(x, y, s=1.0, color=(0.0, 0.0, 0.0, 0.5))
    #
    # for poly in med_poly.geoms:
    #     plt.plot(*poly.exterior.xy, color=(0.0, 0.0, 0.0, 1.0))

    # plt.xticks(r['lons'])
    # plt.yticks(r['lats'])

    plt.grid()
    # plt.tight_layout()
    plt.axis('equal')
    plt.show()









if __name__ == '__main__':
    # check_all()
    # exit()

    prepare()
    exit()
    # prepare_mini_temps()
    # exit()
    # data_source = prepare_sst_and_wind()

    concrete = new_data_from_sector(prepare_sst_and_wind(), (5.5, 43), 1.0, 1)
    concrete_scaled = new_data_from_sector(prepare_sst_and_wind(), (5.5, 43), 1.0, 1)

    display(concrete_scaled[0, :, :], concrete[0, :, :], (np.nanmin(concrete), np.nanmax(concrete)))
    exit()

    print(concrete.size, concrete.shape)

    # mangle(r)
    # concrete[np.isnan(concrete)] = meantemp
    dat = np.vectorize(k_to_deg_c)(concrete)

    show(concrete[0, :, :])

    #print(dat)
    exit()


    mt = np.mean(dat)  #.data)
    print(dat.shape, k_to_deg_c(0))

    dat[dat.data == ds['sst']._FillValue] = mt
    dat = transform.pyramid_expand(dat, upscale=2, sigma=0, preserve_range=True)
    dat = np.vectorize(k_to_deg_c)(dat)

    #dat[dat == -273.15] = np.nan  #-100  #ds['sst']._FillValue  #mt

    r = {
        "data": dat,
        "lons": np.array(ds['longitude']),
        "lats": np.array(ds['latitude'])
    }

    r['lats'] = np.linspace(r['lats'][0], r['lats'][-1], dat.shape[0]-1)
    r['lons'] = np.linspace(r['lons'][0], r['lons'][-1], dat.shape[1]-1)

    concrete = new_data_from_sector(r, (7, 44), 1.0, 1)
    print(concrete, concrete.size, concrete.shape)

    original_dat = np.vectorize(k_to_deg_c)(original_dat)
    ro = {
        "data": original_dat,
        "lons": np.array(ds['longitude']),
        "lats": np.array(ds['latitude'])
    }
    original = new_data_from_sector(ro, (7, 44), 1.0, 1)

    print(original)
    l_max = concrete[1::2, 1::2]
    print(l_max, l_max.size, l_max.shape)

    l_half = original[1::2, 1::2]
    print(l_half, l_half.size, l_half.shape)

    l_min = original[2::2, 2::2]
    print(l_min, l_min.size, l_min.shape)
    #display(concrete, original)

