#!/usr/bin/env python3
import json
from typing import List
from json_encoder import JsonSafeEncoder
import config

import os
import shutil
import copy
import pickle

import netCDF4 as nc
from datetime import datetime
from datetime import timedelta

# import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, box, LineString
from shapely.affinity import translate, scale

from skimage import measure
from scipy.ndimage.filters import gaussian_filter

import json
import math
import pathlib
import numpy as np

import marine_regions

import obspkg_util as utl

from wudi_module import wudi_prepare, wudi_get_sector, wudi_filter_dict


def get_time_as_dict(t, interval='h', time_origin='1900-01-01 00:00:00'):

    date_time_str = time_origin
    date_time_origin = datetime.strptime(date_time_str, '%Y-%m-%d %H:%M:%S')

    xt = float(t)
    delta = timedelta(hours=xt)
    if interval == 'm':
        delta = timedelta(minutes=xt)
    if interval == 's':
        delta = timedelta(seconds=xt)

    delta_readable = date_time_origin + delta

    return {'t': "%.2f" % xt, 'ts': str(delta_readable)}


def deliver(data, prefix, name):
    filename = f'{prefix}/{name}.json'
    with open(filename, 'w') as fp:
        json.dump(data, fp, cls=JsonSafeEncoder)  # , indent=2,
    return filename


def normalize_val(val, mi, ma):
    return (val - mi) / (ma - mi)


class Sector:
    def __init__(self, s_id, bounding_box, map_tuple, levels):
        self.box = bounding_box
        self.bounds_tuple = None
        self.tuple = map_tuple
        self.id = s_id
        self.map_json_temp = {"lines": [], "fills": []}
        self.map_levels = [None] * levels
        self.data_levels = [None] * levels
        self.data_records = {'id': s_id}
        self.data_static = {}
        self.path = f"Sector-{s_id}"
        self.root_path = None
        self.levels = levels

    def __repr__(self):
        return f"Sector-{self.id}"

    def digest(self):
        return self.data_records
        #  f"Sector-{self.id}",

    def validate(self, data_name, level, blob, spec):
        if blob is not None:
            if level is not None:
                if data_name not in self.data_records:
                    self.data_records[data_name] = {"lv": [0] * self.levels}
                    self.data_records[data_name]["lvmax"] = 0

                self.data_records[data_name]["lv"][level] = 1
                self.data_records[data_name]["lvmax"] += 1
            else:
                # no level callout means static data
                self.data_records[data_name] = {"lv": [1], "static": True, "lv_max":3}

            if spec is not None:
                self.data_records[data_name]['spc'] = spec

    def add_data(self, data_name, level=None, blob=None, spec=None):
        self.validate(data_name, level, blob, spec)
        if level is None:
            self.data_static[data_name] = blob
        else:
            # #//FOR static DATA WITH ZOOM LEVELS: ie: depths and countours
            if self.data_levels[level] is None:
                self.data_levels[level] = {}
            if data_name not in self.data_levels[level]:
                self.data_levels[level][data_name] = {}

            self.data_levels[level][data_name] = blob

    def save(self, item=None, data_type=None):
        path_prefix = self.root_path + '/' + self.path
        pathlib.Path(path_prefix).mkdir(exist_ok=True)

        if item == "map_levels":
            for n, lv in enumerate(self.map_levels):
                if lv is not None:
                    self.validate("map", n, 1, None)
                    deliver(lv, path_prefix, f"map-{n}")

        if data_type == "static":
            deliver(self.data_static[item]['data'], path_prefix, f"{item}-0")

        if data_type == "data_levels":
            for n, lv in enumerate(self.data_levels):
                if lv is not None and item in lv:
                    blob = lv[item]
                    rf = len(blob['indx'])
                    if rf:
                        for gn, grp in enumerate(blob['indx']):
                            deliver(lv[item]['data'][gn], path_prefix, f"{item}-{n}-{gn}")
                    else:
                        deliver(lv[item]['data'], path_prefix, f"{item}-{n}")


def init_sector(num, m_wid, m_deg, m_bnd, m_lev) -> Sector:
    y = math.floor(num / m_wid)
    x = num % m_wid
    sector_box = box(
        m_bnd[0][0] + x * m_deg,
        m_bnd[0][3] - y * m_deg,
        m_bnd[0][0] + x * m_deg + m_deg,
        m_bnd[0][3] - y * m_deg - m_deg)
    sector_tuple = (m_bnd[0][0] + (x * m_deg), m_bnd[0][3] - (y * m_deg),)
    return Sector(num, sector_box, sector_tuple, m_lev)


def build_sectors(M_bounds, M_deg, M_levels) -> List:
    width = math.ceil(M_bounds[1] * (1 / M_deg))
    height = math.ceil(M_bounds[2] * (1 / M_deg))
    sector_count = width * height
    print("sector instance count", sector_count)
    return [init_sector(n, width, M_deg, M_bounds, M_levels) for n in range(int(sector_count))]
    pass


def build_contour_levels(source, strata, depth_range, poly_origin):
    contours_all = [[]] * depth_range

    for ra in range(0, depth_range):
        contours_all[ra] = []
        # if ra > len(strata["filter"])-1:
        #     break

        g_data = gaussian_filter(source, sigma=strata["filter"][ra])
        g_range = np.arange(0, strata["depth_max"], strata["depth_interval"][ra])
        for i, lv in enumerate(g_range):
            f_contours = measure.find_contours(g_data, -lv)
            clutch = []
            for ep in f_contours:
                ep = LineString(np.flip(ep))
                ep = scale(ep, xfact=1/60, yfact=-1/60, origin=(0, 0))
                ep = translate(ep, xoff=poly_origin[0], yoff=poly_origin[1])
                clutch.append(ep)

            contours_all[ra].append({'d': float(lv), 'contours': clutch})
            utl.show_progress(f"generate contours {ra} {lv}m", i, len(g_range))

    return contours_all


def get_data_scale_extents(data_packet):
    d = data_packet['data']
    lo = data_packet['lons']
    la = data_packet['lats']
    return lo[0], la[0]


# #// hand the config to the builder
def builder(cfg):
    area_limits = cfg['area_limits']
    simp_limits = cfg['simp_limits']
    areas = np.linspace(area_limits[0], area_limits[1], 4)
    simps = np.linspace(simp_limits[0], simp_limits[1], 4)
    cfg['criteria'] = {'areas': areas, 'simplify': simps}

    session_path = f'{cfg["data_destination"]}/{cfg["name"]}'
    if not os.path.exists(session_path):
        pathlib.Path(session_path).mkdir(exist_ok=True)

    if os.path.exists(session_path):
        try:
            shutil.rmtree(session_path)
            print(f"deleted all files in {session_path}")
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))
            exit()

    pathlib.Path(session_path).mkdir(exist_ok=True)

    # #//LOAD MAP DATA from map_author
    with open(cfg['map_data_source'], "rb") as poly_file:
        med_poly = pickle.load(poly_file)

    cb = cfg['box']
    w = cb[2] - cb[0]
    h = cb[3] - cb[1]
    bounds = [cb, w, h]

    cfg['shape'] = {'w': w, 'h': h}
    cfg['rect'] = {"min_X": cb[0],
                   "min_Y": cb[1],
                   "max_X": cb[2],
                   "max_Y": cb[3]}

    minx, miny, maxx, maxy = bounds[0]
    bounds_box = box(minx, miny, maxx, maxy)
    print('Essential Bounds (deg):', bounds)

    fig, ax = plt.subplots()

    map_collection = utl.parse_map_shapes(med_poly, bounds_box, cfg, ax)

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "guides" in cfg["includes"]:
        guides = []
        [guides.extend(g['guides']) for g in map_collection]
        # guides = get_map_shapes_guides(map_collection)
        out = json.dumps(guides, indent=2, cls=JsonSafeEncoder)
        print(out)

        for g in guides:
            g['geom'] = utl.coords_to_flat(g['geom'].coords)

        print("\rGuides Counted:", len(guides))
        d_path = deliver(guides, session_path, 'guides')
        print(f"Saved guides data to: '{d_path}'")

        cfg['includes_attributes']['guides'] = {"valid": True}
        cfg['includes_data']['guides'] = guides
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "places" in cfg["includes"]:
        print('places')
        places_smash = []
        populated_places = utl.gpd_scope_shp_to_bounds(
            f"{cfg['data_source_path']}/ne_10m_populated_places/ne_10m_populated_places.shp", bounds_box)

        p_max = max(populated_places['SCALERANK'])
        p_min = min(populated_places['SCALERANK'])
        print(p_max, p_min)
        for index, loc in enumerate(populated_places.geometry):

            # if loc.within(bounds_box):
            #print(normalize_val(populated_places.POP_MAX[index], p_max, p_min))
            sca = int(math.ceil(normalize_val(populated_places.SCALERANK[index], p_min, p_max) * 4))
            name = populated_places.NAME[index]
            #print(sca)
            utl.show_progress(f'{name} {sca}', index, len(populated_places.geometry))

            places_smash.append({
                'name': populated_places.NAME[index],
                'country': populated_places.ADM0NAME[index],
                'locale': populated_places.ADM1NAME[index],
                'scale': sca,  # int(math.floor((populated_places.SCALERANK[index] * 4) / 10)) + 1,
                'population': int(populated_places.POP_MAX[index]),
                'class': populated_places.FEATURECLA[index],
                'tz': populated_places.TIMEZONE[index],
                'loc': loc
            })
        cfg['includes_attributes']['places'] = {"valid": True, "levels": 4, "count": len(places_smash), "dimension": 0}
        cfg['includes_data']['places'] = places_smash
        print('')
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "urban" in cfg["includes"]:
        print('urban areas')
        urban_smash = []
        urban_planning = utl.gpd_scope_shp_to_bounds(
            f"{cfg['data_source_path']}/ne_10m_urban_areas/ne_10m_urban_areas.shp", bounds_box)

        u_max = max(urban_planning['scalerank'])
        u_min = min(urban_planning['scalerank'])
        for index, loc in enumerate(urban_planning.geometry):
            utl.show_progress('urban areas', index, len(urban_planning.geometry))
            # if loc.intersects(bounds_box):
            sca = int(math.ceil(normalize_val(urban_planning.scalerank[index], u_min, u_max) * 4))
            area = urban_planning.area_sqkm[index]
            urban_smash.append({
                'scale': sca,
                'area': area,
                'poly': loc
            })
        cfg['includes_attributes']['urban'] = {"valid": True, "levels": 4, "count": len(urban_smash), "dimension": 0}
        cfg['includes_data']['urban'] = urban_smash
        print('')
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "protected_regions" in cfg["includes"]:
        print('protected areas')
        regions = marine_regions.load()
        cfg['includes_attributes']['protected_regions'] = {"valid": True, "levels": 4, "count": len(regions), "dimension": 0}
        cfg['includes_data']['protected_regions'] = regions

        # print(json.dumps(regions, indent=1, cls=JsonSafeEncoder))
        print('')
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "depth_points" in cfg["includes"]:
        print('depth_points')
        depth_points = {
            "data": np.loadtxt(open(f"{cfg['data_source_path']}/bathy_med_tt.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
            "lons": np.loadtxt(open(f"{cfg['data_source_path']}/lon_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
            "lats": np.loadtxt(open(f"{cfg['data_source_path']}/lat_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
            "reso": [10, 6, 2, 1],
            "indx": [[0], [0]],
            "degs": 60,
            "origin": None
        }
        # this is necessary for the contour mapping.
        depth_points['origin'] = get_data_scale_extents(depth_points)
        cfg['includes_attributes']['depth_points'] = {"valid": True, "levels": 4, "dimension": 1}
        cfg['includes_data']['depth_points'] = depth_points
        print('')
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "sst" in cfg["includes"] or "wind" in cfg["includes"]:
        #sst and wind intervals cap at 1/4 degreee.

        def k_to_deg_c(z):
            return z - 273.15

        if "sst" in cfg["includes"]:
            native_degs = 20
            resolution = [10, 5, 4, 2]
            deg = cfg["map_degrees"]
            levels = [lc for lc in resolution if lc <= native_degs * deg]
            print(levels)

            # fn = f"{cfg['data_source_path']}/download.nc"
            new_fn = f"{cfg['data_source_path']}/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1645013561954.nc"
            # #//news in seconds for time
            # #//have to restructure this for it to work: maybe reverse since dimensions:
            #     time = 31;
            #     lat = 318;
            #     lon = 867;
            ds = nc.Dataset(new_fn)
            # nds = nc.Dataset(new_fn)


            # dat = nds['analysed_sst'][:, :, :]
            dat = ds['analysed_sst']  #[:, :, :]
            dat = np.flip(dat, axis=1)
            dat = np.vectorize(k_to_deg_c)(dat)

            surface_temps = {
                "data": dat,
                "lons": ds['lon'][:],
                "lats": np.flip(ds['lat'], axis=0)[:],
                "indx": [ds['time'][:], [0]],
                "reso": levels,
                "degs": 20,
                "origin": None
            }

            times_index = [get_time_as_dict(t, 's', '1981-01-01 00:00:00') for t in surface_temps['indx'][0]]
            d_path = deliver(times_index, session_path, 'sst-indices')
            print(f"Saved sst index data to: '{d_path}'")

            surface_temps['origin'] = get_data_scale_extents(surface_temps)
            cfg['includes_attributes']['sst'] = {
                "valid": True,
                "levels": len(levels),
                "time": True,
                "time-indices": len(surface_temps['indx'][0]),
                "dimension": 1
            }
            cfg['includes_data']['sst'] = surface_temps

        if "wind" in cfg["includes"]:
            native_degs = 4
            resolution = [4, 2, 1]
            deg = cfg["map_degrees"]
            levels = [lc for lc in resolution if lc <= native_degs * deg]
            print(levels)

            fn = f"{cfg['data_source_path']}/download.nc"
            ds = nc.Dataset(fn)

            U = ds['u10']
            V = ds['v10']
            a = np.array([U], dtype=np.str)
            b = np.array([V], dtype=np.str)
            c = np.char.add(a, ' ').astype(str)
            dat = np.char.add(c, b).astype(str)[0]

            surface_wind = {
                "data": dat,
                "lons": ds['longitude'][:],
                "lats": ds['latitude'][:],
                "indx": [ds['time'][:], [0]],
                "reso": levels,
                "degs": 4,
                "origin": None
            }

            times_index = [get_time_as_dict(t, 'h') for t in surface_wind['indx'][0]]
            d_path = deliver(times_index, session_path, 'wind-indices')
            print(f"Saved wind index data to: '{d_path}'")

            surface_wind['origin'] = get_data_scale_extents(surface_wind)
            cfg['includes_attributes']['wind'] = {
                "valid": True,
                "levels": len(levels),
                "time": True,
                "time-indices": len(times_index),
                "dimension": 1
            }
            cfg['includes_data']['wind'] = surface_wind
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "mini-temps" in cfg["includes"]:
        fn = '1999-6-10_med-cmcc-tem-rean-d_1638647808820.nc'
        ds = nc.Dataset(fn)

        dat = ds['thetao'][:, :, :, :]
        #dat = np.vectorize(k_to_deg_c)(dat)

        mini_temps = {
            "data": dat,
            "lons": ds['lon'][:],
            "lats": ds['lat'][:],
            "indx": [ds['time'][50:], ds['depth'][:]],
            "reso": [8, 6, 2, 1],
            "degs": 24,
            "origin": None
        }

        times_index = [get_time_as_dict(t, 'm') for t in mini_temps['indx'][0]]
        d_path = deliver(times_index, session_path, 'mini-temps-indices')
        print(f"Saved mini-temps index data to: '{d_path}'")

        # print(surface_temps['indx'])
        # this is necessary for the contour mapping.
        mini_temps['origin'] = get_data_scale_extents(mini_temps)
        cfg['includes_attributes']['mini-temps'] = {
            "valid": True,
            "levels": 4,
            "time": True,
            "time-indices": len(times_index),
            "dimension": 2
        }
        # exit()
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "contours" in cfg["includes"]:
        print('depth_contours')
        depth_contours = build_contour_levels(depth_points['data'],
                                              cfg["contour_ranges"],
                                              cfg["range"],
                                              depth_points['origin'])
        cfg['includes_attributes']['contours'] = {"valid": True, "levels": 4, "dimension": 1}
        cfg['includes_data']['contours'] = depth_contours
        print('')
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "wudi" in cfg["includes"]:
        wudi_data, w_tree, w_points, w_index, w_maxes = wudi_prepare(cb)
        print("WUDI", len(wudi_data), 'entries')
        cfg['includes_attributes']['wudi'] = {
            "valid": True,
            "levels": 3,
            "count": len(wudi_data),
            "dimension": 0,
            "limits": w_maxes
        }
        cfg['includes_data']['wudi'] = wudi_data
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    all_sector_digests = []
    sst_poach_flat = []
    sector_group = build_sectors(bounds, cfg['map_degrees'], cfg['range'])

    for i, sector_check in enumerate(sector_group):
        utl.show_progress('Sectors', i, len(sector_group))

        for r in range(sector_check.levels):
            collect = [ref for ref in map_collection if ref['geom'][r] and ref['geom'][r].intersects(sector_check.box)]
            fills = []
            lines = []
            for ne, element in enumerate(collect):

                p_geom = sector_check.box.intersection(element['geom'][r])
                fills.extend([utl.package_poly(pl) for pl in utl.poly_s_to_list(p_geom)])

                element_lines = [line for line in element['lines'] if line['lev'] == r and line['shape'] == element['id']]

                for line in element_lines:
                    c_line = copy.deepcopy(line)
                    p_line = sector_check.box.intersection(c_line['geom'])
                    parts = [utl.coords_to_flat(li.coords) for li in utl.line_s_to_list(p_line)]
                    c_line['geom'] = parts
                    lines.append(c_line)

            datum = {"indx": [], "data": fills}
            sector_check.add_data("fills", r, datum)

            datum = {"indx": [], "data": lines}
            sector_check.add_data("lines", r, datum)

        # #//UPDATE ROOT PATH FOR SECTOR INSTANCE
        sector_check.root_path = session_path

        # #//ADD URBAN AREAS IF SO DOING
        if "urban" in cfg["includes"]:

            urban_areas = [ub for ub in cfg["includes_data"]["urban"]
                           if type(ub['poly']) != list
                           and sector_check.box.intersects(ub['poly'])]
            for ub in urban_areas:
                ub['poly'] = utl.coords_to_flat(ub['poly'].exterior.coords)

            if len(urban_areas):
                datum = {"indx": [], "data": urban_areas}
                sector_check.add_data("urban", None, datum)
                sector_check.save("urban", "static")

        # #//ADD PROTECTED AREAS IF SO DOING
        if "protected_regions" in cfg["includes"]:

            protected_areas = [g for g in cfg["includes_data"]["protected_regions"] if Point(g['CENTROID']).within(sector_check.box)]
            if len(protected_areas):
                for ub in protected_areas:
                    all_paths = []
                    mash = ub['geometry']
                    if mash.type == 'Polygon':
                        all_paths.append(utl.package_poly(mash))

                    if mash.type == 'MultiPolygon':
                        for n_shape in mash.geoms:
                            all_paths.append(utl.package_poly(n_shape))

                    ub['geometry'] = all_paths

                datum = {"indx": [], "data": protected_areas}
                sector_check.add_data("protected_regions", None, datum)
                sector_check.save("protected_regions", "static")


                # //exit()

            # protected_areas = [
            #     ub for ub in cfg["includes_data"]["protected_regions"]
            #     if type(ub['geometry']) != list
            #     and sector_check.box.intersects(ub['geometry'])
            # ]
            #
            # for ub in protected_areas:
            #     all_paths = []
            #     mash = ub['geometry']
            #     if mash.type == 'Polygon':
            #         all_paths.append(package_poly(mash))
            #
            #     if mash.type == 'MultiPolygon':
            #         for n_shape in mash.geoms:
            #             all_paths.append(package_poly(n_shape))
            #
            #     ub['geometry'] = all_paths  #coords_to_flat(ub['poly'].exterior.coords)
            #
            # print(protected_areas)
            # exit()
            #
            # if len(protected_areas):
            #     datum = {"indx": [], "data": protected_areas}
            #     sector_check.add_data("protected_areas", None, datum)
            #     sector_check.save("protected_areas", "static")

        # #//ADD PLACES IF SO DOING
        if "places" in cfg["includes"]:
            places = [pl for pl in cfg["includes_data"]["places"]
                      if type(pl['loc']) != list
                      and sector_check.box.intersects(pl['loc'])]
            for pl in places:
                pl['loc'] = utl.coords_to_flat(pl['loc'].coords)
            if len(places):
                datum = {"indx": [], "data": places}
                sector_check.add_data("places", None, datum)
                sector_check.save("places", "static")

        # #//ADD DEPTHS at every level IF SO DOING
        if "depth_points" in cfg["includes"]:
            # generate and save each depth-points level here
            for lv in range(0, len(depth_points['reso'])):
                layer = cfg["includes_data"]["depth_points"]
                get_data = utl.get_data_from_sector(
                    sector_check.tuple,
                    cfg['map_degrees'],
                    depth_points,
                    lv
                )

                if get_data is not None:
                    d_values = get_data['data'].copy()
                    ty = np.where(np.isnan(d_values), None, d_values)
                    datum = {"indx": [], "data": ty.tolist()}
                    sector_check.add_data("depth_points", lv, datum)

            sector_check.save("depth_points", "data_levels")

        # #//ADD CONTOURS at every level IF SO DOING
        if "contours" in cfg["includes"]:
            for lv in range(0, cfg["range"]):  # len(depth_points['reso'])):
                batch = cfg['includes_data']['contours'][lv]
                sector_contours = []
                for conty in batch:
                    relevance = [c for c in conty['contours'] if sector_check.box.intersects(c)]
                    paths = []
                    for lap in relevance:
                        path = sector_check.box.intersection(lap)
                        if path.type == 'MultiLineString':
                            for pg in path.geoms:
                                paths.append(utl.coords_to_flat(pg.coords))
                        else:
                            paths.append(utl.coords_to_flat(path.coords))

                    if len(paths):
                        sector_contours.append({"m": conty['d'], "path": paths})

                datum = {"indx": [], "data": sector_contours}
                sector_check.add_data("contours", lv, datum)

            sector_check.save("contours", "data_levels")

        # #//ADD SST at every level IF SO DOING
        if "sst" in cfg["includes"]:
            # generate and save each depth-points level here
            # cfg['includes_data']['sst'] = surface_temps
            for lv in range(0, len(cfg['includes_data']['sst']['reso'])):
                get_data = utl.get_data_from_sector(
                    sector_check.tuple,
                    cfg['map_degrees'],
                    cfg['includes_data']['sst'],
                    lv
                )

                if get_data is not None:
                    d_values = get_data['data'].copy()
                    ty = np.where(np.isnan(d_values), None, d_values)
                    datum = {"indx": cfg['includes_data']['sst']['indx'][0], "data": ty.tolist()}
                    sector_check.add_data("sst", lv, datum, {"d": 2, "indx": len(cfg['includes_data']['sst']['indx'][0])})

            sector_check.save("sst", "data_levels")
            pass

        # #//ADD WIND at every level IF SO DOING
        if "wind" in cfg["includes"]:
            # generate and save each depth-points level here special case string data!
            for lv in range(0, len(surface_wind['reso'])):
                get_data = utl.get_data_from_sector(
                    sector_check.tuple,
                    cfg['map_degrees'],
                    surface_wind,
                    lv
                )

                #print(get_data[1])

                if get_data is not None:

                    d_values = get_data['data'].copy()
                    sst_poach_flat.append(d_values)

                    datum = {"indx": surface_wind['indx'][0], "data": d_values}
                    sector_check.add_data("wind", lv, datum, {"d": 2, "indx": len(surface_wind['indx'][0])})

            sector_check.save("wind", "data_levels")
            pass

        # #//ADD SST at every level IF SO DOING
        if "mini-temps" in cfg["includes"]:

            # generate and save each depth-points level here
            for lv in range(0, len(mini_temps['reso'])):  #cfg['range']):
                get_data = utl.get_data_from_sector(
                    sector_check.tuple,
                    cfg['map_degrees'],
                    mini_temps,
                    lv
                )

                if get_data is not None:
                    d_values = get_data['data'].copy()
                    datum = {"indx": mini_temps['indx'][0], "data": d_values}
                    sector_check.add_data("mini_temps", lv, datum, {"d": 3, "indx": len(mini_temps['indx'][0])})

            sector_check.save("mini_temps", "data_levels")
            pass

        # #//ADD WUDI DATA!
        if "wudi" in cfg["includes"]:
            selection = wudi_get_sector(sector_check.box.bounds, w_tree, w_points, w_index)
            if selection is not None:
                wudi_data_points = [wudi_filter_dict(wudi_data[index]) for index in selection]
                datum = {"indx": [], "data": wudi_data_points}
                sector_check.add_data("wudi", None, datum)
                sector_check.save("wudi", "static")

        sector_check.save("fills", "data_levels")
        sector_check.save("lines", "data_levels")
        digest = sector_check.digest()
        all_sector_digests.append(digest)

    # #//WRAP-UP
    cfg.pop('includes_data')
    d_path = deliver(cfg, session_path, 'map-spec')
    print(f"\rSaved map data to: '{d_path}'")

    d_path = deliver(all_sector_digests, session_path, 'map-digest')
    print(f"\rSaved digest data to: '{d_path}'")

    d_path = deliver(sst_poach_flat, session_path, 'map-digest-sst')
    print(f"\rSaved flat sst data to: '{d_path}'")

if __name__ == '__main__':
    with open('./configuration_maps.json') as main_config:
        maps = json.load(main_config)

    maps_to_make = []
    if config.call_all:
        maps_to_make.extend([k['name'] for k in maps["maps"]])
    else:
        maps_to_make.extend(config.make_map)

    for the_ref in maps_to_make:
        print(f"\n\nstarting to make map for '{the_ref}'")

        map_main = [m for m in maps["maps"] if m['name'] == the_ref][0]
        with open('./configuration.json') as config_template:
            blank_config = json.load(config_template)

        if blank_config:
            for k, v in map_main.items():
                blank_config[k] = v

        blank_config['data_destination'] = config.data_target
        blank_config['data_source_path'] = config.data_path
        blank_config['map_data_source'] = f'{config.data_path}/{map_main["map_source"]}'
        # blank_config['data_destination'] = maps['data_destination']

        if not config.debugging:
            builder(blank_config)

        out = json.dumps(blank_config, indent=2, cls=JsonSafeEncoder)
        print(out)
        print('debugging', config.debugging)
        print(f"done with '{the_ref}'\n\n")

    maps['data_destination'] = config.data_target
    maps['data_source_path'] = config.data_path
    deliver(maps, config.data_target, 'obspkg-maps')

    # #// exit()
