#!/usr/bin/env python3

import copy
import os
import shutil

import netCDF4 as nc
from datetime import datetime
from datetime import timedelta

import geopandas as gpd
import matplotlib.pyplot as plt
from shapely import ops
from shapely.geometry import Polygon, MultiPolygon, Point, box, LineString
from shapely.validation import make_valid
from shapely.affinity import translate, scale

from skimage import measure
from scipy.ndimage.filters import gaussian_filter

import json
import math
import pathlib
import numpy as np


area_limits = [0.025, 0.000025]
simp_limits = [0.025, 0.0005]
areas = np.linspace(area_limits[0], area_limits[1], 4)
simps = np.linspace(simp_limits[0], simp_limits[1], 4)

bounds_keys = [-1, -1, 1, 1]
save_path = 'chonk'

indices_aggregate = []


def get_time_as_dict(t, interval='h'):
    date_time_str = '1900-01-01 00:00:00'
    date_time_origin = datetime.strptime(date_time_str, '%Y-%m-%d %H:%M:%S')

    xt = float(t)
    delta = timedelta(hours=xt)
    if interval == 'm':
        delta = timedelta(minutes=xt)

    delta_readable = date_time_origin + delta

    return {'t': "%.2f" % xt, 'ts': str(delta_readable)}


def refine_bounds(bounds, no_pad=False):
    tb = []

    for i, b in enumerate(bounds):
        if no_pad:
            tb.append(round(b))
        else:
            tb.append(round(b) + bounds_keys[i])

    w = abs(tb[0] - tb[2])
    h = abs(tb[1] - tb[3])
    return [tb, w, h]


def deliver(data, prefix, name):
    filename = f'{prefix}/{name}.json'
    with open(filename, 'w') as fp:
        json.dump(data, fp, indent=2)
    return filename


def prepare_map_points(the_poly, simplification, min_area):
    shapes = []
    interiors = []
    exterior = []
    multipoly = []
    group = None

    def qualify(ct, sub_shape, min_a):
        area_test = Polygon(sub_shape.coords)
        shape_simple = sub_shape.simplify(simplification, True)
        if area_test.area >= min_a:
            return [ct, shape_simple]
        return None

    if the_poly.type == 'Polygon':
        group = [the_poly]
    if the_poly.type == 'MultiPolygon':
        group = the_poly.geoms

    for the_poly in group:
        sub_interiors = []
        for c, n in enumerate(the_poly.interiors):
            e = qualify(c, n, min_area)
            if e is not None:
                sub_interiors.append(e[1])
                interiors.append(e[1])
                e.append('interior')
                shapes.append(e)

        e = qualify(0, the_poly.exterior, min_area)
        if e is not None:
            e.append('exterior')
            shapes.append(e)
            exterior.append(e[1])

            if len(the_poly.interiors):
                multipoly.append(Polygon(e[1], sub_interiors))
            else:
                multipoly.append(Polygon(e[1]))

    if len(group) == 1:
        print('prepare_map_points poly')
        return shapes, Polygon(exterior[0], interiors)
    else:
        print('prepare_map_points multipoly')
        return shapes, MultiPolygon(multipoly)


def coords_to_flat(coords_list):
    coords_points = []
    list(coords_points.extend((round(cx, 4), round(cy, 4))) for cx, cy in coords_list)
    return coords_points


def package_poly(poly):
    poly_lex = {'out': coords_to_flat(poly.exterior.coords)}
    if len(poly.interiors):
        poly_lex['ins'] = [coords_to_flat(p.coords) for p in poly.interiors]
    return poly_lex


def normalize_val(val, mi, ma):
    return (val - mi) / (ma - mi)


def get_data_from_sector(sector_tuple, sector_width, data_packet, level):
    # data, lons_ax, lats_ax, resolution, el_w)

    interval = data_packet['degs']
    data = data_packet['data'] # #//HERE: IS THE DATA MORE THAN 2 DIMENSIONS?
    resolution = data_packet['reso'][level]
    lons_ax = data_packet['lons']
    lats_ax = data_packet['lats']
    indices = data_packet['indx']

    # print(sector_tuple, sector_width, end=" > ")
    # print("data.shape", len(data.shape), len(indices))

    # #//DATA NEEDS TO BE DIMENSIONAL!!!
    # #//this gets the positional part all done.
    # #//float32 thetao(time, depth, lat, lon)
    # #// so sst_0(time)_0(depth)
    # #// TODO: should the data have a "start-position" index instead or alongside of zeros?
    # #// in an effort to save on empty data being transferred?

    mask_lon = (lons_ax >= sector_tuple[0]) & (lons_ax < sector_tuple[0] + sector_width)
    mask_lat = (lats_ax <= sector_tuple[1]) & (lats_ax > sector_tuple[1] - sector_width)

    xe = [np.where(lons_ax >= i)[0][0] for n, i in enumerate(lons_ax[mask_lon]) if n % resolution == 0]
    ye = [np.where(lats_ax <= i)[0][0] for n, i in enumerate(lats_ax[mask_lat]) if n % resolution == 0]

    if len(xe) == 0 or len(ye) == 0:
        return None

    xm = (np.amin(lons_ax[mask_lon])-sector_tuple[0])
    ym = sector_tuple[1]-(np.amax(lats_ax[mask_lat]))

    xxs = np.sign(xm)
    xxi = int(math.floor(abs(xm) * (interval/resolution))*xxs)
    yyi = int(math.floor(ym * (interval/resolution)))

    try:
        data_min = np.nanmin(data)
        data_max = np.nanmax(data)
    except TypeError:
        data_min = np.nan
        data_max = np.nan
        pass

    data_master = []
    for n_i in range(0, len(indices[0])):
        if len(indices[1]) > 1:
            data_master = [[]] * len(indices[0])

        for n_j in range(0, len(indices[1])):
            #print(data.dtype)

            if data.dtype != '<U17':
                sector_data = np.zeros([int((interval*sector_width)/resolution), int((interval*sector_width)/resolution)])
                sector_data[:] = np.nan
            else:
                sector_data = np.empty([int((interval*sector_width)/resolution), int((interval*sector_width)/resolution)], dtype='<U17')

            if len(indices[0]) > 1 and len(indices[1]) == 1:
                # assume time data here
                for iy, dy in enumerate(ye):
                    for ix, dx in enumerate(xe):
                        try:
                            sector_data[iy+yyi, ix+xxi] = data[n_i, dy, dx]
                        except IndexError:
                            pass
                        #print(sector_data)
                data_master.append(sector_data.tolist() if data.dtype == '<U17' else sector_data)

            elif len(indices[0]) > 1 and len(indices[1]) > 1:
                if n_j % 4 == 0:
                    # assume time data here
                    for iy, dy in enumerate(ye):
                        for ix, dx in enumerate(xe):
                            try:
                                sector_data[iy+yyi, ix+xxi] = data[n_i, n_j, dy, dx]
                            except IndexError:
                                pass

                    if not np.all(np.isnan(sector_data)):
                        ty = np.where(np.isnan(sector_data), None, sector_data)
                        data_master[n_i].append({str(indices[1][n_j]): ty.tolist()})

            else:
                for iy, dy in enumerate(ye):
                    for ix, dx in enumerate(xe):
                        try:
                            sector_data[iy+yyi, ix+xxi] = data[dy, dx]
                        except IndexError:
                            pass
                data_master = sector_data


    extents = [sector_tuple[0], sector_tuple[0]+sector_width, sector_tuple[1]-sector_width, sector_tuple[1]]
    return {'data': data_master, 'extents': extents, 'limits': [data_min, data_max]}


class Sector:
    def __init__(self, s_id, bounding_box, map_tuple):
        self.box = bounding_box
        self.tuple = map_tuple
        self.id = s_id
        self.map_json_temp = {"lines": [], "fills": []}
        self.map_levels = [None, None, None, None]
        self.data_levels = [None, None, None, None]
        self.data_records = {'id': s_id}
        self.data_static = {}
        self.path = f"Sector-{s_id}"
        self.root_path = None

    def __repr__(self):
        return f"Sector-{self.id}"

    def digest(self):
        return self.data_records
        #  f"Sector-{self.id}",

    def validate(self, data_name, level, blob, spec):
        if blob is not None:
            if level is not None:
                if data_name not in self.data_records:
                    self.data_records[data_name] = {"lv": [0] * 4}
                self.data_records[data_name]["lv"][level] = 1
            else:
                # no level callout means static data
                self.data_records[data_name] = {"lv": [1], "static": True}

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


def get_basic_tiles_set(map_poly, map_bounds, deg_per_sector, map_level_range):
    # #//tile object only contains lines and fills.
    # #//that said, lines are associated with pre-existing guides.

    width = math.ceil(map_bounds[1] * (1 / deg_per_sector))
    height = math.ceil(map_bounds[2] * (1 / deg_per_sector))
    sector_count = width * height
    print("sector_count", sector_count)

    def make_sector_box(num):
        y = math.floor(num / width)
        x = num % width
        sector_box = box(
            map_bounds[0][0] + x * deg_per_sector,
            map_bounds[0][3] - y * deg_per_sector,
            map_bounds[0][0] + x * deg_per_sector + deg_per_sector,
            map_bounds[0][3] - y * deg_per_sector - deg_per_sector)
        test_fill = sector_box.intersects(map_poly)

        sector_tuple = (map_bounds[0][0] + (x * deg_per_sector), map_bounds[0][3] - (y * deg_per_sector),)

        if test_fill:
            return [num, sector_box, sector_tuple]

    def get_basic_map(sector, map_shapes, map_at_level, level):
        # if sector.map_levels[level] is None:
        #     sector.map_levels[level] = copy.deepcopy(sector.map_json_temp)
        # #//FILLS/SHAPES/SOLIDS
        blobby = []
        #shape = sector.box.intersection(map_at_level)

        # shape_d = sector.box.intersection(map_at_level)  #test_bounds.intersection(Polygon(n[1]))
        shape = sector.box.difference(map_at_level)
        # #// TODO: this breaks because it's an empty result and everything goes fucky
        # print(shape)
        if shape.type == 'Polygon':
            blobby.append(package_poly(shape))
            #sector.map_levels[level]['fills'].append(package_poly(shape))
        if shape.type == 'MultiPolygon':
            for n_shape in shape.geoms:
                blobby.append(package_poly(n_shape))
                #sector.map_levels[level]['fills'].append(package_poly(n_shape))

        if len(blobby):
            datum = {"indx": [], "data": blobby}
            sector.add_data("fills", level, datum)

        # #//LINES/ETC
        blobby = []
        for element in map_shapes:
            test_line = sector.box.intersects(element[1])
            if test_line:
                lk = sector.box.intersection(element[1])
                if lk.type == 'LineString':
                    tf = [element[0], coords_to_flat(lk.coords)]
                    blobby.append(tf)
                    #sector.map_levels[level]['lines'].append(tf)
                if lk.type == 'MultiLineString':
                    for line in lk.geoms:
                        tfm = [element[0], coords_to_flat(line.coords)]
                        blobby.append(tfm)
                        #sector.map_levels[level]['lines'].append(tfm)

        if len(blobby):
            datum = {"indx": [], "data": blobby}
            sector.add_data("lines", level, datum)

    spf = [Sector(sector[0], sector[1], sector[2])
           for sector in map(make_sector_box, range(int(sector_count))) if sector is not None]

    for level_acquire in range(0, map_level_range):
        reduc_reso = simps[level_acquire]
        reduc_area = areas[level_acquire]
        print(f"level: {level_acquire}, simplification: {reduc_reso}, max area: {reduc_area}")
        s_map_shapes, ms_poly = prepare_map_points(map_poly, reduc_reso, reduc_area)
        s_map_at_level = make_valid(ms_poly)

        for sector in spf:
            get_basic_map(sector, s_map_shapes, s_map_at_level, level_acquire)

    return spf
    pass


def echo_poly(poly, level):
    if level == 4:
        if len(poly.interiors):
            plt.plot(*poly.exterior.xy, color="pink")
            for sub_poly in poly.interiors:
                plt.plot(*sub_poly.xy, color="pink")
        else:
            plt.plot(*poly.exterior.xy, color="pink")


def plot(shape, r=1):
    print(shape.type, '!')

    def plot_poly(polygon):
        plt.plot(*polygon.exterior.xy, color=(0.5, 1/r, 1/r), linewidth="1")
        if len(polygon.interiors):
            for i in polygon.interiors:
                plt.plot(*i.xy, color="blue", linewidth="1")

    if shape.type == 'MultiPolygon':
        for s in shape.geoms:
            plot_poly(s)
    else:
        plot_poly(shape)


def get_map_shape_guides(master_poly, test_bounds, buffer_constant=0.1, simplification_const=0.005):
    # #//TODO, make canonical shapes guide
    master_shapes, simplified_poly = prepare_map_points(master_poly, simps[3], areas[0])
    guides = []

    print("guide master_shapes", len(master_shapes))

    for c, n in enumerate(master_shapes):
        poly_c = test_bounds.intersection(Polygon(n[1]))
        poly_shape = test_bounds.difference(poly_c)

        if poly_shape.type == 'MultiPolygon':
            kbf = poly_shape.buffer(buffer_constant)
            kbf = kbf.simplify(simplification_const)
            print(kbf)

            if kbf.type == 'Polygon':
                poly_shape = kbf
            else:
                for w, pl in enumerate(kbf.geoms):
                    rtl = test_bounds.intersection(pl.exterior)
                    if rtl.type == 'MultiLineString':
                        dz = ops.linemerge(rtl)
                        print(dz.type)
                        if dz.type == 'LineString':
                            plt.plot(*dz.coords.xy, color="red", linewidth="4")
                            guides.append([n[0], len(guides), dz, None, dz.is_ring, n[1].length])
                        else:
                            for dzi in dz.geoms:
                                plt.plot(*dzi.coords.xy, color="red", linewidth="4")
                                guides.append([n[0], len(guides), dzi, None, dzi.is_ring, n[1].length])
                    else:
                        plt.plot(*rtl.coords.xy, color="blue", linewidth="4")
                        guides.append([n[0], len(guides), rtl, None, rtl.is_ring, n[1].length])

        if poly_shape.type == 'Polygon':
            if len(poly_shape.interiors):
                for i in poly_shape.interiors:
                    bf = buffer_constant
                    if n[2] == 'exterior':
                        bf *= -1

                    i_poly = Polygon(i)
                    kbf = i_poly.buffer(bf)
                    kbf = kbf.simplify(simplification_const)

                    if kbf.type == 'MultiPolygon' and n[2] == 'exterior':
                        s_poly = [sp for sp in kbf.geoms if sp.area > 1.0][0]
                        plt.plot(*s_poly.exterior.xy, color="green", linewidth="2")
                        guides.append([n[0], len(guides), s_poly.exterior, n[2],  True, n[1].length])
                    else:
                        plt.plot(*kbf.exterior.xy, color="black", linewidth="2")
                        guides.append([n[0], len(guides), kbf.exterior, n[2], True, n[1].length])

    fresh_guides_json = []
    for gr in guides:
        print(gr)
        gh = gr.copy()
        gh[2] = coords_to_flat(gh[2].coords)
        fresh_guides_json.append(gh)

    used_shapes = set([s[0] for s in guides])

    def validate_shape(shp):
        return {"id": shp[0], "length": shp[1].length, "area": Polygon(shp[1].coords).area, "type": shp[2]}

    og_poly_data = [validate_shape(pn) for pn in master_shapes if pn[0] in used_shapes]

    return guides, fresh_guides_json, og_poly_data


def bind_guides_to_sector_lines(sector, N_guides):
    for lv in sector.data_levels:
        if lv is not None:
            if 'lines' in lv:
                # print(lv['lines'])
                # exit()
                for li in lv['lines']['data']:
                    li.append([g[1] for g in N_guides if g[0] == li[0]])


def build_contour_levels(source, strata, depth_range, poly_origin):
    contours_all = [[]] * depth_range
    for ra in range(0, depth_range):
        contours_all[ra] = []
        g_data = gaussian_filter(source, sigma=strata["filter"][ra])
        g_range = np.arange(0, strata["depth_max"], strata["depth_interval"][ra])
        #print(np.amin(g_data), np.amax(g_data))
        for lv in g_range:
            f_contours = measure.find_contours(g_data, -lv)
            clutch = []
            for ep in f_contours:
                ep = LineString(np.flip(ep))
                ep = scale(ep, xfact=1/60, yfact=-1/60, origin=(0, 0))
                ep = translate(ep, xoff=poly_origin[0], yoff=poly_origin[1])
                clutch.append(ep)

            contours_all[ra].append({'d': float(lv), 'contours': clutch})

    return contours_all


def get_data_scale_extents(data_packet):
    d = data_packet['data']
    print(d.shape)
    lo = data_packet['lons']
    print(lo[0], lo[-1], lo.shape)
    la = data_packet['lats']
    print(la[0], la[-1], la.shape)

    return lo[0], la[0]


# #//original étape
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
if __name__ == '__main__':
    # # Opening JSON file
    # f = open('data.json')
    #
    # # returns JSON object as
    # # a dictionary
    # data = json.load(f)

    main_config = {
        "root_path": "/Users/sac/Sites/chonk/static",
        "session_path_name": "med_mini_halfdeg_feb_2",
        "rect_ignore": "default",
        "rect":
            {
                "min_X": 8,
                "min_Y": 41,
                "max_X": 10,
                "max_Y": 43.5
            },
        "name": "med_full",
        "bounds": [],
        "map_degrees": 2,
        "debug": 2,
        "guides": None,
        "guide_simp": 0.005,
        "guide_buffer": 0.1,
        "includes": ["depth_points", "contours", "urban", "places", "guides"],  #["sst", "urban", "places", "guides", "depth-points", "contour_levels"],  #["sst", "wind", "no mini-temps", "depth-points", "urban", "places", "guides", "contour_levels"],
        "includes_attributes": {},
        "range": 4,
        "contour_ranges":
            {
                "filter": [2.5, 1.5, 0.5, 0.1],
                "depth_interval": [1000, 500, 250, 125],
                "depth_max": 5000
            }
    }

    urban_smash = []
    places_smash = []
    guides = []
    guides_json = []
    depths_points = []
    depth_contours = []
    mini_temps = []
    surface_temps = []
    surface_wind = []

    # #//SHOULD ALL DATA HAVE LEVELS PER SECTOR?
    # #//TODO: make collection of original paths and areas data
    # #//TODO: make the 'reso' variables scale with map_degrees
    # #//TODO: fix indexing!!!

    path = f'{main_config["root_path"]}/data/'
    if not os.path.exists(path):
        pathlib.Path(path).mkdir(exist_ok=True)

    session_path = path+f'{main_config["session_path_name"]}'

    if os.path.exists(session_path):
        try:
            shutil.rmtree(session_path)
            print(f"deleted all files in {session_path}")
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))
            exit()

    pathlib.Path(session_path).mkdir(exist_ok=True)

    # #//LOAD MAP DATA SHAPEFILE
    med_marine_regions = gpd.read_file("data/goas_med/med_local.shp")
    med_poly = med_marine_regions['geometry'][0]
    med_poly_name = med_marine_regions['name'][0]
    print(med_marine_regions.iloc[0])
    print(med_poly.bounds)

    # my_data = np.genfromtxt('corsica/contour20m_Corsica.csv', delimiter=',')
    # print(my_data.shape)
    # med_poly = Polygon(my_data)

    if main_config["rect_ignore"] == "default":
        bounds = refine_bounds(med_poly.bounds)
        main_config['rect'] = {"min_X": bounds[0][0],
                               "min_Y": bounds[0][1],
                               "max_X": bounds[0][2],
                               "max_Y": bounds[0][3]}
    else:
        rect = main_config['rect']
        #bounds = refine_bounds((rect["min_X"], rect["min_Y"], rect["max_X"], rect["max_Y"]), no_pad=True)
        bounds = [[rect["min_X"], rect["min_Y"], rect["max_X"], rect["max_Y"]], 2, 2.5]



    minx, miny, maxx, maxy = bounds[0]
    bounds_box = box(minx, miny, maxx, maxy)
    print(minx, miny, maxx, maxy)



    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # #//BUILD GUIDES FOR MAP
    if "guides" in main_config["includes"]:
        print("making guides", bounds_box)
        guides, guides_json, shape_data_record = get_map_shape_guides(med_poly,
                                                                      bounds_box,
                                                                      main_config["guide_buffer"],
                                                                      main_config["guide_simp"])
        print("guides", len(guides))
        for lshape in shape_data_record:
            print(lshape)

        d_path = deliver(guides_json, session_path, 'guides')
        print(f"Saved guides data to: '{d_path}'")

        main_config['includes_attributes']['guides'] = {"valid": True}
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "places" in main_config["includes"]:
        populated_places = gpd.read_file("data/ne_10m_populated_places/ne_10m_populated_places.shp")
        for index, loc in enumerate(populated_places.geometry):
            if loc.within(bounds_box):
                places_smash.append({
                    'scale': int(math.floor((populated_places.SCALERANK[index] * 4) / 10)) + 1,
                    'pop': int(populated_places.POP_MAX[index]),
                    'class': populated_places.FEATURECLA[index],
                    'tz': populated_places.TIMEZONE[index],
                    'name': populated_places.NAME[index],
                    'country': populated_places.ADM0NAME[index],
                    'locale': populated_places.ADM1NAME[index],
                    'loc': loc
                })
        main_config['includes_attributes']['places'] = {"valid": True, "levels": 4, "count": len(places_smash), "dimension": 0}
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "urban" in main_config["includes"]:
        urban_planning = gpd.read_file("data/ne_10m_urban_areas/ne_10m_urban_areas.shp")
        u_max = max(urban_planning['scalerank'])
        u_min = min(urban_planning['scalerank'])
        for index, loc in enumerate(urban_planning.geometry):
            if loc.intersects(bounds_box):
                sca = int(math.ceil(normalize_val(urban_planning.scalerank[index], u_min, u_max) * 4))
                area = urban_planning.area_sqkm[index]
                urban_smash.append({
                    'scale': sca,
                    'area': area,
                    'poly': loc
                })
        main_config['includes_attributes']['urban'] = {"valid": True, "levels": 4, "count": len(urban_smash), "dimension": 0}
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "depth_points" in main_config["includes"]:
        depth_points = {
            "data": np.loadtxt(open("./data/bathy_med_tt.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
            "lons": np.loadtxt(open("./data/lon_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
            "lats": np.loadtxt(open("./data/lat_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
            "reso": [10, 6, 2, 1],
            "indx": [[0], [0]],
            "degs": 60,
            "origin": None
        }
        # this is necessary for the contour mapping.
        depth_points['origin'] = get_data_scale_extents(depth_points)
        main_config['includes_attributes']['depth_points'] = {"valid": True, "levels": 4, "dimension": 1}
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "sst" in main_config["includes"] or "wind" in main_config["includes"]:
        #sst and wind intervals cap at 1/4 degreee.
        native_degs = 4
        resolution = [4, 2, 1]
        deg = main_config["map_degrees"]
        levels = [lc for lc in resolution if lc <= native_degs*deg]
        print(levels)

        def k_to_deg_c(z):
            return z - 273.15

        fn = 'download.nc'
        ds = nc.Dataset(fn)

        if "sst" in main_config["includes"]:
            dat = ds['sst'][:, :, :]
            dat = np.vectorize(k_to_deg_c)(dat)

            surface_temps = {
                "data": dat,
                "lons": ds['longitude'][:],
                "lats": ds['latitude'][:],
                "indx": [ds['time'][:], [0]],
                "reso": levels,
                "degs": 4,
                "origin": None
            }

            times_index = [get_time_as_dict(t, 'h') for t in surface_temps['indx'][0]]
            d_path = deliver(times_index, session_path, 'sst-indices')
            print(f"Saved sst index data to: '{d_path}'")

            surface_temps['origin'] = get_data_scale_extents(surface_temps)
            main_config['includes_attributes']['sst'] = {
                "valid": True,
                "levels": len(levels),
                "time": True,
                "time-indices": len(times_index),
                "dimension": 1
            }

        if "wind" in main_config["includes"]:
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
            main_config['includes_attributes']['wind'] = {
                "valid": True,
                "levels": len(levels),
                "time": True,
                "time-indices": len(times_index),
                "dimension": 1
            }

        #exit()
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "mini-temps" in main_config["includes"]:
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
        main_config['includes_attributes']['mini-temps'] = {
            "valid": True,
            "levels": 4,
            "time": True,
            "time-indices": len(times_index),
            "dimension": 2
        }
        # exit()
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "contours" in main_config["includes"]:
        depth_contours = build_contour_levels(depth_points['data'],
                                              main_config["contour_ranges"],
                                              main_config["range"],
                                              depth_points['origin'])
        main_config['includes_attributes']['contours'] = {"valid": True, "levels": 4, "dimension": 1}

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    all_sector_digests = []
    sector_group = get_basic_tiles_set(med_poly, bounds, main_config['map_degrees'], main_config['range'])

    for sector_check in sector_group:
        # #//UPDATE ROOT PATH FOR SECTOR INSTANCE
        sector_check.root_path = session_path

        # #//GUIDE MIXIN HERE:
        if "guides" in main_config["includes"]:
            bind_guides_to_sector_lines(sector_check, guides)

        # #//ADD URBAN AREAS IF SO DOING
        if "urban" in main_config["includes"] and len(urban_smash):
            urban_areas = [ub for ub in urban_smash
                           if type(ub['poly']) != list
                           and sector_check.box.intersects(ub['poly'])]
            for ub in urban_areas:
                ub['poly'] = coords_to_flat(ub['poly'].exterior.coords)

            if len(urban_areas):
                datum = {"indx": [], "data": urban_areas}
                sector_check.add_data("urban", None, datum)
                sector_check.save("urban", "static")

        # #//ADD PLACES IF SO DOING
        if "places" in main_config["includes"] and len(places_smash):
            places = [pl for pl in places_smash
                      if type(pl['loc']) != list
                      and sector_check.box.intersects(pl['loc'])]
            for pl in places:
                pl['loc'] = coords_to_flat(pl['loc'].coords)
            if len(places):
                datum = {"indx": [], "data": places}
                sector_check.add_data("places", None, datum)
                sector_check.save("places", "static")

        # #//ADD DEPTHS at every level IF SO DOING
        if "depth_points" in main_config["includes"] and len(depth_points):
            # generate and save each depth-points level here
            for lv in range(0, main_config['range']):
                layer = depth_points
                get_data = get_data_from_sector(
                    sector_check.tuple,
                    main_config['map_degrees'],
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
        if "contours" in main_config["includes"]:
            for lv in range(0, main_config['range']):
                batch = depth_contours[lv]
                sector_contours = []
                for conty in batch:
                    relevance = [c for c in conty['contours'] if sector_check.box.intersects(c)]
                    paths = []
                    for lap in relevance:
                        path = sector_check.box.intersection(lap)
                        if path.type == 'MultiLineString':
                            for pg in path.geoms:
                                paths.append(coords_to_flat(pg.coords))
                        else:
                            paths.append(coords_to_flat(path.coords))

                    if len(paths):
                        sector_contours.append({"m": conty['d'], "path": paths})

                datum = {"indx": [], "data": sector_contours}
                sector_check.add_data("contours", lv, datum)

            sector_check.save("contours", "data_levels")

        # #//ADD SST at every level IF SO DOING
        if "sst" in main_config["includes"]:
            # generate and save each depth-points level here
            for lv in range(0, len(surface_temps['reso'])):
                get_data = get_data_from_sector(
                    sector_check.tuple,
                    main_config['map_degrees'],
                    surface_temps,
                    lv
                )

                if get_data is not None:
                    d_values = get_data['data'].copy()
                    ty = np.where(np.isnan(d_values), None, d_values)
                    datum = {"indx": surface_temps['indx'][0], "data": ty.tolist()}
                    sector_check.add_data("sst", lv, datum, {"d": 2, "indx": len(surface_temps['indx'][0])})

            sector_check.save("sst", "data_levels")
            pass

        # #//ADD WIND at every level IF SO DOING
        if "wind" in main_config["includes"]:
            # generate and save each depth-points level here special case string data!
            for lv in range(0, len(surface_wind['reso'])):
                get_data = get_data_from_sector(
                    sector_check.tuple,
                    main_config['map_degrees'],
                    surface_wind,
                    lv
                )

                if get_data is not None:
                    d_values = get_data['data'].copy()
                    datum = {"indx": surface_wind['indx'][0], "data": d_values}
                    sector_check.add_data("wind", lv, datum, {"d": 2, "indx": len(surface_wind['indx'][0])})

            sector_check.save("wind", "data_levels")
            pass

        # #//ADD SST at every level IF SO DOING
        if "mini-temps" in main_config["includes"]:

            # generate and save each depth-points level here
            for lv in range(0, len(mini_temps['reso'])):  #main_config['range']):
                get_data = get_data_from_sector(
                    sector_check.tuple,
                    main_config['map_degrees'],
                    mini_temps,
                    lv
                )

                if get_data is not None:
                    d_values = get_data['data'].copy()
                    datum = {"indx": mini_temps['indx'][0], "data": d_values}
                    sector_check.add_data("mini_temps", lv, datum, {"d": 3, "indx": len(mini_temps['indx'][0])})

            sector_check.save("mini_temps", "data_levels")
            pass

        sector_check.save("fills", "data_levels")
        sector_check.save("lines", "data_levels")

        digest = sector_check.digest()
        all_sector_digests.append(digest)
        # print(digest)

    d_path = deliver(main_config, session_path, 'map-spec')
    print(f"Saved map data to: '{d_path}'")

    d_path = deliver(all_sector_digests, session_path, 'map-digest')
    print(f"Saved digest data to: '{d_path}'")