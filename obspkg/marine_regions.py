#!/usr/bin/env python3

import geopandas as gpd
import matplotlib.pyplot as plt
import csv
import math
import fiona
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection

# #// https://mapamed.org/index.php?language=en

area_filter = [
    "geometry",
    "NAME",
    "STATUS_YR",
    "STATUS_ENG",
    "MAPAMED_ID",
    "PARENT_ID",
    "REP_AREA",
    "SITE_TYPE_ENG",
    "IUCN_CAT_ENG",
    "WEBSITE"
]

criteria_filter = ['Marine N2000 Proposed Site', 'MPA of National Statute']


# Plots a Polygon to pyplot `ax`
def plot_polygon(ax, poly, **kwargs):
    path = Path.make_compound_path(
        Path(np.asarray(poly.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])

    patch = PathPatch(path, **kwargs)
    collection = PatchCollection([patch], **kwargs)

    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection


def normalize_val(val, mi, ma):
    return (val - mi) / (ma - mi)


def load():
    country_file = '../data/MAPAMED_2019_edition/mapamed_2019_dataset_definition_countries.tsv'
    with open(country_file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        country_data = [row for row in reader]

    shapes_file = '../data/MAPAMED_2019_edition/MAPAMED_2019_spatial_data_epsg3035.gpkg'
    gpkg_tables = {}

    for layer_name in fiona.listlayers(shapes_file):
        # https://gis.stackexchange.com/questions/24340/converting-coordinates-from-meters-to-decimal-degrees-in-qgis
        print(layer_name)
        g_layer = gpd.read_file(shapes_file, layer=layer_name)
        print('converting crs to 4326')
        g_layer = g_layer.to_crs(crs=4326)  # 'EPSG:3857')  #, allow_override=True)
        gpkg_tables[layer_name] = g_layer

    # print(json.dumps(country_data, indent=4))
    return gpkg_parse(gpkg_tables, country_data)


def gpkg_parse(main_tables, country_data):
    regions_all = []
    elements_all = []
    for k, v in main_tables.items():

        print('loading regions', k)
        if k == 'Scope of the Barcelona Convention (IHO-MSFD)':
            regions_all = v.to_dict('records')
            for d in regions_all:
                d.pop('geometry')

        print('loading elements', k)
        if k == 'MAPAMED 2019 edition - EPSG:3035':
            indx = v.shape[0]
            v = v.fillna(np.nan).replace([np.nan], [None])
            # #// HIFIVE TO SELF ^
            mct = [v.iloc[i] for i in range(0, indx) if v.iloc[i]['DESIG_CAT_ENG'] in criteria_filter]
            N = len(mct)
            max_area = 0
            avg_area = 0
            for mv in mct:
                compare = mv['REP_AREA'] if mv['REP_AREA'] is not None else 0.0
                avg_area += compare

                if compare > max_area:
                    max_area = mv['REP_AREA']

            avg_area /= len(mct)

            print(avg_area, 'km')
            print(max_area, 'km')

            for i, elem in enumerate(mct):
                area = {}
                for f in area_filter:
                    area[f] = elem[f]  # if not math.isnan(elem[f]) else 'null'

                if type(area['REP_AREA']) == float:
                    nsca = normalize_val(area['REP_AREA'], 0.1, avg_area)
                    area['scale'] = 1 if math.isnan(nsca) or nsca < 0 else math.ceil(nsca) if nsca < 4 else 4
                else:
                    area['scale'] = 1

                area['CENTROID'] = np.array(elem['geometry'].centroid.coords[0])
                area['COUNTRY'] = [
                    p['COUNTRY_ENG'] for p in country_data if p['ISO3'] in elem['ISO3'][1:-1]
                ]
                area['MED_REGION'] = [
                    p['NAME_ENG'] for p in regions_all if p['MSFD_REGION'] in elem['MSFD_REGION'][1:-1]
                ]

                elements_all.append(area)

    return elements_all





def plot():
    fig, ax = plt.subplots()

    print(area['REP_AREA'], area['scale'])
    pcol = (np.random.rand(1)[0], np.random.rand(1)[0], np.random.rand(1)[0], 1.0)
    # print(json.dumps(area, indent=2, cls=NpEncoder))

    pcol = (0, 0, 0, area['scale'] / 4)

    sx, sy = elem['geometry'].centroid.coords.xy
    plt.plot(sx, sy, 'ro', ms=2.0, color=(0.0, 0.0, 0.0, 0.35))

    for poly in elem['geometry'].geoms:
        # print(poly.type, poly.is_simple)
        if area['scale'] == 2:
            plot_polygon(ax, poly, color=pcol, edgecolor='red')

    plt.show()
