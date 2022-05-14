#!/usr/bin/env python3

import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon, Point, box, LineString, LinearRing
import csv
import json
import math
import fiona
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection

# #// https://mapamed.org/index.php?language=en


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

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


country_file = '../data/MAPAMED_2019_edition/mapamed_2019_dataset_definition_countries.tsv'
with open(country_file) as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    country_data = [row for row in reader]

shapes_file = '../data/MAPAMED_2019_edition/MAPAMED_2019_spatial_data_epsg3035.gpkg'
tables = {}


for layer_name in fiona.listlayers(shapes_file):
    print(layer_name)
    g_layer = gpd.read_file(shapes_file, layer=layer_name)
    # #// https://gis.stackexchange.com/questions/24340/converting-coordinates-from-meters-to-decimal-degrees-in-qgis
    g_layer = g_layer.to_crs(crs=4326)  # 'EPSG:3857')  #, allow_override=True)
    tables[layer_name] = g_layer



print(json.dumps(country_data, indent=4))

fig, ax = plt.subplots()

def normalize_val(val, mi, ma):
    return (val - mi) / (ma - mi)


area_filter = ["NAME", "STATUS_YR", "STATUS_ENG", "MAPAMED_ID", "PARENT_ID", "REP_AREA", "SITE_TYPE_ENG", "IUCN_CAT_ENG", "WEBSITE"]
criteria_filter = ['Marine N2000 Proposed Site', 'MPA of National Statute']

for k, v in tables.items():

    if k == 'Scope of the Barcelona Convention (IHO-MSFD)':
        regions_all = v.to_dict('records')
        for d in regions_all:
            d.pop('geometry')

        print(json.dumps(regions_all, indent=4))

    if k == 'MAPAMED 2019 edition - EPSG:3035':
        indx = v.shape[0]
        counted = 0
        mct = [v.iloc[i] for i in range(0, indx) if v.iloc[i]['DESIG_CAT_ENG'] in criteria_filter]
        N = len(mct)

        max_area = 0
        avg_area = 0
        # #// Total marine extent of the site as officially declared (in km2).
        for mv in mct:
            avg_area += mv['REP_AREA'] if not math.isnan(mv['REP_AREA']) else 0.0
            if mv['REP_AREA'] > max_area:
                max_area = mv['REP_AREA']

        avg_area /= len(mct)

        print(avg_area, 'km')
        print(max_area, 'km')

        # exit()

        for i, elem in enumerate(mct):
            pcol = (np.random.rand(1)[0], np.random.rand(1)[0], np.random.rand(1)[0], 1.0)
            area = {}
            for f in area_filter:
                area[f] = elem[f]

            cx, cy = elem['geometry'].centroid.coords.xy

            nsca = normalize_val(elem['REP_AREA'], 0.1, avg_area)
            area['scale'] = 1 if math.isnan(nsca) or nsca < 0 else math.ceil(nsca) if nsca < 4 else 4

            area['CENTROID'] = np.array(elem['geometry'].centroid.coords[0])
            area['COUNTRY'] = [p['COUNTRY_ENG'] for p in country_data if p['ISO3'] in elem['ISO3'][1:-1]]
            area['MED_REGION'] = [p['NAME_ENG'] for p in regions_all if p['MSFD_REGION'] in elem['MSFD_REGION'][1:-1]]

            print(area['REP_AREA'], area['scale'])

            #print(json.dumps(area, indent=2, cls=NpEncoder))

            pcol = (0, 0, 0, area['scale']/4)

            sx, sy = elem['geometry'].centroid.coords.xy
            plt.plot(sx, sy, 'ro', ms=2.0, color=(0.0, 0.0, 0.0, 0.35))

            for poly in elem['geometry'].geoms:
                # print(poly.type, poly.is_simple)
                if area['scale'] == 2:
                    plot_polygon(ax, poly, color=pcol, edgecolor='red')

        # for sk, sv in v.items():
        #
        #     for
        #     print(len(sv))
        #     #print(v['MSFD_REGION'])

        # for sk, sv in v.items():
        #
        #
        #     print(sv)



#         for i, region in enumerate(v['geometry']):
#             # nat_geom = data[region]['geometry']
#             print(region.type, len(region.geoms), i, '\n')
#             for poly in region.geoms:
#                 # plt.plot(*poly.exterior.xy, linewidth=1, color=(0.0, 0.0, 0.0, 0.35))
#                 plot_polygon(ax, poly, facecolor=(0.0, 0.0, 0.0, 0.25), edgecolor='red')
#
#
plt.show()
exit()



# data = gpd.read_file("./MAPAMED_2019_edition/MAPAMED_2019_spatial_data_epsg3035.gpkg", rows=10)
#print(len(data), data.head(50))

for k, v in tables.items():
    print(k, type(v))
    print(type(v), v.shape)  #info(verbose=True))

    if k == 'MAPAMED 2019 edition - EPSG:3035':
        indx = v.shape[0]
        counted = 0

        # mct = [[k, vn] for k, vn in v.items() ]
        mct = [v.iloc[i] for i in range(0, indx) if v.iloc[i]['DESIG_CAT_ENG'] == 'MPA of National Statute']
        N = len(mct)

        # define the colormap
        cmap = plt.cm.jet
        # extract all colors from the .jet map
        cmaplist = [cmap(i) for i in range(cmap.N)]
        # create the new map
        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

        print(cmap)
        # exit()

        for i, elem in enumerate(mct):

            # print(counted,
            #       elem['MAPAMED_ID'],
            #       i,
            #       elem['NAME'].title(),
            #       elem['ISO3'][1:-1])

            pcol = (np.random.rand(1)[0], np.random.rand(1)[0], np.random.rand(1)[0], 1.0)

            # print("country:", [p for p in country_data if p['ISO3'] == elem['ISO3'][1:-1]])
            counted += 1
            # poly = elem['geometry']
            # print(elem['geometry'].type)
            # print(i, type(i))
            it = math.floor(((i+1)/N)*100)
            print("\r{}% complete".format(it), end='')

            sx, sy = elem['geometry'].centroid.coords.xy
            plt.plot(sx, sy, 'ro', ms=2.0, color=(0.0, 0.0, 0.0, 0.35))

            for poly in elem['geometry'].geoms:
                # print(poly.type, poly.is_simple)
                plot_polygon(ax, poly, color=pcol, edgecolor='red')  # , facecolor=(0.0, 0.0, 0.0, 0.35), edgecolor='red')

                #plt.plot(*poly.exterior.xy, label=elem['DESIG'].title(), c=cmap(i))  # , linewidth=1, color=(0.0, 0.0, 0.0, 0.35))  #, label=elem['SITE_TYPE_ENG'].title())

        # exit()
        #
        #
        #
        # for i in range(0, indx):
        #     elem = v.iloc[i]
        #
        #     if elem['DESIG_CAT_ENG'] == 'MPA of National Statute':
        #         print(counted,
        #               elem['MAPAMED_ID'],
        #               i,
        #               elem['NAME'].title(),
        #               elem['ISO3'][1:-1])
        #
        #         print("country:", [p for p in country_data if p['ISO3'] == elem['ISO3'][1:-1]])
        #         counted += 1
        #         # poly = elem['geometry']
        #         print(elem['geometry'].type)
        #         for poly in elem['geometry'].geoms:
        #             print(poly.type)
        #             plt.plot(*poly.exterior.xy, label=elem['DESIG'].title(), cmap=cmap)  #, linewidth=1, color=(0.0, 0.0, 0.0, 0.35))  #, label=elem['SITE_TYPE_ENG'].title())

plt.show()

        # if 'NAME' in v.iloc[i]:
        #
        # print(v.iloc[i]['NAME_ENG'])

    # break

    # ks = v.keys()[0]
    # print(v[ks], len(v[ks]))
    # for en in v:
    #     print(en, type(v[en]), len(v[en]))

    # print(v, len(v))

    # for i, region in enumerate(v['geometry']):
    #     # nat_geom = data[region]['geometry']
    #     print(region.type, len(region.geoms), i, '\n')
    #     for poly in region.geoms:
    #         plt.plot(*poly.exterior.xy, linewidth=1, color=(0.0, 0.0, 0.0, 0.35))
    #         #plot_polygon(ax, poly)  #, facecolor=(0.0, 0.0, 0.0, 0.35), edgecolor='red')

    # # sx, sy = my_data[:, 0], my_data[:, 1]
    # # plt.plot(sx, sy, 'ro', ms=2.0, color=(0.0, 0.0, 0.0, 0.35))
    # plt.plot(*region..xy, linewidth=1, color=(0.0, 0.0, 0.0, 0.35))

# plt.legend()
#
#
#

