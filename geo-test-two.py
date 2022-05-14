#!/usr/bin/env python3

import os
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import mapping, Polygon, MultiPolygon
from shapely.validation import make_valid
import fiona
from collections import OrderedDict

xmin = -7.0
xmax = 42.5
ymin = 30.0
ymax = 48.0


def update_goas():
    world_marine_regions = gpd.read_file("data/GOaS_v1_20211214/goas_v01.shp")
    print(world_marine_regions.head(50))
    #                  name  latitude   longitude  min_Y      min_X      max_Y      max_X      area_km2
    #6 Mediterranean Region  38.13065   19.70067   30.06809   -6.03255   47.37640   42.35496   2988248
    med_coast = world_marine_regions.cx[xmin:xmax, ymin:ymax]

    polygon_coast = Polygon([(xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin), (xmin, ymax)])
    poly_gdf = gpd.GeoDataFrame([1], geometry=[polygon_coast], crs=med_coast.crs)
    polygon_coast = make_valid(polygon_coast)

    med_clipped = med_coast.clip(polygon_coast, True)
    print(med_clipped.head(50))

    gpd.GeoDataFrame.from_features

    print('done')


    #med_clipped.to_file("data/goas_med/club_med.shp")
    #print('saved')


def save_geometry(s_schema, s_object, fpath):
    # Define a polygon feature geometry with one attribute
    # s_schema = {
    #     'geometry': 'Polygon',
    #     'properties': {'name': 'str'}
    # }

    # Write a new Shapefile
    with fiona.open(fpath, 'w', 'ESRI Shapefile', s_schema) as c:
        ## If there are multiple geometries, put the "for" loop here
        #s_object['geometry'] = mapping(s_object['geometry'])
        c.write(s_object)

        print('poly saved to:', fpath)


def prepare(index, f_schema, f_object):
    sk = {}
    ok = {}
    for k in f_schema:
        sk[k] = str(type(f_object[k][index]).__class)
        ok[k] = f_object[k][index]

    sk['geometry'] = 'Polygon'
    return sk, ok

    #obj_schema = [(sk[k]= v) for k,v in object[index]]

if __name__ == '__main__':
    # update_goas()
    # exit()
    med_marine_regions = gpd.read_file("data/goas_med/med_only.shp")
    # med_marine_regions = gpd.read_file("data/goas_med/club_med.shp", rows=slice(1,2))
    # med_marine_regions.to_file("data/goas_med/med_only.shp")

    print(med_marine_regions.head())
    print(med_marine_regions.iloc[0])
    exit()



    # envgdf = gpd.GeoDataFrame([1])
    #
    # envgdf = gpd.GeoDataFrame(med_only)
    # envgdf = envgdf.rename(columns={0: 'geometry'}).set_geometry('geometry')
    # envgdf.to_file("data/goas_med/med_only.shp")



    #print()



    exit()



    #schema = {}
    with fiona.open('data/goas_med/club_med.shp', 'r') as f_input:
        schema = f_input.schema.copy()
        print(schema)

    index = 1
    g_object = {}
    pobj = {}
    for n in schema['properties']:
        pobj[n] = med_marine_regions[n][index]

    g_object['properties'] = OrderedDict(pobj)
    g_object['geometry'] = mapping(med_marine_regions['geometry'][index])

    save_geometry(schema, g_object, 'data/goas_med/med-only.shp')

    print(g_object)
    exit()


    print(med_marine_regions.head(10))
    print(med_marine_regions.printSchema())
    schema = list(med_marine_regions)
    #
    p_schema, p_object = prepare(1, schema, med_marine_regions)
    print(p_schema)
    print(p_object)
    #//med_marine_regions = med_marine_regions.simplify(1, True)
    exit()

    #med_poly = med_marine_regions["geometry"][1]

    save_geometry(p_schema, p_object, 'data/goas_med/med-only.shp')
    exit()

    fig, ax = plt.subplots(figsize=(10, 6))

    print('has', med_poly.geom_type, med_poly.area, med_poly.exterior.geom_type, med_poly.exterior.length)

    for n in med_poly.interiors:
        #print(n.geom_type, n.length)
        ns = n.simplify(0.008, True)
        ns = ns.buffer(0.005)
        nf = ns.exterior.simplify(0.002, True)
        plt.plot(*nf.xy, color="blue")

    #mext = med_poly.exterior.simplify(1, True)
    exterior_coords = med_poly.exterior.coords[:]
    print(len(exterior_coords))

    plt.plot(*med_poly.exterior.xy, color="red")

    med_marine_regions.plot(
        column='name',
        categorical=True,
        legend=True,
        figsize=(10, 6),
        markersize=25,
        cmap="Set2", ax=ax)

    plt.show()

    #med_marine_regions.to_file('data/goas_med/dataframe.geojson', driver='GeoJSON')
    #med_clipped.to_file("data/goas_med/club_med.shp")
