#!/usr/bin/env python3

import os
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.validation import make_valid


xmin = -7.0
xmax = 37.0
ymin = 29.0
ymax = 47.0


if __name__ == '__main__':
    #world_countries = gpd.read_file("data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
    #world_countries = gpd.read_file("data/ne_10m_ocean/ne_10m_ocean.shp")


    # world_countries = gpd.read_file("data/GOaS_v1_20211214/goas_v01.shp")
    # world_coastlines = gpd.read_file("data/ne_10m_coastline/ne_10m_coastline.shp")
    # //print(shapefile)
    # bounds = world_coastlines.total_bounds
    # print(bounds)
    #polygon = Polygon([(xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin), (xmin, ymax)])


    #world_marine_regions = gpd.read_file("data/GOaS_v1_20211214/goas_v01.shp")  #;//, mask=polygon_v)
    med_marine_regions = gpd.read_file("data/goas_med/club_med.shp")  #;//, mask=polygon_v)

    # print(world_marine_regions["name"][2])
    # med = world_marine_regions["geometry"][2]

    med_coast = med_marine_regions.cx[xmin:xmax, ymin:ymax]

    #polygon_coast = Polygon([(xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin), (xmin, ymax)])
    #poly_gdf = gpd.GeoDataFrame([1], geometry=[polygon_coast], crs=med_coast.crs)
    #polygon_coast = make_valid(polygon_coast)


    #world_marine_regions = gpd.read_file("data/GOaS_v1_20211214/goas_v01.shp")
    print(med_marine_regions.head(10))
    #
    # for n in world_marine_regions:
    #     print(n)

    #med_clipped = med_coast.clip(polygon_coast, True)

    fig, ax = plt.subplots(figsize=(10, 10))
    #poly_gdf.boundary.plot(ax=ax, color="blue")

    med_coast.plot(
        column='name',
        categorical=True,
        legend=True,
        figsize=(10, 10),
        markersize=25,
        cmap="Set2", ax=ax)

    plt.show()

    #//med_clipped.to_file("data/goas_med/club_med.shp")
    exit()

    med_countries = world_countries.cx[xmin:xmax, ymin:ymax]
    med_coast = world_coastlines.cx[xmin:xmax, ymin:ymax]

    polygon_coast = Polygon([(xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin), (xmin, ymax)])
    poly_gdf = gpd.GeoDataFrame([1], geometry=[polygon_coast], crs=med_coast.crs)
    polygon_countries = Polygon([(xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin), (xmin, ymax)])
    poly_2_gdf = gpd.GeoDataFrame([1], geometry=[polygon_countries], crs=med_countries.crs)

    polygon_countries = make_valid(polygon_countries)

    #poly_gdf_2 = gpd.GeoDataFrame([1], geometry=[polygon], crs=med_coast.crs)

    fig, ax = plt.subplots(figsize=(10, 10))
    poly_gdf.boundary.plot(ax=ax, color="red")

    med_clipped = med_coast.clip(polygon_coast, True)
    med_countries = med_countries.clip(polygon_countries, True)

    med_countries.plot(
        #column='NAME_LONG',
        categorical=True,
        legend=True,
        figsize=(10, 10),
        markersize=25,
        cmap="Set2", ax=ax)

    med_clipped.plot(
        column='scalerank',
        categorical=True,
        legend=False,
        figsize=(10, 10),
        markersize=25,
        cmap="Set2", ax=ax)

    #shapefile.plot(ax=ax)
    plt.show()
