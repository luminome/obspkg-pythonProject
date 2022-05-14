#!/usr/bin/env python3
#/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
import sys
import nctoolkit as nc
#from nctoolkit.api import nctoolkit as nc


# if __name__ == '__main__':

#freeze_support()
data_path_fn = '~/Sites/pythonProject/data/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1649679937092.nc'

ds = nc.open_thredds(data_path_fn)
ds.select(timestep = 0)
ds.to_latlon(lon = [-79.5, 79.5], lat = [0.75, 89.75], res = [1, 0.5])
ds.plot()
