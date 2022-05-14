import netCDF4 as nc
from datetime import datetime
from datetime import timedelta
# :units = "minutes since 1900-01-01 00:00:00";

import matplotlib.pyplot as plt
import numpy as np
import math

date_time_str = '1900-01-01 00:00:00'
date_time_origin = datetime.strptime(date_time_str, '%Y-%m-%d %H:%M:%S')

fn = 'download.nc'
# fn = '../1999-6-10_med-cmcc-tem-rean-d_1638647808820.nc'
ds = nc.Dataset(fn)

# print(ds.__dict__)
# print(ds)

for dim in ds.dimensions.values():
    print(dim)

for var in ds.variables.values():
    print(var)

U = ds['u10']
V = ds['v10']
#dat = np.hstack((U, V))
dat = np.concatenate((U[:, :, 0:], V[:, :, 0:]))
print(dat.shape)
print(dat[0,0,0])
print(U[0,0,0])
print(V[0,0,0])

a = np.array([U], dtype=np.str)
b = np.array([V], dtype=np.str)

c = np.char.add(a,' ').astype(str)
d = np.char.add(c,b).astype(str)[0]
print(d.shape)
# data_min = np.nanmin(ds['sst'][4])
# data_max = np.nanmax(ds['sst'][4])
#
# print(data_min, data_max)
#
# t = ds['depth'][0:]
# def d_to_int(z):
#     return int(z)   #math.ceil(z/10)  #;//-273.15
#
# dat = np.vectorize(d_to_int)(t)
#
# print(dat)
# exit()
exit()


xt = float(t[0])
delta_readable = date_time_origin + timedelta(minutes=xt)
d = {'t': "%.2f" % xt, 'ts': str(delta_readable), 'data': []}
print("data timestamp: ", d['ts'])

exit()



def k_to_deg_c(z):
    return z-273.15


dat = ds['sst'][:, :, :]
lat = ds['latitude'][:]
lon = ds['longitude'][:]

dat = np.vectorize(k_to_deg_c)(dat)

print(dat)

# print(dat)
# print(lat)
# print(lon)

plt.imshow(dat, interpolation="none")
plt.show()
