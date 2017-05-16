#!/usr/bin/python3

import numpy as np
from netCDF4 import Dataset
import matplotlib
#matplotlib.use("Agg")
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as  plt 
import sys
import os

vmin = -2
vmax = 30

### figure 1
nc = Dataset(sys.argv[1])
lons = nc["lono"][:]
lats = nc["lato"][:]
values = nc["pot_temp"][:][0]
title = sys.argv[1].split("/")[-1].replace(".nc","")
plt.figure()
lons, lats = np.meshgrid(lons, lats)
m = Basemap(projection='robin', lon_0=0, resolution='c')
x, y = m(lons, lats)
m.drawcoastlines()
plt.contourf(x, y, values,vmax=vmax,vmin=vmin)
plt.colorbar(orientation="horizontal")
plt.title(title)
plt.savefig("map1.ps")

#### figure 2
nc = Dataset(sys.argv[2])
lons = nc["lon"][:]
lats = nc["lat"][:]
values = nc["temp_ann"][:][0]
title = sys.argv[2].split("/")[-1].replace(".nc","")
plt.figure()
lons, lats = np.meshgrid(lons, lats)
m1 = Basemap(projection='robin', lon_0=0, resolution='c')
x, y = m1(lons, lats)
m1.drawcoastlines()
plt.contourf(x, y, values,vmax=vmax,vmin=vmin)
plt.colorbar(orientation="horizontal")
plt.title(title)
plt.savefig("map2.ps")

#### figure 3
nc = Dataset(sys.argv[3])
lons = nc["lon"][:]
lats = nc["lat"][:]
values = nc["diff"][:][0]
title = sys.argv[3].split("/")[-1].replace(".nc","")
plt.figure()
lons, lats = np.meshgrid(lons, lats)
m1 = Basemap(projection='robin', lon_0=0, resolution='c')
x, y = m1(lons, lats)
m1.drawcoastlines()
plt.contourf(x, y, values)
plt.colorbar(orientation="horizontal")
plt.title(title)
plt.savefig("map3.ps")


plt.show()

