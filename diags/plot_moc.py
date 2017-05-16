#!/usr/bin/python3

import numpy as np
from netCDF4 import Dataset
import matplotlib
import numpy.ma as ma
#matplotlib.use("Agg")
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as  plt 
import sys
import os

#vmin = -2
#vmax = 30

### figure 1
nc = Dataset(sys.argv[1])
xs = nc["lato2"][:]
ys = nc["zoce"][:]
values = nc["sf_Atl"][:]
error_value = -6.849315e+26
values = ma.masked_values(values, error_value)
ma.set_fill_value(values, 0)
title = sys.argv[1].split("/")[-1].replace(".nc","")
plt.figure()
X, Y = np.meshgrid(xs, ys)
plt.contourf(X, Y, values)
plt.gca().invert_yaxis()
plt.colorbar()
#title = sys.argv[1].split("/")[-1].replace(".nc","")
plt.title(title)
plt.savefig("mocsection1.ps")

### figure 2
nc = Dataset(sys.argv[2])
xs = nc["lato2"][:]
ys = nc["zoce"][:]
values = nc["sf_Pac"][:]
error_value = -6.849315e+26
values = ma.masked_values(values, error_value)
ma.set_fill_value(values, 0)
title = sys.argv[2].split("/")[-1].replace(".nc","")
plt.figure()
X, Y = np.meshgrid(xs, ys)
plt.contourf(X, Y, values)
plt.gca().invert_yaxis()
plt.colorbar()
#title = sys.argv[1].split("/")[-1].replace(".nc","")
plt.title(title)
plt.savefig("mocsection2.ps")

plt.show()


