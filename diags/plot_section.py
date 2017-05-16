#!/usr/bin/python3

import numpy as np
from netCDF4 import Dataset
import matplotlib
#matplotlib.use("Agg")
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as  plt 
import sys
import os

#vmin = -2
#vmax = 30

### figure 1
nc = Dataset(sys.argv[1])
xs = nc["lat"][:]
ys = nc["zoc"][:]
values = nc["pot_tempAtl"][:]
error_value = -1.e30
values = ma.masked_values(values, error_value)
ma.set_fill_value(values, 0)
title = sys.argv[1].split("/")[-1].replace(".nc","")
plt.figure()
X, Y = np.meshgrid(xs, ys)
plt.contourf(X, Y, values)
plt.gca().invert_yaxis()
plt.colorbar()
plt.title(title)
plt.savefig("section1.ps")

### figure 2
nc = Dataset(sys.argv[1])
xs = nc["lat"][:]
ys = nc["zoc"][:]
values = nc["pot_tempPac"][:]
error_value = -1.e30
values = ma.masked_values(values, error_value)
ma.set_fill_value(values, 0)
title = sys.argv[1].split("/")[-1].replace(".nc","")
plt.figure()
X, Y = np.meshgrid(xs, ys)
plt.contourf(X, Y, values)
plt.gca().invert_yaxis()
plt.colorbar()
plt.title(title)
plt.savefig("section2.ps")

plt.show()


