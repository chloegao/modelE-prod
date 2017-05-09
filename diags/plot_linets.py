#!/usr/bin/python3

import numpy as np
from netCDF4 import Dataset
import matplotlib
#matplotlib.use("Agg")
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as  plt 
import sys
import os

#filename =  "/discover/nobackup/aromanou/TESTNCO/pot_temp.ANN2191-2200.map_lev1.E190F40oQ32.nc"
#nc = Dataset("/discover/nobackup/aromanou/TESTNCO/temp_mon.seasonalCycle_lev1.WOA13_AnnMon_onEgrid.nc")


nc = Dataset(sys.argv[2])
record = [1,2,3,4,5,6,7,8,9,10]
values = nc[sys.argv[1]][:]
title = sys.argv[2].split("/")[-1].replace(".nc","")
plt.figure()
plt.plot(record,values,label="model")
plt.title(title)

plt.savefig("glbavg_ts.ps")

plt.show()

