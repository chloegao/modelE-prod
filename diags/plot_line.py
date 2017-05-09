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

nc = Dataset(sys.argv[1])
record = [1,2,3,4,5,6,7,8,9,10,11,12]
values = nc["seasonalCycle_ts"][:]
title = sys.argv[1].split("/")[-1].replace(".nc","")
plt.figure()
plt.plot(record,values,label="model")
plt.title(title)

nc = Dataset(sys.argv[2])
values = nc["seasonalCycle_ts"][:]
plt.plot(record,values,label="obs")
plt.legend()
plt.savefig("test3.ps")



plt.show()

