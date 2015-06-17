import netCDF4
import os
from modelexe.const import *


nhp = 1		# LISNOW_MODE_LEGACY

def copy_dim(inc, onc, dim_name):
	dim = inc.dimensions[dim_name]
	onc.createDimension(dim.name, len(dim))

input_dir = os.path.join(os.environ['HOME'], 'cmrun')
gic_fname = os.path.join(input_dir, 'GIC.144X90.DEC01.1.ext_1.nc')
landice_ic_fname = os.path.join(input_dir, 'LANDICE_IC.144X90.DEC01.1.ext_1.nc')

print 'Reading ', gic_fname
print 'Writing ', landice_ic_fname
inc = netCDF4.Dataset(gic_fname, 'r')
onc = netCDF4.Dataset(landice_ic_fname, 'w')

copy_dim(inc, onc, 'jm')
copy_dim(inc, onc, 'im')
onc.createDimension('nhp', nhp)
onc.createDimension('nl_max', 4)

dims3 = ('nhp', 'jm', 'im')
dims4 = ('nhp', 'jm', 'im', 'nl_max')

onc.createVariable('nl_snow', 'i', dims3)[:] = 0
onc.createVariable('nl_ice', 'i', dims3)[:] = 4
# No borrowed state for initial condition
onc.createVariable('nl_borrowed', 'i', dims3)[:] = 0

dz_v = onc.createVariable('dz', 'd', dims4)
dz_v.target_depth = 3.0
dz_v[:,:,:,0] = .1
dz_v[:,:,:,1] = 2.9
dz_v[:,:,:,2:] = 4
print 'dz_v', dz_v.shape

wsn_v = onc.createVariable('wsn', 'd', dims4)
wsn_v[:] = dz_v[:] * RHO_ICE
print 'wsn_v', wsn_v.shape

tlandi = inc.variables['tlandi'][:]
print tlandi
hsn_v = onc.createVariable('hsn', 'd', dims4)
print wsn_v.shape
print tlandi.shape
for ihp in xrange(0, nhp):
	# hsn_v[ihp,:] = wsn_v[:] * (tlandi[0,:] * SHW_KG - LHM_KG)
	hsn_v[ihp,:,:,0:2] = wsn_v[ihp,:,:,0:2] * (tlandi[0,:,:,0:2] * SHW_KG - LHM_KG)
	# Copy lower layer to two more layers down.
	hsn_v[ihp,:,:,2] = wsn_v[ihp,:,:,2] * (tlandi[0,:,:,1] * SHW_KG - LHM_KG)
	hsn_v[ihp,:,:,3] = wsn_v[ihp,:,:,3] * (tlandi[0,:,:,1] * SHW_KG - LHM_KG)



onc.close()
inc.close()
