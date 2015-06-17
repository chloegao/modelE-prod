#!/usr/bin/env python
#
# This script add height classes to exiting input files.

import sys
import netCDF4
import argparse
from giss.util import *
import giss.io

# --------------------------------------------------------
# ncin = open netCDF file handle (reading)
def add_required_dims(required_dims, ncin, var_name) :
	var = ncin.variables[var_name]
	required_dims.update(set(var.dimensions))
# --------------------------------------------------------
# Copies a netCDF variable from input to output netCDF file.  Dimensions
# should have been defined already.  Adds a height-classification
# dimension to the beginning.  Copies source variable n times, once
# per height class.
def heightclassify_var(ncin, ncout, var_name) :
	ivar = ncin.variables[var_name]
	dims = ('li_maxheights',) + ivar.dimensions
	ovar = ncout.createVariable(var_name, ivar.dtype, dims)
	for h in range(li_maxheights) :
		ovar[h,:] = ivar[:]

# --------------------------------------------------------
# Check that a variable (var) has dimensions equal to those
# listed in a netCDF file.
# nc = open netCDF file
# var = Numpy (Python) array
# dim_names = Names of dimensions in the netCDF file
def check_dimensions(nc, var, dim_names) :
	mismatch = False
	for i in range(0,len(var.shape)) :
		ncdim = len(nc.dimensions[dim_names[i]])
		vardim = var.shape[i]
		if ncdim != vardim :
			raise Exception('Dimension mismatch, dimension #%d is %d, should be %d (%s)' % (i,vardim,ndim,dim_names[i]))
# --------------------------------------------------------



# Parse Command Line

if len(sys.argv) == 4 :
		gic_fname = sys.argv[1]
		topo_fname = sys.argv[2]
		lic_fname = sys.argv[3]
#elif len(sys.argv) == 2 :
#		gic_fname = sys.argv[1]
#		if (gic_fname[-3:] != '.nc') :
#			print "Input filename '%s' must end in .nc, or output file name must be specified" % (gic_fname)
#			sys.exit()
#		lic_fname = gic_fname[0:-3] + '_LANDICE.nc'
else :
		print 'Usage: %s <GIC(in)> <TOPO(in)> <LIC(out)>'  % (sys.argv[0])
#			'(Output file generated from input if not specified)' % (sys.argv[0])
		sys.exit()

print 'Reading GIC file: %s' % (gic_fname)
print 'Writing LIC file: %s' % (lic_fname)

# =========================================================

# I must write the following variables (indices in C order)
#
# Number of height classes used in each grid cell
#    integer (jm,im) nheight
# Top of each height class (bottom of first height class is sea level)
#    double (li_maxheights,jm,im) height_max
# Portion of each height class in a grid cell. sum pheight(:,j,i) == 1
#    double (li_maxheights,jm,im) pheight
# Mean elevation within each height class, across the grid cell
#    double (li_maxheights,jm,im) mean_elevation
# Height-classified SNOWLI variable
#    double (li_maxheights,jm,im) snowli
# Height-classfied TLANDI variable
#    double (li_maxheights,jm,im,nlayer) tlandi



# Set up variables
li_maxheights=3
#hc_vars_2d = ['snowli', 'tlandi']    # Variables to copy and height-classify
#hc_vars = hc_vars_2d

# Open input and output files
fgic = netCDF4.Dataset(gic_fname, 'r')
flic_out = netCDF4.Dataset(lic_fname, 'w', format='NETCDF3_CLASSIC')

# Query dimensions out of input file
required_dims = set()
add_required_dims(required_dims, fgic, 'snowli')
add_required_dims(required_dims, fgic, 'tlandi')



# Copy over dimensions
flic_out.createDimension('li_maxheights', li_maxheights)
for dim in required_dims :
	flic_out.createDimension(dim, len(fgic.dimensions[dim]));

# Height classify variables
heightclassify_var(fgic, flic_out, 'snowli')
heightclassify_var(fgic, flic_out, 'tlandi')
#snowli.setncattr('long_name', 'Height-classified SNOWLI variable')
#tlandi.setncattr('long_name', 'Height-classfied TLANDI variable')

# -------------------------
# Read ZATMO from GISS-format file, and write it into LIC file
topovars = giss.io.read_vars(topo_fname, ('ZATMO', 'FGICE'))
ZATMO = topovars['ZATMO']
FGICE = topovars['FGICE']
ice_mask = (FGICE > 0)


# Check that dimensions match
check_dimensions(fgic, ZATMO, ('jm','im'))


# ==================================================
# Now we make some dummy height classes, etc.
# ..and store in the netCDF file

nheight = flic_out.createVariable('nheight', 'i4', ('jm', 'im'))
nheight.setncattr('long_name', 'Number of height classes used in each grid cell')
nheight[:,:] = 3 * ice_mask

height_max = flic_out.createVariable('height_max', 'd', ('li_maxheights','jm','im'))
height_max.setncattr('long_name', 'Top of each height class (bottom of first height class is sea level)')
mean_elevation = flic_out.createVariable('mean_elevation', 'd', ('li_maxheights','jm','im'))
mean_elevation.setncattr('long_name', 'Mean elevation within each height class, across the grid cell')

mean_elevation[0,:,:] = .5 * ZATMO[:,:] * ice_mask
height_max[0,:,:] = .75 * ZATMO[:,:] * ice_mask
mean_elevation[1,:,:] = ZATMO[:,:] * ice_mask
height_max[1,:,:] = 1.25 * ZATMO[:,:] * ice_mask
mean_elevation[2,:,:] = 1.5 * ZATMO[:,:] * ice_mask
height_max[2,:,:] = 10000 * ice_mask

pheight = flic_out.createVariable('pheight', 'd', ('li_maxheights','jm','im'))
pheight.setncattr('long_name', 'Portion of each height class in a grid cell. sum pheight(:,j,i) == 1')





pheight[0,:,:] = .25 * ice_mask
pheight[1,:,:] =  .5 * ice_mask
pheight[2,:,:] = .25 * ice_mask

fgic.close()
flic_out.close()

