import netCDF4
import os
this_script_dir = os.path.dirname(os.path.realpath(__file__))

# Retrieve relevant ModelE constants
# Replace this with calls to ExportConstants.F90
nc = netCDF4.Dataset(os.path.join(this_script_dir, 'modele_constants.nc'))
ncv = nc.variables['constants'].__dict__
#RHO_WATER = ncv['constant::rhow']		# [kg m-3]
RHO_ICE = ncv['constant::rhoi']		# [kg m-3]
LHM_KG = ncv['constant::lhm']			# [J kg-1]
LHE_KG = ncv['constant::lhe']			# [J kg-1]
SHI_KG = ncv['constant::shi']			# [J kg-1 K-1]
SHW_KG = ncv['constant::shw']			# [J kg-1 K-1]
TFRZ = ncv['constant::tf']				# [K]
SIGMA = ncv['constant::stbo']			# Stefan-Boltzmann const [W m-2 K-4]
nc.close()

# physical parameters (calculated from constants)
#LAT_FUSION = LHM_KG * RHO_WATER # [J m-3]
#LAT_EVAP = LHE_KG * RHO_WATER # [J m-3]
