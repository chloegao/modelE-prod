import sys
import modele
from modele import f90
from pprint import pprint

# Replace this with calls to ExportConstants.F90
RHOI=916.6
RHOW=1e3

# Set up realistic parmaeters for lisnow
def realistic_params():
	params = modele.Lisnowbase_Mod.Lisnowparams()
	xin = modele.Lisnowbase_Mod.Lisnowin()
	xout = modele.Lisnowbase_Mod.Lisnowout()

	params.max_nl_snow = 0
	params.max_nl_ice = 3

	xin.nl_snow = params.max_nl_snow
	xin.nl_ice = params.max_nl_ice
	params.max_nl_borrowed = 0
	xin.nl_borrowed = params.max_nl_borrowed	# Always

	nl = xin.nl_snow + xin.nl_ice

	# Create layers of ice with no partial melt
	xin.dz[0] = .1
	xin.dz[1:nl] = 15.			# [m]
	for i in xrange(0,nl):
		mass_kg = xin.dz[i] * RHOI		# [kg m-3]
		xin.wsn[i] = mass_kg / RHOW		# [m water equivalent]
		xin.

#TODO: Look up previous all-Ptyhon prototype in exp/ directory


	modele.Lisnowbase_Mod.allocate_snow_adv(params, xin, xout)

	xin.



	return params,xin,xout







# Conversation with Gary Russell, May 7 2015
# 
#         ! Model is:
#         !    nl_snow: 0 or more layers of snow
#         !        (new snow on top of old ice, in lower extent of ice sheet in the winter).
#         !        DO NOT mix snow and ice when regridding
#         !    nl_ice: The usual snow/firn (on upper layers of ice sheet)
#         !        These may be regridded
#         !    nl_borrowed: "Extra" layers below
#         !        Run the heat equation here, even thought we don't "own" this ice.
#         !        And by all means, DO NOT regrid it!
#         integer :: nl_snow
#         integer :: nl_ice
#         integer :: nl_borrowed     ! == params%MAX_NL_BORROWED
# 
#         ! Arrays 1...nl_snow+nl_ice+nl_borrowed
#         double precision, dimension(:), allocatable :: dz
#         double precision, dimension(:), allocatable :: wsn
#         double precision, dimension(:), allocatable :: hsn
# !#ifdef TRACERS_WATER
# !        double precision, dimension(:,:), allocatable :: trsn
# !#endif
# 
#         ! Arrays 1...nl_borrowed
#         ! Ice properties in the borrowed layers
#         double precision, dimension(:), allocatable :: csn_borrowed
#         double precision, dimension(:), allocatable :: ksn_borrowed
# 
#         ! ------------- Input to snow_adv()
#         ! @var srht = solar heat flux [W m-2]
#         ! @var trht = thermal heat flux [W m-2]
#         ! @var snht = sensible heat flux @ beginning of timestep [W m-2]
#         ! @var prht = heat of preciptation [W m-2]
#         double precision :: srht, trht, snht, prht
# 
# sensible heat from surface to atmosphere: 
# latent heat from surface to atmosphere:
# (sign is same as evaporation)
# 
# When evaporation is .5 mm/day, water vapor heat flux is about 20 W/m^2.
# 
# Sensible heat flux in Antarctica and Greenland is always going from air into the ice.  On an annual average, all the ocean and ice things, the flux is the other way.  (Proportional to difference between surface air temp and ground temperature, plus wind speed, stability factor.)  It's also a little interesting in that along Eastern half cost of Antarctica, the heat flux could be 60 W/m^2 into the ice at the edges.  In the interior, it's less than that, like 20 W/m^2 into the ice.
# 
# trht = Thermal radiation (based on blackbody radiation)
# 
# 
#         ! @var t_ground = Temperature of an "extra" layer under our bottom-most layer
#         !      (used to calculate T gradient through bottom)
#         ! @var dz_ground = Thickness of the "extra" layer
#         double precision :: t_ground, dz_ground
# 
#         ! @var evaporation = Evaporation rate [m^3 m-2 s-1 water equiv]
#         double precision :: evaporation
# 
# @ top of ice sheets and center of Antarctica, long-time average evaporation is 0 or possibly negative.  Which means there's dew.  @ bottom of ice sheet, it's more like .5 mm/day (water equivalent).  Interior.. eastern high part of Antarctica, it's negative.  Greenland: in higher parts it's also negative.  Most of it is between 0 and .5 mm/day.  A few edges above .5 mm/day.
# 
# precipitation is about twice that.
# 
#         ! @var pr = precipitation [m water equiv s-1]
#         ! @var dt = timestep [s]
#         double precision :: pr, dt
# 
#         ! @var snsh_dt = d[snht] / dT         (Used for implicit scheme)
#         ! @var evap_dt = d[evaporation] / dT  (Used for implicit scheme)
#         double precision :: snsh_dt, evap_dt
# 
#         ! @var MA1 = mass of lowest atmospheric layer [kg m-2]
#         double precision :: ma1
# 
# gary russell
# 984 at bottom of atm
# 964 mb
# density is just a little above 1
# pressure --> pascals (2000 Pa) / g --> 200 kg/m^2
# 984mb = 98000/g = 10000 kg/m^2
# 
# 
#         ! Q = specific humidity [kg water vapor/kg air]
#         !     (This is the same as mixing ratio)
#         ! Q1 = Specific humidity of the bottom atmosphere layer
#         double precision :: Q1
# 
# QSurf = 2e-3: boundary in plot
# Boundary of that plot goes right along the northern edge of Alaska, across Southern tip of Greenland.  Specific Humidity decreases with pressure.  At top of ice sheet, it will be ~1/4 of that.
# 
#         ! Conductivity and heat capacity of soil below the snow
#         double precision k_ground, c_ground
# 
# 
