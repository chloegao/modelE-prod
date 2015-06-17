import unittest
from fexception import *
import sys
import modele
from modele_exe import f90
from modele_exe.const import *
import traceback
import numpy as np
import gc

#count = 0

def hsn_to_tsnisn(hsn, wsn):
	tsn = np.zeros(len(hsn))
	isn = np.zeros(len(hsn))
	for i in xrange(0,len(hsn)):
		wlhm = wsn[i] * LHM_KG
		if (hsn[i] > 0.):
			tsn[i] = (hsn[i] + wsn[i] * LHM_KG) / (wsn[i] * SHW_KG)
			isn[i] = 0
		elif ( hsn[i] > -wlhm ):
			tsn[i] = 0.
			isn[i] = -hsn[i] / wlhm
		else:
			isn[i] = 1.
			# This assumes the material is 100% frozen (of some density)
			#tsn(n) = &			# [K]
			#   (hsn(n) / wsn(n)  + LHM_KG)   # [J m-2] / [kg m-2] --> [J kg-1] Specific enthalpy
			#   / csn(n)     # / [J kg-1 K-1] = [K]

			# Multiply top and bottom by wsn(n), to remove one floating point division.
			tsn[i] = (hsn[i] + LHM_KG * wsn[i]) / (SHI_KG * wsn[i])

	return tsn,isn

def make_sample_column():

	params = modele_exe.Lisnowbase_Mod.Lisnowparams()
	xin = modele_exe.Lisnowbase_Mod.Lisnowin()
	xout = modele_exe.Lisnowbase_Mod.Lisnowout()

	params.max_nl_snow = 0
	params.max_nl_ice = 3

	xin.nl_snow = params.max_nl_snow
	xin.nl_ice = params.max_nl_ice
	params.max_nl_borrowed = 1			# Simple Dirichlet condition
	xin.nl_borrowed = params.max_nl_borrowed	# Always

	nl_owned = xin.nl_snow + xin.nl_ice

	modele_exe.Lisnowbase_Mod.allocate_snow_adv(params, xin, xout)

	# Create layers of ice with no partial melt
	xin.dz[0] = .1
	xin.dz[1:nl_owned] = 15.			# [m]
	for i in xrange(0,nl_owned):
		xin.wsn[i] = xin.dz[i] * RHO_ICE 		# [kg m-2]
		xin.hsn[i] = -LHM_KG * xin.wsn[i]	# 0C and fully frozen

	# Borrowed layer at bottom
	xin.dz[nl_owned] = 40
	xin.wsn[nl_owned] = xin.dz[nl_owned] * RHO_ICE
	xin.hsn[nl_owned] = -(LHM_KG + SHI_KG) * xin.wsn[nl_owned]	# -1C

	xin.csn_borrowed[:] = SHI_KG		# [J kg-1 K-1]
	xin.ksn_borrowed[:] = 3.22e-6 * RHO_ICE*RHO_ICE

	# Close to equilibrium
	xin.srht = 249.8	# [W m-2]

	# GR: Sensible heat flux in Antarctica and Greenland is always going
	# from air into the ice.  On an annual average, all the ocean and ice
	# things, the flux is the other way.  (Proportional to difference
	# between surface air temp and ground temperature, plus wind speed,
	# stability factor.)  It's also a little interesting in that along
	# Eastern half cost of Antarctica, the heat flux could be 60 W/m^2 into
	# the ice at the edges.  In the interior, it's less than that, like 20
	# W/m^2 into the ice.
	xin.snht = 60	# [W m-2]

	# Decide on how much precipitation and evaporation
	# Aim for a grid cell near the edge of Greenland
	# GR: When evaporation is .5 mm/day, water vapor heat
	#     flux is about 20 W/m^2.
	xin.evaporation_rate = -.005 / 86400.		# [.5 mm/day --> m/s]
	xin.pr_rate = .01 / 86400.				# [1 mm/day --> m/s water equivalent]
	T_precip = -1.
	if T_precip <= 0:
		# Wm-2  = m s-1       K     J kg-1 K-1     kg m-3
		xin.prht = xin.pr_rate * (T_precip * SHI_KG - LHM_KG)
	else:
		raise ValueError('TODO: Add calc for T_precip > 0')

	# Ice-atmosphere interaction runs at 30-minute timesteps
	xin.dt = 1800	# [s]

	# These are too hard to calculate for isolated model
	xin.snsh_dt = 0
	xin.evap_dt = 0


	# The lowest atmospheric layer in ModelE is 20 mbar thick (984 -- 964 mbar)
	# Pressure at bottom of atmosphere = 1000 mbar ~= 100 KPa = 1e5 Pa
	# Next 20 mbar up of air: 20 mbar ~= 2 KPa
	# 2 KPa / g ~= 200 kg m-2      (g = 9.8 [m s-2])
	xin.ma1 = 200		# [kg m-2]

	# QSurf = 2e-3: boundary in plot
	# Boundary of that plot goes right along the northern edge of Alaska,
	# across Southern tip of Greenland.  Specific Humidity decreases with
	# pressure.  At top of ice sheet, it will be ~1/4 of that.
	xin.q1 = 2e-3

	return (params, xin, xout)

class LISnowBaseTestCase(unittest.TestCase):

	def test_tester(self):
		"""Test that we can compare input and output."""

		params, xin, xout = make_sample_column()

		xout.dz[:] = xin.dz[:]
		xout.wsn[:] = xin.wsn[:]
		xout.hsn[:] = xin.hsn[:]
		xout.nl_snow = xin.nl_snow
		xout.nl_ice = xin.nl_ice
		xout.nl_borrowed = xin.nl_borrowed

		inl_owned = xin.nl_snow + xin.nl_ice
		onl_owned = xout.nl_snow + xout.nl_ice

		self.assertAlmostEqual(
			np.sum(xin.dz[:inl_owned]),
			np.sum(xout.dz[:onl_owned]))

		self.assertAlmostEqual(
			np.sum(xin.wsn[:inl_owned]),
			np.sum(xout.wsn[:onl_owned]))

		self.assertAlmostEqual(
			np.sum(xin.hsn[:inl_owned]),
			np.sum(xout.hsn[:onl_owned]))

		# Try to change it up a bit
		xout.dz[0] += .2
		self.assertNotEqual(
			np.sum(xin.dz[:inl_owned]),
			np.sum(xout.dz[:onl_owned]))

		xout.dz[1] -= .2
		self.assertEqual(
			np.sum(xin.dz[:inl_owned]),
			np.sum(xout.dz[:onl_owned]))


	def test_snow_redistr(self):
		"""Test that snow gets redistributed properly without affecting
		any of the conservation properties."""
		params, xin, xout = make_sample_column()

		inl_owned = xin.nl_snow + xin.nl_ice

		xout.dz[:] = xin.dz[:]
#		xout.dz[0] = 3
		xout.wsn[:] = xin.wsn[:]
		xout.hsn[:] = xin.hsn[:]
		xout.nl_snow = xin.nl_snow
		xout.nl_ice = xin.nl_ice
		xout.nl_borrowed = xin.nl_borrowed

	
		fexec(lambda: modele_exe.Lisnowbase_Mod.
			snow_redistr(params, xout.dz, xout.wsn, xout.hsn,
				xout.nl_snow, xout.nl_ice, xout.nl_borrowed, 1.0)
		)

		onl_owned = xout.nl_snow + xout.nl_ice
#		print xout

		# Top layer should always be .1m thick
		self.assertAlmostEqual(.1, xout.dz[0])

		# Check conservation properties
		self.assertAlmostEqual(
			np.sum(xin.dz[:inl_owned]),
			np.sum(xout.dz[:onl_owned]))

		self.assertAlmostEqual(
			np.sum(xin.wsn[:inl_owned]),
			np.sum(xout.wsn[:onl_owned]))

		self.assertAlmostEqual(
			np.sum(xin.hsn[:inl_owned]),
			np.sum(xout.hsn[:onl_owned]))

	def test_heat_eq_melt(self):
		print ' ================= test_heat_eq()'
		params, xin, xout = make_sample_column()
		nl = len(xin.dz)
		xin.dz[:] = 1.
		xin.dz[0] = .1
		xin.dz[-2] = 100.		# Thick layers, prevent heat flow down
		xin.dz[-1] = 100.		# in spite of Dirichlet B.C.
		xin.wsn[:] = xin.dz[:] * RHO_ICE

		tsn = np.zeros(nl)
		tsn[:] = -1.
		tsn[0] = 0.				# Top layer ready to melt a bit
		xin.hsn[:] = (-LHM_KG + tsn[:] * SHI_KG) * xin.wsn[:]

		csn = np.zeros(nl)
		csn[:] = SHI_KG
		ksn = np.zeros(nl)
		ksn[:] = 2.24		# [W m-1 K-1]

		# A bit of energy in at the top
		hfluxes_in = np.array(3.e3)
		hfluxes_in_deriv = np.array(1.e0)
		hfluxes_corr_factor = np.array(0.)
		nfluxes = np.array(1)

		dt = np.array(1000.0)

		hsn0 = np.array(xin.hsn)

		fluxes_in = np.sum(hfluxes_in)
		fluxes_in_deriv = np.sum(hfluxes_in_deriv)
		e = ftry(lambda: modele_exe.Lisnowbase_Mod.
			heat_eq(
				xin.dz, xin.wsn, xin.hsn, csn, ksn, nl,
				fluxes_in, fluxes_in_deriv, hfluxes_corr_factor, dt))
		if e:
			print e
		hfluxes_corr = hfluxes_corr_factor * hfluxes_in_deriv
		print 'hfluxes_corr_factor', hfluxes_corr_factor
		print 'hfluxes_corr', hfluxes_corr

		hsn1 = np.array(xin.hsn)

		self.assertAlmostEqual(1.0, 		# Relative y decimal places
			(sum(hsn0) + (fluxes_in + np.sum(hfluxes_corr)) * dt) /
			(np.sum(hsn1)))

	def check_conservation_energy(self, hsn0, hsn1, hfluxes_in, hfluxes_deriv, flux_corr_factor, dt):
		hfluxes_out = hfluxes_in + flux_corr_factor * flux_corr_factor
		self.assertAlmostEqual(1.0, 		# Relative y decimal places
			(sum(hsn0) + np.sum(hfluxes_out) * dt) /
			(np.sum(hsn1)), places=6)



	def test_pass_all_water(self):
		# Set up column with pure ice and no tolerance for
		# partial meltwater.  Then add water at 0C to it.
		# All the water added should end up coming out the bottom.
		params, xin, xout = make_sample_column()
		params.max_fract_water = 0
		xin.hsn[:] = -LHM_KG * xin.wsn[:]	# 0C and fully frozen

		xout.dz[:] = xin.dz[:]
		xout.wsn[:] = xin.wsn[:]
		xout.hsn[:] = xin.hsn[:]
		xout.nl_snow = xin.nl_snow
		xout.nl_ice = xin.nl_ice
		xout.nl_borrowed = xin.nl_borrowed

		water_in = 1000.	# [kg m-2] Water added to top
		water_T = 0		# [degC]
		heat_in = water_in * SHI_KG * water_T
		water_out = np.array(-1.)	# Allow for intent(inout) parameters
		heat_out = np.array(-1.)

#		print xout
#		print 'BEGIN pass_water', water_in, heat_in, water_out, heat_out
		fexec(lambda: modele_exe.Lisnowbase_Mod.
			pass_water(params, xout.wsn, xout.hsn, xout.dz,
				xout.nl_snow + xout.nl_ice,
				water_in, heat_in, water_out, heat_out)
		)

		# The result: all water should have passed through,
		# nothing else changes

		# The result: all water should have passed through,
		# nothing else changes
		self.assertAlmostEqual(heat_in, heat_out)
		self.assertAlmostEqual(water_in, water_out)
		np.testing.assert_almost_equal(xin.dz, xout.dz)
		np.testing.assert_almost_equal(xin.wsn, xout.wsn)
		np.testing.assert_almost_equal(xin.hsn, xout.hsn)

	def test_pass_to_layer2(self):
		# Set up column with pure ice and no tolerance for
		# partial meltwater.  Then add water at 0C to it.
		# All the water added should end up coming out the bottom.
		params, xin, xout = make_sample_column()
		params.max_fract_water = 0
		xin.hsn[:] = -LHM_KG * xin.wsn[:]	# 0C and fully frozen
		# -200C layer, catch all the water we throw at it.
		xin.hsn[1] = (-LHM_KG - 200. * SHI_KG) * xin.wsn[1]


		xout.dz[:] = xin.dz[:]
		xout.wsn[:] = xin.wsn[:]
		xout.hsn[:] = xin.hsn[:]
		xout.nl_snow = xin.nl_snow
		xout.nl_ice = xin.nl_ice
		xout.nl_borrowed = xin.nl_borrowed

		water_in = 1.	# [m] Water added to top
		water_T = 0		# [degC]
		heat_in = water_in * SHI_KG * water_T
		water_out = np.array(-1.)	# Allow for intent(inout) parameters
		heat_out = np.array(-1.)

#		print xout
#		print 'T=', (xout.hsn[:] / xin.wsn[:] + LHM_KG) / SHI_KG
#		print 'BEGIN pass_water', water_in, heat_in, water_out, heat_out
		fexec(lambda: modele_exe.Lisnowbase_Mod.
			pass_water(params, xout.wsn, xout.hsn, xout.dz,
				xout.nl_snow + xout.nl_ice,
				water_in, heat_in, water_out, heat_out)
		)

#		print 'T=', (xout.hsn[:] / xin.wsn[:] + LHM_KG) / SHI_KG


		self.assertEqual(0., heat_out)
		self.assertEqual(0., water_out)

		# -------------------------
        # The result: all water input should have frozen in second layer

		# Input has enthalpy of 0, so no change anywhere
		np.testing.assert_almost_equal(xin.hsn[:], xout.hsn[:])

		eql = [0,2,3]		# These layers should remain unchanged.
		np.testing.assert_almost_equal(xin.dz[eql], xout.dz[eql])
		np.testing.assert_almost_equal(xin.wsn[eql], xout.wsn[eql])

		# Check layer 1, where all the water should have frozen.
		self.assertAlmostEqual(xin.wsn[1] + water_in, xout.wsn[1])
		self.assertAlmostEqual(xin.dz[1] + water_in /RHO_ICE, xout.dz[1])
	# ---------------------------------------------
	def test_precip_and_heat(self):
		print '-------------------------- test_precip_and_heat()'
		params, xin, xout = make_sample_column()
		xdebug = modele_exe.Lisnowbase_Mod.Debug_Precip_And_Heat()

		xout.dz[:] = xin.dz[:]
		xout.wsn[:] = xin.wsn[:]
		xout.hsn[:] = xin.hsn[:]

		csn = np.array([ 2060.,  2060.,  2060.,  2060.])
		ksn = np.array([ 2.7053009,  2.7053009,  2.7053009,  2.7053009])

		hfluxes_in = np.array([  2.00000000e+02,  -3.15665370e+02,   6.00000000e+01, -1.44675926e-01])
		hfluxes_in_deriv = np.array([ 0.        , -4.62259374,  0.        ,  0.        ])
		flux_corr_factor = np.array(0.)
		water_out = np.array(0.)
		heat_out = np.array(0.)

		nl_owned = xin.nl_snow + xin.nl_ice
		nl_all = nl_owned + xin.nl_borrowed
		evap_rate = np.array(xin.evaporation_rate)

		xdebug.stop_at=5
		print 'Setting stop_at={}'.format(xdebug.stop_at)
		fexec(lambda: modele_exe.Lisnowbase_Mod.
			precip_and_heat(params, xin.dt, \
				xout.dz, xout.wsn, xout.hsn, csn, ksn, \
				nl_owned, nl_all, \
				hfluxes_in, hfluxes_in_deriv, xin.pr_rate, xin.prht, evap_rate, \
				flux_corr_factor, water_out, heat_out, xdebug))

		# ----- Conservation of mass
		wsn_t0 = np.sum(xin.wsn) + (xin.pr_rate + evap_rate) * xin.dt
		wsn_t1 = np.sum(xout.wsn)
		self.assertAlmostEqual(1.0, wsn_t1 / wsn_t0)

		# ----- Conservation of energy
		print 'flux_corr_factor', flux_corr_factor
		hfluxes_out = hfluxes_in + hfluxes_in_deriv * flux_corr_factor
		hsn_t0 = np.sum(xin.hsn) + xin.dt * (xin.prht + np.sum(hfluxes_out)) + heat_out
		hsn_t1 = np.sum(xout.hsn)
		self.assertAlmostEqual(1.0, hsn_t1 / hsn_t0)


	# ---------------------------------------------
	def test_snow_adv3(self):
		print '=================== test_snow_adv3()'
		params, xin, xout = make_sample_column()
		xdebug = modele_exe.Lisnowbase_Mod.Debug_Snow_Adv()

		# Top two layers less dense
		xin.dz[0] = xin.wsn[0] / (RHO_ICE * .8)
		xin.dz[1] = xin.wsn[1] / (RHO_ICE * .9)

		# ------ Adjust insolation to be in equilibrium
		xdebug.stop_at = 0
		e = ftry(lambda: modele_exe.Lisnowbase_Mod.
			snow_adv(params, xin, xout, xdebug))
		if e:
			self.assertEqual(('stop_if', 0), tuple(e))
		xin.srht -= np.sum(xdebug.hfluxes0)
		srht_equilib = xin.srht

		# ------ Run it
		for equilib_diff in [-5., -1., -.1, 0., .1, 1., 5.]:	# [W m-2]
			xin.srht = srht_equilib + equilib_diff

			xin_tsn, xin_isn = hsn_to_tsnisn(xin.hsn, xin.wsn)
			xdebug.stop_at = -1
			print 'Setting xdebug.stop_at = ', xdebug.stop_at
			e = ftry(lambda: modele_exe.Lisnowbase_Mod.
				snow_adv(params, xin, xout, xdebug))
			if e:
				self.assertEqual(('stop_if', 5), tuple(e))

			ihsn0_total = np.sum(xin.hsn) + \
				xin.dt * np.sum(xdebug.hfluxes0 + xdebug.flux_corr_factor0 * xdebug.hfluxes0_deriv)
			ohsn0_total = np.sum(xout.hsn)
			tsn3, isn3 = hsn_to_tsnisn(xdebug.hsn3, xdebug.wsn3)
			xout_tsn, xout_isn = hsn_to_tsnisn(xout.hsn, xout.wsn)

			self.check_conservation_energy(xin.hsn, xdebug.hsn3, xdebug.hfluxes0, xdebug.hfluxes0_deriv, xdebug.flux_corr_factor0, xin.dt)

			# Check conservation on snow_redistr()
			self.assertAlmostEqual(1.0, np.sum(xdebug.hsn4) / np.sum(xdebug.hsn3))
			self.assertAlmostEqual(1.0, np.sum(xdebug.wsn4) / np.sum(xdebug.wsn3))
			self.assertAlmostEqual(1.0, np.sum(xdebug.dz4) / np.sum(xdebug.dz3))
			self.assertAlmostEqual(.1, xdebug.dz4[0])

			# Make sure repack process is sensible
			for i in xrange(0,4):
				self.assertTrue(xdebug.dz5[i] <= xdebug.dz4[i])



if __name__ == "__main__":
	print modele_exe.Lisnowbase_Mod.Lisnowparams()


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
#         ! @var snht = sensible heat flux @ beginning of timestep [W m-2]
#                      sensible heat = conduction and convection
#         ! @var prht = heat of preciptation [W m-2]
#         double precision :: srht, snht, prht
# 
# sensible heat from surface to atmosphere: 
# latent heat from surface to atmosphere:
# (sign is same as evaporation)
# 
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
