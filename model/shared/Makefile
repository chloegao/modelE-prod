# makefile for the component shared
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF
F90SRCS  = AbstractAttribute.F90 \
	ArrayBundle_mod.F90 \
	AttributeDictionary.F90 \
	AttributeHashMap.F90 \
	AttributeReference.F90 \
	Attributes.F90 \
	BaseTime.F90 \
	Calendar.F90 \
	Constants_mod.F90 \
	CubicEquation_mod.F90 \
	Dictionary_mod.F90 \
	FileManager_mod.F90 \
	GaussianQuadrature.F90 \
	GenericType_mod.F90 \
	Geometry_mod.F90 \
	GetTime_mod.F90 \
	JulianCalendar.F90 \
	PlanetCalendar.F90 \
	KeyValuePair_mod.F90 \
	ModelClock.F90 \
	Month.F90 \
	Parser_mod.F90 \
	PlanetParams_mod.F90 \
	PolynomialInterpolator.F90 \
	Precision_mod.F90 \
	Random_mod.F90 \
	Rational.F90 \
	RootFinding_mod.F90 \
	RunTimeControls_mod.F90 \
	SpecialFunctions.F90 \
	stop_model.F90 \
	StringUtilities_mod.F90 \
	System.F90 \
	SystemTimers_mod.F90 \
	TimeAndCalendar.F90 \
	TimeConstants.F90 \
	Time.F90 \
	Utilities.F90 \
	SystemTools.F90

CSRCS = system_tools.c

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk

ifeq ($(USE_PFUNIT),YES)
$(OBJS): F90FLAGS+=-DUSE_PFUNIT -I$(PFUNIT)/mod
endif

