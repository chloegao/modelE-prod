
ifneq ($(MODELE_CPATH),)
ifneq ($(CPATH),)
  CPATH_HACK = CPATH=$(strip $(MODELE_CPATH)):$(CPATH)
else
  CPATH_HACK = CPATH=$(strip $(MODELE_CPATH))
endif
endif

ifneq ($(MODELE_LIBRARY_PATH),)
ifneq ($(LIBRARY_PATH),)
  LIBRARY_PATH_HACK = LIBRARY_PATH=$(strip $(MODELE_LIBRARY_PATH)):$(LIBRARY_PATH)
else
  LIBRARY_PATH_HACK = LIBRARY_PATH=$(strip $(MODELE_LIBRARY_PATH))
endif
endif

# LIBRARY_PATH is specified for F90 because it is used for linking
F90 = $(LIBRARY_PATH_HACK) gfortran
CC = $(CPATH_HACK) gcc
ifneq ($(CPATH_HACK),)
  CPP := $(CPATH_HACK) $(CPP)
endif
GFORTRAN_RELEASE := $(shell gfortran --version | perl -e '<>=~/ (\d+)/; print "$$1";')
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS += -DCOMPILER_G95
FFLAGS = -g -cpp -fconvert=big-endian -O2 -fno-range-check
F90FLAGS = -g -cpp -fconvert=big-endian -O2 -fno-range-check -ffree-line-length-none
LFLAGS =

F90_VERSION = $(shell $(F90) --version | head -1)

# option to treat default real as real*8
R8 = -fdefault-real-8 -fdefault-double-8
I4 = 
#I4 = -fdefault-integer-4
EXTENDED_SOURCE = -ffixed-line-length-132
CRAYPTRFLAGS = -fcray-pointer

#
# Set the following to ensure that the beginning/end of records
# in sequential-access unformatted files have 4-byte markers
# (this does impose a 2 GB limit on record sizes, however)
#
#FFLAGS += -frecord-marker=4
#F90FLAGS += -frecord-marker=4

# check if ABI was specified explicitly
M64 =
ifneq ($(ABI),)
# leave it unset for arm64 Linux
ifneq ($(IS_AARCH64),YES)
  M64 = -m$(ABI)
endif
endif

# machine-specific options
ifeq ($(MACHINE),IRIX64)
  M64 = -mabi=64 
endif

FFLAGS += $(M64)
F90FLAGS += $(M64)
LFLAGS += $(M64)


# flags needed for particular releases

FFLAGS_RELEASE =
# Tested as version ranges rather than fixed lists, so that releases newer
# than the ones known when this was written still get the flags they need.
GFORTRAN_GE_8  := $(shell test "$(GFORTRAN_RELEASE)" -ge 8  2>/dev/null && echo YES)
GFORTRAN_GE_10 := $(shell test "$(GFORTRAN_RELEASE)" -ge 10 2>/dev/null && echo YES)
ifeq ($(GFORTRAN_GE_10),YES)
FFLAGS_RELEASE += -fallow-argument-mismatch
endif
ifeq ($(GFORTRAN_GE_8),YES)
FFLAGS_RELEASE += -fwrapv
endif

# the following flag is not strictly required, but ot makes the code 
# more reproducible, in a sense that the results are the same for -O2 and -O0
# (for all files but CLOUDS2.F90)
ifeq ($(GFORTRAN_GE_10),YES)
FFLAGS_RELEASE += -fno-expensive-optimizations
endif

FFLAGS += $(FFLAGS_RELEASE)
F90FLAGS += $(FFLAGS_RELEASE)

# uncomment next two lines for extensive debugging
# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -fbounds-check -ffpe-trap=invalid,zero,overflow -fbacktrace
F90FLAGS += -fbounds-check -ffpe-trap=invalid,zero,overflow -fbacktrace
FFLAGS += -finit-real=snan
F90FLAGS += -finit-real=snan
#LFLAGS += -lefence
endif
