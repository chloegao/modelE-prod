
F90 = ifort
IFORT_RELEASE := $(shell ifort --version | perl -e \
  'while(<>){ if(/ifort.* (\d+\.\d+)/) { print "$$1"; } }')
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CMP_MOD = $(SCRIPTS_DIR)/compare_module_file.pl -compiler INTEL-ifort-9-0-on-LINUX
FFLAGS = -fpp -O2 -ftz         -convert big_endian 
F90FLAGS = -fpp -O2 -ftz        -convert big_endian -free 
LFLAGS = -O2 -ftz
CPPFLAGS += -DCOMPILER_Intel8 -DCONVERT_BIGENDIAN
F90_VERSION = $(shell $(F90) -v 2>&1)
ifeq ($(MP),YES)
FFLAGS += -openmp
F90FLAGS += -openmp
LFLAGS += -openmp
endif
R8 = -r8
EXTENDED_SOURCE = -extend_source

# flags needed for particular releases

# default flags for latest releases (work for 12.*, 13.*):
FFLAGS_RELEASE = -assume protect_parens -fp-model strict -warn nousage -assume realloc_lhs

# if some releases require different flags enter them here
ifeq ($(IFORT_RELEASE),11.1)
FFLAGS_RELEASE = -assume protect_parens -fp-model strict -fp-speculationoff
endif

ifeq ($(IFORT_RELEASE),11.0)
FFLAGS_RELEASE = -assume protect_parens -fp-model strict -fp-speculationoff
endif

ifeq ($(IFORT_RELEASE),10.1)
FFLAGS_RELEASE = -assume protect_parens -fp-model strict -diag-disable vec  -fp-speculationoff
endif

ifeq ($(IFORT_RELEASE),9.1)
FFLAGS_RELEASE =
endif

FFLAGS += $(FFLAGS_RELEASE)
F90FLAGS += $(FFLAGS_RELEASE)

ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -CB -fpe0 -check uninit -ftrapuv -traceback
LFLAGS += -CB -fpe0 -check uninit -ftrapuv -traceback
F90FLAGS += -CB -fpe0 -check uninit -ftrapuv -traceback
LFLAGSF += -CB -fpe0 -check uninit -ftrapuv -traceback
endif
