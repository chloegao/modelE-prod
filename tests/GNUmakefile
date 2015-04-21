.PHONY: tests all clean
.DEFAULT_GOAL = tests

MODELERC ?= $(HOMEDIR)/.modelErc
ifneq ($(MPI),YES)
  USE_MPI=NO
else
  USE_MPI=YES
endif
FFLAGS = $(CPPFLAGS) $(F90FLAGS)

ifeq ($(wildcard $(MODELERC)),)
   $(error $(MODELERC) not found.)
endif
include $(MODELERC)
ifeq ($(COMPILER),)
   $(error you need to specify COMPILER in your ~.modelErc)
endif
SOURCE_DIR = ../model

# Unfortunately COMPILER is defined differently in modelE (gfortran)
# than in PFUNIT (GNU). So, we rename here to inherit the correct PFUNIT
# configuration.
ifeq ($(COMPILER),gfortran)
   COMPILER=GNU
endif

COMPLIBS = $(shell perl -e 'print reverse <>;' < $(SOURCE_DIR)/.liblist)
LIBS := $(foreach ITEM, $(COMPLIBS), $(SOURCE_DIR)/$(ITEM))

# Determine modelE components to be tested based on availability
ifeq ($(COMPONENTS),)
  TESTDIRS = shared profiler
else
  TESTDIRS := $(foreach ITEM, $(COMPONENTS), $(wildcard $(ITEM)))
endif

# Exclude MPI_Support tests if not running MPI tests
ifeq ($(USE_MPI),NO)
  TESTDIRS := $(filter-out MPI_Support,$(TESTDIRS))
endif

TESTLIBS := $(foreach ITEM, $(TESTDIRS), $(ITEM)/lib$(ITEM).a)
# The include paths are different when compiling tests and executable
TESTINCS := $(foreach ITEM, $(TESTDIRS), -I../$(SOURCE_DIR)/$(ITEM))
FFLAGS += $(TESTINCS) -I$(PFUNIT)/mod
EXEFLAGS += $(FFLAGS) -I.

# Required libraries
PFLIB = -L$(PFUNIT)/lib -lpfunit
IOLIB += -L$(NETCDFHOME)/lib -lnetcdf 

ifeq ($(COMPILER),intel)
  FFLAGS += -DPROTECT_PARENS
endif

# For MPI, use compiler wrapper and link to parallel-netcdf library
ifeq ($(USE_MPI),YES)
  FC = mpif90
  IOLIB += -L$(PNETCDFHOME)/lib -lpnetcdf
  EXEFLAGS += -DUSE_MPI
endif

# Targets:

EXE = tests.x

MYLIBS = $(TESTLIBS) $(LIBS) $(PFLIB) $(IOLIB)
$(EXE): testSuites.inc TESTS
	$(FC) -o $@ $(EXEFLAGS) $(PFUNIT)/include/driver.F90 $(MYLIBS)

all: $(EXE)
ifeq ($(USE_MPI),YES)
	mpirun -np 12 ./$(EXE)
else
	./$(EXE)
endif

TESTS:
	@for dir in $(TESTDIRS); do \
	make -C $$dir all; done

tests: all

clean:
	@for dir in $(TESTDIRS); do \
	make -C $$dir clean; done
	$(RM) tests.x

export FC
export MPIF90
export FFLAGS
export USE_MPI
