.PHONY: tests all createTestSuites clean

COMPILER ?= gfortran
EXTERNAL_CPP = NO
CONFIG_DIR = ../config
SOURCE_DIR = ../model
DECKS_DIR = ../decks
USE_MPI = NO
ifeq ($(MPI),YES)
  USE_MPI = YES
endif

INCS := $(foreach ITEM, $(COMPONENTS), $(SOURCE_DIR)/$(ITEM))
COMPLIBS = $(shell perl -e 'print reverse <>;' < $(SOURCE_DIR)/.liblist)
LIBS := $(foreach ITEM, $(COMPLIBS), $(SOURCE_DIR)/$(ITEM))

# Determine modelE components to be tested based on availability

ifeq ($(COMPONENTS),)
TESTDIRS = shared profiler MPI_Support tracers
else
TESTDIRS := $(foreach ITEM, $(COMPONENTS), $(wildcard $(ITEM)))
endif
ifeq ($(USE_MPI),NO)
  TESTDIRS := $(filter-out MPI_Support,$(TESTDIRS))
endif

TESTLIBS := $(foreach ITEM, $(TESTDIRS), $(ITEM)/lib$(ITEM).a)
# The inclue paths are different when compiling tests and executable
TESTINCS := $(foreach ITEM, $(TESTDIRS), -I../$(SOURCE_DIR)/$(ITEM))
EXEINCS := $(foreach ITEM, $(TESTDIRS), -I$(ITEM))

PFLIB = -L$(PFUNIT)/lib -lpfunit
IOLIB += -L$(NETCDFHOME)/lib -lnetcdf 

# Defaults

ifeq ($(COMPILER), intel)
  PFUNITMK := Intel.mk
endif
ifeq ($(COMPILER), gfortran)
  F90 = gfortran
  F90FLAGS += -g -O0 -fbacktrace
endif
ifeq ($(COMPILER), nag)
  PFUNITMK := NAG.mk
  F90FLAGS += -fpp
endif

# Sets F90 and additional F90FLAGS
F90FLAGS += -I$(PFUNIT)/mod -I$(PFUNIT)/include -I. $(TESTINCS)
EXEFLAGS += -I$(PFUNIT)/mod -I$(PFUNIT)/include $(EXEINCS)

ifeq ($(USE_MPI),YES)
  F90 := mpif90
  EXEFLAGS += -DUSE_MPI -DUSE_PFUNIT
  IOLIB += -L$(PNETCDFHOME)/lib -lpnetcdf
endif

export F90FLAGS
export F90
export USE_MPI

# Targets:

EXE = tests.x

#$(EXE): createTestSuites testSuites.inc tests
$(EXE): testSuites.inc tests
	$(F90) -I. -o $@ $(EXEFLAGS) $(PFUNIT)/include/driver.F90 $(TESTLIBS) $(LIBS) $(PFLIB) $(IOLIB)

tests:
	@for dir in $(TESTDIRS); do \
	make -C $$dir all; done

createTestSuites:
	rm -f testSuites.inc
	@for f in ${TESTDIRS}; do \
	echo \#include \"$${f}/testSuites.inc\" >> testSuites.inc; \
	done;

clean:
	@for dir in $(TESTDIRS); do \
	make -C $$dir clean; done
	$(RM) tests.x

all: $(EXE)
ifeq ($(MPI),YES)
	mpirun -np 12 ./$(EXE)
else
	./$(EXE)
endif

