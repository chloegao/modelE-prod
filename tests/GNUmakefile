.PHONY: tests all createTestSuites clean
.DEFAULT_GOAL = tests

MODELERC ?= $(HOMEDIR)/.modelErc
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
TESTDIRS = shared profiler MPI_Support tracers
else
TESTDIRS := $(foreach ITEM, $(COMPONENTS), $(wildcard $(ITEM)))
endif

ifeq ($(USEMPI),)
  TESTDIRS := $(filter-out MPI_Support,$(TESTDIRS))
endif

TESTLIBS := $(foreach ITEM, $(TESTDIRS), $(ITEM)/lib$(ITEM).a)
# The inclue paths are different when compiling tests and executable
TESTINCS := $(foreach ITEM, $(TESTDIRS), -I../$(SOURCE_DIR)/$(ITEM))
EXEINCS := $(foreach ITEM, $(TESTDIRS), -I$(ITEM))

PFLIB = -L$(PFUNIT)/lib -lpfunit
IOLIB += -L$(NETCDFHOME)/lib -lnetcdf 

FFLAGS += -I$(PFUNIT)/mod -I$(PFUNIT)/include -I. $(TESTINCS)
EXEFLAGS += -I$(PFUNIT)/mod -I$(PFUNIT)/include $(EXEINCS)
ifeq ($(COMPILER),nag)
  FFLAGS += -kind=byte
endif

ifeq ($(COMPILER),intel)
  FFLAGS += -assume realloc_lhs
  FFLAGS += -O0 -g -traceback -warn nousage -convert big_endian -assume protect_parens -ftz -fpp -DPROTECT_PARENS

endif


ifeq ($(USEMPI),YES)
  FFLAGS += -DUSE_MPI
  IOLIB += -L$(PNETCDFHOME)/lib -lpnetcdf
endif

# Targets:

EXE = tests.x

$(EXE): createTestSuites testSuites.inc SUT
	$(FC) -I. -o $@ $(EXEFLAGS) $(FPPFLAGS) $(PFUNIT)/include/driver.F90 $(TESTLIBS) $(LIBS) $(PFLIB) $(IOLIB)

all: $(EXE)
ifeq ($(MPI),YES)
	mpirun -np 12 ./$(EXE)
else
	./$(EXE)
endif

SUT:
	@for dir in $(TESTDIRS); do \
	make -C $$dir all; done

tests: all

createTestSuites:
	rm -f testSuites.inc
	@for f in ${TESTDIRS}; do \
	echo \#include \"$${f}/testSuites.inc\" >> testSuites.inc; \
	done;

clean:
	@for dir in $(TESTDIRS); do \
	make -C $$dir clean; done
	$(RM) tests.x

export FC
export MPIF90
export FFLAGS
export FPPLAGS
export USEMPI
