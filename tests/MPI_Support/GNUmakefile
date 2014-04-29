.PHONY: clean all

SRCS = $(wildcard *.pf)
OBJS = $(SRCS:.pf=.o)

LIB = libMPI_Support.a

$(LIB): $(OBJS)
	ar -r $@ $?
	-@ranlib $@

all: $(OBJS) $(LIB)

testSuites.inc: $(SRCS)

%.F90: %.pf
	python $(PFUNIT)/bin/pFUnitParser.py $<  $@

%.o: %.F90
	$(FC) -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a

