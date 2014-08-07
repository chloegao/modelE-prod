SHELL=/bin/sh

#FC=    ifort -assume byterecl -convert big_endian
#LD=    ifort -assume byterecl -convert big_endian

NC=/usr/local/other/netcdf/3.6.2_intel-10.1.013

LEVEL_OPT = -g  -traceback -C
LEVEL_OPT = -O2 -traceback 

FC=     ifort $(LEVEL_OPT)   -convert big_endian -openmp
FC=     ifort $(LEVEL_OPT)   -convert big_endian

LD=     ifort $(LEVEL_OPT)   -convert big_endian -I${NC}/include -openmp
LD=     ifort $(LEVEL_OPT)   -convert big_endian -I${NC}/include

LDFLAGS =  
CPPFLAGS =
LIBS = -L${NC}/lib -lnetcdf

avg_obj=const_proc.o hycom_dimen.o hycom_arrays.o coagflx.o avg.o \
	gethyb.o extrct.o prtmsk.o bigrid.o gtdpth.o meshsz.o restep.o findmx.o \
	psmoo.o thruflow.o zebra.o TimeConstants.o
avg_src = $(avg_obj:.o=.f)

ov2_obj=const_proc.o hycom_dimen.o hycom_arrays.o pltov_mon.o gethyb.o \
	extrct.o prtmsk.o bigrid.o gtdpth.o meshsz.o restep.o \
	findmx.o thruflow.o zebra.o

latlonz_obj=const_proc.o hycom_dimen.o hycom_arrays.o cpler.o \
	latlonz.o bigrid.o gethyb.o extrct.o gtdpth.o restep.o \
	grdtrns.o findmx.o zebra.o

latlonz_src = $(latlonz_obj:.o=.f)

get2d_obj=const_proc.o hycom_dimen.o get2d.o 
get2d_src = $(get2d_obj:.o=.f)

avg:    $(avg_obj)
	$(LD) -o avg $(avg_obj) $(LIBS) $(LDFLAGS)

pltov_mon:  $(ov2_obj)
	$(LD) -o pltov_mon $(ov2_obj) $(LIBS) $(LDFLAGS)

latlonz:    $(latlonz_obj)
	$(LD) -o latlonz $(latlonz_obj) $(LIBS) $(LDFLAGS)

get2d:  $(get2d_obj)
	$(LD) -o get2d $(get2d_obj) $(LIBS) $(LDFLAGS)

avg.o: avg.f const_proc.f hycom_dimen.f TimeConstants.o
	$(LD) $(LDFLAGS) -c avg.f
bigrid.o: bigrid.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c bigrid.f
coagflx.o: coagflx.f const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c coagflx.f
const_proc.o: const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c const_proc.f
cpler.o: cpler.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c cpler.f
extrct.o: extrct.f const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c extrct.f
get2d.o: get2d.f const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c get2d.f
gethyb.o: gethyb.f const_proc.f hycom_dimen.f hycom_arrays.f
	$(LD) $(LDFLAGS) -c gethyb.f
gtdpth.o: gtdpth.f const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c gtdpth.f
grdtrns.o: grdtrns.f
	$(LD) $(LDFLAGS) -c grdtrns.f
hycom_dimen.o: hycom_dimen.f
	$(LD) $(LDFLAGS) -c hycom_dimen.f
hycom_arrays.o: hycom_arrays.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c hycom_arrays.f
latlonz.o: latlonz.f const_proc.f hycom_dimen.f hycom_arrays.f
	$(LD) $(LDFLAGS) -c latlonz.f
psmoo.o: psmoo.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c psmoo.f
restep.o: restep.f const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c restep.f
TimeConstants.o: ../../model/shared/TimeConstants.F90
	$(LD) $(LDFLAGS) -c ../../model/shared/TimeConstants.F90


tar: makefile input readme $(latlonz_src) $(get2d_src) $(avg_src)
	tar -cvf proc.tar makefile input readme $(latlonz_src) $(get2d_src) $(avg_src)

clean: 
	\rm *.o *.mod
