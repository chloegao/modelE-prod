#!/bin/ksh
#
# Example script to run oneYearTracersToNINT.ksh in a loop over years. That
# script has usage notes you can see by executing it witout arguments.
# Basically, it's a script to create NINT input from tracer code output.
#
######## USER SETS DATA ##################################################
# the run name of the tracers simulation to be used to create NINT input:
run=E14TomaOCNf10_4av
# the year from which to measure time (from it's January):
referenceYear=1850
# input directory (where you have the monthly acc files, and optionally,
# scaled output; symbolic links will be made if this is not local):
inputDirectory=DATA
# output directory where NINT input will arrive; please make local (in .)
outputDirectory=NEW
# Susanne Bauer likes to know when people are creating NINT aerosol input
# so she can check on it. So, by default, she and Greg Faluvegi will be
# e-mailed when this looping script is run, indicating who is making input
# and where. To override this, set private to non-zero:
private=0
# Allow setting of gravitational constant for airmass calculations in case
# that's helpful for changing planets:
gravity=9.80665E0
# Location of the NCO (netCDF operators):
NCO=/usr/local/other/SLES11.1/nco/4.4.4/intel-12.1.0.233/bin
# To skip doing ozone, aerosols, or BCalbedo, set these to > 0.
# For example useful it you want aerosols and ozone to measure from a
# different reference year or you need to use special years for the
# BC albedo, etc... For ozone only, setting a negative value (-N) will
# skip N levels from the top of the tracer output (useful e.g. to avoid
# having to merge tracer 102L output with TOA ozone climatology for
# overlapping levels; set to -1).
skip_ozone=0
skip_aerosols=0
skip_BCalbedo=0
##########################################################################

# Once you set the above, the script could be run in one line, but
# for this example, let's make it a more useful one:
# The would take 1850-2009 transient output and make the climatological
# NINT input. Decadal averages (YYY0-YYY9) were already made with the acc
# files:

# (Note if these had been 9-year averages, a middle year could have been
# chosen as the represetative year. In this example, we follow how we think
# the SST/Sea ice files do it -- making the timestream representative year
# the YYY4...)

x=185                            # These lines are set to do 1850-1859
while [[ ${x} -le 200 ]] ; do    # through 2000-2009 climatologies
  yearLabel=${x}0-${x}9          # representing the years:
  representativeYear=${x}4       # (1854, 1864,... 2004)
  # Pass arguments and execute the script that does one year's worth of data:
  ./oneYearTracersToNINT.ksh $run $yearLabel $representativeYear $inputDirectory $outputDirectory $referenceYear $gravity $NCO $skip_ozone $skip_aerosols $skip_BCalbedo
  let x+=1
done

# Done. Send out a warning :-)
ozoneEmail=gregory.s.faluvegi@nasa.gov
aerosolsEmail=susanne.e.bauer@nasa.gov
if [[ $private -eq 0 ]] ; then
  if [[ $skip_ozone -eq 0 ]] ; then
    echo "${USER} is preparing NINT ozone input from tracers in ${PWD}." | mail -s "NINT_ozone_from_tracers" "$ozoneEmail"
  fi
  if [[ $skip_aerosols -eq 0 || $skip_BCalbedo -eq 0 ]] ; then
    echo "${USER} is preparing NINT aerosol input from tracers in ${PWD}." | mail -s "NINT_aerosols_from_tracers" "$aerosolsEmail"
  fi
fi

