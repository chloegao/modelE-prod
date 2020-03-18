#!/bin/ksh
#
# Example script to run oneYearTracersToNINT.ksh in a loop over years. That
# script has usage notes you can see by executing it witout arguments.
# Basically, it's a script to create NINT input from tracer code output.
#
# The location of NCO operators used to be defined here and passed to the
# script. Now it is defined within the script, which also loads modules.
#
######## USER SETS DATA ##################################################
# the run name of the tracers simulation to be used to create NINT input:
run=E212TomaSSP126aF40oQ40
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
# To skip doing ozone, aerosols, or BCalbedo, set these to > 0.
# For example, useful if you want aerosols and ozone to measure from a
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
# let's at least show a loop over 2 years to make the example less trivial.
# In the case of yearly (non-climatological) files, the 'representative
# year' could be the same as the model year, as it is here.
# You could check out the version of this example script from before
# 2020.03.18 to see an example operating on climatology files.

y1=2084
y2=2085
y=$y1
while [[ ${y} -le $y2 ]] ; do  # Loop over 2 years.
  yearLabel=${y}               # simply year here; climatologies case might look more like '2080-2089'
  representativeYear=${y}      # for climatologies, might look more like '2084' (for 2080-2089 example)
  # Pass arguments and execute the script that does one year's worth of data:
  ./oneYearTracersToNINT.ksh $run $yearLabel $representativeYear $inputDirectory $outputDirectory $referenceYear $gravity $skip_ozone $skip_aerosols $skip_BCalbedo
  let y+=1
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

