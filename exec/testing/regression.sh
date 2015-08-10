#!/bin/bash
#
# Runs regression.py and performs some basic regression testing
#
# Usage: regression.sh [RUNSRC1 RUNSRC2 ... RUNSRCN]
#
# Where [...] are optional arguments representing rundeck template
# names.
#
# Examples:
#
#   From decks directory, without arguments:
#
#   ../exec/testing/regression.sh
#   
#   will test nonProduction_E_AR5_C12 and will write all results
#   in decks directory.
#
#   ../exec/testing/regression.sh E4F40 E4TcadiF40 Earobio_g6c
#
#   will test E4F40, E4TcadiF40 and Earobio_g6c and all runs will 
#   be done in the decks directory.
#
# Caveats:
#
# 1) User must preload working env modules and set MODELERC
#    - Builds with compiler specified in MODELERC
#    - builds with COMPILE_WITH_TRAPS=YES
# 2) Runs interactively
#    - OK for small rundecks
# 3) No separate scratch space
#    - Else builds and runs proceed in the decks directory
#       - Beware of quotas
#    - Builds and runs proceed sequentially
# 4) No baseline testing is done - just internal consistency
#    - serial vs 4 pes
#    - checkpoint/restart vs continuous
# 5) Needs python version 2.7.x
#

root=`pwd`
scripts=$root/../exec/testing

if [ -z "$MODELERC" ]; then
   echo "Please set MODELERC."
   exit 1
else
   compiler=`grep COMPILER $MODELERC | awk -F= '{print $2}'`
fi

cnt=0
if [ "$#" -gt 1 ]; then
   userArgs=( "$@" )
   rundecks=()
   for arg in "${userArgs[@]}"; do
      rundecks=( "${rundecks[@]}" "$arg" )
   done

   node=`uname -n`
   # We need the right python version on DISCOVER
   if [[ "$node" =~ discover || "$node" =~ dali || "$node" =~ borg ]]; then
      export PATH=/usr/local/other/SSSO_Ana-PyD/2.1.0/bin:$PATH
   fi

   repo=${root%/*}
   for run in "${rundecks[@]}"; do
      cp  $scripts/template.cfg $run.cfg
      sed -i "s|COMPILER|${compiler}|g" $run.cfg
      sed -i "s|REPO|${repo}|g" $run.cfg
      sed -i "s/RUNDECK/${run}/g" $run.cfg
   done
else
   useArgs=()
fi

python $scripts/regression.py ${rundecks[@]}
