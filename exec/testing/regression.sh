#!/bin/bash
#
# Runs regression.py and performs some basic regression testing
#
# Usage: regression.sh [RUNSRC1 RUNSRC2 ... RUNSRCN -r]
#
# Where [...] are optional arguments representing rundeck template
# names and one additional flag, -r, to perform restart regression.
#
# Examples:
#
#   From decks directory, without arguments:
#
#   1)   ../exec/testing/regression.sh
#   
#   will compile-only nonProduction_E_AR5_C12 and will write all results
#   in decks directory.
#
#   2)   ../exec/testing/regression.sh E4F40 E4TcadiF40 Earobio_g6c
#
#   will compile-only E4F40, E4TcadiF40 and Earobio_g6c and all runs will 
#   be done in the decks directory.
#
#   3)   ../exec/testing/regression.sh nonProduction_E_AR5_C12 -r
#   
#   will run restart-regression on nonProduction_E_AR5_C12
#
#   4)   ../exec/testing/regression.sh E4F40 E4TcadiF40 E4TctomasF40 -r
#
#   will run restart-regression on E4F40, E4TcadiF40 and E4TctomasF40.
#
# Errors, if any, will be printed on STDOUT.
#
# Caveats:
#
# 1) User must preload working env modules and set MODELERC
#    - Builds with compiler specified in MODELERC
#    - builds with COMPILE_WITH_TRAPS=YES
# 2) Runs interactively
#    - OK for small rundecks
# 3) No separate scratch space
#    -  Builds and runs proceed in the decks directory
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
verification=compileOnly
if [ "$#" -gt 0 ]; then
   userArgs=( "$@" )
   rundecks=()
   for arg in "${userArgs[@]}"; do
      if [ "$arg" == "-r" ]; then
	 verification=restartRun
      else
	 rundecks=( "${rundecks[@]}" "$arg" )
      fi
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
      sed -i "s/VERIFICATION/${verification}/g" $run.cfg
   done
else
   rundecks=("nonProduction_E_AR5_C12")
fi

python $scripts/regression.py ${rundecks[@]}
