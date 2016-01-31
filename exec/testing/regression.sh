#!/bin/bash
#
# Runs regression.py and performs some basic regression testing
#
# Usage: regression.sh [RUNSRC1 RUNSRC2 ... RUNSRCN -c] clean
#
# Where [...] are optional arguments representing rundeck template
# names and one additional flag, -c, to perform compile-only verification.
#
# Examples:
#
#   From decks directory, without arguments:
#
#   1)   ../exec/testing/regression.sh
#   
#   will regression test nonProduction_E_AR5_C12 and will write all results
#   in decks directory.
#
#   2)   ../exec/testing/regression.sh E4F40 E4TcadiF40 Earobio_g6c -c
#
#   will compile-only E4F40, E4TcadiF40 and Earobio_g6c and all work will 
#   be done in the decks directory.
#
#   3)   ../exec/testing/regression.sh E4F40 E4TcadiF40 E4TctomasF40
#
#   will run restart-regression on E4F40, E4TcadiF40 and E4TctomasF40.
#
#   4) ../exec/testing/regression.sh clean
#
#   will remove ALL the temporary files and directories created by the
#   script.
#
# Errors, if any, will be printed on STDOUT.
#
# Caveats:
#
# 1) User must preload working env modules and set MODELERC
#    - Builds with compiler specified in MODELERC
#    - builds with COMPILE_WITH_TRAPS=YES
# 2) Runs interactively
#    - OK for small rundecks (so, example (3) is not recommended)
# 3) No separate scratch space
#    -  Builds and runs proceed in the decks directory
#       - Beware of quotas
#    - Builds and runs proceed sequentially
# 4) No baseline testing is done - just internal consistency
#    - serial vs 4 pes
#    - checkpoint/restart vs continuous
# 5) Needs python version 2.7.x
#
clean () {
   local dirs=()
   rm -rf *.mk *.R *.diff *.log *cfg* *_bin templ
   dirs=`find . -maxdepth 1  -type l -exec ls -d {} \;`
   for d in "${dirs[@]}"; do
      rm -rf $(readlink $d)
   done
   find . -type l -exec rm {} \;
   exit 1
}

root=`pwd`
scripts=$root/../exec/testing

if [ -z "$MODELERC" ]; then
   echo "Please set MODELERC."
   exit 1
else
   compiler=`grep COMPILER $MODELERC | awk -F= '{print $2}'`
fi

cnt=0
verification=restartRun
if [ "$#" -gt 0 ]; then
   if [ "$1" == "clean" ]; then
      clean
   fi
   userArgs=( "$@" )
   rundecks=()
   for arg in "${userArgs[@]}"; do
      if [ "$arg" == "-c" ]; then
	 verification=compileOnly
      else
	 rundecks=( "${rundecks[@]}" "$arg" )
      fi
   done
else
   rundecks=("nonProduction_E_AR5_C12")
fi

   node=`uname -n`
   # We need the right python version on DISCOVER
   if [[ "$node" =~ discover || "$node" =~ dali || "$node" =~ borg ]]; then
      export PATH=/usr/local/other/SSSO_Ana-PyD/2.1.0/bin:$PATH
   fi

   repo=${root%/*}
   for run in "${rundecks[@]}"; do
      cp  $scripts/template.cfg $run.cfg
      sed -i -e "s|COMPILER|${compiler}|g" $run.cfg
      sed -i -e "s|REPO|${repo}|g" $run.cfg
      sed -i -e "s|MODELERC|${MODELERC}|g" $run.cfg
      sed -i -e "s/RUNDECK/${run}/g" $run.cfg
      sed -i -e "s/VERIFICATION/${verification}/g" $run.cfg
   done

python $scripts/regression.py ${rundecks[@]}
