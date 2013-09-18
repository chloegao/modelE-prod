#!/bin/bash

# Main modelE regression tests script. This is an interface to other scripts:
# regressionScripts.pl, modelEunitTests.sh and diffreport.j

# Set up working variables
needDefault(){
   if [ -z "$1" ]; then
     echo "true"
     exit 0
   else
     echo "false"
     exit 1
   fi
}
setupRunEnvVariables()
{
   if  [ "$(needDefault $MODELROOT)" == "true" ] ; then
     export MODELROOT=$NOBACKUP/devel/modelE.clones/master
   fi

   if  [ "$(needDefault $REGWORK)" == "true" ] ; then
     export REGWORK=$NOBACKUP
   fi

   export REGSCRATCH=$REGWORK/regression_scratch
   # create if necessary
   if [ ! -d "$REGSCRATCH" ]; then
     mkdir $REGSCRATCH
   fi

   export REGRESULTS=$REGWORK/regression_results
   # create if necessary
   if [ ! -d "$REGRESULTS" ]; then
     mkdir $REGRESULTS
   fi

   if  [ "$(needDefault $MODELEBASELINE)" == "true" ] ; then
     export MODELEBASELINE=/discover/nobackup/modele/modelE_baseline
   fi

   if  [ "$(needDefault $GCMSEARCHPATH)" == "true" ] ; then
     export GCMSEARCHPATH=/discover/nobackup/projects/giss/prod_input_files
   fi
	
}

watchJob()
{
# Monitor job
# Input arguments: $1=job id
   local jobID=$1

   maxWait=7200
   seconds=0
   done=0
   while [ $seconds -lt $maxWait ];
     do
     qStatus=`qstat | grep $jobID | awk '{print $5}'`
     if [ -z "$qStatus" ]; then
       done=1
       break
     fi
     sleep 30
     let seconds=$seconds+30
   done
}

# MAIN PROGRAM

   # User may just want the help page
   if [ "$1" = "--help" -o "$1" = "-h" ]
   then
     ./regTestsHelp.sh $0 $1; exit 0 ;
   fi

   if [ $# -ne 1 ]; then
     echo "Usage: `basename $0` {cfgFile}"
     exit 65
   fi
   export CONFIG=$1

   setupRunEnvVariables

   unset PFUNIT

   # Run suite of modelE tests
   if [ "$RUN_TESTS" == "YES" ]; then
     echo "Execute regressionTests.pl with $CONFIG configuration..."
     /usr/bin/perl regressionTests.pl $CONFIG.cfg 1>> $CONFIG.out 2>&1
   else
     echo "Skipped regression tests (RUN_TESTS=$RUN_TESTS)"
     echo "Skipped regression tests (RUN_TESTS=$RUN_TESTS)" >> $CONFIG.out
   fi

   # Optionally run unit tests
   if [ "$RUN_UNIT_TESTS" == "YES" ]; then
     echo "Run unit tests..."
     ./modelEunitTests.sh
     wait
   else
     echo "Skipped unit tests (RUN_UNIT_TESTS=$RUN_UNIT_TESTS)"
     echo "Skipped unit tests (RUN_UNIT_TESTS=$RUN_UNIT_TESTS)" >> $CONFIG.out
   fi

   # Optionally create a DIFF report
   if [ "$CREATE_DIFF" == "YES" ]; then
     echo "Run diffreport..."
     # Submit job to run diffreport.x
     if [ -z "$MOCKMODELE" ]; then
       jobID=`qsub $MODELROOT/exec/testing/diffreport.j`
       jobID=`echo $jobID | sed 's/.[a-z]*$//g'`
       watchJob $jobID
     else
       $MODELROOT/exec/testing/diffreport.j > $CONFIG.diff
     fi
   else
     echo "Skipped diffreport (CREATEDIFF=$CREATE_DIFF)"
     echo "Skipped diffreport (CREATEDIFF=$CREATE_DIFF)" > $CONFIG.diff
   fi

   echo "DONE."
   exit 0
