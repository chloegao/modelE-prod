#!/usr/local/bin/bash

# Script to run modelE unit tests on Linux (DISCOVER)

# -------------------------------------------------------------------
# FUNCTIONS
# -------------------------------------------------------------------

# -------------------------------------------------------------------
watchJob()
# -------------------------------------------------------------------
{
# Monitor job
# Input arguments: $1=job id
  jobID=$1

  maxWait=3600
  seconds=0
  jobSuccess=0
  while [ $seconds -lt $maxWait ];
  do
    qStatus=`qstat | grep $jobID | awk '{print $5}'`
    if [ -z "$qStatus" ]; then
      jobSuccess=1
      break
    fi
    sleep 10
    let seconds=$seconds+10
  done
  eval $2="$jobSuccess"
}
# -------------------------------------------------------------------
submitJob()
# -------------------------------------------------------------------
{

  local compiler=$1
  local jobScript=$2
  local testLog=$3
  local mpi=$4

  #local deck=E4TcadiF40
  local deck=nonProduction_E4TcadC12

  MAKELOG=make.log.${compiler}
  FAILLOG=${testLog}.FAILED
  pfunitSuffix="-mpi"
  if [ "$mpi" == "NO" ]; then
    pfunitSuffix="-serial"
    FAILLOG=${testLog}${pfunitSuffix}.FAILED
  fi

# CREATE JOB SCRIPT

  cat << EOF > $jobScript
#!/bin/bash
#PBS -N modelEut
#PBS -l select=1:mpiprocs=12
#PBS -l walltime=0:05:00
#PBS -W group_list=s1001
#PBS -j oe

# set up the modeling environment
. /usr/share/modules/init/bash
module purge
EOF

  if [ "$compiler" == "intel" ]; then

    cat << EOF >> $jobScript
module load comp/intel-14.0.0.080 mpi/impi-3.2.2.006 other/git-1.7.3.4
EOF
   
  else

    cat << EOF >> $jobScript
module load other/comp/gcc-4.8.1 other/mpi/openmpi/1.7.2-gcc-4.8.1-shared other/git-1.7.3.4
EOF

  fi

  cat << EOF >> $jobScript

export PFUNIT=/discover/nobackup/ccruz/Baselibs/pFUnit/${compiler}${pfunitSuffix}
export MODELERC=$REGSCRATCH/${compiler}/modelErc.${compiler}

cd $REGSCRATCH
rm -rf ${deck}.${compiler}
git clone $MODELROOT ${deck}.${compiler} > /dev/null 2>&1

cd $REGSCRATCH/${deck}.${compiler}/decks
make rundeck RUN=$deck RUNSRC=$deck >> make.log.${compiler} 2>&1
EOF

  if [ "$compiler" == "intel" ]; then

    cat << EOF >> $jobScript
make -j gcm RUN=$deck EXTRA_FFLAGS="-O0 -g -traceback" MPI=$mpi >> $MAKELOG 2>&1
wait
EOF

  else

    cat << EOF >> $jobScript
make -j gcm RUN=$deck EXTRA_FFLAGS="-O0 -g -fbacktrace" MPI=$mpi  >> $MAKELOG 2>&1
wait
EOF

  fi

  cat << EOF >> $jobScript
make tests RUN=$deck MPI=$mpi >> $testLog 2>&1
wait
EOF
  chmod +x $jobScript

# SUBMIT JOB SCRIPT

  echo ""  >> $toEmail
  echo "RESULTS [$compiler MPI=$mpi]:" >> $toEmail
  echo ""  >> $toEmail
  jobID=`qsub $jobScript`
  # Not necessary under SLURM
  #jobID=`echo $jobID | sed 's/.[a-z]*$//g'`
  #echo 'jobID='$jobID
  if [ -z "$jobID" ]; then
    echo "There was a queue submission problem" >> $toEmail
    echo ""  >> $toEmail
    return
  fi

  watchJob $jobID jobRan

  if [ $jobRan -eq 0 ]; then
    cp $testLog $FAILLOG
    echo " ### jobID=$jobID wait time (3600 secs) expired." >> $toEmail
    echo " ### Check $FAILLOG" >> $toEmail
    return
  fi  

# PARSE FOR ERRORS

  local anyError=`grep ' FAIL' $testLog` 
  if [ "$anyError" != "" ]; then
    local errMsg=" ### Error detected during unit tests" >> $toEmail
    echo "SUMMARY:"  >> $toEmail
    echo "..."  >> $toEmail
    tail -10 $testLog >> $toEmail
    echo ""  >> $toEmail
  fi 
   
  if [ "$anyError" != "" ]; then
    totLines=`cat $testLog | wc | awk '{print $1}'`
    if [ "$mpi" == "NO" ]; then
      execLine=`cat $testLog | grep -in './tests.x' | awk -F: '{print $1}'`
    else
      execLine=`cat $testLog | grep -in 'mpirun -np' | awk -F: '{print $1}'`
    fi
    # prune the output a little more...
    blockLines=$((totLines-execLine))
    showLines=$((blockLines-1)) 
    tail -$blockLines $testLog > foo
    head -$showLines foo >> $toEmail
    rm foo
  else
    msg=$(tail -3 $testLog | grep "OK")
    if [ "$msg" == "" ]; then
      cp $testLog $FAILLOG
      echo " ### Tests failed to run." >> $toEmail
      echo " ### Check $FAILLOG" >> $toEmail
    else
      tail -3 $testLog | grep "OK" >> $toEmail
      tail -2 $testLog | grep "(" >> $toEmail
      rm -f $testLog
    fi
    echo ""  >> $toEmail
  fi 
  wait

# DELETE FILES
  rm -f $testLog $jobScript

}

# ---------------------
# MAIN
# ---------------------

ROOT=$MODELROOT/exec/testing/testsOutput/
cd $ROOT
toEmail="$CONFIG.unit"
rm -f $toEmail
compilers=(intel gfortran)
mpiMode=(YES NO)
for mpi in "${mpiMode[@]}"; do
  echo " - MPI=$mpi"
  for compiler in "${compilers[@]}"; do 
    echo " -- COMPILER=$compiler"
    job=modelE.${compiler}.j
    log=${ROOT}${compiler}".log"
    submitJob "$compiler" "$job" "$log" "$mpi"
  done 
done 

exit 0
