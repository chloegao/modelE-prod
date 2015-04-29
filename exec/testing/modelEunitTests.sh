#!/bin/bash

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
  local deck=E4TcadiF40

  MAKELOG=make.log.${compiler}

# CREATE JOB SCRIPT

  cat << EOF > $jobScript
#!/bin/bash
#SBATCH --job-name=unitTest
#SBATCH --ntasks=16
#sbatch --constraint=hasw
#SBATCH --time=0:10:00
#SBATCH --account=s1001

# set up the modeling environment
. /usr/share/modules/init/bash
module purge
EOF

  if [ "$compiler" == "intel" ]; then

    cat << EOF >> $jobScript
module load comp/intel-14.0.3.174 mpi/impi-4.1.3.048 other/git-2.3.1
EOF

  else

    cat << EOF >> $jobScript
module load other/comp/gcc-4.9.1 other/mpi/openmpi/1.8.2-gcc-4.9.1 other/git-2.3.1
EOF

  fi

  cat << EOF >> $jobScript

export MODELERC=$REGSCRATCH/${compiler}/modelErc.${compiler}

cd $REGSCRATCH
rm -rf ${deck}.${compiler}

git clone /discover/nobackup/modele/clones/master ${deck}.${compiler} > /dev/null 2>&1

cd $REGSCRATCH/${deck}.${compiler}/decks
make rundeck RUN=$deck RUNSRC=$deck >> $MAKELOG 2>&1
# SERIAL RUN, MPI=NO
export PFUNIT=/discover/nobackup/ccruz/Baselibs/pFUnit/${compiler}"-serial"
make -j gcm RUN=$deck EXTRA_FFLAGS="-O0 -g" MPI=NO  >> $MAKELOG 2>&1
make tests RUN=$deck MPI=NO > $testLog.NO 2>&1
# MPI RUN, MPI=YES
make --quiet clean
export PFUNIT=/discover/nobackup/ccruz/Baselibs/pFUnit/${compiler}"-mpi"
make -j gcm RUN=$deck EXTRA_FFLAGS="-O0 -g" MPI=YES  >> $MAKELOG 2>&1
make tests RUN=$deck MPI=YES > $testLog.YES 2>&1

EOF
  chmod +x $jobScript

# SUBMIT JOB SCRIPT

  jobID=`sbatch $jobScript | awk '{print $4}'`
  if [ -z "$jobID" ]; then
    echo "There was a queue submission problem" >> $toEmail
    echo ""  >> $toEmail
    return
  fi

  watchJob $jobID jobRan

  if [ $jobRan -eq 0 ]; then
    echo " ### jobID=$jobID wait time (3600 secs) expired." >> $toEmail
    return
  fi  

  # Parse log files to generate results for eMail
  mpiMode=(YES NO)
  for mpi in "${mpiMode[@]}"; do
    echo " - MPI=$mpi"
    echo ""  >> $toEmail
    echo "RESULTS [$compiler MPI=$mpi]:" >> $toEmail
    echo ""  >> $toEmail
    parseLog "$testLog.$mpi"
  done

}

# -------------------------------------------------------------------
parseLog()
# -------------------------------------------------------------------
{
  local testLog=$1

# PARSE FOR SUCCESS

  echo "Parsing $testLog..."
  local lineNo=0
  # Find OK string
  local a=`grep -nw OK $testLog | head -1`
  lineNo=${a%%:*}
  # tests ran and all was OK
  if [ ! -z $lineNo ]; then
    msg=$(head -$(( lineNo+1 )) $testLog | tail -1)
    echo " $OK : " $msg  >> $toEmail
    return
  fi

# PARSE FOR FAILURES / BUILD / RUNTIME ERRORS

  lineNo=0
  local failures='FAILURES!!!'
  a=`grep -n $failures $testLog | head -1`
  lineNo=${a%%:*}
  # tests ran and there were failures
  if [ ! -z $lineNo ]; then
    msg=$(head -$(( lineNo+1 )) $testLog | tail -1)
    echo " $failures : " $msg  >> $toEmail
  # tests did not run
  else
    cp $testLog $FAILLOG
    if [ "$mpi" == "NO" ]; then
      execLine=`cat $testLog | grep -in './tests.x' | awk -F: '{print $1}'`
    else
      execLine=`cat $testLog | grep -in 'mpirun -np' | awk -F: '{print $1}'`
    fi
    if [ "$execLine" == "" ]; then
      echo " ### COMPILATION ERROR." >> $toEmail
    else
      msg=`grep SIGSEGV $testLog | grep SIGSEGV`
      if [ "$msg" == "" ]; then
         echo " ### UNEXPECTED RUNTIME ERROR." >> $toEmail
      else
         echo " ### RUNTIME ERROR : $msg, segmentation fault occurred" >> $toEmail
      fi
    fi
    echo ""  >> $toEmail
  fi
}

# ---------------------
# MAIN
# ---------------------

REGSCRATCH=/discover/nobackup/modele/regression_scratch/master
ROOT=`pwd`
toEmail="master.unit"
rm -f $toEmail slurm*out *.YES *.NO
compilers=(intel gfortran)
for compiler in "${compilers[@]}"; do 
  echo " -- COMPILER=$compiler"
  job=modelE.${compiler}.j
  log=${ROOT}"/"${compiler}".log"
  submitJob "$compiler" "$job" "$log"
  rm -f $job
done

MAILTO="giss-modele-regression@lists.nasa.gov"
/usr/bin/mail -s "modelE_RT (unit tests)" $MAILTO < $toEmail

exit 0
