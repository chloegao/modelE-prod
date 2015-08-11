"""
  This script executes modelE serial and MPI runs and generates results
  used to verify various reproducibility measures.
  The script can be executed from the decks directory using default
  options and without arguments as follows:
 
      python  ../exec/testing/regression.py
 
  In that case the script will run the nonProduction_E_AR5_C12 rundeck 
  using the gfortran compiler and in serial and mpi modes. 
  Alternatively one can use the default options with one argument:
 
      python  ../exec/testing/regression.py <runsource>
 
  and run the <runsource> rundeck.
  Finally one can run a set of rundecks by specifying configuration
  files for each runsource, i.e. from the decks subdirectory issue the command:
 
      python  ../exec/testing/regression.py  <runsource1> [<runsource2> ...]
 
  This requires a configuration file name <runsource>.cfg for each runsource.
  When called from a higher level driver (mainDriver.py), this script's results
  are used to generate a report of regression testing reproducibility checks.
  
"""

import sys
import os
import shlex
import subprocess
import logging
import ConfigParser
from regRuns import *

def compare(rundeck, run):
    if run.verification != 'compileOnly':
        if run.mode == 'serial':
            compareBase(run, '1hr')
            if run.verification != 'run1hr':
                compareBase(run, run.endTime)
                # And compare SERIAL checkpoint-restart
                compareRestart(run)
        else:
            for npes in rundeck.npList:
                # Compare runs with baseline
                compareBase(run, '1hr', npes=npes)
                if run.verification != 'run1hr':
                    compareBase(run, run.endTime, npes=npes)
                    compareRestart(run, npes=npes)
                # Compare 1hr run against serial
                compareNPE(run, run.endTime, npes)

def writeDiff(run, fileH):
    fileH.write('%20s' % (run.results[0]))
    fileH.write('%10s' % (run.results[1]))
    fileH.write('%8s'  % (run.results[2]))
    fileH.write('%4s'  % '    ')
    for s in run.results[3:]:
        fileH.write('{: ^5}'.format(s))
        fileH.write('%3s'  % '   ')
    fileH.write('\n')

      
"""
  MAIN PROGRAM
"""
if __name__ == '__main__':

    setRunUtils()

    # System call return code
    OK = 0

    # List of runSources to verify specified on command line
    runSources = []
    if len(sys.argv) > 0:
        for arg in sys.argv[1:]:
            runSources.append(arg)
    else: # if none specified, use nonProduction_E_AR5_C12
        runSources.append('nonProduction_E_AR5_C12')

    # Loop over each run source in list
    for source in runSources:

        # Create rundeck object with default or config properties
        rundeck = newRundeck(source)
        
        # List of rundeck run configurations for each mode
        runs = []
        for mode in rundeck.modeList:
            runs.append(newRun(rundeck, mode))

        # Setup a logging stream
        setupLogging(rundeck)
        logger = logging.getLogger('MAIN    ')
        # All the work is done from the modelE decks directory
        os.chdir(rundeck.decksDir)

        diffFile = rundeck.resultsDir + '/' + rundeck.name + '.diff'
        fileH = open(diffFile, 'a')

        logger.info('Verifying ' + rundeck.name + ': ' + rundeck.verification)

        for run in runs:
            serBuildResult = OK
            mpiBuildResult = OK
            if run.mode == 'serial':
                serBuildResult = build(run)
                if rundeck.verification != 'compileOnly':
                    if serBuildResult == OK:
                        rc = run1hr(run)
                        if rc != 0:
                            continue
                        if rundeck.verification != 'run1hr':
                            rc = runRestart(run, endtime=run.endtime)
                            if rc != 0:
                                continue
            else:
                mpiBuildResult = build(run)
                if run.verification != 'compileOnly':
                    if mpiBuildResult == OK:
                        for npes in rundeck.npList:
                            rc = run1hr(run, npes=npes)
                            if rc != 0:
                                continue
                            if run.verification != 'run1hr':
                                rc = runRestart(run, npes=npes, endtime=run.endtime)
                                if rc != 0:
                                    continue

                elif rundeck.verification == 'customRun':
                    if mpiBuildResult == OK:
                        for npes in run.npList:
                            rc = runLong(run, npes=npes)
                            if rc != 0:
                                continue

            logger.info(rundeck.name + ' ' + run.mode + ' runs complete.')
            # If any build failed, then go on to the next run
            if serBuildResult != OK or mpiBuildResult != OK:
                continue

            if rundeck.standalone == 'yes':
                compare(rundeck, run)
                writeDiff(run, fileH)

        fileH.close()

        logger.info(rundeck.name + ' verification complete.')


    logger.info('Regression testing is done.')            
                    
                    
