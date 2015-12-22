"""
  This script executes modelE serial and MPI run combinations and generates 
  results used to verify various reproducibility measures.

  The script can be executed from the decks directory as follows:

     ../exec/testing/regression.sh [RUNSRC1 RUNSRC2 ... RUNSRCN -r]

  See the documentation in regression.sh 

  If executed standalone as follows:
 
      python  ../exec/testing/regression.py  <runsource1> [<runsource2> ...]
 
  then this requires a configuration file name <runsource>.cfg for each runsource.
  When called from the higher level driver (reg), this script's results
  are used to generate a report of regression testing reproducibility checks.
  
"""

import sys
import os
import shlex
import subprocess
import logging
import ConfigParser
import regCompare
import regRuns as r

def compare(run):
    if run.mode == 'serial':
        regCompare.base(run, 1) # 1hr run
        if run.verification == 'restartRun':
            regCompare.base(run, run.endTime)
            regCompare.restart(run)
    else:
        for npes in run.getOpt('npes'):
            regCompare.base(run, 1, npes=npes)
            if run.verification == 'restartRun':
                regCompare.base(run, run.endTime, npes=npes)
                regCompare.restart(run, npes=npes)
            # Compare NPE vs serial
            regCompare.nPE(run, run.endTime, npes)


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

    # System call return code
    OK = 0

    # List of runSources to verify specified on command line
    runSources = []
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            runSources.append(arg)
    else: # if none specified, use nonProduction_E_AR5_C12
        runSources.append('nonProduction_E_AR5_C12')

    # Loop over each run source in list
    for source in runSources:

        # Create rundeck object
        rundeck = r.newRundeck(source)

        # Setup a logging stream for this rundeck
        rundeck.setLogging()
        logger = logging.getLogger('MAIN    ')
        
        # Create a list of rundeck run configurations for each mode
        runs = []
        for mode in rundeck.modes:
            runs.append(r.newRun(rundeck, mode))

        # All the work is done from the modelE decks directory
        os.chdir(rundeck.decksDir)

        diffFile = rundeck.resultsDir + '/' + rundeck.name + '.diff'
        fileH = open(diffFile, 'a')

        logger.info('Verifying ' + rundeck.name + ': ' + rundeck.verification)

        for run in runs:

            # initialize build results to OK
            serBuildResult = OK
            mpiBuildResult = OK

            if run.mode == 'serial':
                serBuildResult = r.build(run)
                if rundeck.verification == 'compileOnly':
                    continue
                if serBuildResult == OK:
                    # Always run 1hr
                    rc = r.run1hr(run)
                    if rundeck.verification == 'restartRun':
                        rc = r.runRestart(run, endTime=run.endTime)
            else: # MPI
                mpiBuildResult = r.build(run)
                if run.verification == 'compileOnly':
                    continue
                if mpiBuildResult == OK:
                    if run.verification == 'customRun':
                        if mpiBuildResult == OK:
                            for npes in run.npes:
                                rc = r.runLong(run, npes=npes)
                                if rc != OK:
                                    continue
                    else:    
                        for npes in run.npes:
                            # Always run 1hr
                            rc = r.run1hr(run, npes=npes)
                            if run.verification == 'restartRun':
                                rc = r.runRestart(run, npes=npes, endTime=run.endTime)

            logger.info(rundeck.name + ' ' + run.mode + ' runs complete.')

            # If any build failed, then go on to the next run
            if serBuildResult != OK or mpiBuildResult != OK:
                continue

            if rundeck.standalone == 'yes':
                if run.verification != 'customRun' or run.verification == 'compileOnly':
                    compare(run)
                    writeDiff(run, fileH)

        fileH.close()

        logger.info(rundeck.name + ' verification complete.')


    logger.info('Regression testing is done.')            
                    
                    
