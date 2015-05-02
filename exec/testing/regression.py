"""
  This script verifies that MPI and Serial builds produce identical
  results. The script can be executed from the decks directory using
  default options and without arguments as follows:
 
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
  In this case the script can be called from a higher level driver to execute
  a more complicated combination of experiments as is done with the
  nightly regression tests.
 
  ENV Options:
     * If the environment variable DEBUG is set, then the script will
       display all commands, but not actually execute them.
 
"""

import sys
import os
import shlex
import subprocess
import logging
import ConfigParser
from regFuncs import *
       
"""
  MAIN PROGRAM
"""
if __name__ == '__main__':

    setRunUtils()

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

        # Create rundeck object with default or config properties
        rundeck = newRundeck(source)
        
        # List of rundeck run configurations for each mode
        exps = []
        for mode in rundeck.modeList:
            exps.append(newRun(rundeck, mode))
            nmodes = len(exps)

        # Setup a logging stream
        setupLogging(rundeck)
        logger = logging.getLogger('MAIN    ')
        # All the work is done from the modelE decks directory
        os.chdir(rundeck.decksDir)
        # For each rundeck create a diffFile with verification results
        diffFile = rundeck.resultsDir + '/' + rundeck.name + '.diff'
        fileH = open(diffFile, 'w')

        logger.info('Testing ' + rundeck.name + ', testlevel: ' + rundeck.testLevel)

        for exp in exps:

            serBuildResult = OK
            mpiBuildResult = OK
            if exp.mode == 'serial':
                serBuildResult = build(exp)
                if rundeck.testLevel == 'full':
                    if serBuildResult == OK:
                        run1hr(exp)
                        if rundeck.testLevel != 'run1hr':
                            runRestart(exp, endtime=rundeck.endtime)
            else:
                mpiBuildResult = build(exp)
                if rundeck.testLevel == 'full':
                    if mpiBuildResult == OK:
                        for npes in rundeck.npList:
                            run1hr(exp, npes=npes)
                            if rundeck.testLevel != 'run1hr':
                                runRestart(exp, npes=npes, endtime=rundeck.endtime)
                elif rundeck.testLevel == 'long':
                    if mpiBuildResult == OK:
                        for npes in rundeck.npList:
                            runLong(exp, npes=npes)

            logger.info(rundeck.name + ' ' + exp.mode + ' runs complete.')
            if serBuildResult != OK or mpiBuildResult != OK:
                continue

            if rundeck.testLevel == 'full':
                if exp.mode == 'serial':
                    compareBase(exp, '1hr')
                    compareBase(exp, exp.endTime)
                # And compare SERIAL checkpoint-restart 
                    compareRestart(exp)
                else:
                    for npes in rundeck.npList:
                    # Compare runs with baseline
                        compareBase(exp, '1hr', npes=npes)
                        compareBase(exp, exp.endTime, npes=npes)
                        compareRestart(exp, npes=npes)
                        for npes in rundeck.npList:
                    # Compare 1hr run against serial
                            if nmodes > 1:
                                compareNPE(exps[0], exps[1], exp.endTime, npes)
                logger.info(rundeck.name + ' comparisons complete.')
                
        for exp in exps:
            fileH.write('%20s' % (exp.results[0]))
            fileH.write('%10s' % (exp.results[1]))
            fileH.write('%8s'  % (exp.results[2]))
            for s in exp.results[3:]:
                fileH.write(' '.center(3))
                fileH.write(s.center(3))
            fileH.write('\n')
        fileH.close()

    logger.info('Regression testing is done.')
            
                    
                    
