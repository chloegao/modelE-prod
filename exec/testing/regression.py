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
  a more compilcated combination of experiments as is done with the
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

"""
  This class defines a configuration object for each modelE runSource.
  The properties of this object are specified in the configuration options.
"""
class configuration:
    def __init__(self, runSrc, mode, comp):
        runShort = runSrc
        if 'nonProduction' in runSrc:
            runShort = runSrc[14:]
        run = runShort+'.'+ mode+'.'+ comp
        self.run = run
        self.runShort = runShort
        self.mode = mode
        if mode == 'serial':
            self.modeCmd = 'MPI=NO'
        else:
            self.modeCmd = 'MPI=YES'
        self.runSrc = runSrc
        self.runCmd = 'RUN='+run
        self.runSrcCmd = 'RUNSRC='+runSrc 
        self.opts = ' '
        if compOpts == 'debug':
            flags='"-O0 -g"'
            self.opts += 'EXTRA_FFLAGS='+flags
        elif compOpts == 'traps':
            self.opts += 'COMPILE_WITH_TRAPS=YES'
        self.results =  [runShort, comp, mode, ' - ', ' - ', ' - ', ' - ']

"""
  Read options from a config file for each modelE runSource
  configuration. If no file is available then use some sensible
  defaults.
  TODO: This should be a class, in order to avoid all the globals.
"""
def readConfig(rundeck):
    global debug
    global npList
    global modeList
    global compiler
    global compOpts
    global baseDir
    global resdir
    global deck
    global branch
    global decksDir
    global updBase
    global successMark
    global failMark

    successMark = '+'
    failMark    = 'F'

    # Expect to find a configuration file MYCONFIGDIR directory
    if os.environ.has_key('MYCONFIGDIR'):
        myConfigDir = os.environ['MYCONFIGDIR']
    else:
        myConfigDir = './'

    configfile = myConfigDir + '/' + rundeck + '.cfg'
    if os.path.isfile(configfile):
        config = ConfigParser.RawConfigParser()
        config.read(configfile)
    # If there is no config file create default options
    else:
        config = ConfigParser.RawConfigParser()
        config.add_section('regSettings')
        if rundeck == 'nonProduction_E_AR5_C12':
            config.set('regSettings', 'rundeck'   ,'nonProduction_E_AR5_C12')
        else:
            config.set('regSettings', 'rundeck'   ,rundeck)

        config.set('regSettings', 'modelerc'  , os.environ['HOME']+'/.modelErc')
        config.set('regSettings', 'compiler'  ,'gfortran')
        config.set('regSettings', 'modes'     ,'serial,mpi')
        config.set('regSettings', 'nplist'    ,'1,4')
        config.set('regSettings', 'compflags' ,'debug')
        config.set('regSettings', 'branch'    ,'master')
        config.set('regSettings', 'basedir'   ,'.')
        config.set('regSettings', 'updatebase','no')
        config.set('regSettings', 'resultsdir','.')
        config.set('regSettings', 'decksdir'  ,'.')

    deck      = config.get('regSettings','rundeck')
    compiler  = config.get('regSettings','compiler')

    modes     = config.get('regSettings','modes')
    modeList = []
    for mode in modes.split(','):
        if mode == 'serial' or mode == 'mpi':
            modeList.append(mode)
        # Terminate job if there is a non-permitted mode
        else:
            print ' *** Incorrect mode *** ' + mode
            sys.exit(1)

    npes      = config.get('regSettings','nplist')
    npList = []
    for np in npes.split(','):
        npList.append(np)
    compOpts  = config.get('regSettings','compflags')
    baseDir   = config.get('regSettings','basedir')
    # branch is needed to select correct basedir files
    branch    = config.get('regSettings','branch')
    updBase   = config.get('regSettings','updatebase')
    resdir    = config.get('regSettings','resultsdir')
    decksDir  = config.get('regSettings','decksdir')

    # In order to avoid errors in compareBase() when running with
    # default options:
    if baseDir != '.':
        baseDir =  baseDir + '/' + branch + '/' + compiler

    if os.environ.has_key('MODELERC'):
        modelerc = os.environ['MODELERC']
    else:
        modelerc  = config.get('regSettings','modelerc')
        os.environ['MODELERC'] = modelerc

    if os.environ.has_key('DEBUG'):
        debug = os.environ['DEBUG']
    else:
        debug = False

"""
  Setup a logging object for each run. Note that the output file gets
  all the logging output while STDOUT only gets logging INFO in order
  to minimize verbosity.
"""
def setupLogging():
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=resdir+'/'+deck+'-regression.log',
                    filemode='w')
    stdoutLog = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(name)s : %(message)s')
    stdoutLog.setFormatter(formatter)
    stdoutLog.setLevel(logging.INFO)
    logger = logging.getLogger()
    logger.addHandler(stdoutLog)

"""
  Define a subprocess call where we can capture/evaluate stderr output
"""
def sysCall(cmd):
    logger = logging.getLogger('SYSTEM  ')
    if debug:
        logger.debug(cmd)
    else:
        logger.debug(cmd)
        p = subprocess.Popen(shlex.split(cmd), \
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr =  p.communicate()
        return stdout

"""
  A subprocess call that, upon failure rc<>=0, raises an exception
  Subprocess call with shell argument:
  expands environment variables and file globs 
"""
def sysCmd(commandString):
    logger = logging.getLogger('SYSTEM  ')
    if debug:
        logger.debug(commandString)
    else:
        logger.debug(commandString)
        makeLog = resdir + '/'  + deck + '-make.log'
        with open(makeLog,'a') as f:
                status = subprocess.call(commandString, \
                                         stdout=f, stderr=f, shell=True)
        logger.debug('Return code: ' + str(status))
        if (status != 0):
            raise Exception('unix', commandString)

"""
  TODO !!! This is the function we should be using, rather than SysCmd !!!
  A subprocess call that, upon failure rc<>=0, raises an exception
  Subprocess call with shell argument:
  expands environment variables and file globs 
"""
def sysCmd0(commandString, rc=0, sh=False):
    logger = logging.getLogger('SYSTEM  ')
    if debug:
        logger.debug(commandString)
    else:
        logger.debug(commandString)
        makeLog = resdir + '/'  + deck + '-make.log'
        with open(makeLog,'a') as f:
            if sh:
                status = subprocess.call(commandString, \
                                         stdout=f, stderr=f, shell=sh)
            else:
                status = subprocess.call(shlex.split(commandString), \
                                         stdout=f, stderr=f, shell=sh)
        logger.debug('status: ' + str(status))
        if (status != rc):
            raise Exception('unix', commandString)

"""
  Return a checkpoint file name with various identifiers
"""
def checkpointName(exp, duration, npes):
    if exp.mode == 'serial':
        return exp.run  + '.' + duration
    else:
        return exp.run + '.' + duration + '.np=' + str(npes)

"""
  Build a configuration using GNU make
"""
def build(exp):
    logger = logging.getLogger('BUILD   ')
    logger.info(exp.run + ' ' + exp.modeCmd + ' ' + exp.opts)
    try:
        sysCmd('make --quiet clean')
        sysCmd('make rundeck ' + exp.runCmd + ' ' + exp.runSrcCmd)
        sysCmd('make -j4 gcm ' + exp.runCmd + ' ' + exp.modeCmd + ' ' + exp.opts)
        exp.results[3] = successMark
    except:
        exp.results[3] = failMark
        logger.error(' *** Failed to build ' + exp.run)
        raise

"""
  Sets up and runs a 1hr simulation
"""
def run1hr(exp, npes=1):
    logger = logging.getLogger('RUN1HR  ')
    logger.info(exp.run + ', ' + exp.mode + ', npes=' + str(npes))
    try:
        sysCmd('make setup ' + exp.runCmd + ' ' + exp.modeCmd + ' ' + exp.opts)
        sysCmd('../exec/runE ' + exp.run + ' -np ' + str(npes) + ' -cold-restart')
        sysCmd('cd '+exp.run+'; cp fort.2.nc '+checkpointName(exp, '1hr', npes))
    except:
        exp.results[3] = failMark
        message =  ' *** Failed to run 1 hour test for ' + exp.run
        message += ' on ' + str(npes) + ' processors.'
        logger.error(message)
        raise
    
"""
  Runs (N-M)+M hours AND N continuous hours
  Default is to run 24+1 and 25-hour, i.e. N=25 M=1
"""
def runRestart(exp, npes=1, n=25, m=1):
    logger = logging.getLogger('RUNRST  ')
    expectedRC = 13; # modelE convention for successful runs
    restart = './'+exp.run
    if exp.mode == 'mpi':
        restart += ' -np ' + str(npes)
    
    logger.info(exp.run + ', ' + exp.mode + ', npes=' + str(npes))
    try:
        sysCmd('../exec/editRundeck.sh ' + exp.run + ' 48 2 1')
        sysCmd('make setup ' + exp.runCmd + ' ' + exp.modeCmd + ' ' + exp.opts)
        sysCmd('../exec/runE ' + exp.run + ' -np ' + str(npes) + ' -cold-restart')
        sysCmd('cd '+exp.run+'; cp fort.1.nc '+checkpointName(exp, '1dy', npes))
        sysCmd('cd ' + exp.run + '; cp fort.2.nc fort.1.nc')
        sysCmd('cd ' + exp.run + '; rm -f run_status')
#  Need to investigate why the following causes a NameError exception
#  Looks like there is an issue with variable/function/class names in SysCmd
        sysCmd('cd ' + exp.run + '; ' + restart + '; test `head -1 run_status` -eq ' + str(expectedRC))
        sysCmd('cd '+exp.run+';cp fort.2.nc '+checkpointName(exp, 'restart', npes))
    except:
        exp.results[3] = failMark
        message =  ' *** Failed to run 1 day test for ' + exp.run
        message += ' on ' + str(npes) + ' processors.'
        logger.error(message)
        raise
    
"""
  Compare model results with those in the baseline location.
  If no baseline location is specified then comparison will be skipped.
"""
def compareBase(exp, duration, npes=1):
    # Skip if no location given
    if baseDir == '.':
        return

    logger = logging.getLogger('COMPBAS ')
    logger.info('Compare base run: '+exp.run)

    prefix = exp.run + '/'
    file1 = prefix + checkpointName(exp, duration, npes)
    file2 = baseDir + '/' + checkpointName(exp, duration, npes)
    cmp ='diffreport.x ' + file1 + ' ' + file2
    rc = sysCall(cmp)
    if rc == '':
        exp.results[4] = successMark
    else:
        logger.warning(file1 + ' and ' + file1 + ' differ')
        if updBase == 'yes':
            sysCmd('cp ' + file1 + ' ' + file2)
            logger.info('Updated BASELINE')
        exp.results[4] = failMark

"""
  Compare SERIAL vs MPI
"""
def compareNPE(runA, runB, duration, npes):
    logger = logging.getLogger('COMPNPE ')
    logger.info('Compare NPE runs: '+runA.run+' and '+runB.run)

    prefix1 = runA.run + '/'
    prefix2 = runB.run + '/'
    file1 = prefix1 + checkpointName(runA, duration, npes)
    file2 = prefix2 + checkpointName(runB, duration, npes)
    cmp ='diffreport.x ' + file1 + ' ' + file2
    rc = sysCall(cmp)
    if rc == '':
        runB.results[6] = successMark
    else:
        logger.warning('Files ' + file1 + ' and ' + file2 + ' differ')
        runB.results[6] = failMark
            
"""
  Compare full-run (25hr) vs restart run
"""
def compareRestart(exp, npes=1):
    # Hack to skip SCM rundeck
    if exp.run == 'SGP4TESTS':
        return
    logger = logging.getLogger('COMPRST ')
    logger.info('Compare restart run: '+exp.run)
    prefix = exp.run + '/'
    file1 = prefix + checkpointName(exp, '1dy', npes)
    file2 = prefix + checkpointName(exp, 'restart', npes)
    cmp ='diffreport.x ' + file1 + ' ' + file2
    rc = sysCall(cmp)
    if rc == '':
        exp.results[5] = successMark
    else:
        logger.warning('Files ' + file1 + ' and ' + file2 + ' differ')
        # Hack to differentiate the restart errors in CAD rundecks:
        if 'E4Tcad' in exp.run:
            exp.results[5] = failMark+'*'
        else:
            exp.results[5] = failMark
        
"""
  MAIN PROGRAM
"""
if __name__ == '__main__':
    
    runSources = []
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            runSources.append(arg)
    else:
        runSources.append('nonProduction_E_AR5_C12')
    
    for rundeck in runSources:

        readConfig(rundeck)
        setupLogging()
        os.chdir(decksDir)
        logger = logging.getLogger('MAIN    ')
        diffFile = resdir + '/' + rundeck + '.diff'
        fileH = open(diffFile, 'w')
        logger.info('Testing ' + rundeck)

        exps = []
        for mode in modeList:
            exps.append(configuration(rundeck, mode, compiler))
            nmodes = len(exps)

        for exp in exps:
            try:
                if exp.mode == 'serial':
                    build(exp)
                    run1hr(exp)
                    runRestart(exp)
                else:
                    build(exp)
                    try:
                        for npes in npList:
                            run1hr(exp, npes=npes)
                    except:
                        logger.error('  ... abandoning 1hr mpi run.')

                    try:
                        for npes in npList:
                            runRestart(exp, npes=npes)
                    except:
                        logger.error('  ... abandoning restart mpi run.')

                logger.info(rundeck + ' ' + exp.mode + ' runs complete.')
                
            except:
                logger.error(rundeck + ' ' + exp.mode + ' runs FAILED')

        for exp in exps:
            if exp.mode == 'serial':
                compareBase(exp, '1hr')
                compareBase(exp, '1dy')
                # And compare SERIAL checkpoint-restart 
                compareRestart(exp)
            else:
                for npes in npList:
                    # Compare runs with baseline
                    compareBase(exp, '1hr', npes=npes)
                    compareBase(exp, '1dy', npes=npes)
                    compareRestart(exp, npes=npes)

        for npes in npList:
            # Compare 1hr run against serial
            if nmodes > 1:
                compareNPE(exps[0], exps[1], '1dy', npes)

        logger.info(rundeck + ' comparisons complete.')

            
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
            
                    
                    
