"""
  This script verifies that MPI and Serial builds produce identical
  results. The script can be executed from the decks directory using
  default options and without arguments as follows:
 
      python  ../exec/testing/regression.py
 
  In that case the script will run the nonProduction_E_AR5_C12 rundeck 
  using the gfortran compiler and in serial and mpi modes. 
  Alternatively one can use the default options with one argument:
 
      python  ../exec/testing/regression.py <runsource>
 
  and use the <runsource> rundeck.
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

""" 
  This class assigns settings used to test a given rundeck 
"""
class RunSourceProperties:
    def __init__(self, sourceName='nonProduction_E_AR5_C12'):
        # Initialize with defaults
        if sourceName:
            self.name      = sourceName
        else:
            self.name      = 'nonProduction_E_AR5_C12'
        self.compiler      = 'gfortran'
        self.compilerFlags = 'debug'
        self.modes         = 'serial,mpi'
        self.testLevel     = 'full'
        self.endtime       = 25
        self.npes          = '1,4'
        self.modeList      = []
        self.npList        = []
        self.branch        = 'master'
        self.updateBase    = 'no'
        self.systemTests   = 'no'
        self.modelerc      = os.environ['HOME']+'/.modelErc'
        self.baseDir       = '.'
        self.resultsDir    = '.'
        self.decksDir      = '.'
        # Now override defaults with values specified in config file
        getConfigFile(self)

        
"""
  This class defines, among other things, a rundecks's run time options
  passed to the the makefile.
"""
class Arun():
    def __init__(self, rundeck, mode):
        # Arun has-a rundeck
        self.runsrc = rundeck
        shortName = rundeck.name
        if 'nonProduction' in rundeck.name:
            shortName = rundeck.name[14:]
        self.name = shortName+'.'+mode+'.'+rundeck.compiler
        self.runCmd = 'RUN='+self.name
        self.runSrcCmd = 'RUNSRC='+rundeck.name
        if mode == 'serial':
            self.mode = 'serial'
            self.modeCmd = 'MPI=NO'
        else:
            self.mode = 'mpi'
            self.modeCmd = 'MPI=YES'
        self.xflags = ' '
        if rundeck.compilerFlags == 'debug':
            flags='"-O0 -g"'
            self.xflags += 'EXTRA_FFLAGS='+flags
        elif rundeck.compilerFlags == 'traps':
            self.xflags += 'COMPILE_WITH_TRAPS=YES'
        self.results = [shortName, rundeck.compiler, mode, 
                        ' - ', ' - ', ' - ', ' - ']
        self.successMark   = '+'
        self.failMark      = 'F'
        self.etSuffix      = str(rundeck.endtime) + 'hr'

    # A subprocess call that, upon failure rc<>=0, raises an exception.
    # Class membership for this function is one of convenience: need runSource 
    # and Arun data (resultsDir and name)
    def sysCmd(self, commandString, result, stage):
        logger = logging.getLogger('SYSTEM  ')
        status = 0
        if debug:
            logger.debug(commandString)
        else:
            logger.debug(commandString)
            makeLog = self.runsrc.resultsDir + '/'  + self.name + '-make.log'            
            with open(makeLog,'a') as f:
                status = subprocess.call(commandString, \
                                         stdout=f, stderr=f, shell=True)
            logger.debug('Return code: ' + str(status))
            if (status == 0):
                self.results[result] = self.successMark
            else:
                logger.error(commandString+': FAILED')
                self.results[result] = self.failMark+stage
        return status


"""
  Read options from a config file for each modelE rundeck configuration
  If no file is available then use some reasonable defaults.
"""
def getConfigFile(rundeck):

    # If MYCONFIGDIR is defined get configuration file from there
    if os.environ.has_key('MYCONFIGDIR'):
        myConfigDir = os.environ['MYCONFIGDIR']
    else: # it is in the current directory
        myConfigDir = '.'

    configfile = myConfigDir + '/' + rundeck.name + '.cfg'
    if os.path.isfile(configfile):
        config = ConfigParser.RawConfigParser()
        config.read(configfile)
        rundeck.name      = config.get('regSettings', 'rundeck')
        rundeck.compiler  = config.get('regSettings', 'compiler')

        rundeck.modes     = config.get('regSettings', 'modes')
        rundeck.testLevel = config.get('regSettings', 'testlevel')
        rundeck.endtime   = config.getint('regSettings', 'endtime')
        rundeck.npes      = config.get('regSettings', 'nplist')
        rundeck.compilerFlags  = config.get('regSettings', 'compflags')
        rundeck.baseDir   = config.get('regSettings', 'basedir')
        # branch is needed to select correct basedir files
        rundeck.branch    = config.get('regSettings', 'branch')
        rundeck.updateBase   = config.get('regSettings', 'updatebase')
        rundeck.resultsDir    = config.get('regSettings', 'resultsdir')
        rundeck.decksDir  = config.get('regSettings','decksdir')
        rundeck.systemTests = config.get('regSettings', 'systemtests')
    
    else: # There is no config file, so use default options
        config = ConfigParser.RawConfigParser()
        config.add_section('regSettings')
        # Use defaults
        config.set('regSettings', 'rundeck'   , rundeck.name)
        config.set('regSettings', 'compiler'  , rundeck.compiler)
        config.set('regSettings', 'modes'     , rundeck.modes)
        config.set('regSettings', 'testlevel' , rundeck.testLevel)
        config.set('regSettings', 'endtime'   , rundeck.endtime)
        config.set('regSettings', 'nplist'    , rundeck.npList)
        config.set('regSettings', 'compflags' , rundeck.compilerFlags)
        config.set('regSettings', 'branch'    , rundeck.branch)
        config.set('regSettings', 'basedir'   , rundeck.baseDir)
        config.set('regSettings', 'updatebase', rundeck.updateBase)
        config.set('regSettings', 'resultsdir', rundeck.resultsDir)
        config.set('regSettings', 'decksdir'  , rundeck.decksDir)
        config.set('regSettings', 'systemtests', rundeck.systemTests)

    # If defined, use modelErc from environment
    if os.environ.has_key('MODELERC'):
        rundeck.modelerc = os.environ['MODELERC']
    else:
        os.environ['MODELERC'] = config.get('regSettings', 'modelerc')

    # Additional postprocessing
    for mode in rundeck.modes.split(','):
        if mode == 'serial' or mode == 'mpi':
            rundeck.modeList.append(mode)
        # Terminate job if there is a non-permitted mode
        else:
            print ' *** Incorrect mode *** ' + mode
            sys.exit(1)

    for np in rundeck.npes.split(','):
        rundeck.npList.append(np)

    # In order to avoid errors in compareBase() when running with
    # default options:
    if rundeck.baseDir != '.':
        rundeck.baseDir =  rundeck.baseDir + '/' + rundeck.branch + '/' \
        + rundeck.compiler

        
"""
  Setup a logging object for each rundeck. Note that the output file gets
  all the logging output while STDOUT only gets logging INFO in order to
  minimize verbosity.
"""
def setupLogging(rundeck):
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=rundeck.resultsDir+'/'+rundeck.name
                             +'-regression.log',
                    filemode='w')
    stdoutLog = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(name)s : %(message)s')
    stdoutLog.setFormatter(formatter)
    stdoutLog.setLevel(logging.INFO)
    logger = logging.getLogger()
    logger.addHandler(stdoutLog)


"""
  Test if an executable program exists in the path - like unix's which
"""
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
   
 
"""
  Return a checkpoint file name with various identifiers
"""
def checkpointName(exp, duration, npes):
    if exp.mode == 'serial':
        return exp.name  + '.' + duration
    else:
        return exp.name + '.' + duration + '.np=' + str(npes)

    
"""
  Build a configuration using GNU make
"""
def build(exp):
    logger = logging.getLogger('BUILD   ')
    logger.info(exp.name + ' ' + exp.modeCmd + ' ' + exp.xflags)
    status = exp.sysCmd('make --quiet clean', 3, 'b')
    status = exp.sysCmd('make rundeck ' + exp.runCmd + ' ' + exp.runSrcCmd, 
               3, 'b')
    status = exp.sysCmd('make -j4 gcm ' + exp.runCmd + ' ' + exp.modeCmd
               + ' ' + exp.xflags, 3, 'b')
    return status
    
"""
  Sets up and runs a 1hr simulation
"""
def run1hr(exp, npes=1):
    logger = logging.getLogger('RUN1HR  ')
    logger.info(exp.name + ', ' + exp.mode + ', npes=' + str(npes))
    status = exp.sysCmd('make setup ' + exp.runCmd + ' ' + exp.modeCmd + ' '
               + exp.xflags, 3, '1')
    status = exp.sysCmd('../exec/runE ' + exp.name + ' -np ' + str(npes)
               + ' -cold-restart', 3, '1')
    status = exp.sysCmd('cd ' + exp.name + '; cp fort.2.nc ' +
               checkpointName(exp, '1hr', npes), 3, '1')

    
"""
  Run up to ENDTIME hrs with checkpoint at ENDTIME-1 hrs
"""
def runRestart(exp, npes=1, endtime=25):
    logger = logging.getLogger('RUNRST  ')
    expectedRC = 13; # modelE convention for successful runs
    restart = './'+exp.name
    if exp.mode == 'mpi':
        restart += ' -np ' + str(npes)

    logger.info(exp.name + ', ' + exp.mode + ', npes=' + str(npes) \
                    + ', endtime=' + exp.etSuffix)

    checkPt = endtime - 1
    ndisk = checkPt * 2
    if endtime > 24:
        newTime = ' ' + str(ndisk) + ' 2 1'
    else:
        newTime = ' ' +  str(ndisk) + ' 1 ' + str(endtime)

    status = exp.sysCmd('../exec/editRundeck.sh ' + exp.name + newTime,
               3, 'r')
    status = exp.sysCmd('make setup ' + exp.runCmd + ' ' + exp.modeCmd + ' '
               + exp.xflags, 3, 'r')
    status = exp.sysCmd('../exec/runE ' + exp.name + ' -np ' + str(npes)
               + ' -cold-restart', 3, 'r')
    status = exp.sysCmd('cd ' + exp.name + '; cp fort.1.nc '
               + checkpointName(exp, exp.etSuffix, npes), 3, 'r')
    status = exp.sysCmd('cd ' + exp.name + '; cp fort.2.nc fort.1.nc', 3, 'r')
    status = exp.sysCmd('cd ' + exp.name + '; rm -f run_status lock', 3, 'r')
#  Need to investigate why the following causes a NameError exception
#  Looks like there is an issue with variable/function/class names in SysCmd
    status = exp.sysCmd('cd ' + exp.name + '; ' + restart
               + '; test `head -1 run_status` -eq ' + str(expectedRC),
               3, 'r')
    status = exp.sysCmd('cd ' + exp.name + ';cp fort.2.nc '
               + checkpointName(exp, 'restart', npes), 3, 'r')
# Reset rundeck settings for next MPI run
    if npes > 1:
        status = exp.sysCmd('make rundeck ' + exp.runCmd + ' ' 
                            + exp.runSrcCmd, 3, 'r')

    
"""
  Compare model results with those in the baseline location.
  If no baseline location is specified then comparison will be skipped.
"""
def compareBase(exp, duration, npes=1):
    logger = logging.getLogger('COMPBAS ')
    # Skip comparison if no baseDir location is given
    if exp.runsrc.baseDir == '.':
        logger.info('No baseline directory - nothing to do')
        return
    logger.info('Compare base run: '+exp.name)
    prefix = exp.name + '/'
    file1 = prefix + checkpointName(exp, duration, npes)
    file2 = exp.runsrc.baseDir + '/' + checkpointName(exp, duration, npes)
    logger.debug(diffreportExe+' '+file1+' '+file2)
    rc = subprocess.check_output([diffreportExe, file1, file2])
    if rc == '':
        exp.results[4] = exp.successMark
    else:
        exp.results[4] = exp.failMark   
        logger.warning('Baseline reproducibility failed')
        if exp.runsrc.updateBase == 'yes':
            if subprocess.call(['cp', file1, file2]) == 0:
                logger.info('Updated BASELINE')
            else:
                logger.error('Error in: cp '+file1+' ' +file2)
        else:
            logger.info('Consider updating BASELINE')

        
"""
  Compare SERIAL vs MPI
"""
def compareNPE(runA, runB, duration, npes):
    logger = logging.getLogger('COMPNPE ')
    logger.info('Compare NPE runs: '+runA.mode+' and '+runB.mode)
    prefix1 = runA.name + '/'
    prefix2 = runB.name + '/'
    file1 = prefix1 + checkpointName(runA, duration, npes)
    file2 = prefix2 + checkpointName(runB, duration, npes)
    logger.debug(diffreportExe+' '+file1+' '+file2)
    rc = subprocess.check_output([diffreportExe, file1, file2])
    if rc == '':
        runB.results[6] = runB.successMark
    else:
        runB.results[6] = runB.failMark
        logger.warning('NPE reproducibility failed')

        
"""
  Compare full-run (25hr) vs restart run
"""
def compareRestart(exp, npes=1):
    logger = logging.getLogger('COMPRST ')
    logger.info('Compare restart run: '+exp.name)
    prefix = exp.name + '/'
    file1 = prefix + checkpointName(exp, '1dy', npes)
    file2 = prefix + checkpointName(exp, 'restart', npes)
    logger.debug(diffreportExe+' '+file1+' '+file2)
    rc = subprocess.check_output([diffreportExe, file1, file2])
    if rc == '':
        exp.results[5] = exp.successMark
    else:
        # Hack to differentiate the restart errors in CAD/SCM rundecks:
        if 'E4Tcad' in exp.name or 'SGP' in exp.name:
            exp.results[5] = exp.failMark+'*'
        else:
            exp.results[5] = exp.failMark
            logger.warning('Restart reproducibility failed')

        
"""
  MAIN PROGRAM
"""
if __name__ == '__main__':

    global debug
    if os.environ.has_key('DEBUG'):
        debug = os.environ['DEBUG']
    else:
        debug = False

    # System call return code
    OK = 0

    # This is needed to find diffreport.x, assumed to be in $HOME/bin
    os.environ["PATH"] += os.pathsep + os.environ["HOME"] \
      + '/bin'
    diffreportExe = which('diffreport.x')
    if diffreportExe is None:
        print 'No available diffreport.x. Will use diff'
        diffreportExe = 'diff'

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
        rundeck = RunSourceProperties(source)
        
        # List of rundeck run configurations for each mode
        exps = []
        for mode in rundeck.modeList:
            exps.append(Arun(rundeck,mode))
            nmodes = len(exps)

        # Setup a logging stream
        setupLogging(rundeck)
        logger = logging.getLogger('MAIN    ')
        # All the work is done from the modelE decks directory
        os.chdir(rundeck.decksDir)
        # For each rundeck create a diffFile with verification results
        diffFile = rundeck.resultsDir + '/' + rundeck.name + '.diff'
        fileH = open(diffFile, 'w')

        logger.info('Testing ' + rundeck.name)

        for exp in exps:

            serBuildResult = OK
            mpiBuildResult = OK
            if exp.mode == 'serial':
                serBuildResult = build(exp)
                if rundeck.testLevel != 'compileOnly':
                    if serBuildResult == OK:
                        run1hr(exp)
                        if rundeck.testLevel != 'run1hr':
                            runRestart(exp, endtime=rundeck.endtime)
            else:
                mpiBuildResult = build(exp)
                if rundeck.testLevel != 'compileOnly':
                    if mpiBuildResult == OK:
                        for npes in rundeck.npList:
                            run1hr(exp, npes=npes)
                            if rundeck.testLevel != 'run1hr':
                                runRestart(exp, npes=npes, endtime=rundeck.endtime)
            logger.info(rundeck.name + ' ' + exp.mode + ' runs complete.')
            if serBuildResult != OK or mpiBuildResult != OK:
                continue

            if rundeck.testLevel != 'compileOnly':
                if exp.mode == 'serial':
                    compareBase(exp, '1hr')
                    compareBase(exp, '1dy')
                # And compare SERIAL checkpoint-restart 
                    compareRestart(exp)
                else:
                    for npes in rundeck.npList:
                    # Compare runs with baseline
                        compareBase(exp, '1hr', npes=npes)
                        compareBase(exp, '1dy', npes=npes)
                        compareRestart(exp, npes=npes)
                        for npes in rundeck.npList:
                    # Compare 1hr run against serial
                            if nmodes > 1:
                                compareNPE(exps[0], exps[1], '1dy', npes)
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
            
                    
                    
