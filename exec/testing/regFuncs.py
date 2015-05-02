import sys
import os
import shlex
import subprocess
import logging
import ConfigParser
import regUtils

""" 
  This class assigns settings used to test a given rundeck 
"""
class newRundeck:
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
  passed to the the makefile. Note that a "run" has a "rundeck".
"""
class newRun():
    def __init__(self, rundeck, mode):
        # newRun has-a rundeck
        self.runsrc = rundeck
        # shortName is introduced to remove the nonProduction part of rundeck name
        shortName = rundeck.name
        if 'nonProduction' in rundeck.name:
            shortName = rundeck.name[14:]
        # A RUN name is composed of rundeck name, mode and compiler
        self.name = shortName+'.'+mode+'.'+rundeck.compiler
        # The following are (cmd)strings used by gnu-make
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
        # Run results are stored in an array
        self.results = [shortName, rundeck.compiler, mode, 
                        ' - ', ' - ', ' - ', ' - ']
        self.successMark   = '+'
        self.failMark      = 'F'
        # suffix used to identify length of runs
        self.endTime      = str(rundeck.endtime) + 'hr'

    # A subprocess call that, upon failure rc<>=0, raises an exception.
    # Class membership for this function is one of convenience: need runSource 
    # and newRun data (resultsDir and name)
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
   Set debug ENV variable and diffreport executable
"""
def setRunUtils():
    global debug
    global diffreportExe

    if os.environ.has_key('DEBUG'):
        debug = os.environ['DEBUG']
    else:
        debug = False

    # This is needed to find diffreport.x, assumed to be in $HOME/bin
    os.environ["PATH"] += os.pathsep + os.environ["HOME"] \
      + '/bin'
    diffreportExe = regUtils.which('diffreport.x')
    if diffreportExe is None:
        print 'No available diffreport.x. Will use diff'
        diffreportExe = 'diff'


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
        # rundeck settings
        rundeck.name      = config.get('regSettings', 'rundeck')
        rundeck.compiler  = config.get('regSettings', 'compiler')
        rundeck.modes     = config.get('regSettings', 'modes')
        rundeck.testLevel = config.get('regSettings', 'testlevel')
        rundeck.endtime   = config.getint('regSettings', 'endtime')
        rundeck.npes      = config.get('regSettings', 'nplist')
        rundeck.compilerFlags  = config.get('regSettings', 'compflags')
        # system settings
        rundeck.baseDir   = config.get('regSettings', 'basedir')
        rundeck.branch    = config.get('regSettings', 'branch')
        rundeck.updateBase   = config.get('regSettings', 'updatebase')
        rundeck.resultsDir    = config.get('regSettings', 'resultsdir')
        rundeck.decksDir  = config.get('regSettings','decksdir')
        rundeck.systemTests = config.get('regSettings', 'systemtests')
    
    else: # There is no config file, so use default options
        config = ConfigParser.RawConfigParser()
        config.add_section('regSettings')
        # Use defaults
        # rundeck settings
        config.set('regSettings', 'rundeck'   , rundeck.name)
        config.set('regSettings', 'compiler'  , rundeck.compiler)
        config.set('regSettings', 'modes'     , rundeck.modes)
        config.set('regSettings', 'testlevel' , rundeck.testLevel)
        config.set('regSettings', 'endtime'   , rundeck.endtime)
        config.set('regSettings', 'nplist'    , rundeck.npList)
        config.set('regSettings', 'compflags' , rundeck.compilerFlags)
        # system settings
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
  Return a checkpoint file name with various identifiers
"""
def checkpointName(exp, duration, npes):
    if exp.mode == 'serial':
        return exp.name  + '.' + duration
    else:
        return exp.name + '.' + duration + '.np=' + str(npes)

    
"""
  Build model using GNU make
"""
def build(exp):
    logger = logging.getLogger('BUILD   ')
    logger.info(exp.name + ' ' + exp.modeCmd + ' ' + exp.xflags)
    status = exp.sysCmd('make --quiet clean', 3, 'b')
    status = exp.sysCmd('make rundeck ' + exp.runCmd + ' ' + exp.runSrcCmd, 
               3, 'b')
    status = exp.sysCmd('make -j gcm ' + exp.runCmd + ' ' + exp.modeCmd
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
                    + ', endtime=' + exp.endTime)

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
               + checkpointName(exp, exp.endTime, npes), 3, 'r')
    status = exp.sysCmd('cd ' + exp.name + '; cp fort.2.nc fort.1.nc', 3, 'r')
    status = exp.sysCmd('cd ' + exp.name + '; rm -f run_status lock', 3, 'r')
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
  Run 2 month run
"""
def runLong(exp, npes=1):
    logger = logging.getLogger('RUNLONG ')
    logger.info(exp.name + ', ' + exp.mode + ', npes=' + str(npes) \
                    + ', 2 month run')

    # Most rundecks, except hycom, have MONTHI=12, so set MONTHE=2. 
    # For hycom, MONTHI=1, so set MONTHE=3. Note ndisk=480.
    if 'h4c' in exp.name:
        newTime = ' 480 1 0 3'
    else:
        newTime = ' 480 1 0 2'

    status = exp.sysCmd('../exec/editRundeck.sh ' + exp.name + newTime,
               3, 'r')
    status = exp.sysCmd('make setup ' + exp.runCmd + ' ' + exp.modeCmd + ' '
               + exp.xflags, 3, 'r')
    status = exp.sysCmd('../exec/runE ' + exp.name + ' -np ' + str(npes)
               + ' -cold-restart', 3, 'r')

    
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
  Compare continuous-run vs restart run
"""
def compareRestart(exp, npes=1):
    logger = logging.getLogger('COMPRST ')
    logger.info('Compare restart run: '+exp.name)
    prefix = exp.name + '/'
    file1 = prefix + checkpointName(exp, exp.endTime, npes)
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

 
