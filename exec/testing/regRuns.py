# This module contains functions that help setup and compare modelE runs
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
        self.standalone    = 'yes'
        self.verification  = 'restartRun'
        self.compiler      = 'gfortran'
        self.buildType     = 'release'
        self.modes         = 'serial,mpi'
        self.endtime       = 25
        self.npes          = '1,4'
        self.modeList      = []
        self.npList        = []
        self.branch        = 'master'
        self.updateBase    = 'no'
        self.modelerc      = os.environ['HOME']+'/.modelErc'
        self.baseDir       = '.'
        self.resultsDir    = '.'
        self.repository    = '..'
        self.scratchDir    = '.'
        self.decksDir      = '.'
        self.makesystem    = 'makeOld'
        self.savedisk      = '.' # Initialized via .modelErc

        # Now override defaults with values specified in config file
        getConfigFile(self)

        
"""
  This class defines, among other things, a rundecks's run time options
  passed to the the makefile. Note that a "run" is a "rundeck".
"""
class newRun():
    def __init__(self, rundeck, mode):
        # newRun has-a rundeck
        self.runsrc = rundeck
        # shortName is introduced to remove the nonProduction part of a
        # rundeck name
        shortName = rundeck.name
        if 'nonProduction' in rundeck.name:
            shortName = rundeck.name[14:]
        # A RUN name is composed of rundeck name, mode and compiler
        self.shortName = shortName
        self.name = shortName+'.'+mode+'.'+rundeck.compiler
        # The following are (cmd)strings used by gnu-make
        if rundeck.makesystem == 'makeOld':
            self.runCmd = 'RUN='+self.name
            self.runSrcCmd = 'RUNSRC='+rundeck.name
        else:
            self.runCmd = self.name
            self.runSrcCmd = rundeck.name
        if mode == 'serial':
            self.mode = 'serial'
            self.modeCmd = 'MPI=NO'
        else:
            self.mode = 'mpi'
            self.modeCmd = 'MPI=YES'
        self.xflags = ' '
        if rundeck.buildType == 'debug':
            flags='"-O0 -g"'
            self.xflags += 'EXTRA_FFLAGS='+flags
        elif rundeck.buildType == 'traps':
            self.xflags += 'COMPILE_WITH_TRAPS=YES'
        # Run results are stored in an array
        self.results = [shortName, rundeck.compiler, mode, 
                        '  -  ', '  -  ', '  -  ', '  -  ']
        self.successMark   = '+'
        self.failMark      = 'F'
        # suffix used to identify length of runs
        self.endTime      = str(rundeck.endtime) + 'hr'
# use inheritance?
        self.repository = rundeck.repository
        self.makesystem = rundeck.makesystem
        self.compiler = rundeck.compiler
        self.standalone = rundeck.standalone
        self.verification = rundeck.verification
        self.npList = rundeck.npList
        self.endtime = rundeck.endtime


    # System call that records result of subprocess call
    def sysCmd(self, commandString, result, stage):
        logger = logging.getLogger('SYSTEM  ')
        status = 0
        if debug:
            logger.info(commandString)
        else:
            logger.debug(commandString)
            makeLog = self.runsrc.resultsDir + '/'  + self.name + '-make.log'            
            with open(makeLog,'a') as f:
               status = subprocess.call(commandString, \
                                         stdout=f, stderr=f, shell=True)
            if (status == 0):
                self.results[result] = self.successMark
            else:
                logger.error(commandString+': FAILED')
                self.results[result] = self.failMark+stage
            logger.debug('Return code: ' + str(status))
            logger.debug('Result: ' + self.results[result])
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
        rundeck.standalone = config.get('regSettings', 'standalone')
        rundeck.verification = config.get('regSettings', 'verification')
        rundeck.endtime   = config.getint('regSettings', 'endtime')
        rundeck.npes      = config.get('regSettings', 'nplist')
        rundeck.buildType  = config.get('regSettings', 'buildtype')
        # system settings
        rundeck.baseDir   = config.get('regSettings', 'basedir')
        rundeck.branch    = config.get('regSettings', 'branch')
        rundeck.updateBase   = config.get('regSettings', 'updatebase')
        rundeck.resultsDir    = config.get('regSettings', 'resultsdir')
        rundeck.scratchDir    = config.get('regSettings', 'scratchdir')
        rundeck.decksDir  = config.get('regSettings','decksdir')
        rundeck.repository = config.get('regSettings', 'repository')
        rundeck.makesystem = config.get('regSettings', 'makesystem')
    
    else: # There is no config file, so use default options
        config = ConfigParser.RawConfigParser()
        config.add_section('regSettings')
        # Use defaults in newRundeck
        # rundeck settings
        config.set('regSettings', 'rundeck'   , rundeck.name)
        config.set('regSettings', 'compiler'  , rundeck.compiler)
        config.set('regSettings', 'modes'     , rundeck.modes)
        config.set('regSettings', 'standalone' , rundeck.standalone)
        config.set('regSettings', 'verification' , rundeck.verification)
        config.set('regSettings', 'endtime'   , rundeck.endtime)
        config.set('regSettings', 'nplist'    , rundeck.npList)
        config.set('regSettings', 'buildtype' , rundeck.buildType)
        # system settings
        config.set('regSettings', 'branch'    , rundeck.branch)
        config.set('regSettings', 'basedir'   , rundeck.baseDir)
        config.set('regSettings', 'updatebase', rundeck.updateBase)
        config.set('regSettings', 'resultsdir', rundeck.resultsDir)
        config.set('regSettings', 'scratchdir', rundeck.scratchDir)
        config.set('regSettings', 'decksdir'  , rundeck.decksDir)
        config.set('regSettings', 'repository', rundeck.repository)
        config.set('regSettings', 'makesystem', rundeck.makesystem)

    # If defined, use modelErc from environment
    try:
        rundeck.modelerc = os.environ['MODELERC']
    except:
        try:
            os.environ['MODELERC'] = config.get('regSettings', 'modelerc')
        except:
            print ' *** MODELERC is not defined.'
            print ' *** No option modelerc in section [regSettings].'
            sys.exit(1)

    # Extract SAVEDISK from modelErc file - needed in compareNPE            
    cmd = "grep SAVEDISK " + os.environ['MODELERC'] + "| awk -F= '{print $2}'"
    out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,\
                                stderr=subprocess.STDOUT)
    rundeck.savedisk = out.communicate()[0].rstrip()

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
  logging DEBUG while STDOUT only gets logging INFO in order to minimize 
  verbosity.
"""
def setupLogging(rundeck):
# Note filemode is 'append' because we want one logger per rundeck
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=rundeck.resultsDir+'/'+rundeck.name
                             +'-regression.log',
                    filemode='a')
    stdoutLog = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(name)s : %(message)s')
    stdoutLog.setFormatter(formatter)
    stdoutLog.setLevel(logging.INFO)
    logger = logging.getLogger()
    logger.addHandler(stdoutLog)


"""
  Return a checkpoint file name with various identifiers
"""
def checkpointName(run, endTime, npes):
    if run.mode == 'serial':
        return run.name  + '.' + endTime
    else:
        return run.name + '.' + endTime + '.np=' + str(npes)

    
"""
  Build model
"""
def build(run):
    logger = logging.getLogger('BUILD   ')
    logger.info(run.name + ' ' + run.modeCmd + ' ' + run.xflags)

    if run.makesystem == 'makeOld':
        rc = run.sysCmd('make --quiet clean', 3, 'b')
        if rc != 0:
            return 1
        cmd = 'make rundeck '+run.runCmd+' '+run.runSrcCmd
        rc = run.sysCmd(cmd, 3, 'b')
        if rc != 0:
            return 1
        cmd = 'make -j gcm '+run.runCmd+' '+run.modeCmd+' '+run.xflags
        rc = run.sysCmd(cmd, 3, 'b')
        if rc != 0:
            return 1
    else:
        rc = run.sysCmd(run.repository + '/exec/configure ' + run.runCmd + \
                            ' ' + run.runSrcCmd + ' ' + run.modeCmd + \
                            ' ' + run.xflags, 3, 'b')
        if rc != 0:
            return 1
        rc = run.sysCmd('make -j setup', 3, 'b')
        if rc != 0:
            return 1
    return 0
           
"""
  Sets up and runs a 1hr simulation
"""
def run1hr(run, npes=1):
    logger = logging.getLogger('RUN1HR  ')
    mErc = '13'; # modelE convention for successful runs
    logger.info(run.name + ', ' + run.mode + ', npes=' + str(npes))

    if run.makesystem == 'makeOld':
        cmd =  'make setup '+run.runCmd+' '+run.modeCmd+' '+run.xflags
        rc = run.sysCmd(cmd, 3, '1')
    else:
        rc = run.sysCmd('make setup ', 3, '1')
    if rc != 0:
        return 1

    rune = run.repository+'/exec/runE '
    cmd = rune+run.name+' -np '+str(npes)+' -cold-restart'  
    rc = run.sysCmd(cmd, 3, '1')
    if rc != 0:
        return 1

    cmd = 'cd '+run.name+ '; test `head -1 run_status` -eq ' + mErc
    rc = run.sysCmd(cmd, 3, '1')
    if rc != 0:
        return 1

    cmd = 'cd '+run.name+'; cp fort.2.nc '+checkpointName(run, '1hr', npes)
    rc = run.sysCmd(cmd, 3, '1')
    if rc != 0:
        return 1

    logger.info(run.name + ' is DONE')
    return 0

"""
  Run up to ENDTIME hrs with checkpoint at ENDTIME-1 hrs
"""
def runRestart(run, npes=1, endtime=25):
    logger = logging.getLogger('RUNRST  ')
    mErc = '13'; # modelE convention for successful runs
    restart = './'+run.name
    if run.mode == 'mpi':
        restart += ' -np ' + str(npes)

    logger.info(run.name + ', ' + run.mode + ', npes=' + str(npes) + \
                    ', endtime=' + run.endTime)

    # Always checkpoint 1hr prior to enndtime
    # TODO: make this an optional argument?
    checkPt = endtime - 1
    ndisk = checkPt * 2
    if endtime > 24:
        newTime = ' ' + str(ndisk) + ' 2 1'
    else:
        newTime = ' ' +  str(ndisk) + ' 1 ' + str(endtime)

    cmd = run.repository+'/exec/editRundeck.sh '+run.name+newTime 
    rc = run.sysCmd(cmd, 3, 'r')
    # Run make setup again to update the edited rundeck file
    if run.makesystem == 'makeOld':
        cmd =  'make -j setup '+run.runCmd+' '+run.modeCmd+' '+run.xflags
        rc = run.sysCmd(cmd, 3, 'r')
    else:
        rc = run.sysCmd('make setup ', 3, 'r')
    if rc != 0:
        return 1

    rune = run.repository+'/exec/runE '
    cmd = rune+run.name+' -np '+str(npes)+' -cold-restart'  
    rc = run.sysCmd(cmd, 3, 'r')
    if rc != 0:
        return 1

    cmd = 'cd '+run.name+'; cp fort.1.nc '+checkpointName(run, run.endTime, npes)
    rc = run.sysCmd(cmd, 3, 'r')
    if rc != 0:
        return 1

    cmd = 'cd ' + run.name + '; cp fort.2.nc fort.1.nc; rm -f run_status'
    rc = run.sysCmd(cmd, 3, 'r')
    if rc != 0:
        return 1 

    cmd = 'cd '+run.name+'; '+restart+'; test `head -1 run_status` -eq '+mErc
    rc = run.sysCmd(cmd, 3, 'r')
    if rc != 0:
        return 1

    cmd = 'cd '+run.name+';cp fort.2.nc '+checkpointName(run, 'restart', npes)
    rc = run.sysCmd(cmd, 3, 'r')
    if rc != 0:
        return 1

# Reset rundeck settings for next MPI run
    if len(run.npList) > 1:
        if run.makesystem == 'makeOld':
            makecmd = 'make rundeck '+run.runCmd+' '+run.runSrcCmd
            rc = run.sysCmd(makecmd, 3, 'r')
        else:
            rc = run.sysCmd('make rundeck RUN=' + run.runCmd 
                              + ' RUNSRC='+ run.runSrcCmd, 3, 'r')
    if rc != 0:
        return 1

    logger.info(run.name + ' is DONE')
    return 0

   
"""
  Run 2 month run
"""
def runLong(run, npes=1):
    logger = logging.getLogger('RUNLONG ')
    logger.info(run.name + ', ' + run.mode + ', npes=' + str(npes) \
                    + ', 2 month run')

    # Most rundecks, except hycom, have MONTHI=12, so set MONTHE=2. 
    # For hycom, MONTHI=1, so set MONTHE=3. Note ndisk=480.
    if 'h4c' in run.name:
        newTime = ' 480 1 0 3'
    else:
        newTime = ' 480 1 0 2'

    rc = run.sysCmd(run.repository + '/exec/editRundeck.sh ' + run.name + \
                        newTime, 3, 'r')
    if rc != 0:
        return 1
    if run.makesystem == 'makeOld':
        rc = run.sysCmd('make -j setup ' + run.runCmd + ' ' + run.modeCmd + \
                            ' ' + run.xflags, 3, 'r')
    else:
        rc = run.sysCmd('make -j setup', 3, 'r')
    if rc != 0:
        return 1

    rc = run.sysCmd(run.repository + '/exec/runE ' + run.name + ' -np ' + \
                            str(npes) + ' -cold-restart', 3, 'r')
    if rc != 0:
        return 1

    logger.info(run.name + ' is DONE')
    return 0


"""
  Get number of -diffs- when running diffreport
"""
def getNumDiffs(file1,file2):
    cmd = diffreportExe+' '+file1+' '+file2+' | grep diffs | wc -l'
    diff = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,\
                                stderr=subprocess.STDOUT)
    numDiffs = diff.communicate()[0]
    diff.wait()
    return ''.join(numDiffs.split())

"""
  Compare model results with those in the baseline location.
  If no baseline location is specified then comparison will be skipped.
"""
def compareBase(run, endTime, npes=1):
    logger = logging.getLogger('COMPBAS ')

    # Skip comparison if no baseDir location is given
    if run.runsrc.baseDir == '.':
        logger.info('No baseline directory - nothing to do')
        return

    logger.info('Compare '+run.name+' '+endTime+' base run')
    prefix = run.name + '/'

    # Check if run result exists:
    fileTST = prefix + checkpointName(run, endTime, npes)
    if not os.path.exists(fileTST):
        logger.error('---CHECKPOINT file does not exist.')
        return

    fileBAS = run.runsrc.baseDir + '/' + checkpointName(run, endTime, npes)
    logger.debug(diffreportExe+' '+fileTST+' '+fileBAS)

    n = getNumDiffs(fileTST, fileBAS)
    if n == '0':
#    rc = subprocess.check_output([diffreportExe, fileTST, fileBAS])
#    if rc == '':
        run.results[4] = run.successMark
    else:
        run.results[4] = '{: ^5}'.format(n)
#        run.results[4] = run.failMark   
        logger.warning('---Baseline reproducibility failed')
        if run.runsrc.updateBase == 'yes':
            if subprocess.call(['cp', fileTST, fileBAS]) == 0:
                logger.info('Updated BASELINE')
            else:
                logger.error('Error in: cp '+fileTST+' ' +fileBAS)
        else:
            logger.info('BASELINE not updated')

"""
  Compare continuous-run vs restart run
"""
def compareRestart(run, npes=1):
    logger = logging.getLogger('COMPRST ')
    logger.info('Compare '+run.name+' '+run.endTime+' and restart run')
    prefix = run.name + '/'

    # Check if continuous run result exists:
    fileCON = prefix + checkpointName(run, run.endTime, npes)
    if not os.path.exists(fileCON):
        logger.error('---CHECKPOINT file does not exist.')
        return

    # Check if restart run result exists:
    fileRST = prefix + checkpointName(run, 'restart', npes)
    if not os.path.exists(fileRST):
        logger.error('---RESTART file does not exist.')
        return

    logger.debug(diffreportExe+' '+fileCON+' '+fileRST)
    n = getNumDiffs(fileCON, fileRST)
    if n == '0':
#    rc = subprocess.check_output([diffreportExe, fileCON, fileRST])
#    if rc == '':
        run.results[5] = run.successMark
    else:
        # SCM rundeck is not restart reproducible
        if 'SGP' in run.name:
            run.results[5] = run.failMark+'*'
        else:
            run.results[5] = '{: ^5}'.format(n)
#            run.results[5] = run.failMark
            logger.warning('---Restart reproducibility failed')

 
"""
  Compare SERIAL vs MPI
"""
def compareNPE(runMPI, endTime, npes):
    logger = logging.getLogger('COMPNPE ')
    logger.info('Compare serial and '+runMPI.mode+ \
                    ' (' +str(npes)+' npes) runs')

    # Check if SERIAL result exists:
    if runMPI.makesystem == 'makeOld':
        ddir = '/decks/'
    else:
        ddir = '/'
    fileSER = runMPI.runsrc.savedisk + '/' + \
        runMPI.shortName+'.serial.' + \
        runMPI.runsrc.compiler + '/' + \
        runMPI.shortName+'.serial.' + \
        runMPI.runsrc.compiler + '.' + str(endTime)
    logger.debug('serial file: '+fileSER)
    if not os.path.exists(fileSER):
        logger.warning('---SERIAL file does not exist.')
        return

    # Check if MPI result exists:
    fileMPI = runMPI.name + '/' + checkpointName(runMPI, endTime, npes)
    if not os.path.exists(fileMPI):
        logger.warning('---MPI file does not exist.')
        return

    logger.debug(diffreportExe+' '+fileSER+' '+fileMPI)
    n = getNumDiffs(fileSER, fileMPI)
    if n == '0':
#    rc = subprocess.check_output([diffreportExe, fileSER, fileMPI])
#    if rc == '':
        runMPI.results[6] = runMPI.successMark
    else:
        runMPI.results[6] = '{: ^5}'.format(n)
#        runMPI.results[6] = runMPI.failMark
        logger.warning('---NPE reproducibility failed')

        
