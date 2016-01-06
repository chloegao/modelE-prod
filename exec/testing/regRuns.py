# This module contains functions that help setup and compare modelE runs
import sys
import os
import shutil
import shlex
import time
import subprocess
import logging
import ConfigParser
import regUtils as utils
import regTest as t

""" 
  This class assigns settings used to test a given rundeck 
"""
class newRundeck(t.regTest):

    def __init__(self, sourceName='nonProduction_E_AR5_C12'):
        # Overide name
        if sourceName:
            self.name      = sourceName
        else:
            self.name      = 'nonProduction_E_AR5_C12'
        super(newRundeck, self).__init__(sourceName)

        # modelE rundeck defaults    
        self.standalone    = 'yes'
        self.branch        = ''
        self.updateBase    = 'no'
        self.baseDir       = '.'
        self.repository    = '..'
        self.decksDir      = '.'
        self.makesystem    = 'makeOld'
        self.compiler      = '' 
        self.savedisk      = '.'
        # Now override defaults with values specified in config file
        getConfigSettings(self)

    #  Setup a logging object for each rundeck. Note that the output file gets
    #  logging DEBUG while STDOUT only gets logging INFO in order to minimize 
    #  verbosity.
    def setLogging(self):
        # Note filemode is 'append' because we want one logger per rundeck
        logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=self.resultsDir+'/'+self.name
                             +'-regression.log',
                    filemode='a')
        stdoutLog = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(name)s : %(message)s')
        stdoutLog.setFormatter(formatter)
        stdoutLog.setLevel(logging.INFO)
        logger = logging.getLogger()
        logger.addHandler(stdoutLog)
        

"""
  This class defines, among other things, a rundecks's run time options
  passed to the the makefile. Note that a "run" is a "rundeck".
"""
class newRun(newRundeck):

    def __init__(self, rundeck, mode):
        super(newRun, self).__init__(rundeck.name)
        self.runsrc = rundeck.name
        # shortName is introduced to remove the nonProduction part of a
        # rundeck name
        shortName = rundeck.name
        if 'nonProduction' in rundeck.name:
            shortName = rundeck.name[14:]
        # A RUN name is composed of rundeck name, mode and compiler
        self.shortName = shortName
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
        if rundeck.buildType == 'debug':
            flags='"-O0 -g"'
            self.xflags += 'EXTRA_FFLAGS='+flags
        elif rundeck.buildType == 'traps':
            self.xflags += 'COMPILE_WITH_TRAPS=YES'
        # Run results are stored in an array
        self.results = [shortName, rundeck.compiler, mode, 
                        '  -  ', '  -  ', '  -  ', '  -  ', '  -  ']
        self.successMark   = '+'
        self.failMark      = 'F'
        self.createMark    = 'C'
        self.naMark        = '-'
        # In case we want to debug this run:
        if os.environ.has_key('DEBUG'):
            self.debug = os.environ['DEBUG']
        else:
            self.debug = False

    # System call that records result of subprocess call
    def sysCmd(self, commandString, resultIndex, stageID, makeLog='yes'):

        logger = logging.getLogger('SYSTEM  ')
        status = 1
        grepResult = 'OK'
        if self.debug:
            logger.info(commandString)
        else:
            logger.debug(commandString)
            if makeLog == 'yes':
                makeLog = self.resultsDir + '/'  + self.name + '-make.log'
                with open(makeLog,'a') as f:
                    status = subprocess.call(commandString, \
                                             stdout=f, stderr=f, shell=True)
            else: # we do not log the unit tests, we just capture their output
                proc = subprocess.Popen(commandString,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        )
                proc.wait()
                status = proc.returncode
                if status != 0:
                    # Are there unit test "Failues"
                    grep = subprocess.Popen(shlex.split('grep Failures'),
                                            stdin=proc.stdout,
                                                        stdout=subprocess.PIPE,
                    )
                    # If so, how many?
                    cut = subprocess.Popen(shlex.split('cut -f 2 -d,'),
                                            stdin=grep.stdout,
                                            stdout=subprocess.PIPE,
                    )
                    awk =subprocess.Popen(shlex.split("awk '{print $2}'"),
                                            stdin=cut.stdout,
                                            stdout=subprocess.PIPE,
                    ) 
                    grep.stdout.close()
                    out,err = awk.communicate()
                    grepResult = out.strip()
                    # If grep was empty then it was a build error:
                    if grepResult == '':
                        grepResult = 'Fb'
                    # Write result to small file for diffreport
                    f = open(".unit", "w+")
                    f.write(grepResult)
                    f.close()
                    
            logger.debug('Return code: ' + str(status))
            logger.debug('Unit tests result: ' + grepResult)

            if (status == 0):
                self.results[resultIndex] = self.successMark
            else:
                if grepResult == 'OK':
                    logger.error(commandString+': FAILED')
                    self.results[resultIndex] = self.failMark+stageID
                else:
                    self.results[resultIndex] = grepResult
                    
 
""" 
  Set configuration for a rundeck
"""
def getConfigSettings(rundeck):

    # If MYCONFIGDIR is defined get configuration file from there
    if os.environ.has_key('MYCONFIGDIR'):
        myConfigDir = os.environ['MYCONFIGDIR']
    else: # it is in the current directory
        myConfigDir = '.'
    configfile = myConfigDir + '/' + rundeck.name + '.cfg'
    
    if os.path.isfile(configfile):
        config = ConfigParser.ConfigParser()
        config.read(configfile)
        # rundeck settings
        rundeck.name      = config.get('regSettings', 'rundeck')
        rundeck.modelerc  = config.get('regSettings', 'modelerc')
        rundeck.compiler  = config.get('regSettings', 'compiler')
        # rundeck.modes is a list - config setting is a string, so...
        rundeck.modes     = config.get('regSettings', 'modes').split(',')
        # npes is a string in the config file but a list of ints 
        # in the data structure...
        npes = config.get('regSettings', 'npes')
        nint = []
        for n in npes.split():
            nint.append(int(n))
        rundeck.npes = nint
        rundeck.standalone = config.get('regSettings', 'standalone')
        rundeck.verification = config.get('regSettings', 'verification')
        rundeck.unitTest = config.get('regSettings', 'unittest')
        rundeck.endTime   = int(config.get('regSettings', 'endtime'))
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
    else:
        print ' *** Configuration file does not exist *** ' + configfile
        sys.exit(1)

    # Extract SAVEDISK from modelErc file - needed in compareNPE            
    cmd = "cat "+rundeck.modelerc+"| grep SAVEDISK"+"| awk -F= '{print $2}'"
    out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,\
                               stderr=subprocess.STDOUT)
    rundeck.savedisk = out.communicate()[0].rstrip()

    for mode in rundeck.modes:
        if mode != 'serial' and mode != 'mpi':
            print ' *** Incorrect mode *** ' + mode

    # In order to avoid errors in compareBase() when running with
    # default options:
    if rundeck.baseDir != '.':
        rundeck.baseDir =  rundeck.baseDir + '/' + rundeck.branch + '/' \
            + rundeck.compiler

    
"""
  Build model
"""
def build(run):
    logger = logging.getLogger('BUILD   ')
    logger.info(run.name + ' ' + run.modeCmd + ' ' + run.xflags)

    if run.makesystem == 'makeOld':
        try:
            cmd = 'make rundeck '+run.runCmd+' '+run.runSrcCmd
            run.sysCmd(cmd, 3, 'b')
        except Exception, e:
            logger.exception(str(e))
            return 1
        try:
            cmd = 'make -j gcm '+run.runCmd+' '+run.modeCmd+' '+run.xflags
            run.sysCmd(cmd, 3, 'b')
        except Exception, e:
            logger.exception(str(e))
            return 1
    else:
        try:
            cmd = 'make rundeck '+run.runCmd+' '+run.runSrcCmd
            run.sysCmd(cmd, 3, 'b')
        except Exception, e:
            logger.exception(str(e))
            return 1
        try:
            os.chdir(run.decksDir+'/../'+run.name)
            cmd = '../configme/discover_' + run.compiler + ' ../decks/' \
                  + run.name+'.R ..'
            run.sysCmd(cmd, 3, 'b')
        except Exception, e:
            logger.exception(str(e))
            return 1
        try:
            cmd = 'make -j'
            run.sysCmd(cmd, 3, 'b')
        except Exception, e:
            logger.exception(str(e))
            return 1

    if run.unitTest == 'yes':
        logger.info('Run unit tests...')
        try:
            cmd = 'make tests '+run.runCmd+' '+run.modeCmd
            run.sysCmd(cmd, 4, 't', 'no')
        except Exception, e:
            logger.exception(str(e))
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
        try:
            cmd =  'make setup '+run.runCmd+' '+run.modeCmd+' '+run.xflags
            run.sysCmd(cmd, 3, '1')
        except Exception, e:
            logger.exception(str(e))
            return 1
    else:
        try:
            os.chdir(run.decksDir+'/..')
            cmd = 'python python/rune/make_rundir.py decks/' +run.name \
                  + '.R ' + run.name + '-scratch'
            run.sysCmd(cmd, 3, '1')
        except Exception, e:
            logger.exception(str(e))
            return 1

    if run.makesystem == 'makeOld':
        try:
            rune = run.repository+'/exec/runE '
            cmd = rune+run.name+' -np '+str(npes)+' -cold-restart'  
            run.sysCmd(cmd, 3, '1')
        except Exception, e:
            logger.exception(str(e))
            cmd = 'cd '+run.name+'; touch '+run.name+'.1hr.FAILED'
            rc = run.sysCmd(cmd, 3, '1')
            return 1

        try:
            cmd = 'cd '+run.name+ '; test `head -1 run_status` -eq ' + mErc
            run.sysCmd(cmd, 3, '1')
        except Exception, e:
            logger.exception(str(e))
            return 1

        try:
            cmd = 'cd ' + run.name + '; cp fort.2.nc ' \
                  + utils.checkpointName(run.name, run.mode, '1hr', npes)
            run.sysCmd(cmd, 3, '1')
        except Exception, e:
            logger.exception(str(e))
            return 1
    else:
        os.chdir(run.decksDir+'/../'+run.name+'-scratch')
        rune = ''
        if run.mode == 'mpi':
            rune = 'mpirun -np '+str(npes)
        try:
            cmd = rune+' ../'+run.name+'/model/modelexe -i I -cold-restart'  
            run.sysCmd(cmd, 3, '1')
        except Exception, e:
            logger.exception(str(e))
            return 1
        try:
            cmd = 'test `head -1 run_status` -eq ' + mErc
            run.sysCmd(cmd, 3, '1')
        except Exception, e:
            logger.exception(str(e))
            return 1
        try:
            cmd = 'cp fort.2.nc ' \
                  + utils.checkpointName(run.name, run.mode, '1hr', npes)
            run.sysCmd(cmd, 3, '1')
        except Exception, e:
            logger.exception(str(e))
            return 1

    logger.info(run.name + ' is DONE')
    return 0

"""
  Run up to ENDTIME hrs with checkpoint at ENDTIME-1 hrs
"""
def runRestart(run, npes=1, endTime=25):
    logger = logging.getLogger('RUNRST  ')
    mErc = '13'; # modelE convention for successful runs
    restart = './'+run.name
    if run.mode == 'mpi':
        restart += ' -np ' + str(npes)

    logger.info(run.name + ', ' + run.mode + ', npes=' + str(npes) + \
                    ', endTime=' + str(endTime))

    checkPt = endTime - 1
    ndisk = checkPt * 2
    if endTime > 24:
        newTime = ' ' + str(ndisk) + ' 2 1'
    else:
        newTime = ' ' +  str(ndisk) + ' 1 ' + str(endTime)

    try:
        cmd = run.repository+'/exec/editRundeck.sh '+run.name+newTime 
        run.sysCmd(cmd, 3, 'e')
    except Exception, e:
        logger.exception(str(e))
        return 1

    # Run make setup again to update the edited rundeck file
    if run.makesystem == 'makeOld':
        try:
            cmd =  'make -j setup '+run.runCmd+' '+run.modeCmd+' '+run.xflags
            run.sysCmd(cmd, 3, 's')
        except Exception, e:
            logger.exception(str(e))
            return 1
    else:
        try:
            run.sysCmd('make setup ', 3, 's')
        except Exception, e:
            logger.exception(str(e))
            return 1

    try:
        rune = run.repository+'/exec/runE '
        cmd = rune+run.name+' -np '+str(npes)+' -cold-restart'  
        run.sysCmd(cmd, 3, 'r')
    except Exception, e:
        logger.exception(str(e))
        cmd = 'cd '+run.name+'; touch '+run.name+'.'+str(endTime)+'.FAILED'
        run.sysCmd(cmd, 3, 'a')
        return 1

    try:
        cmd = 'cd ' + run.name + '; cp fort.1.nc ' \
        +utils.checkpointName(run.name, run.mode, str(endTime)+'hr', npes)
        run.sysCmd(cmd, 3, 'b')
    except Exception, e:
        logger.exception(str(e))
        return 1

    try:
        cmd = 'cd ' + run.name + '; cp fort.2.nc fort.1.nc; rm -f run_status'
        run.sysCmd(cmd, 3, 'c')
    except Exception, e:
        logger.exception(str(e))
        return 1 

    try:
        cmd = 'cd ' + run.name + '; ' + restart \
            + '; test `head -1 run_status` -eq ' + mErc
        run.sysCmd(cmd, 3, 'R')
    except Exception, e:
        logger.exception(str(e))
        cmd = 'cd '+run.name+'; touch '+run.name+'.restart.FAILED'
        run.sysCmd(cmd, 3, 'd')
        return 1

    try:
        cmd = 'cd ' + run.name + ';cp fort.2.nc ' \
            + utils.checkpointName(run.name, run.mode, 'restart', npes)
        run.sysCmd(cmd, 3, 're')
    except Exception,e:
        logger.exception(str(e))
        return 1

# Reset rundeck settings for next MPI run
    if len(run.npes) > 1:
        if run.makesystem == 'makeOld':
            try:
                makecmd = 'make rundeck '+run.runCmd+' '+run.runSrcCmd
                run.sysCmd(makecmd, 3, 'f')
            except Exception, e:
                logger.exception(str(e))
                return 1
        else:
            try:
                run.sysCmd('make rundeck RUN=' + run.runCmd 
                              + ' RUNSRC='+ run.runSrcCmd, 3, 'f')
            except Exception, e:
                logger.exception(str(e))
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

    try:
        run.sysCmd(run.repository + '/exec/editRundeck.sh ' + run.name + \
                        newTime, 3, 'l')
    except Exception, e:
        logger.exception(str(e))
        return 1

    if run.makesystem == 'makeOld':
        try:
            run.sysCmd('make -j setup ' + run.runCmd + ' ' + run.modeCmd + \
                            ' ' + run.xflags, 3, 'l')
        except Exception, e:
            logger.exception(str(e))
            return 1
    else:
        try:
            run.sysCmd('make -j setup', 3, 'l')
        except Exception, e:
            logger.exception(str(e))
            return 1

    try:
        run.sysCmd(run.repository + '/exec/runE ' + run.name + ' -np ' + \
                            str(npes) + ' -cold-restart', 3, 'l')
    except Exception, e:
        logger.exception(str(e))
        return 1

    logger.info(run.name + ' is DONE')
    return 0

