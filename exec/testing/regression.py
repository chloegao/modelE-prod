#-------------------------------------------------------------------------------
# This script verifies that MPI and Serial builds produce identical results
# for a specified set of rundecks.
# Usage:
#    From the decks subdirectory issue the command:
#      ../exec/testing/regression.py  <runsource1> [<rundeck2> ...]
#
#    Requires a configuration file name <rundeckName>.cfg
#
#
# ENV Options:
#    * If the environment variable DEBUG is set, then the script will
#      display all commands, but not actually execute them.
#
#-------------------------------------------------------------------------------

import sys
import os
import shlex
import subprocess
import logging
import ConfigParser

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
        self.results =  [runShort, comp, mode, '---', '---', '---', '---']

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
    global modeleDir

    if os.environ.has_key('DECKSDIR'):
        decksdir = os.environ['DECKSDIR']
    else:
        decksdir = './'

    subprocess.call(['pwd'])
    configfile = decksdir + '/' + rundeck + '.cfg'
    if os.path.isfile(configfile):
        config = ConfigParser.RawConfigParser()
        config.read(configfile)
    else:
        print 'Using default config options.'
        config = ConfigParser.RawConfigParser(\
            {'rundeck'   :'nonProduction_E_AR5_C12', \
             'modelerc'  :'~/.modelErc', \
             'compiler'  :'gfortran',    \
             'modes'     :'serial,mpi',  \
             'nplist'    :'1,4',         \
             'compflags' :'default',     \
             'branch':'master',      \
             'basedir'   :'.',           \
             'updatebase':'no',          \
             'resultsdir':'.',           \
             'modeledir' :'.'            \
            })

    deck      = config.get('regSettings','rundeck')
    modelerc  = config.get('regSettings','modelerc')
    compiler  = config.get('regSettings','compiler')
    modes     = config.get('regSettings','modes')
    modeList = []
    for mode in modes.split(','):
        modeList.append(mode)
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
    modeleDir = config.get('regSettings','modeledir')

    baseDir =  baseDir + '/' + compiler + '/' + branch
    os.environ['MODELERC'] = modelerc

    if os.environ.has_key('DEBUG'):
        debug = os.environ['DEBUG']
    else:
        debug = False

#-------------------------------------------------------------------------------
def setuplogging():
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

#-------------------------------------------------------------------------------
# Define a subprocess call where we can capture/evaluate stderr output
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

#-------------------------------------------------------------------------------
# A subprocess call that, upon failure rc<>=0, raises an exception
# Subprocess call with shell argument:
# expands environment variables and file globs 
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
        logger.debug('status: ' + str(status))
        if (status != 0):
            raise Exception('unix', commandString)

#-------------------------------------------------------------------------------
# !!! This is the function we should be using, rather than SysCmd !!!
# A subprocess call that, upon failure rc<>=0, raises an exception
# Subprocess call with shell argument:
# expands environment variables and file globs 
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

#-------------------------------------------------------------------------------
# Retrun a checkpoint file name with various identifiers
def checkpointName(exp, duration, npes):
    if exp.mode == 'serial':
        return exp.run  + '.' + duration
    else:
        return exp.run + '.' + duration + '.np=' + str(npes)

#-------------------------------------------------------------------------------
# Build a configuration.
def build(exp):
    logger = logging.getLogger('BUILD   ')
    logger.info(exp.run + ' ' + exp.modeCmd + ' ' + exp.opts)
    try:
        sysCmd('make --quiet clean')
        sysCmd('make rundeck ' + exp.runCmd + ' ' + exp.runSrcCmd)
        sysCmd('make -j4 gcm ' + exp.runCmd + ' ' + exp.modeCmd + ' ' + exp.opts)
        exp.results[3] = 'OK'
    except:
        exp.results[3] = 'bld'
        logger.error('   Failed to build ' + exp.run)
        raise

#-------------------------------------------------------------------------------
# Runs a 1hr simulation
def run1hr(exp, npes=1):
    logger = logging.getLogger('RUN1HR  ')
    logger.info(exp.mode + ', npes=' + str(npes))
    try:
        sysCmd('make setup ' + exp.runCmd + ' ' + exp.modeCmd + ' ' + exp.opts)
        sysCmd('../exec/runE ' + exp.run + ' -np ' + str(npes) + ' -cold-restart')
        sysCmd('cd '+exp.run+'; cp fort.2.nc '+checkpointName(exp, '1hr', npes))
    except:
        exp.results[3] = '1hr'
        message = '   Failed to run 1 hour test for ' + exp.run
        message += ' on ' + str(npes) + ' processors.'
        logger.error(message)
        raise
    
#-------------------------------------------------------------------------------
# Runs a 1dy (25hr) AND a 24hr restart simulation
def run1dy(exp, npes=1):
    logger = logging.getLogger('RUN1DY  ')
    expectedRC = 13; # modelE convention
    restart = './'+exp.run
    if exp.mode == 'mpi':
        restart += ' -np ' + str(npes)
    
    logger.info(exp.mode + ', npes=' + str(npes))
    try:
        sysCmd('../exec/editRundeck.sh ' + exp.run + ' 48 2 1')
        sysCmd('make setup ' + exp.runCmd + ' ' + exp.modeCmd + ' ' + exp.opts)
        sysCmd('../exec/runE ' + exp.run + ' -np ' + str(npes) + ' -cold-restart')
        sysCmd('cd '+exp.run+'; cp fort.1.nc '+checkpointName(exp, '1dy', npes))
        sysCmd('cd ' + exp.run + '; cp fort.2.nc fort.1.nc')
        sysCmd('cd ' + exp.run + '; rm -f run_status')
# Need to investigate why the following causes a NameError exception
# Looks like there is an issue with variable/function/class names in SysCmd
        sysCmd('cd ' + exp.run + '; ' + restart + '; test `head -1 run_status` -eq ' + str(expectedRC))
        sysCmd('cd '+exp.run+';cp fort.2.nc '+checkpointName(exp, 'restart', npes))
    except:
        exp.results[3] = '1dy'
        message = '   Failed to run 1 day test for ' + exp.run
        message += ' on ' + str(npes) + ' processors.'
        logger.error(message)
        raise
    
#-------------------------------------------------------------------------------
# Compare model results with baseline
def compareBase(exp, duration, npes=1):
    logger = logging.getLogger('COMPBAS ')
    prefix = exp.run + '/'
    file1 = prefix + checkpointName(exp, duration, npes)
    file2 = baseDir + '/' + checkpointName(exp, duration, npes)
    logger.info(file1 + ' ' + file2)
    cmp ='diffreport.x ' + file1 + ' ' + file2
    rc = sysCall(cmp)
    if rc == '':
        exp.results[4] = 'OK'
    else:
        logger.debug('   BASE and ' + file1 + ' differ')
        exp.results[4] = 'bas'

#-------------------------------------------------------------------------------
# Compare SERIAL vs MPI
def compareNPE(runA, runB, duration, npes):
    logger = logging.getLogger('COMPNPE ')
    prefix1 = runA.run + '/'
    prefix2 = runB.run + '/'
    file1 = prefix1 + checkpointName(runA, duration, npes)
    file2 = prefix2 + checkpointName(runB, duration, npes)
    logger.info(file1 + ' ' + file2)
    cmp ='diffreport.x ' + file1 + ' ' + file2
    rc = sysCall(cmp)
    if rc == '':
        runB.results[6] = 'OK'
    else:
        logger.debug('   Files ' + file1 + ' and ' + file2 + ' differ')
        runB.results[6] = 'npe'
            
#-------------------------------------------------------------------------------
# Compare full-run (25hr) vs restart run
def compareRestart(exp, npes=1):
    logger = logging.getLogger('COMPRST ')
    numLinesExpected = '0'
#    numLinesExpected = '12' # E4TcadC12
    prefix = exp.run + '/'
    file1 = prefix + checkpointName(exp, '1dy', npes)
    file2 = prefix + checkpointName(exp, 'restart', npes)
    logger.info(file1 + ' ' + file2)
    cmp ='diffreport.x ' + file1 + ' ' + file2
    rc = sysCall(cmp)
    if rc == '':
        exp.results[5] = 'OK'
    else:
        logger.debug('   Files ' + file1 + ' and ' + file2 + ' differ')
        exp.results[5] = 'rst'
        
#-------------------------------------------------------------------------------
# MAIN PROGRAM
if __name__ == '__main__':
    global resultTemplateS
    global resultTemplateM
    
    for rundeck in sys.argv[1:]:

        readConfig(rundeck)
        setuplogging()
        os.chdir(modeleDir)
        logger = logging.getLogger('MAIN    ')
        diffFile = resdir + '/' + rundeck + '.diff'
        fileH = open(diffFile, 'w')
        logger.info('Testing ' + rundeck)

        exps = []
        for mode in modeList:
            exps.append(configuration(rundeck, mode, compiler))
            
        try:
            build(exps[0])
            run1hr(exps[0])
            run1dy(exps[0])
            # Done with SERIAL runs - compare runs with baseline
            compareBase(exps[0], '1hr')
            compareBase(exps[0], '1dy')
            # And compare SERIAL checkpoint-restart 
            compareRestart(exps[0])

            try:
                build(exps[1])
                for npes in npList:
                    run1hr(exps[1], npes=npes)
                    # Compare 1hr run against serial
                    compareNPE(exps[0], exps[1], '1hr', npes)
                    # Compare runs with baseline
                    compareBase(exps[1], '1hr', npes=npes)
            except:
                logger.error('  ... abandoning 1hr mpi test.')

            try:
                build(exps[1])
                for npes in npList:
                    run1dy(exps[1], npes=npes)
                    # Compare 1hr run against serial
                    compareNPE(exps[0], exps[1], '1dy', npes)
                    # Compare runs with baseline
                    compareBase(exps[1], '1dy', npes=npes)
            except:
                logger.error('  ... abandoning 1dy mpi test.')

            logger.info(rundeck + ' Testing complete.')

        except:
            logger.error(rundeck + ' verification FAILED')

        fileH.write('%20s' % (exps[0].results[0]))
        fileH.write('%10s' % (exps[0].results[1]))
        fileH.write('%8s' % (exps[0].results[2]))
        for s in exps[0].results[3:]:
            fileH.write('%6s' % (s))
        fileH.write('\n')
        fileH.write('%20s' % (exps[1].results[0]))
        fileH.write('%10s' % (exps[1].results[1]))
        fileH.write('%8s' % (exps[1].results[2]))
        for s in exps[1].results[3:]:
            fileH.write('%6s' % (s))
        fileH.write('\n')
        fileH.close()

    logger.info('Done')
            
                    
                    
