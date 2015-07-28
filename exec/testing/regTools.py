# This module contains tools to help setup the modelE regression tests
import string
import ConfigParser
import re
import os
import sys
import errno
import shutil
import subprocess
import glob
import logging
import time
import regUtils
from regRuns import *

logger = logging.getLogger('tools')

#-------------------------------------------------------------------------------
# Setup modelE testing environment:
def setupEnv(config, compconfig):
    logger.info('Setup testing environment')
    userconfig = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    branch =  userconfig['repobranch']
    resultsDir = userconfig['scratchdir'] + '/results/' + userconfig['repobranch']
    scratchDir = userconfig['scratchdir'] + '/scratch/' + userconfig['repobranch']

    if userconfig['cleanscratch'] == 'yes':
        if not os.path.exists(resultsDir):
            regUtils.mkdir_p(resultsDir)    
            regUtils.mkdir_p(scratchDir)
        else:
            regUtils.cleanDir(scratchDir)
            regUtils.cleanDir(resultsDir)

        setupModelEenv(config, compconfig)
        gitCloneRepository(config)


#-------------------------------------------------------------------------------
# Clone the model from the user-specified git repository
def gitCloneRepository(config):
    logger.info('Clone user-specified git repository')
    userconfig = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    scratch = userconfig['scratchdir']
    repo = userconfig['repository']
    branch =  userconfig['repobranch']
    clone = scratch + '/scratch/' + branch + '/' + branch
    logger.debug('Clone repository %s',clone)

    cwd = os.getcwd()
    logger.debug('Cloning %s into %s', repo, clone)
    cmd = 'git clone -b ' + branch + ' ' + repo + ' ' + clone \
        + '> /dev/null 2>&1'
    subprocess.check_call(cmd, shell=True)
    os.chdir(cwd)

#-------------------------------------------------------------------------------
# ModelE specific setup
def setupModelEenv(config, compconfig):
    userconfig =regUtils. ConfigSectionMap(config, 'USERCONFIG')
    branch =  userconfig['repobranch']
    makesystem =  userconfig['makesystem']
    resultsDir = userconfig['scratchdir'] + '/results/' + branch
    scratchDir = userconfig['scratchdir'] + '/scratch/' + branch

# the following directories are modelE specific:
    regUtils.mkdir_p(scratchDir+'/decks_repository')
    regUtils.mkdir_p(scratchDir+'/cmrun')
    regUtils.mkdir_p(scratchDir+'/exec')
    regUtils.mkdir_p(scratchDir+'/savedisk')

# We need to get a list of compilers...
    compilers = regUtils.getCompilers(compconfig)
    libsconfig = regUtils.ConfigSectionMap(compconfig, 'COMPCONFIG')

# ... to create modelErc file(s)
    for comp in compilers:
       if not os.path.exists(scratchDir + comp):
          regUtils.mkdir_p(resultsDir + '/' + comp)
          regUtils.mkdir_p(scratchDir + '/' + comp)
       writeModelErc(libsconfig, scratchDir, comp)

#-------------------------------------------------------------------------------
# Write a compiler-specific modelErc file
def writeModelErc(cfg, scratchDir, compiler):

   s = string.Template('\
   DECKS_REPOSITORY=$scr/decks_repository\n\
   CMRUNDIR=$scr/cmrun\n\
   EXECDIR=$scr/exec\n\
   SAVEDISK=$scr/savedisk\n\
   GCMSEARCHPATH=$datadir\n\
   COMPILER=$cm\n\
   MPIDISTR=$mn\n\
   MPIDIR=$md\n\
   NETCDFHOME=$nd\n\
   PNETCDFHOME=$pd\n\
   BASELIBDIR5=$bd\n\
   BUILD_OUT_OF_SOURCE=NO\n\
   OVERWRITE=YES\n\
   OUTPUT_TO_FILES=NO\n\
   VERBOSE_OUTPUT=NO')

   if compiler == 'gfortran':
      modelErc = s.substitute(cm=compiler,\
                              scr=scratchDir,\
                              datadir=cfg['modeldatadir'],\
                              mn=cfg['gccmpi'],\
                              md=cfg['gccmpidir'],\
                              nd=cfg['gccnetcdf'],\
                              pd=cfg['gccpnetcdf'],\
                              bd=cfg['gccesmf'])
   elif compiler == 'intel':
      modelErc = s.substitute(cm=compiler,\
                              scr=scratchDir,\
                              datadir=cfg['modeldatadir'],\
                              mn=cfg['intelmpi'],\
                              md=cfg['intelmpidir'],\
                              nd=cfg['intelnetcdf'],\
                              pd=cfg['intelpnetcdf'],\
                              bd=cfg['intelesmf'])
   elif compiler == 'nag':
      modelErc = s.substitute(cm=compiler,\
                              scr=scratchDir,\
                              datadir=cfg['modeldatadir'],\
                              mn=cfg['nagmpi'],\
                              md=cfg['nagmpidir'],\
                              nd=cfg['nagnetcdf'],\
                              pd=cfg['nagpnetcdf'],\
                              bd=cfg['nagesmf'])
   else:
      modelErc = s.substitute(cm=compiler,\
                              scr=scratchDir,\
                              datadir=cfg['modeldatadir'],\
                              mn=cfg['gccmpi'],\
                              md=cfg['gccmpidir'],\
                              nd=cfg['gccnetcdf'],\
                              pd=cfg['gccpnetcdf'],\
                              bd=cfg['gccesmf'])

   rcfile = open(scratchDir + '/' + compiler + '/modelErc.' + compiler, "w")
   rcfile.write(modelErc)
   rcfile.write('\n')
   rcfile.close()
   logger.debug('Created modelErc file for compiler %s', compiler)

#-------------------------------------------------------------------------------
def setupCloneTasks(config, compconfig, decklist):
    userconfig =regUtils.ConfigSectionMap(config, 'USERCONFIG')
    compilers = regUtils.getCompilers(compconfig)

    cloneTasks = []
    for deck in decklist:
      # since nonProduction rundeck names can be quite long, extract the
      # nonProduction_ part...
        dName = deck.name
        if re.search('nonProduction', deck.name):
            start = deck.name.find('nonProduction') + 14
            dName = deck.name[start:]
         
        for comp in deck.getOpt('compilers').split(','):
            for mode in deck.getOpt('modes').split(','):
                cmode = '.' + mode
                if comp in compilers:
                   commandString = regUtils.gitCloneCommand(config, dName, comp, cmode)
                   cloneTasks.append(commandString)
                else:
                   logger.error(comp+' is not defined in COMPCONFIG')
            
    for t in cloneTasks:
        logger.debug('CLONE TASK %s', t)
    return cloneTasks

#-------------------------------------------------------------------------------
# For out-of source builds, create directory for each rundeck/compiler/mode combo
def setupRuns(config, compconfig, decklist):
    userconfig =regUtils.ConfigSectionMap(config, 'USERCONFIG')
    compilers = regUtils.getCompilers(compconfig)

    tasks = []
    for deck in decklist:
      # since nonProduction rundeck names can be quite long, extract the
      # nonProduction_ part...
        dName = deck.name
        if re.search('nonProduction', deck.name):
            start = deck.name.find('nonProduction') + 14
            dName = deck.name[start:]

        for comp in deck.getOpt('compilers').split(','):
            for mode in deck.getOpt('modes').split(','):
                cmode = '.' + mode
                if comp in compilers:
                    regUtils.mkdirCommand(config, dName, comp, cmode)
                else:
                    logger.error(comp+' is not defined in COMPCONFIG')

#-------------------------------------------------------------------------------
# Return a command to submit/execute a [batch] job
def setupScriptTasks(config, compconfig, decklist):
    logger.info('Prepare and execute tasks...')
    compilers = regUtils.getCompilers(compconfig)

    scriptTasks = []
    for deck in decklist:
        for comp in deck.getOpt('compilers').split(','):
            for mode in deck.getOpt('modes').split(','):
                if comp in compilers:
                    commandString = \
                        createScriptTask(config, compconfig, deck, comp, mode)
                    scriptTasks.append(commandString)
                else:
                    logger.error(comp+' is not defined in COMPCONFIG')

    for t in scriptTasks:
        logger.debug('SCRIPT TASK %s', t)
    return scriptTasks

#-------------------------------------------------------------------------------
# Creates script to be submitted to batch system OR to be executed interactively
# Batch system is assumed to be the one on NCCS-DISCOVER machines
def createScriptTask(config, compconfig, deck, comp, mode):
    userconfig  = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    modules    = userconfig['modules']
    useBatch   = userconfig['usebatch']
    branch     = userconfig['repobranch']
    scriptsDir = userconfig['scriptsdir'] + '/exec/testing/'
    useMods    = userconfig['modules']
    resultsDir = userconfig['scratchdir'] + '/results/' + \
                 branch + '/' + comp
    scratchDir = userconfig['scratchdir'] + '/scratch/' + \
                 branch + '/' + comp
    sponsorID  = userconfig['sponsorid']

    deckName = deck.name
    jobName = deckName
    if re.search('nonProduction', deckName):
       start = deckName.find('nonProduction') + 14
       jobName = deckName[start:]
    filename = resultsDir + '/' + jobName + '.' + mode + '.bash'
    fileHandle = open ( filename, 'w' ) 

    if useBatch == 'yes':
        # If we are just compiling this rundeck
        if deck.getOpt('verification') == 'compileOnly':
            cores = 4
            walltime = '00:30:00'

        # customRun is a 2-month run
        elif deck.getOpt('verification') == 'customRun':
            cores = 88                
            if re.search('tomas', deckName):
                walltime = '8:00:00'
            elif re.search('amp', deckName):
                walltime = '2:00:00'
            elif re.search('cadi', deckName):
                walltime = '2:00:00'
            elif re.search('obio', deckName):
                walltime = '1:00:00'
            elif re.search('C12', deckName):
                cores = 22
                walltime = '0:30:00'
            elif re.search('M20', deckName):
                cores = 44
                walltime = '0:30:00'
            else:
                cores = 44
                walltime = '1:00:00'
            
        # regular runs (1hr and/or restart)
        else:               
            if 'mpi' in mode: 
                cores = 8
                walltime = '01:00:00'
                if re.search('tomas', deckName):
                    cores = 88                
                elif re.search('amp', deckName):
                    cores = 44        
                elif re.search('E_AR5_V2', deckName):
                    if re.search('NINT', deckName):
                        cores = 8
                    else: # CADI and CAMP
                        cores = 44

            # serial
            else:
                cores = 1
                walltime = '00:30:00'
                if re.search('obio', deckName):
                    walltime = '01:30:00'
                elif re.search('cadi', deckName):
                    walltime = '04:00:00'

            # Adjust the walltime for some rundecks
            if re.search('C12', deckName):
                walltime = '00:30:00'
            elif re.search('Mars', deckName):
                walltime = '00:30:00'
            elif re.search('SGP', deckName):
                walltime = '00:10:00'
            elif re.search('M20', deckName):
                walltime = '00:30:00'

        outname = resultsDir + '/' + jobName + '.' + comp + '.out'
        errname = resultsDir + '/' + jobName + '.' + comp + '.err'
        fileHandle.write ('#!/bin/bash' + '\n')
        fileHandle.write ('#SBATCH -J ' + jobName + '\n')
        fileHandle.write ('#SBATCH -o ' + outname + '\n')
        fileHandle.write ('#SBATCH -e ' + errname + '\n')
        fileHandle.write ('#SBATCH --account='  + sponsorID + '\n')
        fileHandle.write ('#SBATCH --time='     + walltime + '\n')
        fileHandle.write ('#SBATCH --ntasks=' + str(cores) + '\n')
        # Use Haswell NODES
        fileHandle.write ('#SBATCH --constraint=hasw' + '\n')
        #if walltime == '00:30:00':
        #    fileHandle.write ('#SBATCH --qos=debug' + '\n')
           
    # Create rest of script used in batch OR interactive jobs:

    # Do we have modules to 'load'?
    if modules == 'yes':
        machine = subprocess.check_output(['uname','-n'])
        if 'borg' in machine or 'discover' in machine or 'dali' in machine:
            fileHandle.write ('. /usr/share/modules/init/bash' + '\n')
            fileHandle.write ('module purge' + '\n')
        # Need the following module on DISCOVER to get python 2.7.x
            fileHandle.write ('module load other/SSSO_Ana-PyD/SApd_2.1.0' + '\n')
        # If not on DISCOVER
        else:
            fileHandle.write ('#!/bin/bash' + '\n')
            # This is not portable...just my MAC so far
            fileHandle.write ('. /opt/local/share/Modules/3.2.10/init/bash' + '\n')
            fileHandle.write ('module purge' + '\n')

        # Using different naming convention for module names
        compvendor = comp
        if comp == 'gfortran':
            compvendor = 'gcc'

        modsconfig = regUtils.ConfigSectionMap(compconfig, 'COMPCONFIG')
        for mod in modsconfig['modulelist'].split(','):
            if re.search(compvendor, mod):
                for mm in modsconfig[mod].split(','):
                    cmd = 'module load ' + mm +'\n'
                    fileHandle.write (cmd)

    makesystem =  userconfig['makesystem']
    if makesystem == 'makeOld':
        decksDir = scratchDir + '/' + jobName +  '.' + mode + '/decks/'
    else:
        decksDir = scratchDir + '/' + jobName +  '.' + mode
        fileHandle.write ('export BUILD_OUT_OF_SOURCE=YES\n')

    fileHandle.write ('export DECKSDIR=' + decksDir + '\n')

    # cd to the working dir and run the script
    fileHandle.write ('cd ' + decksDir + '\n')
    fileHandle.write ('python ' + scriptsDir + '/' + 'regression.py ' + deckName + '\n')
    fileHandle.write (' ' + '\n')
    fileHandle.close()

    modelErc = scratchDir + '/modelErc.' + comp
    createRegConfig(config, deck, modelErc, comp, jobName, mode)

    if useBatch == 'yes':
        commandString = 'sbatch ' + filename
    else:
        commandString = 'chmod +x ' + filename + ';' + filename
               
    return commandString

#-------------------------------------------------------------------------------
# Create a config file for regression.py script. 
# Note: there is one config file for each rundeck/compiler combination
def createRegConfig(config, deck, modelerc, comp, jobName, mode):
    cfg  = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    branch     = cfg['repobranch']
    resultsDir = cfg['scratchdir'] + '/results/' + \
            branch + '/' + comp
    scratch = cfg['scratchdir'] + '/scratch/' + \
            branch + '/' + comp
    makesystem = cfg['makesystem']
    if makesystem == 'makeOld':
        decksDir = scratch + '/' + jobName + '.' + mode + '/decks/'
    else:
        decksDir = scratch + '/' + jobName + '.' + mode

# Regression tests do not run the regression script in standalone mode    
    standalone = 'no'
    regconfig  = ConfigParser.RawConfigParser()
    regconfig.add_section('regSettings')
    regconfig.set('regSettings', 'rundeck', deck.name)
    regconfig.set('regSettings', 'modelerc', modelerc)
    regconfig.set('regSettings', 'compiler', comp)
    regconfig.set('regSettings', 'modes', mode)
    regconfig.set('regSettings', 'standalone', standalone)
    regconfig.set('regSettings', 'verification', deck.getOpt('verification'))
    regconfig.set('regSettings', 'endtime', deck.getOpt('endtime'))
    regconfig.set('regSettings', 'nplist', deck.getOpt('npes'))
    regconfig.set('regSettings', 'buildtype', cfg['buildtype'])
    regconfig.set('regSettings', 'repository', cfg['repository'])
    regconfig.set('regSettings', 'branch', branch)
    regconfig.set('regSettings', 'basedir', cfg['basedir'])
    regconfig.set('regSettings', 'updatebase', cfg['updatebase'])
    regconfig.set('regSettings', 'makesystem', makesystem)
    regconfig.set('regSettings', 'resultsdir', resultsDir)
    regconfig.set('regSettings', 'scratchdir', scratch)
    regconfig.set('regSettings', 'decksdir', decksDir)

    filename = decksDir + '/' + deck.name + '.cfg'
    logger.debug(filename)
    with open(filename, 'w') as configfile:
        regconfig.write(configfile)


#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
def writeDiff(run, fileH):
    fileH.write('%20s' % (run.results[0]))
    fileH.write('%10s' % (run.results[1]))
    fileH.write('%8s'  % (run.results[2]))
    fileH.write('%4s'  % '    ')
    for s in run.results[3:]:
        fileH.write('{: ^5}'.format(s))
        fileH.write('%3s'  % '   ')
    fileH.write('\n')

#-------------------------------------------------------------------------------
# This function performs a verification of the run output produced by the
# regression tests. 
def verifyRuns(config, runSources):
    logger.info('Verifying...')
    setRunUtils()

    userconfig = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    makesystem =  userconfig['makesystem']
    scratchDir = userconfig['scratchdir'] + '/scratch/' + \
        userconfig['repobranch'] + '/' 

    # Loop over each run source in list
    for source in runSources:

        dirName = source.name 
        if re.search('nonProduction', source.name):
            start = source.name.find('nonProduction') + 14
            dirName = source.name[start:]

        # List of rundeck run configurations for each mode
        runs = []
        for mode in source.getOpt('modes').split(','):
           for comp in source.getOpt('compilers').split(','):
              if makesystem == 'makeOld':
                 decksDir = scratchDir+comp+'/'+dirName+'.'+mode+'/decks'
              else:
                 decksDir = scratchDir+comp+'/'+dirName+'.'+mode
              os.chdir(decksDir)
              os.environ['MYCONFIGDIR'] = decksDir
              # Create rundeck object with default or config properties
              rundeck = newRundeck(source.name)
              runs.append(newRun(rundeck, mode))

        logger.info('Verifying ' + rundeck.name + ': ' + rundeck.verification)

        for run in runs:
           if makesystem == 'makeOld':
              decksDir = scratchDir+run.compiler+'/'+dirName+'.'+run.mode+'/decks'
           else:
              decksDir = scratchDir+run.compiler+'/'+dirName+'.'+run.mode
           os.chdir(decksDir)

           # For each rundeck/compiler/mode combination
           # create a diffFile with verification results
           diffFile = rundeck.resultsDir + '/' + run.name + '.diff'
           fileH = open(diffFile, 'w')

           # Did executable build?
           exe = dirName+'.'+run.mode+'.'+run.compiler
           cmd = 'ls '+exe+'_bin/'+(exe+'.exe')
           status = run.sysCmd(cmd, 3, 'b')
           # If not, then go on to next experiment
           if status != 0:
              continue
           else:
               compare(rundeck, run)
               writeDiff(run, fileH)

           fileH.close()

    logger.info(rundeck.name + ' verification complete.')

#-------------------------------------------------------------------------------
# Create a diff report and notify via email
def sendDiffreport(config, compconfig, eTime):
    userconfig  = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    mailto     = userconfig['mailto']
    branch     = userconfig['repobranch']
    resultsDir = userconfig['scratchdir'] + '/results/' + branch
    buildtype  = userconfig['buildtype']
    message    = userconfig['message']
    sortdiff   = userconfig['sortdiff']
    compilers  = regUtils.getCompilers(compconfig)

    diffFile = resultsDir + '/' + 'diffreport.txt'
    fp = open(diffFile, 'w')
    fp.write(message + ' \n')
    fp.write('-'*72+'\n')
    fp.write('Branch: ' + branch)
    fp.write('  --  Build type: ' + buildtype +  '\n')
    fp.write('-'*72+'\n')
    fp.write('%70s\n' % ('    -REPRODUCIBILITY   '))
    fp.write('%20s%10s%8s%8s%8s%8s%8s\n' % \
        ('RUNDECK', 'COMPILER', 'MODE', 'RUN', 'BAS', 'RST', 'NPE'))
    fp.write('-'*72+'\n')


    # Look at run diffs and check for build failures (Fb)
    subprocess.call('find '+resultsDir+' -name \*.diff -exec cat {} \; >' \
                        +resultsDir + '/' + 'alldiffs', shell=True)
    # If failures exist, signal it with a file
    subprocess.call('cat '+resultsDir + '/' + 'alldiffs | ' \
                        + "awk '{print $4}' | grep Fb > " \
                        + resultsDir + '/' + 'compileFail', shell=True)

    if sortdiff == 'yes':
        # sort compiler column
        subprocess.call('cat '+resultsDir + '/' + 'alldiffs | sort -k 2,2 >' \
            +resultsDir + '/' + 'sorteddiffs', shell=True)
        with open(resultsDir + '/' + 'sorteddiffs','r') as inf:
            fp.write(inf.read())
    else:
        for comp in compilers:
            diffs = glob.glob(resultsDir + '/' + comp + '/*.diff')
            for f in diffs:
                with open(f,'r') as inf:
                    fp.write(inf.read())

    fp.write('-'*72+'\n')
    hhmmss = time.strftime('%H:%M:%S', time.gmtime(eTime))
    fp.write('Time taken = %s \n' %(hhmmss))
    fp.write('-'*72+'\n')
    fp.write('Legend:\n')
    fp.write('-'*7+'\n')
    fp.write('+   : success\n')
    fp.write('NUM : number of reproducibility differences\n')
    fp.write('Fb  : build failure\n')
    fp.write('F1  : 1hr run-time failure\n')
    fp.write('Fr  : restart run-time failure\n')
    fp.write('F*  : expected failure\n')
    fp.write('-   : not available\n')
    fp.write('Notes:\n')
    fp.write('-'*6+'\n')
    compconfig = regUtils.ConfigSectionMap(compconfig, 'COMPCONFIG')
    compVers =  compconfig['compiler_versions'].split(",")
    i=0
    for comp in compilers:
        fp.write(comp+' compiler version: '+compVers[i]+'\n')
        i+=1
    fp.write('Results in: ' + resultsDir +  '\n')
    fp.close()

    subject = '"[modelE-regression]" '
    cmd = '/usr/bin/mail -s ' + subject + mailto + ' < ' + diffFile
    subprocess.check_call(cmd, shell=True)

    # If there were no compilation failures then signal it with a file
    # Longer regression tests will run only if this file exists
    if os.path.exists(resultsDir + '/' + 'compileFail'):
        if os.stat(resultsDir + '/' + 'compileFail').st_size == 0:
            subprocess.call('touch ' + os.environ['HOME'] \
                            + '/.CompileModelEOK', shell=True)


