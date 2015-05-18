# This module contains tools to help setup the regression tests
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
from regTest import *
from regFuncs import *

logger = logging.getLogger('modelE')

#-------------------------------------------------------------------------------
# From config file create/return a list of model run configurations
def getModelConfigurations(config):
    sections = config.sections()
    modelConfig = {}

# For convenience divide the sections in the configuration file into
# two types: CONFIG and RUNDECKS. The former have CONFIG in their names.
# Thus, if a section name does NOT have CONFIG in its name then it is
# a rundeck configuration:
    for sect in sections:
        match = not re.search("CONFIG",sect)
        # get all rundeck sections from config file
        if (match):
            modelConfig[sect] = regUtils.ConfigSectionMap(config, sect)

# map modelConfig to a more manageable list
    runList = []
    for name,options in modelConfig.items():
       runList.append(regTest(name))
# Each item in runList (each rundeck) also needs user-defined options
    userconfig = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    for d in runList:
       d.setOpts(userconfig, modelConfig)
        
    return runList

#-------------------------------------------------------------------------------
# Setup testing environment:
# 1) Create working directories
# 2) Clone model from git repository and...
# 3) Perform additional model "specific" setup
def setupEnv(config, bpconfig):
   userconfig = regUtils.ConfigSectionMap(config, 'USERCONFIG')
   branch =  userconfig['repobranch']
   resultsDir = userconfig['scratchdir'] + '/regression_results/' + branch
   scratchDir = userconfig['scratchdir'] + '/regression_scratch/' + branch

   if not os.path.exists(resultsDir):
      regUtils.mkdir_p(resultsDir)    
      regUtils.mkdir_p(scratchDir)
   else:
      if userconfig['cleanscratch'] == 'yes':
         regUtils.cleanDir(scratchDir)
         regUtils.cleanDir(resultsDir)

   setupModelEenv(config, bpconfig)
   gitCloneRepository(config)


#-------------------------------------------------------------------------------
# Clone the model from the user-specified git repository
def gitCloneRepository(config):
   userconfig = regUtils.ConfigSectionMap(config, 'USERCONFIG')
   scratch = userconfig['scratchdir']
   repo = userconfig['repository']
   branch =  userconfig['repobranch']
   clone = scratch + '/regression_scratch/' + branch + '/' + branch
   logger.debug('Clone repository %s',clone)

# if we have a valid clone
   cwd = os.getcwd()
   if os.path.isdir(clone):
       os.chdir(clone)
       if os.path.isdir('.git'):
           logger.warning('%s is already a git repository', clone)
##   if os.path.isdir(clone) and git("rev-parse" "--is-inside-work-tree"):
##   if os.path.isdir(clone) and os.path.isdir('.git'):
##      ##os.chdir(clone)
##      logger.warning('%s is already a git repository', clone)
# else clone it
   else:
      logger.debug('Cloning %s into %s', repo, clone)
      cmd = 'git clone -b ' + branch + ' ' + repo + ' ' + clone \
          + '> /dev/null 2>&1'
      subprocess.check_call(cmd, shell=True)

   os.chdir(cwd)

#-------------------------------------------------------------------------------
# ModelE specific setup
def setupModelEenv(config, bpconfig):
   userconfig =regUtils. ConfigSectionMap(config, 'USERCONFIG')
   branch =  userconfig['repobranch']
   resultsDir = userconfig['scratchdir'] + '/regression_results/' + branch
   scratchDir = userconfig['scratchdir'] + '/regression_scratch/' + branch

# the following directories are modelE specific:
   regUtils.mkdir_p(scratchDir+'/decks_repository')
   regUtils.mkdir_p(scratchDir+'/cmrun')
   regUtils.mkdir_p(scratchDir+'/exec')
   regUtils.mkdir_p(scratchDir+'/savedisk')

# We need to get a list of compilers...
   compilers = regUtils.getCompilers(bpconfig)
   libsconfig = regUtils.ConfigSectionMap(bpconfig, 'COMPCONFIG')

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
   OVERWRITE=YES\n\
   OUTPUT_TO_FILES=NO\n\
   VERBOSE_OUTPUT=YES\n\
   COMPILER=$cm\n\
   MPIDISTR=$mn\n\
   MPIDIR=$md\n\
   NETCDFHOME=$nd\n\
   PNETCDFHOME=$pd\n\
   BASELIBDIR5=$bd')

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
def createDiffreport(config, runSources):
    setRunUtils()

    sysconfig  = regUtils.ConfigSectionMap(config, 'SYSCONFIG')
    scratchDir = sysconfig['scratchdir'] + '/regression_scratch/' + \
        sysconfig['repobranch'] + '/' 

    # Loop over each run source in list
    for source in runSources:
        dirName = source.name 
        if re.search('nonProduction', source.name):
            start = source.name.find('nonProduction') + 14
            dirName = source.name[start:]

        for comp in source.getOpt('compilers').split(','):
            nmodes=0
            exps = []
            for mode in source.getOpt('modes').split(','):
                nmodes = nmodes + 1
                decksDir = scratchDir+comp+'/'+dirName+'.' + mode.upper()+'/decks'
                os.chdir(decksDir)
                os.environ['MYCONFIGDIR'] = decksDir
        # Create rundeck object with default or config properties
                rundeck = RunSourceProperties(source.name)
                exps.append(Arun(rundeck, mode))

        # List of rundeck run configurations for each mode
 #               for mode in source.getOpt('modes').split(','):
                     
        # For each rundeck create a diffFile with verification results
                diffFile = rundeck.resultsDir + '/' + rundeck.name + '.diff'
                fileH = open(diffFile, 'w')
                
                for exp in exps:
                    # Did executable build?
                    exe = dirName+'.'+mode+'.'+comp
                    cmd = 'ls '+exe+'_bin/'+(exe+'.exe')
                    status = exp.sysCmd(cmd, 3, 'b')
                    # If not, then go on to next experiment
                    if status != 0:
                        continue

                    # Else, do comparisons
                    if source.getOpt('testlevel') == 'full':
                        if exp.mode == 'serial':
                            compareBase(exp, '1hr')
                            compareBase(exp, exp.etSuffix)
                            # And compare SERIAL checkpoint-restart 
                            compareRestart(exp)
                        else:
                            for npes in source.getOpt('npes').split(','):
                                # Compare runs with baseline
                                compareBase(exp, '1hr', npes=npes)
                                compareBase(exp, exp.etSuffix, npes=npes)
                                compareRestart(exp, npes=npes)
                                for npes in source.getOpt('npes').split(','):
                                    # Compare 1hr run against serial
                                    if nmodes > 1:
                                        compareNPE(exps[0], exps[1], exp.etSuffix, npes)
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

#-------------------------------------------------------------------------------
# Create a diff report and notify via email
def sendDiffreport(config, bpconfig, eTime):
    userconfig  = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    mailto     = userconfig['mailto']
    branch     = userconfig['repobranch']
    resultsDir = userconfig['scratchdir'] + '/regression_results/' + branch
    compflags  = userconfig['compflags']
    sortdiff   = userconfig['sortdiff']
    testType   = userconfig['testtype']
    compilers  = regUtils.getCompilers(bpconfig)

    diffFile = resultsDir + '/' + 'diffreport.txt'
    fp = open(diffFile, 'w')
    fp.write('ModelE test results \n')
    fp.write('-'*20+'\n')
    fp.write('Branch:' + branch + ',  Compiler flags:' + compflags + \
        ',  Test type:' + testType + '\n')
    fp.write('-'*62+'\n')
    fp.write('%62s\n' % ('-REPRODUCIBILITY'))
    fp.write('%20s%10s%8s%6s%6s%6s%6s\n' % \
        ('RUNDECK', 'COMPILER', 'MODE', 'RUN', 'BAS', 'RST', 'NPE'))
    fp.write('-'*62+'\n')

    if sortdiff == 'yes':
	# sort mode column
        subprocess.call('find '+resultsDir+' -name \*.diff -exec cat {} \; | sort -k 3,3 >' \
            +resultsDir + '/' + 'alldiffs', shell=True)
        with open(resultsDir + '/' + 'alldiffs','r') as inf:
            fp.write(inf.read())
    else:
        for comp in compilers:
            diffs = glob.glob(resultsDir + '/' + comp + '/*.diff')
            for f in diffs:
                with open(f,'r') as inf:
                    fp.write(inf.read())
    
    fp.write('-'*62+'\n')
    hhmmss = time.strftime('%H:%M:%S', time.gmtime(eTime))
    fp.write('Time taken = %s \n' %(hhmmss))
    fp.write('-'*62+'\n')
    fp.write('Legend:\n')
    fp.write('-'*7+'\n')
    fp.write('+  : success\n')
    fp.write('F  : reproducibility failure\n')
    fp.write('Fb : build failure\n')
    fp.write('F1 : 1hr run failure\n')
    fp.write('Fr : restart run failure\n')
    fp.write('F* : expected failure\n')
    fp.write('-  : not available\n')
    fp.write('Notes:\n')
    fp.write('-'*6+'\n')
    compconfig = regUtils.ConfigSectionMap(bpconfig, 'COMPCONFIG')
    compVers =  compconfig['compiler_versions'].split(",")
    i=0
    for comp in compilers:
        fp.write(comp+' compiler version: '+compVers[i]+'\n')
        i+=1
    fp.close()

    subject = '"modelE_RT (' + branch + ')" '
    cmd = '/usr/bin/mail -s ' + subject + mailto + ' < ' + diffFile
    subprocess.check_call(cmd, shell=True)

#-------------------------------------------------------------------------------
# Each task is roughly comprised of these steps:
# 0) create clones of reference clone for each rundeck/compiler/mode combination
# 1) create a rundeck instance from the modelE templates
# 2) build/setup the model
# 3) run the model for the specified duration
# 4) verification of model results
# Only (0) is done here. The rest is done by regression.py
def setupCloneTasks(config, bpconfig, decklist):
    userconfig =regUtils.ConfigSectionMap(config, 'USERCONFIG')
    compilers = regUtils.getCompilers(bpconfig)

    cloneTasks = []
    for deck in decklist:
      # since nonProduction rundeck names can be quite long, extract the
      # nonProduction_ part...
        dName = deck.name
        if re.search('nonProduction', deck.name):
            start = deck.name.find('nonProduction') + 14
            dName = deck.name[start:]
         
        for comp in deck.getOpt('compilers').split(','):
            if comp in compilers:
                commandString = regUtils.gitCloneCommand(config, dName, comp)
                cloneTasks.append(commandString)
            else:
                logger.error(comp+' is not defined in COMPCONFIG')
            
    for t in cloneTasks:
        logger.debug('CLONE TASK %s', t)
    return cloneTasks

#-------------------------------------------------------------------------------
# Return a command to submit/execute a [batch] job
def setupScriptTasks(config, bpconfig, decklist):
    compilers = regUtils.getCompilers(bpconfig)

    scriptTasks = []
    for deck in decklist:
        for comp in deck.getOpt('compilers').split(','):
            if comp in compilers:
                commandString = createScriptTask(config, bpconfig, deck, comp)
                scriptTasks.append(commandString)
            else:
                logger.error(comp+' is not defined in COMPCONFIG')
            
    for t in scriptTasks:
        logger.debug('SCRIPT TASK %s', t)
    return scriptTasks

#-------------------------------------------------------------------------------
# Creates script to be submitted to batch system OR to be executed interactively
# Note that there are (still)) several hardwired batch parameters.
def createScriptTask(config, bpconfig, deck, comp):
    userconfig  = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    modules    = userconfig['modules']
    useBatch   = userconfig['usebatch']
    branch     = userconfig['repobranch']
    compopts   = userconfig['compflags']
    debugReg   = userconfig['debugscript']
    scriptsDir = userconfig['scriptsdir'] + '/exec/testing/'
    useMods    = userconfig['modules']
    resultsDir = userconfig['scratchdir'] + '/regression_results/' + \
                 branch + '/' + comp
    scratchDir = userconfig['scratchdir'] + '/regression_scratch/' + \
                 branch + '/' + comp
    sponsorID  = userconfig['sponsorid']

    deckName = deck.name
    jobName = deckName
    if re.search('nonProduction', deckName):
       start = deckName.find('nonProduction') + 14
       jobName = deckName[start:]
    filename = resultsDir + '/' + jobName + '.bash'
    fileHandle = open ( filename, 'w' ) 

    if useBatch == 'yes':
        # If we are just compiling this rundeck
        if deck.getOpt('testlevel') == 'compileOnly':
            cores = 1
            walltime = '00:30:00'

        elif deck.getOpt('testlevel') == 'long':
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
            
        else:               
            if 'mpi' in deck.modes: 
            # Set number of cores (tasks)
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
                elif re.search('cadi', deckName):
                    cores = 8              
                    walltime = '06:00:00'

            else:
                cores = 1
                walltime = '03:00:00'
                if re.search('obio', deckName):
                    walltime = '01:30:00'
                elif re.search('cadi', deckName):
                    walltime = '06:00:00'
                elif re.search('tomas', deckName):
                    logger.error('This rundeck will not run in SERIAL:', deckName)
                    sys.exit()
                elif re.search('amp', deckName):
                    logger.error('This rundeck will not run in SERIAL:', deckName)
                    sys.exit()
                elif re.search('E_AR5_V2', deckName):
                    if re.search('NINT', deckName):
                        walltime = '01:00:00'
                    else: # CADI and CAMP
                        logger.error('This rundeck will not run in SERIAL:', deckName)
                        sys.exit()

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
        fileHandle.write ('#SBATCH --output='   + outname + '\n')
        fileHandle.write ('#SBATCH --error='    + errname + '\n')
        fileHandle.write ('#SBATCH --account='  + sponsorID + '\n')
        fileHandle.write ('#SBATCH --job-name=' + jobName + '\n')
        fileHandle.write ('#SBATCH --time='     + walltime + '\n')
        fileHandle.write ('#SBATCH --ntasks=' + str(cores) + '\n')
        # Use Haswell NODES
        fileHandle.write ('#SBATCH --constraint=hasw\n')

    # ELSE create rest of script for batch AND interactive job:

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
            fileHandle.write ('. /opt/local/share/Modules/3.2.10/init' + '\n')
            fileHandle.write ('module purge' + '\n')

        # Using different naming convention for module names
        compvendor = comp
        if comp == 'gfortran':
            compvendor = 'gcc'

        modsconfig = regUtils.ConfigSectionMap(bpconfig, 'COMPCONFIG')
        for mod in modsconfig['modulelist'].split(','):
            if re.search(compvendor, mod):
                for mm in modsconfig[mod].split(','):
                    cmd = 'module load ' + mm +'\n'
                    fileHandle.write (cmd)

    decksDir = scratchDir + '/' + jobName + '/decks/'
    fileHandle.write ('export DECKSDIR=' + decksDir + '\n')
    # The following variable is (optionally) exported to regression.py
    if debugReg == 'yes':
        fileHandle.write ('export DEBUG=1' + '\n')

    # cd to the working dir and run the script
    fileHandle.write ('cd ' + decksDir + '\n')
    fileHandle.write ('python ' + scriptsDir + '/' + 'regression.py ' + deckName + '\n')
    fileHandle.write (' ' + '\n')
    fileHandle.close()

    modelErc = scratchDir + '/modelErc.' + comp
    createRegConfig(config, deck, modelErc, comp, jobName)

    if useBatch == 'yes':
        commandString = 'sbatch ' + filename
    else:
        commandString = 'chmod +x ' + filename + ';' + filename
               
    return commandString

#-------------------------------------------------------------------------------
# Create a config file for regression.py script. 
# Note: there is one config file for each rundeck/compiler/mode combination
def createRegConfig(config, deck, modelerc, comp, jobName):
    cfg  = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    branch     = cfg['repobranch']
    regconfig = ConfigParser.RawConfigParser()
    decksDir = cfg['scratchdir'] + '/regression_scratch/' + \
        branch + '/' + comp + '/' + jobName + '/decks/'
    resultsDir = cfg['scratchdir'] + '/regression_results/' + branch + '/' + comp

    regconfig.add_section('regSettings')
    regconfig.set('regSettings', 'rundeck', deck.name)
    regconfig.set('regSettings', 'modelerc', modelerc)
    regconfig.set('regSettings', 'compiler', comp)
    regconfig.set('regSettings', 'modes', deck.getOpt('modes'))
    regconfig.set('regSettings', 'testlevel', deck.getOpt('testlevel'))
    regconfig.set('regSettings', 'endtime', deck.getOpt('endtime'))
    regconfig.set('regSettings', 'nplist', deck.getOpt('npes'))
    regconfig.set('regSettings', 'compflags', cfg['compflags'])
    regconfig.set('regSettings', 'repository', cfg['repository'])
    regconfig.set('regSettings', 'branch', cfg['repobranch'])
    regconfig.set('regSettings', 'basedir', cfg['basedir'])
    regconfig.set('regSettings', 'updatebase', cfg['updatebase'])
    regconfig.set('regSettings', 'systemtests', cfg['systemtests'])
    regconfig.set('regSettings', 'resultsdir', resultsDir)
    regconfig.set('regSettings', 'decksdir', decksDir)

    filename = decksDir + deck.name + '.cfg'
    logger.debug(filename)
    with open(filename, 'w') as configfile:
        regconfig.write(configfile)


