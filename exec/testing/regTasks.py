# This module creates commands to run the modelE regression tests
import string
import ConfigParser
import re
import os
import sys
import regTools
import regTest
import logging
import subprocess
import ConfigParser

logger = logging.getLogger('regTasks')

#-------------------------------------------------------------------------------
# Return a command that creates a clone of the reference clone
def gitCloneCommand(config, deckname, compiler):
    sysconfig = regTools.ConfigSectionMap(config, 'SYSCONFIG')
    branch    = sysconfig['repobranch']
    scratch   = sysconfig['scratchdir'] + '/regression_scratch/' + branch + '/'
    reference = scratch + '/' + branch
    clone     = scratch + '/' + compiler + '/' + deckname
# if clone does not exist then create it
    if not os.path.isdir(clone):
        s = string.Template('git clone -b $b $r $t > /dev/null 2>&1')
        return s.substitute(b=branch, r=reference, t=clone)
   
#-------------------------------------------------------------------------------
# Each task is roughly comprised of these steps:
# 0) create clones of reference clone for each rundeck/compiler/mode combination
# 1) create a rundeck instance from the modelE templates
# 2) build/setup the model
# 3) run the model for the specified duration
# 4) verification of model results
# Only (0) is done here. The rest is done by regression.py
def setupCloneTasks(config, decklist):
    sysconfig = regTools.ConfigSectionMap(config, 'SYSCONFIG')
    compilers = regTools.getCompilers(config)

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
                commandString = gitCloneCommand(config, dName, comp)
                cloneTasks.append(commandString)
            else:
                logger.error(comp+' is not defined in COMPCONFIG')
            
    for t in cloneTasks:
        logger.debug('CLONE TASK %s', t)
    return cloneTasks

#-------------------------------------------------------------------------------
# Return a command to submit/execute a [batch] job
def setupScriptTasks(config, decklist):
    compilers = regTools.getCompilers(config)

    scriptTasks = []
    for deck in decklist:
        for comp in deck.getOpt('compilers').split(','):
            if comp in compilers:
                commandString = createScriptTask(config, deck, comp)
                scriptTasks.append(commandString)
            else:
                logger.error(comp+' is not defined in COMPCONFIG')
            
    for t in scriptTasks:
        logger.debug('SCRIPT TASK %s', t)
    return scriptTasks

#-------------------------------------------------------------------------------
# Creates script to be submitted to batch system OR to be executed interactively
# Note that there are (still)) several hardwired batch parameters.
def createScriptTask(config, deck, comp):
    sysconfig  = regTools.ConfigSectionMap(config, 'SYSCONFIG')
    modules    = sysconfig['modules']
    useBatch   = sysconfig['usebatch']
    branch     = sysconfig['repobranch']
    compopts   = sysconfig['compflags']
    debugReg   = sysconfig['debugscript']
    scriptsDir = sysconfig['scriptsdir']
    useMods    = sysconfig['modules']
    resultsDir = sysconfig['scratchdir'] + '/regression_results/' + \
                 branch + '/' + comp
    scratchDir = sysconfig['scratchdir'] + '/regression_scratch/' + \
                 branch + '/' + comp
    sponsorID  = sysconfig['sponsorid']

    deckName = deck.name
    jobName = deckName
    if re.search('nonProduction', deckName):
       start = deckName.find('nonProduction') + 14
       jobName = deckName[start:]
    filename = resultsDir + '/' + jobName + '.bash'
    fileHandle = open ( filename, 'w' ) 

    if useBatch == 'yes':
        # Set number of cores (tasks)
        if 'MPI' in deckName:
            cores = 8
            if re.search('tomas', deckName):
                cores = 88
            elif re.search('amp', deckName):
                cores = 44
            elif re.search('E_AR5_V2', deckName):
                if re.search('NINT', deckName):
                    cores = 8
                else: # CADI and CAMP
                    cores = 44
        else:
            cores = 1
        # Set the walltime
        walltime = '03:00:00'
        if re.search('C12', deckName):
            walltime = '00:30:00'
        elif re.search('Mars', deckName):
            walltime = '00:30:00'
        elif re.search('SGP', deckName):
            walltime = '00:30:00'
        elif re.search('M20', deckName):
            walltime = '00:30:00'
        elif re.search('obio', deckName):
            walltime = '02:00:00'
        elif re.search('cadi', deckName):
            walltime = '03:00:00'
        elif re.search('tomas', deckName):
            logger.error('This rundeck will not run in SERIAL:', deckName)
            sys.exit()
        elif re.search('amp', deckName):
            logger.error('This rundeck will not run in SERIAL:', deckName)
            sys.exit()
        elif re.search('E_AR5_V2', deckName):
            if re.search('NINT', deckName):
                walltime = '00:30:00'
            else: # CADI and CAMP
                logger.error('This rundeck will not run in SERIAL:', deckName)
                sys.exit()

        # Redine wall-time if compileOnly test
        if deck.getOpt('testlevel')=='compileOnly':
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

    # ELSE create rest of scriptfor batch AND interactive job:

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

        modsconfig = regTools.ConfigSectionMap(config, 'COMPCONFIG')
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
    cfg  = regTools.ConfigSectionMap(config, 'SYSCONFIG')
    branch     = cfg['repobranch']
    regconfig = ConfigParser.RawConfigParser()
    decksDir = cfg['scratchdir'] + '/regression_scratch/' + \
        branch + '/' + comp + '/' + jobName + '/decks/'
    resultsDir = cfg['scratchdir'] + '/regression_results/' + branch + '/' + comp

    if 'MPI' in jobName:
        mode = 'mpi'
        npes = deck.npes
    else:
        mode = 'serial'
        npes = 1

    regconfig.add_section('regSettings')
    regconfig.set('regSettings', 'rundeck', deck.name)
    regconfig.set('regSettings', 'modelerc', modelerc)
    regconfig.set('regSettings', 'compiler', comp)
    regconfig.set('regSettings', 'modes', mode)
    regconfig.set('regSettings', 'testlevel', deck.getOpt('testlevel'))
    regconfig.set('regSettings', 'endtime', deck.getOpt('endtime'))
    regconfig.set('regSettings', 'nplist', npes)
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

 
    



