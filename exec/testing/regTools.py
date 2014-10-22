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
from regTest import *
import logging

logger = logging.getLogger('regTools')

#-------------------------------------------------------------------------------
# Workaround for "mkdir -p" command
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

#-------------------------------------------------------------------------------
# Read a config file, return config object
def readConfig(cfgfile):
    config = ConfigParser.ConfigParser()
    config.read(cfgfile)

    logger.debug('Read configuration file %s',cfgfile)
    return config

#-------------------------------------------------------------------------------
# Get a list of compilers used
def getCompilers(config):
   compconfig = ConfigSectionMap(config, 'COMPCONFIG')
   compilers = compconfig['compilers'].split(",")
   return compilers

#-------------------------------------------------------------------------------
# "Safe" way to clean the contents of a directory
def cleanDir(adir):
   if(adir == '/' or adir == "\\"): 
       logger.error('Cannot clean %s',adir)
       return
   else:
       for file_object in os.listdir(adir):
           logger.debug('Will clean up %s',adir)
           file_object_path = os.path.join(adir, file_object)
           if os.path.isfile(file_object_path):
               os.unlink(file_object_path)
           else:
               shutil.rmtree(file_object_path) 

#-------------------------------------------------------------------------------
# Return a dict (i.e. a key,value pair) from each section in config
def ConfigSectionMap(config, section):
    adict = {}
    options = config.options(section)

    for option in options:
        try:
            adict[option] = config.get(section, option)
            if adict[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            logger.error("exception on %s!" % option)
            adict[option] = None

    logger.debug('Read section %s',section)

    return adict

#-------------------------------------------------------------------------------
# From config file create/return a list of model run configurations
def getModelConfigurations(config):
    sections = config.sections()
    modelConfig = {}

    for sect in sections:
        match = not re.search("CONFIG",sect)
        # get all rundeck sections from config file
        if (match):
            modelConfig[sect] = ConfigSectionMap(config, sect)

# map modelConfig to a more manageable list
# first extract names
    runList = []
    for name,options in modelConfig.items():
       runList.append(regTest(name))
# then specified options
    sysconfig = ConfigSectionMap(config, 'SYSCONFIG')
    for d in runList:
       d.setOpts(sysconfig, modelConfig)
        
    return runList

#-------------------------------------------------------------------------------
# Setup testing environment:
# 1) Create working directories
# 2) Clone model from git repository and...
# 3) Perform additional model "specific" setup
def setupEnv(config):
   sysconfig = ConfigSectionMap(config, 'SYSCONFIG')
   branch =  sysconfig['repobranch']
   resultsDir = sysconfig['scratchdir'] + '/regression_results/' + branch
   scratchDir = sysconfig['scratchdir'] + '/regression_scratch/' + branch

   if not os.path.exists(resultsDir):
      mkdir_p(resultsDir)    
      mkdir_p(scratchDir)
   else:
      if sysconfig['cleanscratch'] == 'yes':
         cleanDir(scratchDir)
         cleanDir(resultsDir)

   setupModelEenv(config)
   gitCloneRepository(config)


#-------------------------------------------------------------------------------
# Clone the model from the user-specified git repository
def gitCloneRepository(config):
   sysconfig = ConfigSectionMap(config, 'SYSCONFIG')
   scratch = sysconfig['scratchdir']
   repo = sysconfig['repository']
   branch =  sysconfig['repobranch']
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
def setupModelEenv(config):
   sysconfig = ConfigSectionMap(config, 'SYSCONFIG')
   branch =  sysconfig['repobranch']
   resultsDir = sysconfig['scratchdir'] + '/regression_results/' + branch
   scratchDir = sysconfig['scratchdir'] + '/regression_scratch/' + branch

# the following directories are modelE specific:
   mkdir_p(scratchDir+'/decks_repository')
   mkdir_p(scratchDir+'/cmrun')
   mkdir_p(scratchDir+'/exec')
   mkdir_p(scratchDir+'/savedisk')

# We need to get a list of compilers...
   compilers = getCompilers(config)
   libsconfig = ConfigSectionMap(config, 'LIBSCONFIG')

# ... to create modelErc file(s)
   for comp in compilers:
      if not os.path.exists(scratchDir + comp):
         mkdir_p(resultsDir + '/' + comp)
         mkdir_p(scratchDir + '/' + comp)
      writeModelErc(libsconfig, scratchDir, comp)

# DISCOVER hack:
# while discover git installation is too old...
# or just load the git module
   machine = os.getenv('HOST')
   if machine == 'discover' or machine == 'borg':
      os.environ["PATH"] += os.pathsep + \
      '/usr/local/other/SLES11.1/git/1.8.5.2/libexec/git-core/git'

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
# Create a diff report and notify via email
def sendDiffreport(config):
    sysconfig = ConfigSectionMap(config, 'SYSCONFIG')
    branch =  sysconfig['repobranch']
    resultsDir = sysconfig['scratchdir'] + '/regression_results/' + branch
    mailto =  sysconfig['mailto']
    compilers = getCompilers(config)

    diffFile = resultsDir + '/' + 'diffreport.txt'
    fp = open(diffFile, 'w')
    fp.write('ModelE test results, branch=' + branch + '\n')
    fp.write('-'*62+'\n')
    fp.write('%62s\n' % ('-REPRODUCIBILITY'))
    fp.write('%20s%10s%8s%6s%6s%6s%6s\n' % \
        ('RUNDECK', 'COMPILER', 'MODE', 'RUN', 'BAS', 'RST', 'NPE'))
    fp.write('-'*62+'\n')
    fp.close()
    for comp in compilers:
        diffs = glob.glob(resultsDir + '/' + comp + '/*.diff')
        with open(diffFile, 'a') as out:
            for f in diffs:
                with open(f,'r') as inf:
                    out.write(inf.read())

    subject = '"modelE_RT (' + branch + ')" '
    cmd = '/usr/bin/mail -s ' + subject + mailto + ' < ' + diffFile
    subprocess.check_call(cmd, shell=True)

