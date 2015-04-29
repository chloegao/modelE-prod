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
def getCompilers(cfg):
   compconfig = ConfigSectionMap(cfg, 'COMPCONFIG')
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

# For convenience divide the sections in the configuration file into
# two types: CONFIG and RUNDECKS. The former have CONFIG in their names.
# Thus, if a section name does NOT have CONFIG in its name then it is
# a rundeck configuration:
    for sect in sections:
        match = not re.search("CONFIG",sect)
        # get all rundeck sections from config file
        if (match):
            modelConfig[sect] = ConfigSectionMap(config, sect)

# map modelConfig to a more manageable list
    runList = []
    for name,options in modelConfig.items():
       runList.append(regTest(name))
# Each item in runList (each rundeck) also needs user-defined options
    userconfig = ConfigSectionMap(config, 'USERCONFIG')
    for d in runList:
       d.setOpts(userconfig, modelConfig)
        
    return runList

#-------------------------------------------------------------------------------
# Setup testing environment:
# 1) Create working directories
# 2) Clone model from git repository and...
# 3) Perform additional model "specific" setup
def setupEnv(config, bpconfig):
   userconfig = ConfigSectionMap(config, 'USERCONFIG')
   branch =  userconfig['repobranch']
   resultsDir = userconfig['scratchdir'] + '/regression_results/' + branch
   scratchDir = userconfig['scratchdir'] + '/regression_scratch/' + branch

   if not os.path.exists(resultsDir):
      mkdir_p(resultsDir)    
      mkdir_p(scratchDir)
   else:
      if userconfig['cleanscratch'] == 'yes':
         cleanDir(scratchDir)
         cleanDir(resultsDir)

   setupModelEenv(config, bpconfig)
   gitCloneRepository(config)


#-------------------------------------------------------------------------------
# Clone the model from the user-specified git repository
def gitCloneRepository(config):
   userconfig = ConfigSectionMap(config, 'USERCONFIG')
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
   userconfig = ConfigSectionMap(config, 'USERCONFIG')
   branch =  userconfig['repobranch']
   resultsDir = userconfig['scratchdir'] + '/regression_results/' + branch
   scratchDir = userconfig['scratchdir'] + '/regression_scratch/' + branch

# the following directories are modelE specific:
   mkdir_p(scratchDir+'/decks_repository')
   mkdir_p(scratchDir+'/cmrun')
   mkdir_p(scratchDir+'/exec')
   mkdir_p(scratchDir+'/savedisk')

# We need to get a list of compilers...
   compilers = getCompilers(bpconfig)
   libsconfig = ConfigSectionMap(bpconfig, 'COMPCONFIG')

# ... to create modelErc file(s)
   for comp in compilers:
      if not os.path.exists(scratchDir + comp):
         mkdir_p(resultsDir + '/' + comp)
         mkdir_p(scratchDir + '/' + comp)
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
# Create a diff report and notify via email
def sendDiffreport(config, bpconfig):
    userconfig  = ConfigSectionMap(config, 'USERCONFIG')
    mailto     = userconfig['mailto']
    branch     = userconfig['repobranch']
    resultsDir = userconfig['scratchdir'] + '/regression_results/' + branch
    compflags  = userconfig['compflags']
    sortdiff   = userconfig['sortdiff']
    compilers  = getCompilers(bpconfig)

    diffFile = resultsDir + '/' + 'diffreport.txt'
    fp = open(diffFile, 'w')
    fp.write('ModelE test results, branch=' + branch + \
        ', compiler flags=' + compflags + '\n')
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
    compconfig = ConfigSectionMap(bpconfig, 'COMPCONFIG')
    compVers =  compconfig['compiler_versions'].split(",")
    i=0
    for comp in compilers:
        fp.write(comp+' compiler version: '+compVers[i]+'\n')
        i+=1
    fp.close()

    subject = '"modelE_RT (' + branch + ')" '
    cmd = '/usr/bin/mail -s ' + subject + mailto + ' < ' + diffFile
    subprocess.check_call(cmd, shell=True)

