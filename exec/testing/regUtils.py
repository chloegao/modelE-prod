# This module contains uyilities to help setup the regression tests
import string
import ConfigParser
import os
import sys
import errno
import shutil
import subprocess
import logging
from regTest import *

logger = logging.getLogger('regUtils')

#-------------------------------------------------------------------------------
# Return a command that creates a clone of the reference clone
def gitCloneCommand(config, deckname, compiler):
    userconfig = ConfigSectionMap(config, 'USERCONFIG')
    branch    = userconfig['repobranch']
    scratch   = userconfig['scratchdir'] + '/regression_scratch/' + branch + '/'
    reference = scratch + '/' + branch
    clone     = scratch + '/' + compiler + '/' + deckname
# if clone does not exist then create it
    if not os.path.isdir(clone):
        s = string.Template('git clone -b $b $r $t > /dev/null 2>&1')
        return s.substitute(b=branch, r=reference, t=clone)

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
#  Test if an executable program exists in the path - like unix's which
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
   
#-------------------------------------------------------------------------------
# Read a config file, return config object
def readConfig(cfgfile):
    config = ConfigParser.ConfigParser()
    config.read(cfgfile)

    logger.debug('Read configuration file %s',cfgfile)
    return config

#-------------------------------------------------------------------------------
# Setup testing environment:
# 1) Create working directories
# 2) Clone model from git repository and...
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

   getCodeFromRepository(config)

#-------------------------------------------------------------------------------
# Get a list of compilers used
def getCompilers(config):
   compconfig = ConfigSectionMap(config, 'COMPCONFIG')
   compilers = compconfig['compilers'].split(",")
   return compilers


#-------------------------------------------------------------------------------
# Clone the model from the user-specified git repository
def getCodeFromRepository(config):
   sysconfig = ConfigSectionMap(config, 'SYSCONFIG')
   scratch = sysconfig['scratchdir']
   repo = sysconfig['repository']
   branch =  sysconfig['repobranch']
# Create a  clone named 'branch' in a  directory named 'branch'
   clone = scratch + '/regression_scratch/' + branch + '/' + branch
   logger.debug('Will access repository %s',repo)

# Check if remote is a valid repository
   p = subprocess.Popen(["git", "ls-remote", repo], stdout=subprocess.PIPE)
   output = p.communicate()[0]
   if p.returncode != 0:
      logger.error("Specified git repository does not exists. RC=["+str(p.returncode)+"]")
      # If we have no code, exit :-(
      sys.exit(0)

   cwd = os.getcwd()
# Always fresh-clone      
   logger.debug('Cloning %s into %s', repo, clone)
   # If scratch space is not "cleaned" then there may be a repository
   if os.path.isdir(clone):
       os.chdir(clone)
       if os.path.isdir('.git'):
           logger.warning('%s is already a git repository', clone)
           subprocess.check_call('git pull', shell=True)
   else:
       cmd = 'git clone -b ' + branch + ' ' + repo + ' ' + clone \
          + '> /dev/null 2>&1'
       subprocess.check_call(cmd, shell=True)

   os.chdir(cwd)
