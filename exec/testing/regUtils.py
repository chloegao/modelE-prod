# This module contains utilities to help setup the regression tests
import string
import ConfigParser
import os
import re
import sys
import errno
import shutil
import subprocess
import logging
from regTest import *

logger = logging.getLogger('utils')

import os, datetime

#-------------------------------------------------------------------------------
# Create directory wit timestamp
def mkdirTimeSTamp(list, filename):
    mydir = os.path.join(os.getcwd(), datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    try:
        os.makedirs(mydir)
    except OSError, e:
        if e.errno != 17:
            raise # This was not a "directory exist" error..
    with open(os.path.join(mydir, filename), 'w') as d:
        d.writelines(list)

#-------------------------------------------------------------------------------
# Create/return a list of model run configurations specified in config file
def getModelConfigurations(config):
    sections = config.sections()
    modelConfig = {}

# Retrieve all model run configurations. For convenience divide the sections
# in the configuration file into two types: CONFIG and others. The former 
# have CONFIG in their names. Thus, if a section name does NOT have CONFIG 
# in its name then it is a rundeck configuration.
    for sect in sections:
        match = not re.search("CONFIG",sect)
        # get all rundeck sections from config file
        if (match):
            modelConfig[sect] = ConfigSectionMap(config, sect)

# Store each model run configurations in a list and let each item in the list
# have access to the user-defined options
    runList = []
    for name,options in modelConfig.items():
        # Each item is a regression test (regTest) instance
        runList.append(regTest(name))
    userconfig = ConfigSectionMap(config, 'USERCONFIG')
    for d in runList:
        d.setOpts(userconfig, modelConfig)
        
    return runList

#-------------------------------------------------------------------------------
# Create a directory composed of various user-defined attributes
def mkdirCommand(config, deckname, compiler, cmode):
    userconfig = ConfigSectionMap(config, 'USERCONFIG')
    branch    = userconfig['repobranch']
    scratch   = userconfig['scratchdir'] + '/scratch/' + branch + '/'
    reference = scratch + '/' + branch
    adir      = scratch + '/' + compiler + '/' + deckname + cmode
    if not os.path.isdir(adir):
        mkdir_p(adir)

#-------------------------------------------------------------------------------
# Return a command that creates a clone of the reference clone
def gitCloneCommand(config, deckname, compiler, cmode):
    userconfig = ConfigSectionMap(config, 'USERCONFIG')
    branch    = userconfig['repobranch']
    scratch   = userconfig['scratchdir'] + '/scratch/' + branch + '/'
    reference = scratch + '/' + branch
    clone     = scratch + '/' + compiler + '/' + deckname + cmode
# if clone does not exist then create it
    if not os.path.isdir(clone):
        s = string.Template('git clone -b $b $r $t > /dev/null 2>&1')
        return s.substitute(b=branch, r=reference, t=clone)
    else:
        logger.debug('Git clone %s exists', clone)
        return clone

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
    logger.info('Cleaning up scratch space...')
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

    logger.info('Read configuration file %s',cfgfile)
    return config

#-------------------------------------------------------------------------------
# Get a list of compilers used
def getCompilers(config):
   compconfig = ConfigSectionMap(config, 'COMPCONFIG')
   compilers = compconfig['compilers'].split(",")
   return compilers

