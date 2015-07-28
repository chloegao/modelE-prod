# This is the main driver for the modelE regression tests. To run the scripts:
#
#   $  python mainDriver.py [configuration file name]
#
import time
import sys
import os.path
import regUtils
import regPool
import regTools as tools
import logging

logger = logging.getLogger('main')

#-------------------------------------------------------------------------------
# MAIN DRIVER
def main():
    starttime = time.time()
    useMessage = 'Usage: python ' + sys.argv[0] + ' <configFileName> # no file extension'
    if len(sys.argv)==1:
        print useMessage
        sys.exit()
    else:
        cfgfile = str(sys.argv[1]) + '.cfg'
        if not os.path.isfile(cfgfile):
            print 'Error: ' + cfgfile + ' : file does not exist'
            print useMessage
            sys.exit()

# Logger setup       
    logging.basicConfig(
        filename = str(sys.argv[1]) + '.LOG',
        format = "%(levelname) -10s %(module)s:%(lineno)s %(funcName)s %(message)s",
        level = logging.DEBUG,
        filemode = 'w'
    )
    stdoutLog = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(name)s : %(message)s')
    stdoutLog.setFormatter(formatter)
    if os.environ.has_key('DEBUG'):
        stdoutLog.setLevel(logging.DEBUG)
    else:
        stdoutLog.setLevel(logging.INFO)
    logger = logging.getLogger()
    logger.addHandler(stdoutLog)

# Read user-defined config file and store in a config object
    config = regUtils.readConfig(cfgfile)

# COMPCONFIG section contains computational configuration information (compilers,
# libraries, etc). There are COMPCONFIG defaults in file comp.cfg but those can be
# overridden in the user-defined config file by re-defining the defaults.
    if not config.has_section("COMPCONFIG"):
        compconfig = regUtils.readConfig('comp.cfg')
    else:
        compconfig = config

    userconfig = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    makesystem = userconfig['makesystem']

# --- modelE specific workflow ---
# Config file contains a list of model configurations pertinent to rundecks.
# For convenience store that list separately:
    runList = regUtils.getModelConfigurations(config)

# Let's setup the testing environment, specific to modelE
    tools.setupEnv(config, compconfig)

    if makesystem == 'makeOld':
        # Create gitTasks
        gitTasks = tools.setupCloneTasks(config, compconfig, runList)
        # ... and execute them (if NOT debugging)
        if not os.environ.has_key('DEBUG'):
            regPool.runCommands(gitTasks, 'no')
    else:
        # Setup run directories for out of source builds
        tools.setupRuns(config, compconfig, runList)

# Create scripts
    scriptTasks = tools.setupScriptTasks(config, compconfig, runList)

# ... and run them
    useBatch = userconfig['usebatch']
    regPool.runCommands(scriptTasks, useBatch)

    eTime =  time.time()-starttime

# Verify runs and notify (only if "mailto" field is not empty)
    tools.verifyRuns(config, runList)
    if userconfig['mailto']:
        tools.sendDiffreport(config, compconfig, eTime)
# -------------------------------

    logger.info('Time taken = %f' %(eTime))

#-------------------------------------------------------------------------------
# MAIN PROGRAM 
if __name__ == "__main__":
    main()
    logger.info('Regression tests are done!')
