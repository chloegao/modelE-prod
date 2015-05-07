# This is the main driver for the modelE regression tests. To run the scripts:
# python mainDriver.py [configuration file]
import time
import sys
import os.path
import regUtils
import regPool
import modelE
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
       
    logging.basicConfig(
        filename = str(sys.argv[1]) + '.LOG',
        format = "%(levelname) -10s %(module)s:%(lineno)s %(funcName)s %(message)s",
        level = logging.DEBUG,
        filemode = 'w'
    )

# Read user-defined config file and store in a config object
    config = regUtils.readConfig(cfgfile)
# Get default COMP options from a config file
    if not config.has_section("COMPCONFIG"):
        bpconfig = regUtils.readConfig('system.cfg')
    else:
        bpconfig = config

# config file contains a list of model configurations. In modelE these
# configurations are called rundecks. For convenience store that information in 
# a separate list:
    runList = modelE.getModelConfigurations(config)

# Let's setup the testing environment:
    modelE.setupEnv(config, bpconfig)

# Create gitTasks
    gitTasks = modelE.setupCloneTasks(config, bpconfig, runList)
# ... and execute them
    regPool.runCommands(gitTasks, 'no')

# Create scripts
    scriptTasks = modelE.setupScriptTasks(config, bpconfig, runList)
# ... and run them
    userconfig = regUtils.ConfigSectionMap(config, 'USERCONFIG')
    useBatch = userconfig['usebatch']
    regPool.runCommands(scriptTasks, useBatch)

    eTime =  time.time()-starttime
# Gather results and notify
    if userconfig['diffreport'] == 'yes':
#        modelE.createDiffreport(config, runList)
        modelE.sendDiffreport(config, bpconfig, eTime)
    logger.info('Time taken = %f' %(eTime))

#-------------------------------------------------------------------------------
# MAIN PROGRAM 
if __name__ == "__main__":
    main()
    logger.info('Regression tests are done!')
