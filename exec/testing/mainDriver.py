# This is the main driver for the modelE regression tests. To run the scripts:
#    python mainDriver.py [configuration file]
import time
import sys
import regTools
import regTasks
import regPool
import logging
import os.path

logger = logging.getLogger('main')

#-------------------------------------------------------------------------------
# MAIN DRIVER
def main():
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
    config = regTools.readConfig(cfgfile)
# Get default COMP options from a config file
    if not config.has_section("COMPCONFIG"):
        bpconfig = regTools.readConfig('system.cfg')
    else:
        bpconfig = config

# config file contains a list of model configurations. In modelE these
# configurations are called rundecks. For convenience store that information in 
# a separate list:
    runList = regTools.getModelConfigurations(config)

# Let's setup the testing environment:
    regTools.setupEnv(config, bpconfig)

# Create gitTasks
    gitTasks = regTasks.setupCloneTasks(config, bpconfig, runList)
# ... and execute them
    regPool.runCommands(gitTasks, 'no')

# Create scripts
    scriptTasks = regTasks.setupScriptTasks(config, bpconfig, runList)
# ... and run them
    userconfig = regTools.ConfigSectionMap(config, 'USERCONFIG')
    useBatch = userconfig['usebatch']
    regPool.runCommands(scriptTasks, useBatch)

# Gather results and notify
    if userconfig['diffreport'] == 'yes':
        regTools.sendDiffreport(config, bpconfig)

#-------------------------------------------------------------------------------
# MAIN PROGRAM 
if __name__ == "__main__":
    starttime = time.time()
    main()
    logger.info('Regression tests are done!')
    logger.info('Time taken = %f' %(time.time()-starttime))
