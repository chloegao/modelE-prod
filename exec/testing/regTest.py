#  Base class for regression tests
class regTest:
    def __init__(self, name):
        self.name = name
        self.modes = ['mpi']
        self.compilers = ['gfortran']
        self.compflags = 'default'
        self.npes = [8]
        self.testLevel = 'full'
        self.duration = 1
        self.endtime = 25
        # SYSCONFIG options 
        self.useBatch = 'no'
        self.modules = 'no'
        self.scrDir = '/tmp'
        self.resDir = '/tmp'
        self.nsteps = 2

    def setOpts(self, sysconfig, deckconfig):
        for name,options in deckconfig.items():
            if self.name == name:
                for kk,vv in options.items():
                    if kk=='modes':
                        self.modes = vv
                    elif kk=='compilers':
                        self.compilers = vv
                    elif kk=='endtime':
                        self.endtime = vv
                    elif kk=='testlevel':
                        self.testLevel = vv
                    elif kk=='npes':
                        self.npes = vv
                    elif kk=='compflags':
                        self.compile_only = vv
                    elif kk=='duration':
                        self.duration = vv
        self.useBatch = sysconfig['usebatch']
        self.modules = sysconfig['modules']
        self.resDir = sysconfig['scratchdir'] + '/regression_results/'
        self.scrDir = sysconfig['scratchdir'] + '/regression_scratch/'
        self.nsteps = self.duration * 2
            
    def getOpt(self, opt):
        if opt=='modes':
            return self.modes
        elif opt=='compilers':
            return self.compilers
        elif opt=='compflags':
            return self.compile_only
        elif opt=='npes':
            return self.npes
        elif opt=='duration':
            return self.duration
        elif opt=='testlevel':
            return self.testLevel
        elif opt=='endtime':
            return self.endtime
        elif opt=='nsteps':
            return self.nsteps
        elif opt=='usebatch':
            return self.useBatch
        elif opt=='modules':
            return self.modules
        elif opt=='resDir':
            return self.resDir
        elif opt=='scrDir':
            return self.scrDir

    def dump(self):
        print self.__dict__
