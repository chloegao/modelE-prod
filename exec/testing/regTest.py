#  Base class for modelE regression tests
class regTest:
    def __init__(self, name):
        self.name = name
        self.modes = ['mpi']
        self.compilers = ['gfortran']
        self.buildtype = 'release'
        self.npes = [1,4]
        self.verification = 'restartRun'
        # for restartRun, 25hr endtime
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
                    elif kk=='verification':
                        self.verification = vv
                    elif kk=='npes':
                        self.npes = vv
                    elif kk=='buildtype':
                        self.compile_only = vv
        self.useBatch = sysconfig['usebatch']
        self.modules = sysconfig['modules']
        self.resDir = sysconfig['scratchdir'] + '/results/'
        self.scrDir = sysconfig['scratchdir'] + '/scratch/'
            
    def getOpt(self, opt):
        if opt=='modes':
            return self.modes
        elif opt=='compilers':
            return self.compilers
        elif opt=='buildtype':
            return self.compile_only
        elif opt=='npes':
            return self.npes
        elif opt=='verification':
            return self.verification
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
