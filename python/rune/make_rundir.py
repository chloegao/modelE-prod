import os
srcdir = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append(os.path.join(srcdir, '..', 'lib'))

import modele.rundir

import string
import tempfile
import filecmp
import shutil
import modele.modelerc


rundeck = sys.argv[1]
rundir = sys.argv[2]

modele.rundir.make_rundir(rundeck, rundir)
