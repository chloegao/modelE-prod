#!/usr/bin/env python
#
# Executes Python with a PYTHONPATH specific for one build (rundeck)

import sys
import os

# Get directory of this script
this_script_dir = os.path.dirname(os.path.realpath(__file__))

# Extract name of script to run from the arguments
build_name = sys.argv[1]
target_script_dir = os.path.dirname(os.path.realpath(sys.argv[2]))

# Locate the build directory based on its name.  Right now,
# just look in ${MODELE}/<build-name>.  But in the future,
# we might scan down a search path, look in cmrun, look
# at env variables, etc.
modele_root = os.path.join(this_script_dir, '..')
build_dir = os.path.join(modele_root, build_name)

# Add build-specific directory to Python path
sys.path = [
	os.path.join(build_dir, 'pyext'),
	os.path.join(modele_root, 'python', 'lib'),
	target_script_dir,
	] + sys.path

# Strip our stuff off of argv and execute the resulting file
sys.argv = sys.argv[2:]
execfile(sys.argv[0])
