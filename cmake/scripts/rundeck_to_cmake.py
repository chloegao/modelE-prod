import os
srcdir = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append(os.path.join(srcdir, '..', '..', 'python', 'lib'))
import modele.rundeck
import re
import tempfile
import filecmp
import shutil

# ------------------------------------------------------
class WriteIfDifferent(object):
    """Allows user to write to a temporary file, then move it
    to the destination only if it is different from the destination."""
    def __init__(self, ofname, **kwargs):
        """ofname: Name we ultimately want to write to."""
        self.ofname = ofname
        self.file = tempfile.NamedTemporaryFile(delete=False, **kwargs)
        self.tfname = self.file.name

    def close(self):
        self.file.close()
        try:
            if filecmp.cmp(self.tfname, self.ofname):
                # Files are equal, we are done!
                os.remove(self.tfname)
                return
        except: pass # Error means the files were NOT equal.

        # Files are not equal, so copy the temporary file over.
        shutil.copyfile(self.tfname, self.ofname)
        os.remove(self.tfname)

# ------------------------------------------------------
print('==================== BEGIN rundeck_to_cmake.py')
print('sys.argv', sys.argv)

# Get command line arguments
source_root = sys.argv[1]       # ${PROJECT_SOURCE_ROOT}
build_root = sys.argv[2]        # ${PROJECT_BUILD_ROOT}
rundeck_fname = sys.argv[3]
defines = sys.argv[4:]          # Add these to rundeck_opts.h

rundeck = modele.rundeck.load_rundeck(rundeck_fname, auto_download=False)
build = rundeck.build    # We only care about the build part of the rundeck
build.components['profiler'] = None     # Include profiler in all builds
#build.components['landice'] = None     # Include landice in all builds


src_files = []

# Find the source files (Object Modules)
for obj_module in build.sources:
    src_file = None
    for ext in ('.f', '.F', '.f90', '.F90'):
        fname = os.path.join(source_root, 'model', obj_module + ext)
        if os.path.exists(fname):
            src_file = obj_module+ext
            break

    if src_file is None:
        raise ValueError('Cannot find source file for object module {0}'.format(obj_module))
    src_files.append(src_file)

# Write out our source files for top-level ModelE build
with open(os.path.join(build_root, 'model', 'modele_SOURCES.cmake'), 'w') as out:
    out.write('# Machine-generated, DO NOT EDIT.\n')

    # Write out component options
    out.write('\n# Component Options\n')
    for component,options in build.components.items():
        if options is None:
            continue
        for name,value in options.items():
            out.write('set({0} {1})\n'.format(name,value))

    # Now include component lists of files
    # (which could depend on component options)
    out.write('\n# Get lists of files in components we will use\n')
    for component in build.components.keys():
        out.write('include(${CMAKE_CURRENT_SOURCE_DIR}/' + '{0}/{0}_SOURCES.cmake)\n'.format(component))

    out.write('\n#Set main list of files for the ModelE library\n')
    out.write('set(modele_SOURCES')

    for component in build.components.keys():
        out.write('\n\t${%s_SOURCES}' % component)

    for src_file in src_files:
        out.write('\n\t${CMAKE_CURRENT_SOURCE_DIR}/'+src_file)

    out.write(')\n')

# Write preprocessor #define's into the rundeck_opts.h file
rundeck_opts = WriteIfDifferent(os.path.join(build_root, 'model', 'rundeck_opts.h'), mode='w')
out = rundeck_opts.file

for symbol,definition in build.defines.items():
    if definition is None:
        out.write('#define {0}\n'.format(symbol))
    else:
        out.write('#define {0} {1}\n'.format(symbol, definition))

# From command line
for deff in defines:
    out.write('#define {0}\n'.format(deff))
rundeck_opts.close()

print('==================== END rundeck_to_cmake.py')
