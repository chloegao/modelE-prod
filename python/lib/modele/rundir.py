from modele import rundeck

import sys
import os
import string
import tempfile
import filecmp
import shutil
from modele import pathutil

# TODO: Be careful not to leave around zero-length files when downloading

# --------------------------------------

def namelist_time(suffix, dt):
    return 'YEAR{0}={1},MONTH{0}={2},DATE{0}={3},HOUR{0}={4},' \
        .format(suffix,dt.year,dt.month,dt.day,dt.hour)

def make_rundir(rd, rundir):
    ret = True

    # output line sections
    parameters = []
    data_files = []
    data_lines = []
    inputz = []
    inputz_cold = []

    # Organize parameters into ModelE sections
    for param in sorted(list(rd.params.values())):
        pname = param.pname
        if len(pname) == 1:
            if (id(param.type) == id(rundeck.FILE)):
                if param.value is not None:
                    data_lines.append(" _file_{}='{}'".format(pname[0], param.value))
                    data_files.append((pname[0], param.value))
#                else:
#                    parameters.append("! Not Found: {}={}".format(pname[0], param.sval))

            elif (id(param.type) == id(rundeck.GENERAL)):
                parameters.append(' %s=%s' % (param.pname[0], param.value))
            elif (id(param.type) == id(rundeck.DATETIME)):
                raise ValueError('Cannot put DATETIME values into parameters section of rundeck.')
            else:
                raise ValueError('Unknown parameter type %s' % param.type)
        elif len(pname) == 2:
            if pname[0].lower() == 'inputz':
                iz = inputz
            elif pname[0].lower() == 'inputz_cold':
                iz = inputz_cold
            else:
                raise ValueError('Unknown compund name: {}'.format(pname))

            if pname[1].upper() == 'END_TIME':
                iz.append(namelist_time('E', param.value))
            elif pname[1].upper() == 'START_TIME':
                iz.append(namelist_time('I', param.value))
            else:
                iz.append('{}={},'.format(pname[1],param.value))



    # ------- Make the rundir
    try:
        os.makedirs(rundir)
    except OSError:
        pass
    try:
        os.remove(os.path.join(rundir, 'I'))
    except OSError:
        pass

    # -------- Remove old symlinks
    for label, fname in data_files:
        try:
            os.remove(os.path.join(rundir, label))
        except OSError:
            pass

    # -------- Link data files
    for label, fname in data_files:
        os.symlink(fname, os.path.join(rundir, label))

    # Write them out to the I file
    fname = os.path.join(rundir, 'I')
    out = open(fname, 'w')
    out.write(rd.preamble[0])    # First line of preamble
    out.write('\n')

    out.write('&&PARAMETERS\n')
    out.write('\n'.join(parameters))
    out.write('\n')
    out.write('\n'.join(data_lines))
    out.write('\n&&END_PARAMETERS\n')

    out.write('\n&INPUTZ\n')
    out.write('\n'.join(inputz))
    out.write('\n/\n\n')

    out.write('&INPUTZ_cold\n')
    out.write('\n'.join(inputz_cold))
    out.write('\n/\n')
