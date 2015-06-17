from modele import rundeckio

import sys
import os
import string
import tempfile
import filecmp
import shutil
import modele.modelerc


# http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
def search_file(filename, search_path):
   """Given a search path, find file
   """
   file_found = 0
   paths = string.split(search_path, os.pathsep)
   for path in paths:
	  if os.path.exists(os.path.join(path, filename)):
		  file_found = 1
		  break
   if file_found:
	  return os.path.abspath(os.path.join(path, filename))
   else:
	  return None
#	   raise IOError('File not fond in $GCMSEARCHPATH: {}'.format(filename))

# --------------------------------------
def make_rundir(rundeck, rundir):
    ret = True

    # ------- Read input
    rc = modele.modelerc.read_modelerc()
    sections = rundeckio.read_rundeck(rundeck)

    sys.stderr.write('GCMSEARCHPATH={}\n'.format(rc['GCMSEARCHPATH']))

    # ------- Make the rundir
    try: os.makedirs(rundir)
    except: pass
    try: os.remove(os.path.join(rundir, 'I'))
    except: pass

    # -------- Find the data files
    data_files = []
    for label,fname0 in sections['Data input files']:
    	fname = search_file(fname0, rc['GCMSEARCHPATH'])
    	data_files.append((label, fname0, fname))

    # -------- Remove old symlinks
    for label, fname0, fname in data_files:
    	if fname is None: continue
    	try: os.remove(os.path.join(rundir, label))
    	except: pass

    # -------- Link data files
    for label, fname0, fname in data_files:
    	if fname is None:
    		sys.stderr.write('Cannot find file: {}\n'.format(fname0))
    		ret = False
    	else:
    		os.symlink(fname, os.path.join(rundir, label))

    # -------- Write I file
    out = open(os.path.join(rundir, 'I'), 'w')
    out.write(sections['preamble'][0])
    out.write('\n')

    out.write('&&PARAMETERS\n')
    for key,val in sections['Parameters']:
    	out.write(' {}={}\n'.format(key,val))
    for label,fname0,fname in data_files:
    	if fname is not None:
    		out.write(" _file_{}='{}'\n".format(label, fname))
    	else:
    		out.write(" # Not Found: {}={}\n".format(label, fname0))
    out.write('&&END_PARAMETERS\n')

    out.write('\n&INPUTZ\n')
    for line in sections['InputZ'][:-2]:
    	out.write(line)
    	out.write('\n')
    out.write('/\n\n')
    out.write('&INPUTZ_cold\n')
    for line in sections['InputZ'][-2:]:
    	out.write(line)
    	out.write('\n')

    return ret
