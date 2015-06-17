#!/usr/bin/python
# Copyright (C) 2006-2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import re, sys, os, getopt

def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]

def usage():
	print """f90_mod_deps.py [options] <file ...>

Generates dependencies for the given Fortran 90 source files based on
module and use statements in them. Options are:

  -h, --help	 This help output.

  -o, --output <file>
				 Specify the output file. If unspecified then output
				 is to stdout.

  -d, --dep-re <regexp>
				 Regular expression to match against each module name
				 within use statements. Defaults to matching
				 everything.

  -D, --dep-template <template>
				 Template expression for the dependency produced from
				 each module matched by the dep-re regular expression.

  -m, --mod-re <regexp>
				 Regular expression to match against each module name
				 within module definition statements. Defaults to
				 matching everything.

  -M, --mod-template <template>
				 Template expression for the dependency target
				 produced from each module matched by the mod-re
				 regular expression.

  -v, --verbose	 Turn on verbose debugging output.

For a discussion of managing Fortran 90 dependencies see:
http://tableau.stanford.edu/~mwest/group/Fortran_90_Module_Dependencies

Example:
f90_mod_deps.py --output src/myfile.deps --dep-re "(pmc_.*)" \\
	  --dep-template "src/\1.mod" --mod-re "(.*)" \\
	  --mod-template "src/\1.mod" src/myfile.f90
"""

# default options
class Opts:
	output = None
	dep_re = "(.*)"
	dep_template = "\\1.mod"
	mod_re = "(.*)"
	mod_template = "\\1.mod"
	verbose = False

def get_deps_and_mods(filename, opts, base='.'):
	if opts.verbose:
		sys.stderr.write("Processing %s\n" % filename)
	deps = []
	mods = []
	f = open(os.path.join(base, filename))
	if not f:
		print "ERROR: unable to open %s%s" % filename
		sys.exit(1)
	use_line_re = re.compile("^\s*use\s+(\S.+)\s*$")
	cont_line_re = re.compile("^(.*)&\s*$")
	mod_line_re = re.compile("^\s*module\s+(\S+)\s*$")
	split_re = re.compile("\s*,\s*")
	dep_re = re.compile(opts.dep_re)
	mod_re = re.compile(opts.mod_re)
	within_use_statement = False
	line_with_use = False
	for line in f:
		line = line.lower()
		match = use_line_re.search(line)
		if match:
			within_use_statement = True
			rest_line = match.group(1)
		else:
			rest_line = line
		if within_use_statement:
			match = cont_line_re.search(rest_line)
			if match:
				rest_line = match.group(1)
			else:
				within_use_statement = False
			line_items = split_re.split(rest_line.strip())
			for item in line_items:
				if item:
					if opts.verbose:
						sys.stderr.write("use: %s\n" % item)
					match = dep_re.match(item)
					if match:
						dep = match.expand(opts.dep_template)
						if opts.verbose:
							sys.stderr.write("matched to: %s\n" % dep)
						if dep not in deps:
							deps.append(dep.strip())
		else:
			# not within_use_statement
			match = mod_line_re.search(line)
			if match:
				mod_name = match.group(1)
				if opts.verbose:
					sys.stderr.write("module: %s\n" % mod_name)
				match = mod_re.match(mod_name)
				if match:
					mod = match.expand(opts.mod_template).lower()
					if opts.verbose:
						sys.stderr.write("matched to: %s\n" % mod)
					if mod not in mods:
						mods.append(mod.strip())
	f.close()
	return (deps, mods)

def get_alldeps(filenames, base='.'):
	allmods = set()
	alldm = list()

	opts = Opts()		# Fake out bad code we inherited that expects to only run on command line.
	for filename in filenames:
		deps, mods = get_deps_and_mods(filename, opts, base=base)
		for mod in mods :
			allmods.add(mod)
		alldm.append((filename, deps, mods))

	depRE = re.compile(r'\s*(.*(:|=>))?(.*)\.mod\s*')
	for filename, deps, mods in alldm:

		# Fix up parsing of deps
		deps2 = set()
		for dep in deps:
			match = depRE.match(dep)
			dep = match.group(3)
			deps2.add(dep)

		localdeps = set()
		for dep in deps2 :
			if dep in allmods :			# Only depend on things WE generate
				if dep not in mods :	# No circular dependencies
					localdeps.add(dep)

		yield (filename, sorted(list(deps2)),
			sorted([x[:-4] for x in mods]))


# ----------------------------------------
# Remove common path prefix from filenames
fnames = sys.argv[1:]
common_prefix = os.path.commonprefix(fnames)
fnames = [x[len(common_prefix):] for x in fnames]



# ----------------------------------
def get_abstract_fname(fname):
	dir,name = os.path.split(fname)
	if len(dir) > 0:
		abstract_fname = dir	# Sources in directories are included wholesale by the directory name
	else:
		abstract_fname = os.path.splitext(name)[0]
	return abstract_fname
# ----------------------------------
findREs = []


for fname in fnames:
	dir,name = os.path.split(fname)
	if len(dir) > 0:
		abstract_fname = dir	# Sources in directories are included wholesale by the directory name
	else:
		abstract_fname = os.path.splitext(name)[0]
#	print 'abstract_fname', abstract_fname

	findREs.append(re.compile('('+abstract_fname+')', re.MULTILINE))

#fnames_todo = fnames
#
#
## Construct regular expressions to grep through...
#while len(fnames_todo) > 0:
#	fnames_thistime = fnames_todo[:50]
#	fnames_todo = fnames_todo[50:]
#
#	# Find files in the top-level templates/ directory that use the source files identified
#	toplevel_abstract = []
#	for fname in fnames_thistime:
#		dir,name = os.path.split(fname)
#		abstract_fname = os.path.splitext(name)[0]
#		print 'abstract_fname', abstract_fname
#		toplevel_abstract.append(abstract_fname)
#
#	re_str = '|'.join(['('+x+')' for x in toplevel_abstract])
#	print re_str
#	findREs.append(re.compile(re_str, re.MULTILINE))

used_abstract_fnames = set()
# List through templates directory, finding abstract fnames we've used
for template_leaf in os.listdir('templates'):
	# Read the template file
	try:
		with open(os.path.join('templates', template_leaf)) as content_file:
			template = content_file.read()
	except:		# Don't worry if we can't open directories, etc.
		continue

	# Grep through it for mentions of abstract fnames
	for findRE in findREs:
		for x in findRE.findall(template):
			xx = x
			if len(xx) == 0: continue
			used_abstract_fnames.add(xx)

# Remove obsolete files (not used in any rundeck)
good_fnames = set(('LIGlint2.F90','LIGrid.F90','LISNOWBASE.f','LISheet.F90','LISheetDummy.F90','LISnoGli.F90','LISnow.F90'))
bad_fnames = set()

for fname in fnames:
	abstract_fname = get_abstract_fname(fname)

	if fname in good_fnames:
		pass	# Already whitelisted
	elif abstract_fname in used_abstract_fnames:
		good_fnames.add(fname)
	else:
		bad_fnames.add(fname)

#print '--------- Good fnames'
#print '\n'.join(sorted(list(good_fnames)))
#print '--------- Bad (obsolete) fnames'
#print '\n'.join(sorted(list(bad_fnames)))

# ----------------------------------
# Set up table (bymods) of which fnames define each module

bymods = dict()
mods_in_fname = dict()
for fname,deps,mods in get_alldeps(good_fnames, base=common_prefix):
	mods_in_fname[fname] = mods
	for mod in mods:
		if mod in bymods:
			bymods[mod].append(fname)
		else:
			bymods[mod] = [fname]

# -----------------------------------------
# Merge libraries together to create alternative packages

# Libraries we will merge together.
# new : (old libraries)
mergelibs_rev = {
	('icedyn','ICEDYN') : (('icedyn','ICEDYN'), ('icedyn_com','ICEDYN_DRV')),
	('icedyn','ICEDYN_DUM') : (('icedyn','ICEDYN_DUM'), ('icedyn_com','ICEDYN_DUM_DRV'), ('icedyn_com','ICEDYN_DUM')),
	('icedyn','ICEDYNo_DUM') : (('icedyn','ICEDYNo_DUM'), ('icedyn_com','ICEDYNo_DUM_DRV'), ('icedyn_com','ICEDYNo_DUM')),


	('geom_cs_support','') : (('int_ag2og_mod','OCN_Int_CS'), ('int_og2ag_mod','OCN_Int_CS'), ('rad_cosz0','COSZ_2D')),

	# COSZ_2D is repeated here.  That eliminates need to add it if GEOM_CS or GNOM_CS is linked.
	('geom','GEOM_B') : (('rad_cosz0','GEOM_B'),('int_ag2og_mod', 'OCN_Int_LATLON'), ('int_og2ag_mod', 'OCN_Int_LATLON')),
	('geom','GEOM_CS') : (),		# Link ('geom_cs_support','') if you use this library
	('geom','GNOM_CS') : (),		# Link ('geom_cs_support','') if you use this library


	('giss_oturb', 'OCNGISSVM') : (('gissmix_com', 'OCNGISSVM'),),
	('giss_oturb', 'OCNGISS_TURB') : (('gissmix_com', 'OCNGISS_TURB'),),

	# No 'geom', 'domain_decomp_atm' or 'horizontalres' alternative needed if ('soatm', 'yes') is linked.
	('soatm', 'yes') :
		(('diag_com', 'SOATM_COM'), ('fluxes', 'SOATM_COM'), ('geom', 'SOATM_COM'), ('domain_decomp_atm', 'SOATM_COM'),
		('rad_com', 'SOATM_DRV'), ('surf_albedo', 'SOATM_DRV')),

	('soatm', 'no') :
		(('diag_com', 'DIAG_COM'), ('fluxes', 'FLUXES'),
		('rad_com', 'RAD_COM'), ('surf_albedo', 'ALBEDO')),


	# If you choose ('scm','yes') then you don't need 'geom' and 'horizontalres' alternatives.
	('scm', 'yes') :
		(('horizontalres', 'SCM_COM'), ('geom', 'SCM_COM'))

}


mergelibs = {}
for dest, srcs in mergelibs_rev.items():
	for src in srcs:
		if src in mergelibs:
			mergelibs[src].append(dest)
		else:
			mergelibs[src] = [dest]

# ----------------------------------


# Examine files that require alternatives libraries (i.e. can't be linked together)
# We will start here by creating one library per module declaration per source file
libs = {}
for mod,fnames in bymods.items():
	if len(fnames) > 1:
		opt_name = mod
		for fname in fnames:
			dir,fn = os.path.split(fname)
			base,ext = os.path.splitext(fn)
			orig_lib_name = (opt_name, base)

			if orig_lib_name in mergelibs:
				iter = mergelibs[orig_lib_name]
				print 'rename',orig_lib_name,'->',iter
			else:
				iter = (orig_lib_name,)

			# Put this source file in possibly more than one library.
			for lib_name in iter:
				if lib_name in libs:
					libs[lib_name].add(fname)
				else:
					libs[lib_name] = set([fname])


# Find files in more than one library
libs_for_file = {}
for lib_name,fnames in libs.items():
	for fname in fnames:
		if fname in libs_for_file:
			libs_for_file[fname].append(lib_name)
		else:
			libs_for_file[fname] = [lib_name]

print '---------- Files in more than one "original" library (should be none)'
for fname, lib_names in libs_for_file.items():
	abstract_fname = get_abstract_fname(fname)
	llen = 0
	for lib_name in lib_names:
		if lib_name[1] == abstract_fname: llen += 1

	if llen > 1:
		print fname,lib_names


# ----------- Put all files not yet in any library into a "base" library
base_fnames= []
for fname in good_fnames:
	if fname not in libs_for_file:
		base_fnames.append(fname)
libs[('base','')] = base_fnames


print '---------- Libraries with conflicting source files'
llibs = sorted(list(libs.items()))
for lib_name,fnames in llibs:
	mods_this_lib = set()
	err = False
	for fname in fnames:
		for mod in mods_in_fname[fname]:
			if mod in mods_this_lib:
				err = True
			mods_this_lib.add(mod)
	if err:
		print lib_name,fnames,sorted(list(mods_this_lib))

print '---------- Library Definitions'
for lib_name, fnames in llibs:
	if len(lib_name[1]) == 0:
		slib_name = 'lib{}.so'.format(lib_name[0])
	else:
		slib_name = '{}/lib{}.so'.format(*lib_name)

	sfnames = ' '.join(fnames)

	print '{} : {}'.format(slib_name, sfnames)


sys.exit(0)

def write_deps(outf, filename, deps, mods):
	filebase, fileext = os.path.splitext(filename)
	dir,leafbase = os.path.split(filebase)
	outf.write('%s : %s.lo\n' % (' '.join(mods), leafbase))
	outf.write('%s.lo : %s %s\n' % (leafbase, filename, ' '.join(deps)))


# 	outf.write("%s%s%s.lo : %s %s\n" %
# 		(" ".join(mods),
# 		("" if len(mods) == 0 else " "),
# 		leafbase, filename, " ".join(deps)))
# 	outf.write("\t$(LTFCCOMPILE) -c -o %s.lo $<\n\n" %
# 		leafbase)

def process_args():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ho:d:D:m:M:v",
								   ["help", "output=", "dep-re=",
									"dep-template=", "mod-re=",
									"mod-template=", "verbose"])
	except getopt.GetoptError:
		print "ERROR: invalid commandline options"
		usage()
		sys.exit(1)
	myopts = Opts()
	for o, a in opts:
		if o in ("-h", "--help"):
			usage()
			sys.exit()
		if o in ("-o", "--output"):
			myopts.output = a
		if o in ("-d", "--dep-re"):
			myopts.dep_re = a
		if o in ("-D", "--dep-template"):
			myopts.dep_template = a
		if o in ("-m", "--mod-re"):
			myopts.mod_re = a
		if o in ("-M", "--mod-template"):
			myopts.mod_template = a
		if o in ("-v", "--verbose"):
			myopts.verbose = True
			sys.stderr.write("Verbose output on\n")
	if len(args) < 1:
		usage()
		sys.exit(1)
	if myopts.verbose:
		sys.stderr.write("output = %s\n" % myopts.output)
		sys.stderr.write("dep-re = %s\n" % myopts.dep_re)
		sys.stderr.write("dep-template = %s\n" % myopts.dep_template)
		sys.stderr.write("mod-re = %s\n" % myopts.mod_re)
		sys.stderr.write("mod-template = %s\n" % myopts.mod_template)
	return (myopts, args)

def main():
	(opts, filenames) = process_args()
	if opts.output:
		outf = open(opts.output, "w")
		if opts.verbose:
			sys.stderr.write("Output to %s\n" % opts.output)
	else:
		outf = sys.stdout
		if opts.verbose:
			sys.stderr.write("Output to STDOUT\n")
	outf.write("# DO NOT EDIT --- auto-generated file\n")
	alldeps = list()
	allmods = set()
	for filename in filenames:
		(deps, mods) = get_deps_and_mods(filename, opts)
#		print filename, deps, mods
		for mod in mods :
			allmods.add(mod)
		alldeps.append((filename, deps, mods))
#	print '--------------------------'

	for (filename, deps, mods) in alldeps :
		if opts.verbose:
			sys.stderr.write("deps: %s\n" % " ".join(deps))
			sys.stderr.write("mods: %s\n" % " ".join(mods))
#		if deps:
		if True :
			localdeps = list()
			for mod in deps :
				if mod in allmods :			# Only depend on things WE generate
					if mod not in mods :	# No circular dependencies
						localdeps.append(mod)
			write_deps(outf, filename, localdeps, mods)
	if outf != sys.stdout: outf.close()

if __name__ == "__main__":
	main()
								   

