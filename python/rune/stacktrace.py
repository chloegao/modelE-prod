import sys
import re
import subprocess
import StringIO

log_fname = sys.argv[1]
exe_fname = sys.argv[2]


# ------------ Parse the log
ref_addrRE = re.compile(r'^REFERENCE_ADDRESS (.*?) (0x[0-9a-fA-F]+)')
tracelineRE = re.compile(r'#\d*\s*(0x[0-9a-fA-F]+)')

stacktrace_r = []

with open(log_fname, 'r') as fin:
	for line in fin:
		if line.startswith('REFERENCE_ADDRESS'):
			match = ref_addrRE.match(line)
			if match is None: continue
			ref_symbol = match.group(1)
			ref_addr_r = int(match.group(2), 0)
			continue

		if line.startswith('Backtrace for this error'):
			break
		if line.startswith('Program aborted. Backtrace:'):
			break

	# ---------- Read the numeric stacktrace
	for line in fin:
		match = tracelineRE.match(line)
		if match is None: continue
		stacktrace_r.append(int(match.group(1), 0))

# ------------- Get symbols out of the binary file
# - Compute the ref_offset from this, used to convert from
symtab_f = dict()

nmRE = re.compile(r'([0-9a-fA-F]+)\s+T\s+_(.+)?\s*')
cmd = ['nm', exe_fname]
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
out,err = proc.communicate()

for line in StringIO.StringIO(out):
	match = nmRE.match(line)
	if match is None: continue
	#print match.group(1), match.group(2)
	if match.group(2) == ref_symbol:
		ref_offset = int(match.group(1), 16) - ref_addr_r
		break

# ------------- Relocate the stacktrace to file coordinates and stringify
#print '\n'.join(['%x' % (x + 0) for x in stacktrace_r])
#print '--------'
stacktrace_fs = '\n'.join(['%x' % (x + ref_offset) for x in stacktrace_r])
#print stacktrace_fs

# ------------------ Convert numeric coordinates
cmd = ['xcrun', 'atos', '-o', exe_fname]
proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
ret = proc.communicate(stacktrace_fs)[0]
print ret

#print 'ref_offset', ref_offset
