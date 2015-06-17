import re
import sys

# Python package to read and parse rundecks.

# ----------------------------------------------------------

# Finds the division between preprocessor sections
sectionRE = re.compile(r'\s*((Preamble:?)|(Preprocessor\s+Options:?)|(Run\s+Options:?)|(Object\s+modules:?)|(Components:?)|(Component\s+Options:?)|(Data\s+input\s+files:?)|(&&PARAMETERS)|(&INPUTZ))\s*')
num_sections = 9

def parse_section(fin, parse_line):
	"""Rundecks are grouped in sections.  This parses the current section, and
	reports the section type of the next section."""
	while True:
		line = next(fin)	# Pass along EOFException

		match = sectionRE.match(line)
		if match is None:
			parse_line(line)
		else:
			groups = match.groups()[1:]		# Skip the top-level group
			for i in xrange(0,len(groups)):
				if groups[i] is not None:
					return i		# Next section
			raise ValueError('Could not find next section for ' + line)
# ----------------------------------------------------------
# Parsers for the different sections
class Parser(object):
	def __init__(self):
		self.section = list()

class CopyLines(Parser):
	def __call__(self, line):
		self.section.append(line)

class CopyLinesNoComments(Parser):
	def __call__(self, line):
		exp = line.find('!')
		if (exp >= 0): line = line[:exp]
		line = line.strip()
		if len(line) > 0:
			self.section.append(line)



preprocessorRE = re.compile(r'\s*#define\s+.*')
class PreprocessorOptions(Parser):
	def __call__(self, line):
		exp = line.find('!')
		if (exp >= 0): line = line[:exp]

		if preprocessorRE.match(line) is not None:
			self.section.append(line)

key_eq_valueRE = re.compile(r'\s*([^=\s]*)\s*=\s*([^!]*).*')
class KeyEqValue(Parser):
	def __call__(self, line):
		exp = line.find('!')
		if (exp >= 0): line = line[:exp]

		match = key_eq_valueRE.match(line)
		if match is not None:
			self.section.append((match.group(1), match.group(2).strip()))

class ComponentOptions(Parser):
	def __call__(self, line):
		exp = line.find('!')
		if (exp >= 0): line = line[:exp]

		match = key_eq_valueRE.match(line)
		if match is not None:
			scomp = match.group(1)
			options = match.group(2).strip()
			parsed_options = []

			if scomp.startswith('OPTS_'):
				component = scomp[5:]
			else:
				component = scomp

			for opt in options.split(' '):
				if len(opt) == 0: continue
				words = opt.split('=')
				if len(words) != 2:
					raise ValueError('Bad component option {}'.format(opt))
				parsed_options.append((words[0].strip(), words[1].strip()))
			self.section.append((component, tuple(parsed_options)))



class InputZ(Parser):
	def __call__(self, line):
		exp = line.find('!')
		if (exp >= 0): line = line[:exp]

		for ll in line.split(','):
			match = key_eq_valueRE.match(ll)
			if match is not None:
				self.section.append((match.group(1), match.group(2).strip()))

class WordList(Parser):
	def __call__(self, line):
		exp = line.find('!')
		if (exp >= 0): line = line[:exp]

		words = line.split(' ')
		for word in words:
			word = word.strip()
			if len(word) > 0: self.section.append(word)
# ----------------------------------------------------------


def read_rundeck(fname):
	buf = []
	section = 0

	section_parsers = [
		('preamble', CopyLines()),
		('Preprocessor Options', PreprocessorOptions()),
		('Run Options', KeyEqValue()),
		('Object Modules', WordList()),
		('Components', WordList()),
		('Component Options', ComponentOptions()),
		('Data input files', KeyEqValue()),
		('Parameters', KeyEqValue()),
		('InputZ', CopyLinesNoComments())
	]

	with open(fname) as fin:
		try:
			while True:
				section = parse_section(fin, section_parsers[section][1])
		except StopIteration:
			pass

	ret = dict()
	for parser in section_parsers:
		ret[parser[0]] = parser[1].section
	return ret

#with open(sys.argv[1], 'r') as fin:
#	sections = read_rundeck(fin)
#
#print sections


