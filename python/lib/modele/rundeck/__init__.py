from modele.rundeck import legacy
import datetime
import os
import sys
from modele import pathutil
import copy
import urllib2
from modele import xhash

# parameter types
GENERAL = 'GENERAL'
FILE = 'FILE'
DATETIME = 'DATETIME'

class Param(object):
    def __init__(self, pname, type, value, sval=None):
        self.pname = pname
        self.type = type
        self.value = value
        self.sval = sval    # Raw unparsed value, if appropriate

        if id(self.type) == id(DATETIME):
            dt = self.value
            if not isinstance(dt, datetime.datetime):
                raise ValueError('Values of type DATETYPE must have Python type datetime.datetime')
            if (dt.minute!=0) or (dt.second !=0) or (dt.microsecond!=0):
                raise ValueError('Values of type DATETYPE must be on the hour.  Error in: {}'.format(dt))
            if (dt.tzinfo is not None):
                raise ValueError('Values of type DATETYPE cannot have a timezone.  Error in: {}'.format(dt))

    def __lt__(self, other):
        return self.pname < other.pname

    def __repr__(self):
        return repr((self.pname, self.type, self.value))

    def sname(self):
        return '.'.join(pname)

def replace_date(dd, suffix, result):
    try:
        yeari = int(dd['YEAR'+suffix])
        monthi = int(dd['MONTH'+suffix])
        datei = int(dd['DATE'+suffix])
        houri = int(dd['HOUR'+suffix])
        dd[result] = datetime.datetime(yeari, monthi, datei, houri, 0, 0)
        del dd['YEAR'+suffix]
        del dd['MONTH'+suffix]
        del dd['DATE'+suffix]
        del dd['HOUR'+suffix]
    except KeyError:
        pass



class Params(dict):
    def __init__(self, file_path, auto_download=False):
        """auto_download:
            Download files we can't find locally?"""
        self.file_path = file_path
        self.auto_download = auto_download

    def add(self,param):
        self[param.pname] = param

    def __repr__(self):
        return '\n'.join([repr(x) for x in sorted(list(self.values()))])
        

    def download_file(self, sval):
        # Only try to download if it's a raw leafname
        if len(os.path.split(sval)[0]) > 0:
            return None

        # Try to download the file
        # http://stackoverflow.com/questions/22676/how-do-i-download-a-file-over-http-using-python
        file_name = os.path.join(self.file_path[0], sval)
        url = 'http://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/' + sval

        with open(file_name, 'wb') as fout:
            u = urllib2.urlopen(url)
            meta = u.info()
            file_size = int(meta.getheaders("Content-Length")[0])
            print "Downloading: %s Bytes: %s" % (file_name, file_size)

            file_size_dl = 0
            block_sz = 8192
            while True:
                buffer = u.read(block_sz)
                if not buffer:
                    break

                file_size_dl += len(buffer)
                fout.write(buffer)
                status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
                status = status + chr(8)*(len(status)+1)
                print status,
        return file_name


    def set_file(self, symbol, fname):
        ret = True
        try:
            fname_full = pathutil.search_file(fname, self.file_path)
        except IOError as e:
#            try:
#                fname_full = self.download_file(fname)
#            except Exception as e2:
#                print(e2)
                sys.stderr.write('{}: {}\n'.format(symbol, e))
                ret = False    # Error condition
                fname_full = None

        self.add(Param((symbol,), FILE, fname_full, sval=fname))
        return ret


    def add_legacy(self, legacy):
        """Extract rundeck parametesr from a legacy rundeck."""
        ret = True
        for symbol,fname in legacy['Data input files']:
                ret = ret and self.set_file(symbol, fname)


        for symbol,value in legacy['Parameters']:
            param = Param((symbol,), GENERAL, value)
            self.add(param)

        # ------- Deal with the namelists

        # Split into a series of namelists, splitting on 'ISTART=...'
        inputz = legacy['InputZ']
        inputzs = list()
        inputz_cur = list()
        for item in inputz:
            if item[0].upper() == 'ISTART':
                if len(inputz_cur) > 0:
                    inputzs.append(dict(inputz_cur))
                inputz_cur = list()
            inputz_cur.append((item[0].upper(), item[1]))
        inputzs.append(dict(inputz_cur))
    
        for inputz in inputzs:
            replace_date(inputz, 'I', 'START_TIME')
            replace_date(inputz, 'E', 'END_TIME')

        if len(inputzs) > 2:
            raise ValueError('At most one ISTART line is allowed')

        prefixes = ('INPUTZ', 'INPUTZ_cold')
        for prefix,inputz in zip(prefixes,inputzs):
            for symbol,value in inputz.items():
                type = DATETIME if isinstance(value, datetime.datetime) else GENERAL
                self.add(Param((prefix,symbol), type, value))

        return ret
# ------------------------------------------
class Build(object):
    def __init__(self):
        self.sources = set()    # Object Modules
        self.components = dict()    # Directories of sources --> options

        self.defines = dict()    # Preprocessor Options

    def update_hash(self, hash):
        xhash.update(self.sources, hash)
        xhash.update(self.components, hash)
        xhash.update(self.defines, hash)

    def add_legacy(self, legacy):
        for src in legacy['Object Modules']:
            self.sources.add(src)
        for symbol,value in legacy['Preprocessor Options']:
            self.defines[symbol] = value
        for component in legacy['Components']:
            self.components[component] = None

        for component,options in legacy['Component Options']:
            if component not in self.components:
                raise ValueError('Options found for non-existant component %s' % component)
            self.components[component] = dict(options)



# ------------------------------------------
class ChangeSysPath(object):
    def __init__(self, new_path):
        self.new_path = new_path
    def __enter__(self):
        self.old_path = sys.path
        sys.path = self.new_path
    def __exit__(self, type, value, traceback):
        sys.path = self.old_path

def load_rundeck(fname, template_path=None, file_path=None, auto_download=True):
    if template_path is None:
        template_path = copy.copy(default_template_path)

    # Resolve the rundeck filename
    fname = pathutil.search_file(fname, default_template_path)

    # Add the directory of the rundeck to the path
    dirname,leafname = os.path.split(fname)
    template_path = [dirname] + template_path

    # Create a blank rundeck
    rd = Rundeck(template_path=template_path, auto_download=auto_download)

    root,ext = os.path.splitext(leafname)
    if ext == '.py':
        # Load Python-format rundeck
        with ChangeSysPath(rd.template_path + sys.path):
            globals = dict()
            rd_code = __import__(root, globals)
            rd_code.setup(rd)
    else:
        # Load legacy-format rundeck
        rd.load_legacy(fname)

    return rd


# ------------------------------------------
try:
    MODELE_TEMPLATE_PATH = os.environ['MODELE_TEMPLATE_PATH'].split(os.pathsep)
except:
    MODELE_TEMPLATE_PATH = []

default_template_path = MODELE_TEMPLATE_PATH + [os.path.join(pathutil.modele_root(), 'templates')]

# Search for input files
try:
    default_file_path = os.environ['MODELE_FILE_PATH'].split(os.pathsep)
except Exception as e:
    default_file_path = ['.']
# ------------------------------------------

class Rundeck(object):
    def __init__(self, template_path=None, file_path=None, auto_download=False):
        if template_path is None:
            self.template_path = default_template_path
        else:
            self.template_path = template_path
        if file_path is None:
            self.file_path = default_file_path
        else:
            self.file_path = file_path

        self.preamble = None
        self.params = Params(self.file_path, auto_download=auto_download)
        self.build = Build()

    def update_hash(self, hash):
        xhash.update(self.build, hash)

    def set(self, name, type, value):
        if id(type) == id(FILE):
            self.params.set_file(name, value)
        else:
            return self.params.add(Param((name,), type, value))

    def add_legacy(self, legacy):
        self.preamble = legacy['preamble']
        self.params.add_legacy(legacy)
        self.build.add_legacy(legacy)

    def load_legacy(self, fname, template_path=None, file_path=None):
        if template_path is None:
            template_path = default_template_path
        if file_path is None:
            file_path = default_file_path

        fname_full = pathutil.search_file(fname, template_path)

        fin = legacy.preprocessor(fname_full, template_path)
        legacy_rundeck = legacy.read_rundeck(fin)       # Auto-closes
        self.add_legacy(legacy_rundeck)

    def __repr__(self):
        return '\n'.join(('============= Rundeck',
            '--------- Preamble', \
            repr(self.preamble),
            '--------- Params', \
            repr(self.params),
            '-------- Sources', \
            repr(self.build.sources),
            '-------- Components', \
            repr(self.build.components),
            '-------- Defines', \
            repr(self.build.defines)))
