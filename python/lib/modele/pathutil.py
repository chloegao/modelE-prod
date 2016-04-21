import os

# http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
def search_file(filename, search_path):
    """Given a search path, find file
    """
    if os.path.exists(filename):
        return os.path.abspath(filename)

    for path in search_path:
        fname = os.path.abspath(os.path.join(path, filename))
        if os.path.exists(fname):
            return fname
    raise IOError('File not found in search path: {}'.format(filename))


# Returns the root of this ModelE installation
def modele_root():
    dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.realpath(os.path.join(dir, '..', '..', '..'))
