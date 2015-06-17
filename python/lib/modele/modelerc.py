import os
import bashfile

def read_modelerc(fname = None) :
	"""Reads the settings in the user's .modelErc file.

	Returns:	{string : string}
		Dictionary of the name/value pairs found in the file.

	See:
		giss.bashfile.read_env()"""
		
	if fname is None :
		try:
			fname = os.environ['MODELERC']
		except:
			fname = os.path.join(os.environ['HOME'], '.modelErc')

	return bashfile.read_env(fname)
