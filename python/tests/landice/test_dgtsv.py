import unittest
from fexception import *
import sys
import modele
import _f90wrap_modele as modelex
from modele import f90
from modele.const import *
import traceback
import numpy as np

class DgtsvTestCase(unittest.TestCase):

	def test_dgtsv(self):
		# See: http://www.nag.co.uk/lapack-ex/node7.html
		n=np.array(5)
		nrhs=np.array(1)		# Number of right-hand sides
		dl = np.array([3.4,3.6,7.0,-6.0])	# sub-diagonal
		d = np.array([3.0,2.3,-5.0,-0.9,7.1])	# diagonal
		du = np.array([2.1,-1.0,1.9,8.0])		# super-diagonal
		b = np.array([2.7,-0.5,2.6,0.6,2.7])
		ldb = 5		# Leading dimension of b (equal to n)
		info = np.array(289)

		#print modelex.dgtsv.__doc__
		modelex.dgtsv(n, nrhs, dl, d, du, b, info,ldb=ldb)

		np.testing.assert_almost_equal(b, np.array([-4.,7.,3.,-4.,-3.]))
