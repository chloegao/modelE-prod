from modele import xhash,rundeck
import hashlib
import unittest
import os

srcdir = os.path.dirname(os.path.abspath(__file__))

class TestHash(unittest.TestCase):

    def test_simple_hash(self):
        """Just tests that the hash algorithm runs"""
        hash = hashlib.md5()
        xhash.update({'a' : 5, 'b' : 4}, hash)
        xhash.update(17, hash)
        xhash.update(('foo', 'bar'), hash)

    def test_hash_rundeck(self):
        hash1a = xhash.hexdigest(rundeck.load_rundeck(os.path.join(srcdir, 'rundeck1a.R')))
        hash1b = xhash.hexdigest(rundeck.load_rundeck(os.path.join(srcdir, 'rundeck1b.R')))
        hash2 = xhash.hexdigest(rundeck.load_rundeck(os.path.join(srcdir, 'rundeck2.R')))

        self.assertEqual(hash1a, hash1b)
        self.assertNotEqual(hash1a, hash2)


if __name__ == "__main__":
    unittest.main()
