import unittest
import numpy

from svdb.export_module import db_header

#mock argeparse arguments
class args:
   version="9000"

class TestExport(unittest.TestCase):

    #test that the header function do not crash
    def test_header(self):
        assert (db_header(args))
