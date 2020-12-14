import unittest

from svdb.readVCF import readVCFLine


class TestReadVCFLine(unittest.TestCase):

    def test_comment(self):
        line = '# This is a comment'
        self.assertIsNone(readVCFLine(line))
