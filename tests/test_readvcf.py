import unittest

from svdb.readVCF import readVCFLine


class TestReadVCFLine(unittest.TestCase):

    def test_comment(self):
        line = '# This is a comment'
        self.assertIsNone(readVCFLine(line))

    def test_del(self):
        line = 'A\t1\tCNVnator_del_1281\tN\t<DEL>\t.\tPASS\tEND=84000;SVTYPE=DEL;SVLEN=-84000;IMPRECISE;natorRD=0.469202\tGT:CN\t0/1:1'
        assert(readVCFLine(line) == ("A",1,"A",84000,"DEL",{"END":"84000","SVTYPE":"DEL","SVLEN":"-84000","natorRD":"0.469202"},{"GT":["0/1"],"CN":["1"] }) )

    def test_bnd(self):
        line = '18\t34863882\tSV_28046_1\tN\t]X:146727655]ACACACAC\t8\tPASS\tSVTYPE=BND;CIPOS=0,242;CIEND=-143\tGT:CN\t0/1:.'
        assert(readVCFLine(line) == ("18",34863882,"X",146727655,"BND",{"SVTYPE":"BND","CIPOS":"0,242","CIEND":"-143"},{"GT":["0/1"],"CN":["."] }) )

    def test_bnd2(self):
        line = 'X\t146727655\tSV_28046_1\tACACACA\t[18:34863882[N\t8\tPASS\tSVTYPE=BND;CIPOS=0,242;CIEND=-143\tGT:CN\t0/1:.'
        assert(readVCFLine(line) == ("18",34863882,"X",146727655,"BND",{"SVTYPE":"BND","CIPOS":"0,242","CIEND":"-143"},{"GT":["0/1"],"CN":["."] }) )



