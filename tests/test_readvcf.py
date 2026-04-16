import unittest

from svdb.read_vcf import readVCFLine


class TestReadVCFLine(unittest.TestCase):

    def test_comment(self):
        line = '# This is a comment'
        self.assertIsNone(readVCFLine(line))

    def test_del(self):
        line = 'A\t1\tCNVnator_del_1281\tN\t<DEL>\t.\tPASS\tEND=84000;SVTYPE=DEL;SVLEN=-84000;IMPRECISE;natorRD=0.469202\tGT:CN\t0/1:1'
        v = readVCFLine(line)
        assert v.chrA == "A"
        assert v.posA == 1
        assert v.chrB == "A"
        assert v.posB == 84000
        assert v.event_type == "DEL"
        assert v.info == {"END": "84000", "SVTYPE": "DEL", "SVLEN": "-84000", "natorRD": "0.469202"}
        assert v.fmt == {"GT": ["0/1"], "CN": ["1"]}

    def test_bnd(self):
        line = '18\t34863882\tSV_28046_1\tN\t]X:146727655]ACACACAC\t8\tPASS\tSVTYPE=BND;CIPOS=0,242;CIEND=-143\tGT:CN\t0/1:.'
        v = readVCFLine(line)
        assert v.chrA == "18"
        assert v.posA == 34863882
        assert v.chrB == "X"
        assert v.posB == 146727655
        assert v.event_type == "BND"
        assert v.info == {"SVTYPE": "BND", "CIPOS": "0,242", "CIEND": "-143"}
        assert v.fmt == {"GT": ["0/1"], "CN": ["."]}

    def test_bnd2(self):
        line = 'X\t146727655\tSV_28046_1\tACACACA\t[18:34863882[N\t8\tPASS\tSVTYPE=BND;CIPOS=0,242;CIEND=-143\tGT:CN\t0/1:.'
        v = readVCFLine(line)
        assert v.chrA == "18"
        assert v.posA == 34863882
        assert v.chrB == "X"
        assert v.posB == 146727655
        assert v.event_type == "BND"



