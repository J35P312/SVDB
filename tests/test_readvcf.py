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

    # --- Delly TRA format (lines 40-47) ---

    def test_delly_tra_no_swap(self):
        """Delly TRA where chrA < chrB — no chromosome swap needed."""
        line = '1\t1000\t.\tN\t<TRA>\t.\tPASS\tSVTYPE=BND;CHR2=2;END=2000\tGT\t0/1'
        v = readVCFLine(line)
        assert v.chrA == "1"
        assert v.posA == 1000
        assert v.chrB == "2"
        assert v.posB == 2000
        assert v.event_type == "BND"

    def test_delly_tra_swap(self):
        """Delly TRA where chrA > chrB — chromosomes must be reordered."""
        line = '2\t3000\t.\tN\t<TRA>\t.\tPASS\tSVTYPE=BND;CHR2=1;END=1000\tGT\t0/1'
        v = readVCFLine(line)
        # After swap chrA becomes "1", chrB becomes "2"
        assert v.chrA == "1"
        assert v.chrB == "2"
        assert v.event_type == "BND"

    # --- Intrachromosomal paths (lines 55-88) ---

    def test_svlen_fallback(self):
        """When INFO has SVLEN but not END, posB is derived from SVLEN."""
        line = '1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVLEN=-500\tGT\t0/1'
        v = readVCFLine(line)
        assert v.chrA == "1"
        assert v.posB == 1500  # posA + abs(SVLEN)
        assert v.event_type == "DEL"

    def test_fermikit_inverted_coords(self):
        """Fermikit can emit variants where END < POS; posA/posB must be swapped."""
        line = '1\t2000\t.\tN\t<DEL>\t.\tPASS\tEND=1500;SVTYPE=DEL\tGT\t0/1'
        v = readVCFLine(line)
        assert v.posA == 1500
        assert v.posB == 2000

    def test_nucleotide_alt_ins(self):
        """Alt is a raw nucleotide sequence longer than ref → INS, treated as single point."""
        line = '1\t1000\t.\tA\tATGCATGC\t.\tPASS\t.\tGT\t0/1'
        v = readVCFLine(line)
        assert v.event_type == "INS"
        assert v.posA == 1000
        assert v.posB == 1000  # INS is collapsed to single point

    def test_nucleotide_alt_del(self):
        """Alt is a raw nucleotide sequence shorter than ref → DEL."""
        line = '1\t1000\t.\tATGCATGC\tA\t.\tPASS\t.\tGT\t0/1'
        v = readVCFLine(line)
        assert v.event_type == "DEL"
        assert v.posA == 1000
        assert v.posB == 1007  # posA + len(ref) - 1

    def test_dup_normalised_from_subtype(self):
        """DUP:TANDEM and similar subtypes should be normalised to DUP."""
        line = '1\t1000\t.\tN\t<DUP:TANDEM>\t.\tPASS\tEND=5000;SVTYPE=DUP\tGT\t0/1'
        v = readVCFLine(line)
        assert v.event_type == "DUP"

    # --- BND breakpoint notation edge cases (lines 91-115) ---

    def test_bnd_double_bracket_normalisation(self):
        """BND alt starting with [[ is normalised to [ before parsing."""
        # [[chr2:5000[N  →  [chr2:5000[N  after normalisation
        line = '1\t1000\t.\tN\t[[2:5000[N\t.\tPASS\tSVTYPE=BND\tGT\t0/1'
        v = readVCFLine(line)
        assert v.chrA == "1"
        assert v.chrB == "2"
        assert v.posB == 5000
        assert v.event_type == "BND"

    def test_bnd_intrachromosomal_posB_lt_posA(self):
        """Intrachromosomal BND where posB < posA must be swapped."""
        # BND on same chromosome, but breakpoint B is upstream of A
        line = '1\t5000\t.\tN\t]1:1000]N\t.\tPASS\tSVTYPE=BND\tGT\t0/1'
        v = readVCFLine(line)
        assert v.chrA == "1"
        assert v.chrB == "1"
        assert v.posA == 1000
        assert v.posB == 5000

