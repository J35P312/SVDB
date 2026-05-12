import unittest

from svdb.export_module import db_header, make_representing_variant, build_genotype_columns, vcf_line as make_vcf_line

#mock argeparse arguments
class args:
   version="9000"

class TestExport(unittest.TestCase):

    #test that the header function do not crash
    def test_header(self):
        assert (db_header(args))


class TestMakeRepresentingVariant(unittest.TestCase):

    def test_fields_are_set(self):
        v = make_representing_variant("DEL", "1", "1", 100, 90, 110, 200, 190, 210)
        assert v["type"] == "DEL"
        assert v["chrA"] == "1"
        assert v["chrB"] == "1"
        assert v["posA"] == 100
        assert v["ci_A_start"] == 90
        assert v["ci_A_end"] == 110
        assert v["posB"] == 200
        assert v["ci_B_start"] == 190
        assert v["ci_B_end"] == 210

    def test_precise_variant_has_equal_ci_bounds(self):
        # for a unique variant, ci = pos on both sides
        v = make_representing_variant("INS", "2", "2", 500, 500, 500, 500, 500, 500)
        assert v["ci_A_start"] == v["posA"] == v["ci_A_end"]
        assert v["ci_B_start"] == v["posB"] == v["ci_B_end"]


class TestBuildGenotypeColumns(unittest.TestCase):

    def test_absent_samples_get_ref(self):
        cols = build_genotype_columns(["s1", "s2", "s3"], ["s2"])
        assert cols == ["0/0", "./1", "0/0"]

    def test_all_samples_present(self):
        cols = build_genotype_columns(["s1", "s2"], ["s1", "s2"])
        assert cols == ["./1", "./1"]

    def test_no_hits(self):
        cols = build_genotype_columns(["s1", "s2"], [])
        assert cols == ["0/0", "0/0"]

    def test_order_follows_sample_ids(self):
        cols = build_genotype_columns(["b", "a"], ["a"])
        assert cols == ["0/0", "./1"]


class TestVcfLineStripChr(unittest.TestCase):

    def _cluster(self, chrA, chrB, vtype="DEL"):
        rep = make_representing_variant(vtype, chrA, chrB, 100, 100, 100, 200, 200, 200)
        variants = {0: {"posA": 100, "posB": 200, "sample_id": "s1"}}
        return [rep, variants]

    def test_chr_prefix_preserved_by_default(self):
        line = make_vcf_line(self._cluster("chr1", "chr1"), "id1", ["s1"])
        assert line.startswith("chr1\t")

    def test_strip_chr_removes_prefix(self):
        line = make_vcf_line(self._cluster("chr1", "chr1"), "id1", ["s1"], strip_chr=True)
        assert line.startswith("1\t")

    def test_strip_chr_bnd_chrB(self):
        rep = make_representing_variant("BND", "chr1", "chr2", 100, 100, 100, 500, 500, 500)
        cluster = [rep, {0: {"posA": 100, "posB": 500, "sample_id": "s1"}}]
        line = make_vcf_line(cluster, "id1", ["s1"], strip_chr=True)
        fields = line.split("\t")
        assert fields[0] == "1"
        assert "chr2" not in fields[4]
        assert "2:" in fields[4]

    def test_no_chr_prefix_unaffected(self):
        line = make_vcf_line(self._cluster("1", "1"), "id1", ["s1"], strip_chr=True)
        assert line.startswith("1\t")
