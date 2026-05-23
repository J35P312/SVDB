import subprocess
import sys
import unittest
from pathlib import Path

from svdb.export_module import db_header, make_representing_variant, build_genotype_columns, vcf_line as make_vcf_line

#mock argeparse arguments
class args:
   version="9000"

class args_no_samples:
   version="9000"
   samples="off"

FIXTURES = Path(__file__).parent / "fixtures"
MANTA_CHR1  = FIXTURES / "manta_chr1_del.vcf"
HG002_CHR1  = FIXTURES / "dragen_hg002_chr1_del.vcf"
SVDB = [sys.executable, "-m", "svdb"]


def _chrom_header(text: str) -> str:
    for line in text.splitlines():
        if line.startswith("#CHROM"):
            return line
    return ""


def _data_lines(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln and not ln.startswith("#")]

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


class TestNoSamplesFlag(unittest.TestCase):
    """--samples off: sites-only output — no FORMAT or GT columns."""

    def _cluster(self):
        rep = make_representing_variant("DEL", "1", "1", 100, 100, 100, 200, 200, 200)
        variants = {
            0: {"posA": 100, "posB": 200, "sample_id": "s1"},
            1: {"posA": 102, "posB": 198, "sample_id": "s2"},
        }
        return [rep, variants]

    def test_vcf_line_no_samples_has_eight_fields(self):
        line = make_vcf_line(self._cluster(), "id1", ["s1", "s2"], no_samples=True)
        assert len(line.split("\t")) == 8

    def test_vcf_line_no_samples_omits_gt(self):
        line = make_vcf_line(self._cluster(), "id1", ["s1", "s2"], no_samples=True)
        assert "GT" not in line

    def test_vcf_line_with_samples_has_format_and_gt(self):
        line = make_vcf_line(self._cluster(), "id1", ["s1", "s2"], no_samples=False)
        fields = line.split("\t")
        assert fields[8] == "GT"
        assert len(fields) == 11  # 8 fixed + FORMAT + 2 sample GT cols

    def test_header_no_samples_omits_format_line(self):
        hdr = db_header(args_no_samples)
        assert "##FORMAT" not in hdr

    def test_header_with_samples_includes_format_line(self):
        hdr = db_header(args)
        assert "##FORMAT" in hdr


class TestNoSamplesIntegration:
    """End-to-end: build a 2-sample DB, export with --samples off."""

    def test_chrom_header_has_no_sample_columns(self, tmp_path):
        prefix = tmp_path / "pop"
        subprocess.run(SVDB + ["--build", "--files", str(MANTA_CHR1), str(HG002_CHR1),
                               "--prefix", str(prefix)], check=True, capture_output=True)
        subprocess.run(SVDB + ["--export", "--db", str(tmp_path / "pop.db"),
                               "--prefix", str(tmp_path / "out"),
                               "--samples", "off"], check=True, capture_output=True)
        vcf = (tmp_path / "out.vcf").read_text()
        header = _chrom_header(vcf)
        fields = header.split("\t")
        assert fields == ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    def test_data_lines_have_eight_fields(self, tmp_path):
        prefix = tmp_path / "pop"
        subprocess.run(SVDB + ["--build", "--files", str(MANTA_CHR1), str(HG002_CHR1),
                               "--prefix", str(prefix)], check=True, capture_output=True)
        subprocess.run(SVDB + ["--export", "--db", str(tmp_path / "pop.db"),
                               "--prefix", str(tmp_path / "out"),
                               "--samples", "off"], check=True, capture_output=True)
        vcf = (tmp_path / "out.vcf").read_text()
        for line in _data_lines(vcf):
            assert len(line.split("\t")) == 8, f"Expected 8 fields, got: {line}"

    def test_occ_and_frq_still_in_info(self, tmp_path):
        prefix = tmp_path / "pop"
        subprocess.run(SVDB + ["--build", "--files", str(MANTA_CHR1), str(HG002_CHR1),
                               "--prefix", str(prefix)], check=True, capture_output=True)
        subprocess.run(SVDB + ["--export", "--db", str(tmp_path / "pop.db"),
                               "--prefix", str(tmp_path / "out"),
                               "--samples", "off"], check=True, capture_output=True)
        vcf = (tmp_path / "out.vcf").read_text()
        data = _data_lines(vcf)
        assert len(data) > 0
        assert all("OCC=" in ln and "FRQ=" in ln for ln in data)
