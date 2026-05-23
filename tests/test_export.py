import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np

from svdb.export_module import cluster_variants, cluster_variants_union_find, db_header, make_representing_variant, build_genotype_columns, vcf_line as make_vcf_line

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


class TestClusterVariantsNoDuplicateMembership(unittest.TestCase):
    """Each variant must appear in exactly one output cluster.

    Regression test for the greedy-star duplicate-membership bug: a variant
    that is a neighbour of two different representatives gets added to both
    clusters because the inner loop does not check whether the variant has
    already been claimed.

    Scenario
    --------
    Representatives 0 (3 neighbours) and 3 (2 neighbours) both have variant 1
    in their neighbour lists.  Representative 0 is processed first (highest
    degree) and claims variant 1.  Representative 3 is not in 0's list so it
    stays unclaimed; when it is processed it adds variant 1 a second time.

    This test FAILS with the current (buggy) cluster_variants and will PASS
    once the implementation skips already-claimed neighbours.
    """

    def _make_inputs(self):
        variant_dictionary = {
            0: {"posA": 100, "posB": 200, "sample_id": "s"},
            1: {"posA": 200, "posB": 300, "sample_id": "s"},
            2: {"posA": 150, "posB": 250, "sample_id": "s"},
            3: {"posA": 500, "posB": 600, "sample_id": "s"},
        }
        similarity_matrix = {
            0: np.array([0, 1, 2]),  # highest degree → processed first
            1: np.array([0, 1, 2]),
            2: np.array([0, 1, 2]),
            3: np.array([1, 3]),     # shares neighbour 1 with rep 0; not in 0's list
        }
        return variant_dictionary, similarity_matrix

    def test_no_variant_in_two_clusters(self):
        variant_dictionary, similarity_matrix = self._make_inputs()
        clusters = cluster_variants(variant_dictionary, similarity_matrix)

        membership: dict[int, list[int]] = {}
        for cluster_idx, cluster in enumerate(clusters):
            for var_idx in cluster[1]:
                membership.setdefault(var_idx, []).append(cluster_idx)

        duplicates = {k: v for k, v in membership.items() if len(v) > 1}
        self.assertFalse(
            duplicates,
            f"variants appear in multiple clusters (bug): {duplicates}",
        )


class TestClusterVariantsUnionFind(unittest.TestCase):
    """Union-Find clustering: transitivity and no duplicate membership."""

    def _make_chain_inputs(self):
        """A-B and B-C overlap; A-C do not.

        greedy_star produces two clusters ({A,B} and {C} or {A} and {B,C}
        depending on degree); union_find must produce one: {A, B, C}.
        """
        variant_dictionary = {
            0: {"posA": 100, "posB": 200, "sample_id": "s"},  # A
            1: {"posA": 150, "posB": 250, "sample_id": "s"},  # B — overlaps A and C
            2: {"posA": 200, "posB": 300, "sample_id": "s"},  # C
        }
        # A-B edge, B-C edge, no A-C edge
        similarity_matrix = {
            0: np.array([0, 1]),
            1: np.array([0, 1, 2]),
            2: np.array([1, 2]),
        }
        return variant_dictionary, similarity_matrix

    def test_transitivity_merges_chain(self):
        """A-B and B-C must all end up in one cluster."""
        variant_dictionary, similarity_matrix = self._make_chain_inputs()
        clusters = cluster_variants_union_find(variant_dictionary, similarity_matrix)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(set(clusters[0][1].keys()), {0, 1, 2})

    def test_no_duplicate_membership(self):
        """Each variant appears in exactly one cluster."""
        variant_dictionary, similarity_matrix = self._make_chain_inputs()
        clusters = cluster_variants_union_find(variant_dictionary, similarity_matrix)
        membership: dict[int, list[int]] = {}
        for idx, cluster in enumerate(clusters):
            for var in cluster[1]:
                membership.setdefault(var, []).append(idx)
        duplicates = {k: v for k, v in membership.items() if len(v) > 1}
        self.assertFalse(duplicates)

    def test_disjoint_groups_stay_separate(self):
        """Two groups with no edges between them produce two clusters."""
        variant_dictionary = {
            0: {"posA": 100, "posB": 200, "sample_id": "s"},
            1: {"posA": 110, "posB": 210, "sample_id": "s"},
            2: {"posA": 900, "posB": 1000, "sample_id": "s"},
            3: {"posA": 910, "posB": 1010, "sample_id": "s"},
        }
        similarity_matrix = {
            0: np.array([0, 1]),
            1: np.array([0, 1]),
            2: np.array([2, 3]),
            3: np.array([2, 3]),
        }
        clusters = cluster_variants_union_find(variant_dictionary, similarity_matrix)
        self.assertEqual(len(clusters), 2)
        members = [set(c[1].keys()) for c in clusters]
        self.assertIn({0, 1}, members)
        self.assertIn({2, 3}, members)

    def test_representative_is_highest_degree(self):
        """The representative of a component is the variant with most neighbours."""
        variant_dictionary = {
            0: {"posA": 100, "posB": 200, "sample_id": "s"},
            1: {"posA": 150, "posB": 250, "sample_id": "s"},
            2: {"posA": 200, "posB": 300, "sample_id": "s"},
        }
        similarity_matrix = {
            0: np.array([0, 1]),
            1: np.array([0, 1, 2]),  # highest degree
            2: np.array([1, 2]),
        }
        clusters = cluster_variants_union_find(variant_dictionary, similarity_matrix)
        self.assertEqual(len(clusters), 1)
        # representative posA should be variant 1's posA=150
        self.assertEqual(clusters[0][0]["posA"], 150)


class TestClusterMethodIntegration:
    """--cluster_method union_find produces fewer lines and no duplicate membership."""

    def test_union_find_export_fewer_lines_than_star(self, tmp_path):
        """Union-Find merges transitively; star splits a 4-node chain into two clusters.

        Chain: ins0(1000) - ins1(1015) - ins2(1030) - ins3(1045), step=15 ≤ ins_distance=25.
        Non-adjacent pairs (gap=30) do NOT overlap. Degrees: ins0=2, ins1=3, ins2=3, ins3=2.
        Greedy star: ins1 (first highest-degree) claims ins0,ins1,ins2; ins3 is left alone → 2 clusters.
        Union-Find: transitive closure connects all four → 1 cluster.
        """
        import subprocess as _subprocess
        vcf = tmp_path / "sample.vcf"
        seq = "ACGT" * 20   # 80 bp, identical
        lines = [
            "##fileformat=VCFv4.1",
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            f"1\t1000\tins0\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1000",
            f"1\t1015\tins1\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1015",
            f"1\t1030\tins2\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1030",
            f"1\t1045\tins3\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1045",
        ]
        vcf.write_text("\n".join(lines) + "\n")
        _subprocess.run(SVDB + ["--build", "--files", str(vcf),
                                "--prefix", str(tmp_path / "svdb")],
                        check=True, capture_output=True)

        # greedy star: ins1 has degree 3, claims ins0+ins1+ins2; ins3 is a singleton → 2 lines
        _subprocess.run(SVDB + ["--export", "--db", str(tmp_path / "svdb.db"),
                                "--prefix", str(tmp_path / "out_star"),
                                "--no_ins_seq"],
                        check=True, capture_output=True)
        star_lines = [ln for ln in (tmp_path / "out_star.vcf").read_text().splitlines()
                      if ln and not ln.startswith("#")]

        # union_find: transitive closure merges all four into one → 1 line
        _subprocess.run(SVDB + ["--export", "--db", str(tmp_path / "svdb.db"),
                                "--prefix", str(tmp_path / "out_uf"),
                                "--cluster_method", "union_find",
                                "--no_ins_seq"],
                        check=True, capture_output=True)
        uf_lines = [ln for ln in (tmp_path / "out_uf.vcf").read_text().splitlines()
                    if ln and not ln.startswith("#")]

        assert len(star_lines) == 2, f"expected 2 star clusters, got {len(star_lines)}"
        assert len(uf_lines) == 1, f"expected 1 union_find cluster, got {len(uf_lines)}"
