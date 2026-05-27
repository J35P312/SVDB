import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np

from svdb.export_module import (
    cluster_variants, cluster_variants_union_find,
    _expand_and_cluster_worker, _resolve_workers,
    db_header, make_representing_variant, build_genotype_columns,
    vcf_line as make_vcf_line,
    _pick_ins_seq,
)

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

    def test_vcf_line_no_samples_variants_has_no_sample_ids(self):
        line = make_vcf_line(self._cluster(), "id1", ["s1", "s2"], no_samples=True)
        info = line.split("\t")[7]
        variants_field = next(f for f in info.split(";") if f.startswith("VARIANTS="))
        assert "s1" not in variants_field
        assert "s2" not in variants_field

    def test_vcf_line_with_samples_variants_includes_sample_ids(self):
        line = make_vcf_line(self._cluster(), "id1", ["s1", "s2"], no_samples=False)
        info = line.split("\t")[7]
        variants_field = next(f for f in info.split(";") if f.startswith("VARIANTS="))
        assert "s1" in variants_field
        assert "s2" in variants_field

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

    def test_variants_field_has_no_sample_ids(self, tmp_path):
        prefix = tmp_path / "pop"
        subprocess.run(SVDB + ["--build", "--files", str(MANTA_CHR1), str(HG002_CHR1),
                               "--prefix", str(prefix)], check=True, capture_output=True)
        subprocess.run(SVDB + ["--export", "--db", str(tmp_path / "pop.db"),
                               "--prefix", str(tmp_path / "out"),
                               "--samples", "off"], check=True, capture_output=True)
        vcf = (tmp_path / "out.vcf").read_text()
        for line in _data_lines(vcf):
            info = line.split("\t")[7]
            variants = next((f for f in info.split(";") if f.startswith("VARIANTS=")), "")
            assert "manta_chr1" not in variants
            assert "hg002_chr1" not in variants

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


class TestClusterVariantsMultiMembership(unittest.TestCase):
    """Greedy star intentionally allows a variant to appear in multiple clusters.

    A variant that overlaps two different representatives contributes to both
    their clusters so that OCC/FRQ reflect every group it genuinely belongs to,
    rather than being arbitrarily assigned to one.

    Scenario
    --------
    Representatives 0 (3 neighbours) and 3 (2 neighbours) both have variant 1
    in their neighbour lists.  Representative 0 is processed first (highest
    degree) and claims variant 1 (marking it non-representative).  Representative
    3 is not in 0's list so it remains active; when processed it also adds
    variant 1 — intentional multi-cluster membership.
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

    def test_variant_can_appear_in_multiple_clusters(self):
        variant_dictionary, similarity_matrix = self._make_inputs()
        clusters = cluster_variants(variant_dictionary, similarity_matrix)

        membership: dict[int, list[int]] = {}
        for cluster_idx, cluster in enumerate(clusters):
            for var_idx in cluster[1]:
                membership.setdefault(var_idx, []).append(cluster_idx)

        multi = {k: v for k, v in membership.items() if len(v) > 1}
        self.assertTrue(
            multi,
            "expected variant 1 to appear in both clusters (multi-membership), but got exclusive assignment",
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

        Chain: ins0(1000) - ins1(1030) - ins2(1060) - ins3(1090), step=30 ≤ position_only ins_distance=50.
        Non-adjacent pairs (gap=60) exceed ins_distance=50. Degrees: ins0=2, ins1=3, ins2=3, ins3=2.
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
            f"1\t1030\tins1\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1030",
            f"1\t1060\tins2\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1060",
            f"1\t1090\tins3\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1090",
        ]
        vcf.write_text("\n".join(lines) + "\n")
        _subprocess.run(SVDB + ["--build", "--files", str(vcf),
                                "--prefix", str(tmp_path / "svdb")],
                        check=True, capture_output=True)

        # greedy star: ins1 has degree 3, claims ins0+ins1+ins2; ins3 is a singleton → 2 lines
        _subprocess.run(SVDB + ["--export", "--db", str(tmp_path / "svdb.db"),
                                "--prefix", str(tmp_path / "out_star"),
                                "--data_profile", "position_only"],
                        check=True, capture_output=True)
        star_lines = [ln for ln in (tmp_path / "out_star.vcf").read_text().splitlines()
                      if ln and not ln.startswith("#")]

        # union_find: transitive closure merges all four into one → 1 line
        _subprocess.run(SVDB + ["--export", "--db", str(tmp_path / "svdb.db"),
                                "--prefix", str(tmp_path / "out_uf"),
                                "--cluster_method", "union_find",
                                "--data_profile", "position_only"],
                        check=True, capture_output=True)
        uf_lines = [ln for ln in (tmp_path / "out_uf.vcf").read_text().splitlines()
                    if ln and not ln.startswith("#")]

        assert len(star_lines) == 2, f"expected 2 star clusters, got {len(star_lines)}"
        assert len(uf_lines) == 1, f"expected 1 union_find cluster, got {len(uf_lines)}"


class TestWorkers(unittest.TestCase):
    """Unit tests for worker-count resolution and the worker function itself."""

    def test_resolve_workers_auto_bounded(self):
        import os
        n = _resolve_workers(0)
        self.assertGreaterEqual(n, 1)
        self.assertLessEqual(n, max(1, os.cpu_count() or 1))

    def test_resolve_workers_explicit(self):
        self.assertEqual(_resolve_workers(2), 2)
        self.assertEqual(_resolve_workers(8), 8)

    def test_resolve_workers_serial(self):
        self.assertEqual(_resolve_workers(1), 1)

    def test_worker_fn_del(self):
        sub_dict = {
            0: {"posA": 1000, "posB": 2000, "sample_id": "s1", "ins_seq": None, "ins_len": None},
            1: {"posA": 1050, "posB": 2050, "sample_id": "s2", "ins_seq": None, "ins_len": None},
        }
        sub_coords = np.array([[0, 1000, 2000], [1, 1050, 2050]])
        task = (sub_dict, sub_coords, "DEL", "1", "1",
                25, 2500, 0.8, None, None, False, "star")
        result = _expand_and_cluster_worker(task)
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)
        rep, cdict = result[0]
        self.assertEqual(rep["type"], "DEL")
        self.assertEqual(rep["chrA"], "1")
        self.assertEqual(rep["chrB"], "1")

    def test_worker_fn_union_find(self):
        """Worker dispatches correctly to union_find when cluster_method='union_find'."""
        sub_dict = {
            0: {"posA": 100, "posB": 200, "sample_id": "s1", "ins_seq": None, "ins_len": None},
            1: {"posA": 150, "posB": 250, "sample_id": "s2", "ins_seq": None, "ins_len": None},
            2: {"posA": 200, "posB": 300, "sample_id": "s3", "ins_seq": None, "ins_len": None},
        }
        sub_coords = np.array([[0, 100, 200], [1, 150, 250], [2, 200, 300]])
        # chain: 0-1 and 1-2 overlap at bnd_distance=2500, 0-2 also do → one cluster either way
        task = (sub_dict, sub_coords, "DEL", "1", "1",
                25, 2500, 0.6, None, None, False, "union_find")
        result = _expand_and_cluster_worker(task)
        total_members = sum(len(cdict) for _, cdict in result)
        self.assertEqual(total_members, 3)


class TestWorkersIntegration:
    """--workers N produces the same variant count as --workers 1."""

    def test_workers_flag_same_line_count(self, tmp_path):
        """Parallel export (--workers 2) yields the same number of VCF lines as serial."""
        import subprocess as _subprocess
        vcf = tmp_path / "sample.vcf"
        seq = "ACGT" * 20
        lines = [
            "##fileformat=VCFv4.1",
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            f"1\t1000\ti0\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1000",
            f"1\t1015\ti1\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1015",
            f"1\t5000\ti2\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=5000",
            f"1\t5015\ti3\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=5015",
        ]
        vcf.write_text("\n".join(lines) + "\n")
        _subprocess.run(SVDB + ["--build", "--files", str(vcf),
                                "--prefix", str(tmp_path / "svdb")],
                        check=True, capture_output=True)

        for label, w in [("serial", "1"), ("parallel", "2")]:
            _subprocess.run(SVDB + ["--export", "--db", str(tmp_path / "svdb.db"),
                                    "--prefix", str(tmp_path / f"out_{label}"),
                                    "--workers", w, "--data_profile", "position_only"],
                            check=True, capture_output=True)

        def data_lines(p):
            return [ln for ln in p.read_text().splitlines() if ln and not ln.startswith("#")]

        serial = data_lines(tmp_path / "out_serial.vcf")
        parallel = data_lines(tmp_path / "out_parallel.vcf")
        assert len(serial) == len(parallel), (
            f"serial={len(serial)} lines, parallel={len(parallel)} lines"
        )


class TestPickInsSeq:
    """_pick_ins_seq: most-common sequence or None for missing/mixed clusters."""

    def _cluster(self, entries):
        """Build a variant_dict {i: {"ins_seq": ...}} from a list of seq-or-None values."""
        return {i: {"ins_seq": s} for i, s in enumerate(entries)}

    def test_all_same_sequence_returns_it(self):
        d = self._cluster(["ACGT", "ACGT", "ACGT"])
        assert _pick_ins_seq(d) == "ACGT"

    def test_most_common_sequence_wins(self):
        d = self._cluster(["ACGT", "TTTT", "ACGT"])
        assert _pick_ins_seq(d) == "ACGT"

    def test_all_none_returns_none(self):
        d = self._cluster([None, None])
        assert _pick_ins_seq(d) is None

    def test_mixed_some_none_returns_none(self):
        d = self._cluster(["ACGT", None, "ACGT"])
        assert _pick_ins_seq(d) is None

    def test_single_member_with_sequence(self):
        d = self._cluster(["TTTT"])
        assert _pick_ins_seq(d) == "TTTT"

    def test_single_member_none_returns_none(self):
        d = self._cluster([None])
        assert _pick_ins_seq(d) is None

    def test_mixed_cluster_vcf_line_outputs_symbolic_ins(self, tmp_path):
        """vcf_line emits <INS> ALT when cluster has mixed ins_seq membership."""
        import subprocess as _subprocess
        seq = "ACGTACGT" * 10  # 80 bp, under default 1000 cap
        vcf_with = tmp_path / "with_seq.vcf"
        vcf_sym = tmp_path / "symbolic.vcf"
        vcf_with.write_text("\n".join([
            "##fileformat=VCFv4.1",
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            f"1\t1000\ti0\tN\tN{seq}\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1000",
        ]) + "\n")
        vcf_sym.write_text("\n".join([
            "##fileformat=VCFv4.1",
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            f"1\t1005\ti1\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN={len(seq)};END=1005",
        ]) + "\n")
        _subprocess.run(
            SVDB + ["--build", "--files", str(vcf_with), str(vcf_sym),
                    "--prefix", str(tmp_path / "svdb")],
            check=True, capture_output=True,
        )
        _subprocess.run(
            SVDB + ["--export", "--db", str(tmp_path / "svdb.db"),
                    "--prefix", str(tmp_path / "out")],
            check=True, capture_output=True,
        )
        lines = _data_lines((tmp_path / "out.vcf").read_text())
        assert lines, "exported VCF has no data lines"
        alt_col = lines[0].split("\t")[4]
        assert alt_col == "<INS>", f"expected symbolic <INS>, got {alt_col!r}"
