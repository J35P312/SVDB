"""Unit tests for svdb.ins_similarity.

Covers extract_ins_sequence, levenshtein_similarity, and sequence_gate using
both simple synthetic sequences and sequences extracted directly from the
ins_similarity fixture VCFs.
"""
import argparse
import unittest
import unittest.mock
import pytest
from pathlib import Path

from svdb.__main__ import _fraction, _positive_int
from svdb.ins_similarity import (
    apply_ins_profile,
    cap_seq,
    decompress_ins_seq,
    extract_ins_sequence,
    has_excess_n,
    levenshtein_similarity,
    parse_svlen,
    sequence_gate,
)

FIXTURES = Path(__file__).parent / "fixtures" / "ins_similarity"


def _read_vcf_seq(vcf_path: Path) -> str:
    """Extract the insertion sequence from the single data line of a fixture VCF."""
    for line in vcf_path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        return extract_ins_sequence(fields[3], fields[4])
    return ""


class TestExtractInsSequence:

    def test_symbolic_ins_returns_empty(self):
        assert extract_ins_sequence("G", "<INS>") == ""

    def test_symbolic_ins_me_returns_empty(self):
        assert extract_ins_sequence("N", "<INS:ME>") == ""

    def test_any_bracket_alt_returns_empty(self):
        assert extract_ins_sequence("N", "<DUP>") == ""

    def test_anchored_n_ref_strips_first_base(self):
        assert extract_ins_sequence("N", "NATCGATCG") == "ATCGATCG"

    def test_anchored_g_ref_strips_first_base(self):
        assert extract_ins_sequence("G", "GTCGATCGAT") == "TCGATCGAT"

    def test_single_base_alt_same_as_ref_returns_empty(self):
        assert extract_ins_sequence("N", "N") == ""

    def test_unanchored_alt_returned_as_is(self):
        # First base of ALT does not match REF — treat whole ALT as sequence
        assert extract_ins_sequence("N", "ATCGATCG") == "ATCGATCG"

    def test_fixture_grch37_pos0_not_empty(self):
        seq = _read_vcf_seq(FIXTURES / "grch37_pos0_sim1.00" / "caller_a.vcf")
        assert len(seq) == 113

    def test_fixture_symbolic_alt_returns_empty(self):
        seq = _read_vcf_seq(FIXTURES / "symbolic_alt" / "caller_a.vcf")
        assert seq == ""

    def test_fixture_low_svlen_sequences_non_empty(self):
        # low_svlen_ratio: both have sequences (SVLEN=100/220, rejected by SVLEN gate, not seq gate)
        seq_a = _read_vcf_seq(FIXTURES / "low_svlen_ratio" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "low_svlen_ratio" / "caller_b.vcf")
        assert seq_a != ""
        assert seq_b != ""
        assert len(seq_b) > len(seq_a)


class TestLevenshteinSimilarity:

    def test_identical_returns_one(self):
        assert levenshtein_similarity("ATCGATCG", "ATCGATCG") == pytest.approx(1.0)

    def test_one_substitution_in_eight(self):
        sim = levenshtein_similarity("ATCGATCG", "ATCGATCC")
        assert sim == pytest.approx(7 / 8)

    def test_completely_different_is_low(self):
        assert levenshtein_similarity("AAAAAAAA", "TTTTTTTT") < 0.2

    def test_empty_both_returns_one(self):
        assert levenshtein_similarity("", "") == pytest.approx(1.0)

    def test_fixture_pos0_sim_is_one(self):
        seq_a = _read_vcf_seq(FIXTURES / "grch37_pos0_sim1.00" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch37_pos0_sim1.00" / "caller_b.vcf")
        assert levenshtein_similarity(seq_a, seq_b) == pytest.approx(1.0)

    def test_fixture_pos25_high_sim(self):
        # PAIR07327 CuteSV×DYSGU: dist=1, expected ~0.977
        seq_a = _read_vcf_seq(FIXTURES / "grch37_pos25_sim0.985" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch37_pos25_sim0.985" / "caller_b.vcf")
        assert levenshtein_similarity(seq_a, seq_b) > 0.96

    def test_fixture_pos25_medium_sim(self):
        # PAIR08263 CuteSV×GIAB: dist=6, expected ~0.789
        seq_a = _read_vcf_seq(FIXTURES / "grch37_pos25_sim0.789" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch37_pos25_sim0.789" / "caller_b.vcf")
        sim = levenshtein_similarity(seq_a, seq_b)
        assert 0.76 < sim < 0.82

    def test_fixture_neg_c_medium_high_sim(self):
        # PAIR30601 adversarial complexity mismatch: expected ~0.826
        seq_a = _read_vcf_seq(FIXTURES / "grch37_neg_c_sim0.837" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch37_neg_c_sim0.837" / "caller_b.vcf")
        sim = levenshtein_similarity(seq_a, seq_b)
        assert 0.80 < sim < 0.86

    def test_fixture_neg_c_low_sim(self):
        # PAIR30604 adversarial: expected ~0.700
        seq_a = _read_vcf_seq(FIXTURES / "grch37_neg_c_sim0.672" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch37_neg_c_sim0.672" / "caller_b.vcf")
        sim = levenshtein_similarity(seq_a, seq_b)
        assert 0.67 < sim < 0.74

    def test_fixture_grch38_neg_c_very_low_sim(self):
        # PAIR73997 adversarial ONT: expected ~0.573
        seq_a = _read_vcf_seq(FIXTURES / "grch38_neg_c_sim0.573" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch38_neg_c_sim0.573" / "caller_b.vcf")
        sim = levenshtein_similarity(seq_a, seq_b)
        assert sim < 0.65


class TestSequenceGate:

    def test_empty_seq_a_skips_gate(self):
        """Symbolic ALT on one side → no sequence to compare, gate passes (True = allow merge)."""
        assert sequence_gate("", "ATCGATCGATCGATCG", threshold=0.75) is True

    def test_empty_seq_b_skips_gate(self):
        assert sequence_gate("ATCGATCGATCGATCG", "", threshold=0.75) is True

    def test_both_empty_skips_gate(self):
        assert sequence_gate("", "", threshold=0.75) is True

    def test_identical_above_any_threshold(self):
        seq = "ATCGATCGATCGATCGATCGATCGATCG"
        assert sequence_gate(seq, seq, threshold=0.95) is True

    def test_below_threshold_rejects(self):
        # sim("AAAAAAAA", "TTTTTTTT") ≈ 0; threshold=0.75 → reject (False)
        assert sequence_gate("AAAAAAAA", "TTTTTTTT", threshold=0.75) is False

    def test_at_threshold_passes(self):
        # Construct sequences with known similarity just at the threshold.
        # "AAAAAATT" vs "AAAAAAAA": 2 subs in 8 → sim=6/8=0.75
        assert sequence_gate("AAAAAATT", "AAAAAAAA", threshold=0.75) is True

    def test_cohort_threshold_passes_medium_sim(self):
        # grch37_pos25_sim0.789 sequences: sim≈0.789 > 0.75 → passes cohort threshold
        seq_a = _read_vcf_seq(FIXTURES / "grch37_pos25_sim0.789" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch37_pos25_sim0.789" / "caller_b.vcf")
        assert sequence_gate(seq_a, seq_b, threshold=0.75) is True

    def test_sample_threshold_rejects_medium_sim(self):
        # same sequences: sim≈0.789 < 0.85 → rejected at sample threshold
        seq_a = _read_vcf_seq(FIXTURES / "grch37_pos25_sim0.789" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch37_pos25_sim0.789" / "caller_b.vcf")
        assert sequence_gate(seq_a, seq_b, threshold=0.85) is False

    def test_cohort_threshold_rejects_low_sim(self):
        # grch37_neg_c_sim0.672: sim≈0.700 < 0.75 → rejected by cohort
        seq_a = _read_vcf_seq(FIXTURES / "grch37_neg_c_sim0.672" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch37_neg_c_sim0.672" / "caller_b.vcf")
        assert sequence_gate(seq_a, seq_b, threshold=0.75) is False

    def test_sample_threshold_rejects_adversarial_neg_c(self):
        # grch37_neg_c_sim0.837 adversarial pair: sim≈0.826 < 0.85 → rejected by sample
        seq_a = _read_vcf_seq(FIXTURES / "grch37_neg_c_sim0.837" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch37_neg_c_sim0.837" / "caller_b.vcf")
        assert sequence_gate(seq_a, seq_b, threshold=0.85) is False

    def test_cohort_threshold_passes_adversarial_neg_c(self):
        # same pair: sim≈0.826 > 0.75 → passes cohort threshold (known limitation)
        seq_a = _read_vcf_seq(FIXTURES / "grch37_neg_c_sim0.837" / "caller_a.vcf")
        seq_b = _read_vcf_seq(FIXTURES / "grch37_neg_c_sim0.837" / "caller_b.vcf")
        assert sequence_gate(seq_a, seq_b, threshold=0.75) is True


class TestHasExcessN:

    def test_no_n_returns_false(self):
        assert has_excess_n("ATCGATCGATCG") is False

    def test_none_input_returns_false(self):
        assert has_excess_n(None) is False

    def test_empty_string_returns_false(self):
        assert has_excess_n("") is False

    def test_above_default_threshold_returns_true(self):
        # 2/10 N = 0.2 > 0.1 default
        assert has_excess_n("NNAAAAAAAA") is True

    def test_exactly_at_default_threshold_returns_false(self):
        # 1/10 N = 0.1, not > 0.1 -- boundary is not flagged
        assert has_excess_n("NAAAAAAAAA") is False

    def test_all_n_returns_true(self):
        assert has_excess_n("NNNNNNNN") is True

    def test_case_insensitive(self):
        assert has_excess_n("nnnnaaaa") is True
        assert has_excess_n("NnNnAAAA") is True

    def test_custom_threshold(self):
        seq = "NAAAAAAAAA"  # 1/10 = 0.1
        assert has_excess_n(seq, max_n_fraction=0.2) is False
        assert has_excess_n(seq, max_n_fraction=0.05) is True


class TestCapSeq:

    def test_none_input_returns_empty(self):
        assert cap_seq(None, 500) == ""

    def test_empty_string_returns_empty(self):
        assert cap_seq("", 500) == ""

    def test_no_cap_returns_seq_unchanged(self):
        assert cap_seq("ATCG", None) == "ATCG"

    def test_under_cap_returns_seq_unchanged(self):
        assert cap_seq("ATCG", 10) == "ATCG"

    def test_over_cap_returns_empty(self):
        assert cap_seq("A" * 600, 500) == ""

    def test_exactly_at_cap_is_not_capped(self):
        seq = "A" * 500
        assert cap_seq(seq, 500) == seq

    def test_excess_n_returns_empty(self):
        """N-heavy sequence is treated the same as over-cap: deferred to
        position+SVLEN matching, same as a symbolic ALT -- issue #95."""
        assert cap_seq("NNNNNNNNAA", None) == ""

    def test_excess_n_then_sequence_gate_defers(self):
        """End-to-end: an N-heavy sequence run through cap_seq (as merge and
        query do) makes sequence_gate defer (return True) regardless of how
        dissimilar the other side is."""
        n_heavy = cap_seq("NNNNNNNNNN", None)
        other = cap_seq("ATCGATCGAT", None)
        assert sequence_gate(n_heavy, other, threshold=0.99) is True


class TestDecompressInsSeq:

    def test_none_returns_none(self):
        assert decompress_ins_seq(None) is None

    def test_str_passthrough(self):
        assert decompress_ins_seq("ATCG") == "ATCG"

    def test_bytes_decompresses(self):
        import zlib
        compressed = zlib.compress(b"ATCGATCG")
        assert decompress_ins_seq(compressed) == "ATCGATCG"


class TestParseSvlen:

    def test_positive_svlen(self):
        assert parse_svlen("END=1000;SVLEN=150;SVTYPE=INS") == 150

    def test_negative_svlen_returns_abs(self):
        assert parse_svlen("SVLEN=-200;SVTYPE=DEL") == 200

    def test_absent_returns_none(self):
        assert parse_svlen("END=1000;SVTYPE=DEL") is None

    def test_at_start_of_info(self):
        assert parse_svlen("SVLEN=42") == 42


class TestApplyInsProfile:

    def _args(self, **kwargs):
        ns = argparse.Namespace(
            data_profile=None,
            ins_distance=None,
            ins_svlen_ratio=None,
            ins_seq_similarity=None,
        )
        for k, v in kwargs.items():
            setattr(ns, k, v)
        return ns

    def test_no_profile_sets_defaults(self):
        args = self._args()
        apply_ins_profile(args)
        assert args.ins_distance == 25
        assert args.ins_svlen_ratio == 0.90
        assert args.ins_seq_similarity == 0.75
        assert args.no_ins_seq is False

    def test_cohort_profile(self):
        args = self._args(data_profile="cohort")
        apply_ins_profile(args)
        assert args.ins_distance == 50
        assert args.ins_svlen_ratio == 0.80
        assert args.ins_seq_similarity == 0.75
        assert args.no_ins_seq is False

    def test_sample_profile(self):
        args = self._args(data_profile="sample")
        apply_ins_profile(args)
        assert args.ins_distance == 25
        assert args.ins_svlen_ratio == 0.90
        assert args.ins_seq_similarity == 0.85
        assert args.no_ins_seq is False

    def test_position_only_profile(self):
        args = self._args(data_profile="position_only")
        apply_ins_profile(args)
        assert args.ins_distance == 50
        assert args.ins_svlen_ratio == 0.90
        assert args.no_ins_seq is True

    def test_explicit_flag_overrides_profile(self):
        args = self._args(data_profile="cohort", ins_distance=10)
        apply_ins_profile(args)
        assert args.ins_distance == 10
        assert args.ins_svlen_ratio == 0.80  # from cohort

    def test_explicit_sim_overrides_position_only(self):
        args = self._args(data_profile="position_only", ins_seq_similarity=0.90)
        apply_ins_profile(args)
        assert args.ins_seq_similarity == 0.90
        assert args.no_ins_seq is True  # profile still sets this


class TestArgValidators(unittest.TestCase):

    def test_fraction_valid_boundary_values(self):
        f = _fraction("--overlap")
        assert f("0.0") == 0.0
        assert f("1.0") == 1.0
        assert f("0.75") == 0.75

    def test_fraction_above_one_raises(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            _fraction("--overlap")("1.2")

    def test_fraction_negative_raises(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            _fraction("--ins_seq_similarity")("-0.1")

    def test_positive_int_valid(self):
        f = _positive_int("--min_pts")
        assert f("1") == 1
        assert f("5") == 5

    def test_positive_int_zero_raises(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            _positive_int("--min_pts")("0")

    def test_positive_int_negative_raises(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            _positive_int("--min_pts")("-1")

    def test_positive_int_float_string_raises(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            _positive_int("--min_pts")("1.2")


class TestCLIArgValidation(unittest.TestCase):
    """Verify that out-of-range values are rejected by the real argparse wiring,
    not just by the validator functions in isolation."""

    def _assert_cli_error(self, argv, expected_fragment):
        """Run svdb with argv; assert it exits non-zero and stderr contains expected_fragment."""
        from svdb.__main__ import main
        import io
        with unittest.mock.patch("sys.argv", argv), \
             unittest.mock.patch("sys.stderr", new_callable=io.StringIO) as mock_err, \
             self.assertRaises(SystemExit) as cm:
            main()
        self.assertNotEqual(cm.exception.code, 0)
        self.assertIn(expected_fragment, mock_err.getvalue())

    def test_merge_overlap_above_one(self):
        self._assert_cli_error(
            ["svdb", "--merge", "--overlap", "1.5", "--vcf", "/dev/null"],
            "--overlap must be in [0.0, 1.0]",
        )

    def test_merge_ins_svlen_ratio_above_one(self):
        self._assert_cli_error(
            ["svdb", "--merge", "--ins_svlen_ratio", "1.1", "--vcf", "/dev/null"],
            "--ins_svlen_ratio must be in [0.0, 1.0]",
        )

    def test_merge_ins_seq_similarity_above_one(self):
        self._assert_cli_error(
            ["svdb", "--merge", "--ins_seq_similarity", "1.2", "--vcf", "/dev/null"],
            "--ins_seq_similarity must be in [0.0, 1.0]",
        )

    def test_export_overlap_negative(self):
        self._assert_cli_error(
            ["svdb", "--export", "--db", "x.db", "--overlap", "-0.1"],
            "--overlap must be in [0.0, 1.0]",
        )

    def test_export_min_pts_zero(self):
        self._assert_cli_error(
            ["svdb", "--export", "--db", "x.db", "--min_pts", "0"],
            "--min_pts must be ≥ 1",
        )

    def test_query_max_frq_above_one(self):
        self._assert_cli_error(
            ["svdb", "--query", "--query_vcf", "x.vcf", "--db", "x.vcf", "--max_frq", "1.5"],
            "--max_frq must be in [0.0, 1.0]",
        )
