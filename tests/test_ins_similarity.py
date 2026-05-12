"""Unit tests for svdb.ins_similarity.

Covers extract_ins_sequence, levenshtein_similarity, and sequence_gate using
both simple synthetic sequences and sequences extracted directly from the
ins_similarity fixture VCFs.
"""
import pytest
from pathlib import Path

from svdb.ins_similarity import extract_ins_sequence, levenshtein_similarity, sequence_gate

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
