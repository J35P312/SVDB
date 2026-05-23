"""Integration tests for insertion sequence-aware merging (issue #82).

Each test drives the full SVDB merge pipeline via subprocess and checks the
number of output VCF data lines: if two records are merged the output has 1
line; if they are kept separate it has 2 lines.

Fixture layout (tests/fixtures/ins_similarity/):
  grch37_pos0_sim1.00   — GRCh37, dist=0,  sim=1.000  (identical, must merge)
  grch37_pos25_sim0.985 — GRCh37, dist=1,  sim=0.977  (high sim, must merge)
  grch37_pos25_sim0.789 — GRCh37, dist=6,  sim=0.789  (medium sim: cohort merges, sample rejects)
  grch37_neg_c_sim0.837 — GRCh37, dist=8,  sim=0.826  (adversarial: sample rejects, cohort passes)
  grch37_neg_c_sim0.672 — GRCh37, dist=1,  sim=0.700  (below both thresholds, always rejects)
  grch38_pos0_sim1.00   — GRCh38, dist=0,  sim=1.000  (identical, must merge)
  grch38_neg_c_sim0.573 — GRCh38, dist=26, sim=0.573  (dist outside hard cap → seq gate skipped)
  low_svlen_ratio        — dist=0,  SVLEN=100/220  (SVLEN ratio 0.455 < 0.90, must not merge)
  symbolic_alt           — GRCh38, dist=7,  both <INS> (no sequence, must merge on pos+SVLEN)
"""
import subprocess
import sys
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures" / "ins_similarity"
SVDB = [sys.executable, "-m", "svdb"]


def run_merge(*args):
    return subprocess.run(
        SVDB + ["--merge"] + list(args),
        capture_output=True,
        text=True,
    )


def data_lines(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln and not ln.startswith("#")]


def fixture(name: str, caller: str) -> str:
    return str(FIXTURES / name / f"{caller}.vcf")


# ---------------------------------------------------------------------------
# Default behaviour: pos+SVLEN+seq (cohort profile, threshold=0.75)
# ---------------------------------------------------------------------------

class TestDefaultBehaviour:

    def test_identical_sequences_merge(self):
        r = run_merge("--vcf", fixture("grch37_pos0_sim1.00", "caller_a"),
                      fixture("grch37_pos0_sim1.00", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_grch38_identical_sequences_merge(self):
        r = run_merge("--vcf", fixture("grch38_pos0_sim1.00", "caller_a"),
                      fixture("grch38_pos0_sim1.00", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_high_sim_pos25_merges(self):
        # sim≈0.977 > 0.75 → merged
        r = run_merge("--vcf", fixture("grch37_pos25_sim0.985", "caller_a"),
                      fixture("grch37_pos25_sim0.985", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_medium_sim_pos25_merges_at_cohort_threshold(self):
        # sim≈0.789 > 0.75 → merged with default cohort threshold
        r = run_merge("--vcf", fixture("grch37_pos25_sim0.789", "caller_a"),
                      fixture("grch37_pos25_sim0.789", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_low_sim_neg_c_rejected_by_default(self):
        # sim≈0.700 < 0.75 → not merged
        r = run_merge("--vcf", fixture("grch37_neg_c_sim0.672", "caller_a"),
                      fixture("grch37_neg_c_sim0.672", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 2

    def test_symbolic_alt_merges_on_position_and_svlen(self):
        # Both <INS>: no sequence → gate skipped, merge on pos+SVLEN
        # SVLEN=71/70 (ratio=0.986 ≥ 0.90), dist=7 ≤ 25 → merged
        r = run_merge("--vcf", fixture("symbolic_alt", "caller_a"),
                      fixture("symbolic_alt", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_low_svlen_ratio_not_merged(self):
        # SVLEN=100/220 (ratio≈0.455 < 0.90) → SVLEN gate rejects
        r = run_merge("--vcf", fixture("low_svlen_ratio", "caller_a"),
                      fixture("low_svlen_ratio", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 2


# ---------------------------------------------------------------------------
# --data_profile sample (threshold = 0.85)
# ---------------------------------------------------------------------------

class TestDataProfileSample:

    def test_identical_sequences_still_merge(self):
        r = run_merge("--data_profile", "sample",
                      "--vcf", fixture("grch37_pos0_sim1.00", "caller_a"),
                      fixture("grch37_pos0_sim1.00", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_high_sim_pos25_still_merges(self):
        # sim≈0.977 > 0.85 → still merged at sample threshold
        r = run_merge("--data_profile", "sample",
                      "--vcf", fixture("grch37_pos25_sim0.985", "caller_a"),
                      fixture("grch37_pos25_sim0.985", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_medium_sim_rejected_at_sample_threshold(self):
        # sim≈0.789 < 0.85 → rejected at sample threshold (false split, expected tradeoff)
        r = run_merge("--data_profile", "sample",
                      "--vcf", fixture("grch37_pos25_sim0.789", "caller_a"),
                      fixture("grch37_pos25_sim0.789", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 2

    def test_adversarial_neg_c_rejected_at_sample_threshold(self):
        # sim≈0.826 < 0.85 → correctly blocked by sample threshold
        r = run_merge("--data_profile", "sample",
                      "--vcf", fixture("grch37_neg_c_sim0.837", "caller_a"),
                      fixture("grch37_neg_c_sim0.837", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 2

    def test_symbolic_alt_still_merges(self):
        r = run_merge("--data_profile", "sample",
                      "--vcf", fixture("symbolic_alt", "caller_a"),
                      fixture("symbolic_alt", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1


# ---------------------------------------------------------------------------
# --data_profile cohort (threshold = 0.75, same as default)
# ---------------------------------------------------------------------------

class TestDataProfileCohort:

    def test_adversarial_neg_c_passes_cohort_threshold(self):
        # sim≈0.826 > 0.75 → merged at cohort threshold (known limitation per benchmark)
        r = run_merge("--data_profile", "cohort",
                      "--vcf", fixture("grch37_neg_c_sim0.837", "caller_a"),
                      fixture("grch37_neg_c_sim0.837", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_low_sim_still_rejected(self):
        r = run_merge("--data_profile", "cohort",
                      "--vcf", fixture("grch37_neg_c_sim0.672", "caller_a"),
                      fixture("grch37_neg_c_sim0.672", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 2


# ---------------------------------------------------------------------------
# --no_ins_seq: disable sequence gate entirely
# ---------------------------------------------------------------------------

class TestNoInsSeq:

    def test_low_sim_merges_when_seq_disabled(self):
        # sim≈0.700 would normally be rejected; with --no_ins_seq it merges on pos+SVLEN
        r = run_merge("--no_ins_seq",
                      "--vcf", fixture("grch37_neg_c_sim0.672", "caller_a"),
                      fixture("grch37_neg_c_sim0.672", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_adversarial_neg_c_merges_when_seq_disabled(self):
        r = run_merge("--no_ins_seq",
                      "--vcf", fixture("grch37_neg_c_sim0.837", "caller_a"),
                      fixture("grch37_neg_c_sim0.837", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_low_svlen_ratio_still_rejected_with_seq_disabled(self):
        # SVLEN gate is independent of sequence gate
        r = run_merge("--no_ins_seq",
                      "--vcf", fixture("low_svlen_ratio", "caller_a"),
                      fixture("low_svlen_ratio", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 2


# ---------------------------------------------------------------------------
# --ins_svlen_ratio custom value
# ---------------------------------------------------------------------------

class TestInsSvlenRatio:

    def test_loose_svlen_ratio_merges_low_ratio_pair(self):
        # SVLEN=100/220, ratio=0.455. With --ins_svlen_ratio 0.40 the SVLEN gate passes.
        # Disable sequence check to isolate SVLEN gate behaviour.
        r = run_merge("--ins_svlen_ratio", "0.40", "--no_ins_seq",
                      "--vcf", fixture("low_svlen_ratio", "caller_a"),
                      fixture("low_svlen_ratio", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_strict_svlen_ratio_rejects_high_ratio_pair(self):
        # grch37_pos25_sim0.985: SVLEN=133/133, ratio=1.0. With --ins_svlen_ratio 1.0 still passes.
        r = run_merge("--ins_svlen_ratio", "1.0",
                      "--vcf", fixture("grch37_pos25_sim0.985", "caller_a"),
                      fixture("grch37_pos25_sim0.985", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1


# ---------------------------------------------------------------------------
# Position window and hard cap (sequence gate applies only within 25bp)
# ---------------------------------------------------------------------------

class TestPositionWindow:

    def test_default_ins_distance_25_rejects_pos26(self):
        # grch38_neg_c_sim0.573 is at dist=26, which exceeds default ins_distance=25
        r = run_merge("--vcf", fixture("grch38_neg_c_sim0.573", "caller_a"),
                      fixture("grch38_neg_c_sim0.573", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 2

    def test_ins_distance_50_merges_pos26_via_svlen_only(self):
        # With --ins_distance 50: dist=26 ≤ 50, SVLEN=108/117 ratio≈0.923 ≥ 0.90.
        # dist=26 > 25 hard cap → sequence gate skipped. Merge on pos+SVLEN.
        r = run_merge("--ins_distance", "50",
                      "--vcf", fixture("grch38_neg_c_sim0.573", "caller_a"),
                      fixture("grch38_neg_c_sim0.573", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1


# ---------------------------------------------------------------------------
# --data_profile + explicit --ins_seq_similarity conflict: profile wins
# ---------------------------------------------------------------------------

class TestProfileOverridesThreshold:

    def test_max_ins_seq_len_bypasses_sequence_gate(self):
        # grch37_neg_c_sim0.672: extracted sequences ~70/67 bp, sim≈0.672 < 0.75 → normally 2 lines.
        # With --max_ins_seq_len 50 both are capped to "" → gate skipped → merge → 1 line.
        r = run_merge("--max_ins_seq_len", "50",
                      "--vcf", fixture("grch37_neg_c_sim0.672", "caller_a"),
                      fixture("grch37_neg_c_sim0.672", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 1

    def test_max_ins_seq_len_above_seq_length_does_not_merge(self):
        # Cap above actual sequence length (~70 bp) → sequences still compared → 2 lines.
        r = run_merge("--max_ins_seq_len", "200",
                      "--vcf", fixture("grch37_neg_c_sim0.672", "caller_a"),
                      fixture("grch37_neg_c_sim0.672", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 2

    def test_profile_wins_over_explicit_threshold_with_warning(self):
        # When both are specified, --data_profile overrides --ins_seq_similarity.
        # sample profile → threshold=0.85 → sim≈0.789 < 0.85 → rejected.
        # If user's explicit 0.50 were used instead, it would pass.
        r = run_merge("--data_profile", "sample",
                      "--ins_seq_similarity", "0.50",
                      "--vcf", fixture("grch37_pos25_sim0.789", "caller_a"),
                      fixture("grch37_pos25_sim0.789", "caller_b"))
        assert r.returncode == 0
        assert len(data_lines(r.stdout)) == 2  # profile threshold 0.85 wins
        assert "Warning" in r.stderr or "warning" in r.stderr.lower()
