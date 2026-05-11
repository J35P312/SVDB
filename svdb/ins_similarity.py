"""Sequence-based gate for insertion merging (issue #82).

Provides three public functions used by the merge pipeline:

  extract_ins_sequence(ref, alt) -> str
      Normalise a VCF REF/ALT pair to a bare insertion sequence.
      Returns "" for symbolic ALTs (e.g. <INS>) — the caller signal to skip
      the sequence check.

  levenshtein_similarity(seq_a, seq_b) -> float
      Normalised Levenshtein similarity in [0, 1] via rapidfuzz.

  sequence_gate(seq_a, seq_b, threshold) -> bool
      True = allow the merge to proceed; False = reject.
      Returns True immediately when either sequence is empty (symbolic ALT
      or otherwise absent), deferring the merge decision to position+SVLEN.
"""

from rapidfuzz.distance import Levenshtein


def extract_ins_sequence(ref: str, alt: str) -> str:
    """Return the bare insertion sequence from a VCF REF/ALT pair.

    Handles three encodings:
      - Symbolic  (ALT starts with '<'):  return ""
      - Anchored  (ALT[0] == REF[0]):     strip the anchor base, return the rest
      - Single-base ALT matching REF:     return "" (no inserted sequence)
      - Unanchored (Sniffles v1 style):   return ALT as-is
    """
    if alt.startswith("<"):
        return ""
    if len(alt) == 1:
        return ""
    if alt[0] == ref[0]:
        return alt[1:]
    return alt


def levenshtein_similarity(seq_a: str, seq_b: str) -> float:
    """Normalised Levenshtein similarity in [0, 1].

    Uses rapidfuzz for efficiency (~14–33 µs/pair on typical insertion lengths).
    Two empty strings are considered identical (returns 1.0).
    """
    return Levenshtein.normalized_similarity(seq_a, seq_b)


def sequence_gate(seq_a: str, seq_b: str, threshold: float) -> bool:
    """Return True if the merge is allowed, False if it should be rejected.

    Skips the similarity check (returns True) when either sequence is absent —
    this preserves the existing position+SVLEN merge behaviour for callers that
    emit symbolic ALTs (<INS>, <INS:ME>, etc.).
    """
    if not seq_a or not seq_b:
        return True
    return levenshtein_similarity(seq_a, seq_b) >= threshold
