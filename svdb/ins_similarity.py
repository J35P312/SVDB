"""Sequence-based gate for insertion merging (issue #82).

Public API used by merge, query, and export pipelines:

  extract_ins_sequence(ref, alt) -> str
      Normalise a VCF REF/ALT pair to a bare insertion sequence.
      Returns "" for symbolic ALTs (e.g. <INS>) — signal to skip sequence check.

  levenshtein_similarity(seq_a, seq_b) -> float
      Normalised Levenshtein similarity in [0, 1] via rapidfuzz.

  sequence_gate(seq_a, seq_b, threshold) -> bool
      True = allow merge; False = reject.  Returns True when either sequence is
      empty (symbolic ALT or absent), deferring to position+SVLEN.

  parse_svlen(info_str) -> Optional[int]
      Extract |SVLEN| from a VCF INFO string.

  resolve_ins_seq_threshold(args) -> float
      Resolve effective sequence similarity threshold from args, applying
      --data_profile presets and emitting a warning on conflict.
"""

import logging
import re
import zlib
from typing import Optional, Union

from rapidfuzz.distance import Levenshtein

logger = logging.getLogger(__name__)

_DATA_PROFILE_THRESHOLDS = {
    "sample": 0.85,
    "cohort": 0.75,
}
_DEFAULT_INS_SEQ_SIMILARITY = 0.75


def decompress_ins_seq(raw: Union[bytes, str, None]) -> Optional[str]:
    """Recover a stored ins_seq value to a plain string.

    New databases store ins_seq as a zlib-compressed BLOB (bytes); legacy
    databases store it as TEXT (str).  Returns None when raw is None.
    """
    if raw is None:
        return None
    if isinstance(raw, bytes):
        return zlib.decompress(raw).decode()
    return raw


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


def levenshtein_similarity(seq_a: str, seq_b: str, score_cutoff: float = 0.0) -> float:
    """Normalised Levenshtein similarity in [0, 1].

    Uses rapidfuzz for efficiency (~14–33 µs/pair on typical insertion lengths).
    Two empty strings are considered identical (returns 1.0).
    score_cutoff enables early termination: returns 0.0 immediately when the
    similarity cannot reach the cutoff, avoiding the full O(n*m) computation.
    """
    return Levenshtein.normalized_similarity(seq_a, seq_b, score_cutoff=score_cutoff)


def sequence_gate(seq_a: str, seq_b: str, threshold: float) -> bool:
    """Return True if the merge is allowed, False if it should be rejected.

    Skips the similarity check (returns True) when either sequence is absent —
    this preserves the existing position+SVLEN merge behaviour for callers that
    emit symbolic ALTs (<INS>, <INS:ME>, etc.).
    """
    if not seq_a or not seq_b:
        return True
    return levenshtein_similarity(seq_a, seq_b, score_cutoff=threshold) >= threshold


def cap_seq(seq: Optional[str], max_len: Optional[int]) -> str:
    """Return seq unchanged if within max_len; return "" if it exceeds the cap.

    An empty return causes sequence_gate to skip similarity and defer to
    position+SVLEN matching — consistent with symbolic ALT behaviour.
    """
    if seq and max_len is not None and len(seq) > max_len:
        return ""
    return seq or ""


def parse_svlen(info_str: str) -> Optional[int]:
    """Extract |SVLEN| from a VCF INFO string. Returns None if absent or unparseable."""
    m = re.search(r'(?:^|;)SVLEN=(-?\d+)', info_str)
    if m:
        return abs(int(m.group(1)))
    return None


def resolve_ins_seq_threshold(args) -> float:
    """Return the effective insertion sequence similarity threshold.

    --data_profile takes precedence over an explicit --ins_seq_similarity.
    If both are supplied, emits a warning and uses the profile threshold.
    """
    profile = getattr(args, "data_profile", None)
    explicit = getattr(args, "ins_seq_similarity", None)

    if profile is not None and explicit is not None:
        resolved = _DATA_PROFILE_THRESHOLDS[profile]
        logger.warning(
            "Both --data_profile %s and --ins_seq_similarity %.2f were specified. "
            "--data_profile takes precedence; using threshold %.2f.",
            profile, explicit, resolved,
        )
        return resolved
    if profile is not None:
        return _DATA_PROFILE_THRESHOLDS[profile]
    if explicit is not None:
        return explicit
    return _DEFAULT_INS_SEQ_SIMILARITY
