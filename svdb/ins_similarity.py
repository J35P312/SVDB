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

  cap_seq(seq, max_len) -> str
      Cap insertion sequence to max_len before comparison; returns "" if over
      the cap, which causes sequence_gate to defer to position+SVLEN.

  parse_svlen(info_str) -> Optional[int]
      Extract |SVLEN| from a VCF INFO string.

  apply_ins_profile(args) -> None
      Resolve all INS matching parameters onto args, applying --data_profile
      presets.  Profile sets the base; explicit flags override per-parameter.
      All args.ins_* attributes are guaranteed non-None after this call.
"""

import logging
import re
import zlib
from typing import Optional, Union

from rapidfuzz.distance import Levenshtein

logger = logging.getLogger(__name__)

# Effective defaults when no profile and no explicit flag is set.
_INS_DEFAULTS: dict = {
    "ins_distance": 25,
    "ins_svlen_ratio": 0.90,
    "ins_seq_similarity": 0.75,
    "no_ins_seq": False,
}

# Full presets.  Profile values are overridden by any individually specified flag.
_PROFILES: dict = {
    "sample": {
        "ins_distance": 25,
        "ins_svlen_ratio": 0.90,
        "ins_seq_similarity": 0.85,
        "no_ins_seq": False,
    },
    "cohort": {
        "ins_distance": 50,
        "ins_svlen_ratio": 0.80,
        "ins_seq_similarity": 0.75,
        "no_ins_seq": False,
    },
    "position_only": {
        "ins_distance": 50,
        "ins_svlen_ratio": 0.90,
        "ins_seq_similarity": 0.75,  # unused when no_ins_seq=True
        "no_ins_seq": True,
    },
}


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


def apply_ins_profile(args) -> None:
    """Resolve INS matching parameters from --data_profile and explicit flags.

    --data_profile sets the base for all INS parameters; any individually
    specified --ins_* flag overrides the profile for that parameter only.
    All args.ins_* attributes are guaranteed non-None after this call.
    """
    profile_name = getattr(args, "data_profile", None)
    profile = _PROFILES.get(profile_name, {})

    # no_ins_seq: always derived from profile (position_only sets it to True)
    args.no_ins_seq = profile.get("no_ins_seq", False)

    # Numeric params: None = not explicitly set → use profile value then default
    for key in ("ins_distance", "ins_svlen_ratio", "ins_seq_similarity"):
        if getattr(args, key, None) is None:
            setattr(args, key, profile.get(key, _INS_DEFAULTS[key]))
