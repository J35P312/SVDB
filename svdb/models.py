"""Shared data structures for SVDB."""
from typing import NamedTuple


class VCFVariant(NamedTuple):
    """Parsed representation of a single VCF data line."""
    chrA: str
    posA: int
    chrB: str
    posB: int
    event_type: str
    info: dict
    fmt: dict


class MergeVariant(NamedTuple):
    """A variant entry in the merge path variant dictionary."""
    chrB: str
    event_type: str
    posA: int
    posB: int
    source: str    # vcf file path or priority tag
    index: int     # global sort index
    raw_line: str  # original VCF line, unparsed
