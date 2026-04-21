"""Shared data structures for SVDB."""
from dataclasses import dataclass, field
from typing import Dict, List, NamedTuple


@dataclass(frozen=True)
class VCFVariant:
    """Parsed representation of a single VCF data line.

    Immutable: fields are set at parse time and must not change.
    Dict fields (info, fmt) are excluded from hashing because dicts
    are not hashable; equality and hash are based on coordinates and type.
    """
    chrA: str
    posA: int
    chrB: str
    posB: int
    event_type: str
    info: Dict[str, str] = field(hash=False)
    fmt: Dict[str, List[str]] = field(hash=False)

    def is_interchromosomal(self) -> bool:
        """True when the variant spans two different chromosomes (e.g. BND translocation)."""
        return self.chrA != self.chrB

    def is_insertion(self) -> bool:
        """True for insertion-type variants (INS, MEINS, etc.).

        Insertions are treated as single points in the coordinate model
        (posA == posB), so distance-based overlap does not apply.
        """
        return "INS" in self.event_type

    def is_precise(self) -> bool:
        """True for BND breakend variants represented as precise breakpoints."""
        return self.event_type == "BND"


class MergeVariant(NamedTuple):
    """A variant entry in the merge path variant dictionary."""
    chrB: str
    event_type: str
    posA: int
    posB: int
    source: str      # vcf file path or priority tag
    sort_index: int  # global sort index (renamed from index to avoid shadowing tuple.index)
    raw_line: str    # original VCF line, unparsed

    def is_insertion(self) -> bool:
        """True for insertion-type variants — mirrors VCFVariant.is_insertion()."""
        return "INS" in self.event_type
