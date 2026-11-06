"""Shared data structures for SVDB."""
import re
from dataclasses import dataclass, field
from typing import Dict, List, NamedTuple, Optional, Set, Tuple

_GT_SEP = re.compile(r"[/|]")


def combine_cnv_types(types: Set[str]) -> str:
    """One type for a set of DEL/DUP/CNV allele types: itself if unique, else CNV."""
    return next(iter(types)) if len(types) == 1 else "CNV"


def gt_type(alt_types: Tuple[str, ...], gt: str) -> Optional[str]:
    """Type carried by a genotype over alt_types; None if it carries no ALT allele."""
    types = {alt_types[int(a) - 1] for a in _GT_SEP.split(gt)
             if a.isdigit() and 0 < int(a) <= len(alt_types)}
    return combine_cnv_types(types) if types else None


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
    ins_seq: str = field(default="", compare=False, hash=False)
    svlen: Optional[int] = field(default=None, compare=False, hash=False)
    vcf_filter: str = field(default="", compare=False, hash=False)
    # per-ALT DEL/DUP/CNV types of a multi-ALT copy-number record; empty otherwise
    alt_types: Tuple[str, ...] = field(default=(), compare=False, hash=False)

    def type_for_gt(self, gt: str) -> Optional[str]:
        """Type carried by one sample's GT; None if a multi-ALT GT carries no ALT allele."""
        if not self.alt_types:
            return self.event_type
        return gt_type(self.alt_types, gt)

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
    source: str           # vcf file path or priority tag
    sort_index: int       # global sort index (renamed from index to avoid shadowing tuple.index)
    raw_line: str         # original VCF line, unparsed
    ins_seq: str = ""            # bare insertion sequence (empty for non-INS or symbolic ALT)
    svlen: Optional[int] = None  # |SVLEN|: from INFO, or derived from len(ins_seq) when absent
    vcf_filter: str = ""         # VCF FILTER column value

    def is_insertion(self) -> bool:
        """True for insertion-type variants — mirrors VCFVariant.is_insertion()."""
        return "INS" in self.event_type
