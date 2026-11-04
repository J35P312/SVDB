"""SV type matching for merge and query (--svtype_match).

exact: types must be equal.
hierarchical: a type also matches its parent or child, never its siblings.
Parents follow the VCF spec: colon-separated subtypes (INS:ME:ALU is an INS:ME
and an INS), and CNV covers DEL and DUP.
"""
from functools import lru_cache
from typing import FrozenSet, Iterable, List


@lru_cache(maxsize=None)
def lineage(svtype: str) -> FrozenSet[str]:
    """The type and all its parents: DEL:ME:LINE1 -> {DEL:ME:LINE1, DEL:ME, DEL, CNV}."""
    parts = svtype.split(":")
    types = {":".join(parts[:i]) for i in range(1, len(parts) + 1)}
    if parts[0] in ("DEL", "DUP"):
        types.add("CNV")
    return frozenset(types)


def svtypes_match(a: str, b: str, hierarchical: bool) -> bool:
    if a == b:
        return True
    return hierarchical and (a in lineage(b) or b in lineage(a))


def matching_svtypes(svtype: str, candidates: Iterable[str], hierarchical: bool) -> List[str]:
    """The candidate types that svtype matches."""
    return [c for c in candidates if svtypes_match(svtype, c, hierarchical)]
