# check the "overlap" of interchromosomaltranslocations
from typing import Optional, Tuple


def precise_overlap(
    chrApos_query: int, chrBpos_query: int,
    chrApos_db: int, chrBpos_db: int,
    distance: int,
) -> Tuple[Optional[float], bool]:
    """Return (max_breakpoint_distance, True) if both breakpoints are within distance, else (None, False)."""
    Adist = abs(chrApos_query - chrApos_db)
    Bdist = abs(chrBpos_query - chrBpos_db)
    if max([Adist, Bdist]) <= distance:
        return float(max([Adist, Bdist])), True
    return None, False


def isSameVariation(
    chrApos_query: int, chrBpos_query: int,
    chrApos_db: int, chrBpos_db: int,
    ratio: float,
    distance: int,
) -> Tuple[Optional[float], bool]:
    """Return (overlap_ratio, True) if variants overlap sufficiently, else (None, False)."""
    if abs(chrApos_query - chrApos_db) <= distance and abs(chrBpos_query - chrBpos_db) <= distance:
        region_start = min([chrApos_db, chrApos_query])
        overlap_start = max([chrApos_db, chrApos_query])

        region_end = max([chrBpos_db, chrBpos_query])
        overlap_end = min([chrBpos_db, chrBpos_query])

        event_ratio = float(overlap_end - overlap_start + 1) / float(region_end - region_start + 1)

        if event_ratio >= ratio:
            return event_ratio, True
        return None, False
    return None, False


def variant_overlap(
    chrA: str, chrB: str,
    chrApos_query: int, chrBpos_query: int,
    chrApos_db: int, chrBpos_db: int,
    ratio: float,
    distance: int,
) -> Tuple[Optional[float], bool]:
    """Dispatch to precise_overlap (interchromosomal) or isSameVariation (intrachromosomal)."""
    if chrA == chrB:
        return isSameVariation(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance)
    else:
        return precise_overlap(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, distance)
