from __future__ import absolute_import

# check the "overlap" of interchromosomaltranslocations


def precise_overlap(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, distance):
    Adist = abs(chrApos_query - chrApos_db)
    Bdist = abs(chrBpos_query - chrBpos_db)
    if max([Adist, Bdist]) <= distance:
        return max([Adist, Bdist]), True
    return False, False

# check if intrachromosomal vaiants overlap


# event is in the DB, variation is the new variation I want to insert
def isSameVariation(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance):
    if abs(chrApos_query - chrApos_db) <= distance and abs(chrBpos_query - chrBpos_db) <= distance:

        region_start = min([chrApos_db, chrApos_query])
        overlap_start = max([chrApos_db, chrApos_query])

        region_end = max([chrBpos_db, chrBpos_query])
        overlap_end = min([chrBpos_db, chrBpos_query])

        try:
            event_ratio = float(overlap_end - overlap_start + 1) / \
                float(region_end - region_start + 1)
        except Exception:
            event_ratio = 0

        if event_ratio >= ratio:
            return event_ratio, True
        return None, False
    else:
        return None, False


def variant_overlap(chrA, chrB, chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance):
    match = False
    overlap = False
    if chrA == chrB:
        overlap, match = isSameVariation(
            chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance)
    else:
        overlap, match = precise_overlap(
            chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, distance)
    return overlap, match
