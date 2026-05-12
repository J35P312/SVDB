import unittest

from svdb.overlap_module import isSameVariation, precise_overlap, variant_overlap, insertion_svlen_match


class TestOverlapModule(unittest.TestCase):

    #test that disjunct SV return 0 overlap
    def test_disjunct(self):
        chrApos_query=1
        chrBpos_query=10
        chrApos_db=11
        chrBpos_db=20
        ratio=0.5
        distance=1000
        assert( not isSameVariation(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance)[1] )

    #test that identical SV return 1 overlap
    def test_identical(self):
        chrApos_query=1
        chrBpos_query=10
        chrApos_db=1
        chrBpos_db=10
        ratio=0.5
        distance=1000
        assert( 1 == isSameVariation(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance)[0] )

    #test that SV sharing 50% of bases return 0.5 overlap
    def test_semi_overlap(self):
        chrApos_query=1
        chrBpos_query=10
        chrApos_db=6
        chrBpos_db=10
        ratio=0.5
        distance=1000
        assert( 0.5 == isSameVariation(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance)[0] )

    #test that large sv differing more than the set "distance" treshold are not merged
    def test_large_similar(self):
        chrApos_query=1
        chrBpos_query=1000000
        chrApos_db=1002
        chrBpos_db=1000000
        ratio=0.5
        distance=1000
        assert( not isSameVariation(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance)[1] )

    #test that precise variants closer than the set distance are considered the same
    def test_similar_precise(self):
        chrApos_query=1000
        chrBpos_query=10000
        chrApos_db=1100
        chrBpos_db=10100
        distance=200

        assert(precise_overlap(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, distance))

    #test that precise variants more distant than the than the set distance are considered different
    def test_different_precise(self):
        chrApos_query=1000
        chrBpos_query=10000
        chrApos_db=1100
        chrBpos_db=10100
        distance=50

        assert(precise_overlap(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, distance))

    #test the variant_overlap overlap module
    def test_variant_overlap(self):
        chrA="1"
        chrB="1"
        chrApos_query=1
        chrBpos_query=10
        chrApos_db=1
        chrBpos_db=10
        ratio=0.5
        distance=1000
        assert( (1,True) == variant_overlap(chrA,chrB,chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance) )

    def test_variant_overlap_interchromosomal(self):
        """When chrA != chrB, variant_overlap dispatches to precise_overlap (line 52)."""
        overlap, match = variant_overlap("1", "2", 1000, 5000, 1050, 5050, 0.5, 200)
        assert match is True
        assert overlap == 50.0  # max(abs(1000-1050), abs(5000-5050)) = 50

    def test_variant_overlap_interchromosomal_miss(self):
        """Interchromosomal breakpoints beyond bnd_distance do not match."""
        _, match = variant_overlap("1", "2", 1000, 5000, 2000, 5000, 0.5, 500)
        assert match is False


class TestInsertionSvlenMatch(unittest.TestCase):

    def test_identical_svlen_matches(self):
        assert insertion_svlen_match(100, 100, 0.90) is True

    def test_ratio_exactly_at_threshold(self):
        # 90/100 = 0.90 — exactly at the default threshold
        assert insertion_svlen_match(90, 100, 0.90) is True

    def test_ratio_just_below_threshold(self):
        # 89/100 = 0.89 < 0.90
        assert insertion_svlen_match(89, 100, 0.90) is False

    def test_ratio_is_symmetric(self):
        # 100/89 = 0.89 — order of arguments must not matter
        assert insertion_svlen_match(100, 89, 0.90) is False

    def test_loose_threshold_accepts_low_ratio(self):
        assert insertion_svlen_match(50, 100, 0.40) is True

    def test_fixture_low_svlen_ratio_rejected(self):
        # low_svlen_ratio fixture: SVLEN=100/220 → ratio=100/220≈0.455 < 0.90
        assert insertion_svlen_match(100, 220, 0.90) is False

    def test_fixture_grch37_pos25_svlen_passes(self):
        # grch37_pos25_sim0.789: SVLEN=51/52 → ratio≈0.981 ≥ 0.90
        assert insertion_svlen_match(51, 52, 0.90) is True

    def test_fixture_neg_c_svlen_passes(self):
        # grch37_neg_c_sim0.837: SVLEN=180/190 → ratio≈0.947 ≥ 0.90
        assert insertion_svlen_match(180, 190, 0.90) is True

