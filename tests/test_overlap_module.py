import unittest

from svdb.overlap_module import isSameVariation, precise_overlap, variant_overlap


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

