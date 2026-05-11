import unittest

from svdb.merge_vcf_module_cython import collect_info, skip_variant, collect_sample, sanitize_id, format_tag


class TestSanitizeId(unittest.TestCase):

    def test_replaces_semicolons(self):
        assert sanitize_id("a;b") == "a_b"

    def test_replaces_colons(self):
        assert sanitize_id("a:b") == "a_b"

    def test_replaces_pipes(self):
        assert sanitize_id("a|b") == "a_b"

    def test_replaces_all_three_delimiters(self):
        assert sanitize_id("a;b:c|d") == "a_b_c_d"

    def test_plain_string_unchanged(self):
        assert sanitize_id("SV001") == "SV001"


class TestFormatTag(unittest.TestCase):

    def test_basic(self):
        assert format_tag("SV001", "chr1") == "SV001|chr1"

    def test_sanitizes_var_id(self):
        assert format_tag("SV;001:x|y", "chr1") == "SV_001_x_y|chr1"

    def test_value_is_not_sanitized(self):
        # the value side (e.g. a chromosome name) is kept verbatim
        assert format_tag("SV001", "chr1;extra") == "SV001|chr1;extra"


class TestMerge(unittest.TestCase):

    #check that the info field is summarised properly
    def test_collect_info(self):
        info = ["chr1", "1" , "hej" , "." ,"<DEL>", "." , "PASS" ,"END=5;SVTYPE=DEL;TEST=1,2,3,4,5"]
        result=collect_info(info)
        assert (result=="hej|END:5|SVTYPE:DEL|TEST:1:2:3:4:5")

    #check that sample columns are retrieved properly
    def test_collect_collect_sample(self):
        vcf_line = ["chr1", "1" , "hej" , "." ,"<DEL>", "." , "PASS" ,"END=5;SVTYPE=DEL;TEST=1,2,3,4,5","GT:CN","1/1:0"]
        samples=["bob"]
        sample_order={"bob":{"cnvnator_bob":0}}
        f="cnvnator_bob"
        result=collect_sample(vcf_line,samples,sample_order,f)
        assert (result=="hej|bob|GT:1/1|CN:0")


    #test the skip_variant filter
    def test_skip_variant_different_var(self):
        no_var=False
        analysed_variants=set([1,2])
        current_variant=3
        pass_only=False
        vcf_line_A=["chr1", "1" , "hej" , "." ,"<DEL>", "." , "PASS" ,"END=5;SVTYPE=DEL;TEST=1,2,3,4,5"]
        vcf_line_B=["chr1", "1" , "hej" , "." ,"<DEL>", "." , "FAIL" ,"END=5;SVTYPE=DEL;TEST=1,2,3,4,5"]
        type_A="DEL"
        type_B="DEL"
        chrA="chr1"
        chrB="chr1"

        #Do not merge variants of different type
        result=skip_variant(chrA,chrB,"DUP",type_B,vcf_line_A,vcf_line_B,pass_only,current_variant,analysed_variants,no_var)
        assert (result)

        #merge variants of different types if no_var is True
        result=skip_variant(chrA,chrB,"DUP",type_B,vcf_line_A,vcf_line_B,pass_only,current_variant,analysed_variants,True)
        assert (not result)

        #Do not cluster already clustered variants
        result=skip_variant(chrA,chrB,type_A,type_B,vcf_line_A,vcf_line_B,pass_only,2,analysed_variants,True)
        assert (result)

        #Do not cluster variants located on different chromosomes
        result=skip_variant(chrA,"X",type_A,type_B,vcf_line_A,vcf_line_B,pass_only,current_variant,analysed_variants,True)
        assert (result)

        #Skip filtered variants (if pass_only =True)
        result=skip_variant(chrA,"X",type_A,type_B,vcf_line_A,vcf_line_B,True,current_variant,analysed_variants,True)
        assert (result)





