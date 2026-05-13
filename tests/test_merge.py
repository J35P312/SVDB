import argparse
import unittest

from svdb.merge_vcf_module_cython import collect_info, skip_variant, collect_sample, sanitize_id, format_tag, sort_format_field


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


class TestSortFormatField(unittest.TestCase):

    def _args(self):
        return argparse.Namespace(same_order=False)

    def test_format_union_includes_lower_priority_tags(self):
        """Lower-priority caller's unique FORMAT tags must appear in the merged FORMAT."""
        # manta (priority 1): FORMAT GT:GQ
        manta = "chr1\t100\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT:GQ\t0/1:30"
        # lumpy (priority 2): FORMAT GT:SU — SU is lumpy-specific
        lumpy = "chr1\t100\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT:SU\t0/1:15"

        samples = ["sample_A"]
        sample_order = {"sample_A": {"manta.vcf": 0, "lumpy.vcf": 0}}
        priority_order = ["manta.vcf", "lumpy.vcf"]
        files = {"manta.vcf": manta, "lumpy.vcf": lumpy}
        line = manta.split("\t")

        result = sort_format_field(line, samples, sample_order, priority_order, files, self._args())

        fmt_keys = result[8].split(":")
        assert "GT" in fmt_keys
        assert "GQ" in fmt_keys
        assert "SU" in fmt_keys, "lumpy-specific SU tag was dropped (pre-fix regression)"

        # sample values: GT and GQ come from manta; SU not in manta → "."
        vals = result[9].split(":")
        assert vals[fmt_keys.index("GT")] == "0/1"
        assert vals[fmt_keys.index("GQ")] == "30"
        assert vals[fmt_keys.index("SU")] == "."

    def test_format_placeholder_uses_correct_width(self):
        """Missing multi-value FORMAT entry produces the right number of '.' placeholders."""
        # caller1 (priority 1): sample_A with FORMAT GT:AD, AD is 3-value
        c1 = "chr1\t100\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT:AD\t0/1:10,20,5"
        # caller2 (priority 2): sample_B with FORMAT GT:DP
        c2 = "chr1\t100\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT:DP\t0/1:30"

        samples = ["sample_A", "sample_B"]
        sample_order = {"sample_A": {"c1.vcf": 0}, "sample_B": {"c2.vcf": 0}}
        priority_order = ["c1.vcf", "c2.vcf"]
        files = {"c1.vcf": c1, "c2.vcf": c2}
        line = c1.split("\t")

        result = sort_format_field(line, samples, sample_order, priority_order, files, self._args())

        fmt_keys = result[8].split(":")
        assert fmt_keys == ["GT", "AD", "DP"]

        # sample_A: GT=0/1, AD=10,20,5 from c1; DP missing → "."
        a = result[9].split(":")
        assert a[fmt_keys.index("GT")] == "0/1"
        assert a[fmt_keys.index("AD")] == "10,20,5"
        assert a[fmt_keys.index("DP")] == "."

        # sample_B: GT=0/1, DP=30 from c2; AD missing → ".,.,." (3 values, not ".")
        b = result[10].split(":")
        assert b[fmt_keys.index("GT")] == "0/1"
        assert b[fmt_keys.index("AD")] == ".,.,.", "j=0 bug: AD placeholder width used wrong index"
        assert b[fmt_keys.index("DP")] == "30"





