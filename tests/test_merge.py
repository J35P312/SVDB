import unittest
import numpy

from svdb.merge_vcf_module_cython import collect_info, skip_variant


class TestReadVCFLine(unittest.TestCase):

    #check that the info field is summarised properly
    def test_collect_info(self):
        info = ["chr1", "1" , "hej" , "." ,"<DEL>", "." , "PASS" ,"END=5;SVTYPE=DEL;TEST=1,2,3,4,5"]
        result=collect_info(info)
        assert (result=="hej|END:5|SVTYPE:DEL|TEST:1:2:3:4:5")

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




 
