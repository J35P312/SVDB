import unittest
import numpy

from svdb.merge_vcf_module_cython import collect_info, skip_variant, retrieve_key, collect_sample


class TestMerge(unittest.TestCase):

    #check that we find and retrieve correct entry of info column
    def test_retrieve_key(self):
        line="Y\t13799001\tCNVnator_dup_810:concatenated_ACC5821A7_XXXXXX_R_CNVnator|CNVnator_dup_1313:concatenated_ACC5838A1_XXXXXX_R_CNVnator\tN\t<DUP>\t.\tPASS\tEND=13870000;SVTYPE=DUP;SVLEN=71000;IMPRECISE;natorRD=25.6613;natorP1=0.000938744;natorP2=1.33071e-34;natorP3=0.00215104;natorP4=2.21186e-33;natorQ0=1;VARID=CNVnator_dup_1313:concatenated_ACC5838A1_XXXXXX_R_CNVnator;set=Intersection;concatenated_ACC5821A7_XXXXXX_R_CNVnator_FILTERS=CNVnator_dup_810|PASS;concatenated_ACC5838A1_XXXXXX_R_CNVnator_FILTERS=CNVnator_dup_1313|PASS;concatenated_ACC5821A7_XXXXXX_R_CNVnator_SAMPLES=CNVnator_dup_810|concatenated_ACC5821A7_XXXXXX_R_CNVnator|GT:./1|CN:26;concatenated_ACC5838A1_XXXXXX_R_CNVnator_SAMPLES=CNVnator_dup_1313|concatenated_ACC5838A1_XXXXXX_R_CNVnator|GT:./1|CN:35;concatenated_ACC5821A7_XXXXXX_R_CNVnator_INFO=CNVnator_dup_810|END:13870000|SVTYPE:DUP|SVLEN:71000|IMPRECISE|natorRD:25.6613|natorP1:0.000938744|natorP2:1.33071e-34|natorP3:0.00215104|natorP4:2.21186e-33|natorQ0:1;concatenated_ACC5838A1_XXXXXX_R_CNVnator_INFO=CNVnator_dup_1313|END:13870000|SVTYPE:DUP|SVLEN:71000|IMPRECISE|natorRD:35.4121|natorP1:0.000300039|natorP2:0|natorP3:0.00099868|natorP4:0|natorQ0:1\tGT:CN"
        key="SVLEN"
        assert(retrieve_key(line, key) == "71000")

    def test_retrieve_key2(self):
        line="Y\t13799001\tCNVnator_dup_810:concatenated_ACC5821A7_XXXXXX_R_CNVnator|CNVnator_dup_1313:concatenated_ACC5838A1_XXXXXX_R_CNVnator\tN\t<DUP>\t.\tPASS\tEND=13870000;SVTYPE=DUP;SVLEN=71000;IMPRECISE;natorRD=25.6613;natorP1=0.000938744;natorP2=1.33071e-34;natorP3=0.00215104;natorP4=2.21186e-33;natorQ0=1;VARID=CNVnator_dup_1313:concatenated_ACC5838A1_XXXXXX_R_CNVnator;set=Intersection;concatenated_ACC5821A7_XXXXXX_R_CNVnator_FILTERS=CNVnator_dup_810|PASS;concatenated_ACC5838A1_XXXXXX_R_CNVnator_FILTERS=CNVnator_dup_1313|PASS;concatenated_ACC5821A7_XXXXXX_R_CNVnator_SAMPLES=CNVnator_dup_810|concatenated_ACC5821A7_XXXXXX_R_CNVnator|GT:./1|CN:26;concatenated_ACC5838A1_XXXXXX_R_CNVnator_SAMPLES=CNVnator_dup_1313|concatenated_ACC5838A1_XXXXXX_R_CNVnator|GT:./1|CN:35;concatenated_ACC5821A7_XXXXXX_R_CNVnator_INFO=CNVnator_dup_810|END:13870000|SVTYPE:DUP|SVLEN:71000|IMPRECISE|natorRD:25.6613|natorP1:0.000938744|natorP2:1.33071e-34|natorP3:0.00215104|natorP4:2.21186e-33|natorQ0:1;concatenated_ACC5838A1_XXXXXX_R_CNVnator_INFO=CNVnator_dup_1313|END:13870000|SVTYPE:DUP|SVLEN:71000|IMPRECISE|natorRD:35.4121|natorP1:0.000300039|natorP2:0|natorP3:0.00099868|natorP4:0|natorQ0:1\tGT:CN"
        key="VLEN"
        assert(retrieve_key(line, key) == False)
 
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




 
