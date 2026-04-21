import gzip


from svdb.vcf_utils import normalize_chrom, open_vcf, parse_ci, parse_info_field


class TestOpenVcf:

    def test_reads_plain_vcf(self, tmp_path):
        vcf = tmp_path / "test.vcf"
        vcf.write_text("line1\nline2\n")
        with open_vcf(str(vcf)) as f:
            lines = f.readlines()
        assert lines == ["line1\n", "line2\n"]

    def test_reads_gzip_vcf(self, tmp_path):
        vcf = tmp_path / "test.vcf.gz"
        with gzip.open(str(vcf), "wt") as f:
            f.write("line1\nline2\n")
        with open_vcf(str(vcf)) as f:
            lines = f.readlines()
        assert lines == ["line1\n", "line2\n"]

    def test_plain_file_is_context_manager(self, tmp_path):
        vcf = tmp_path / "test.vcf"
        vcf.write_text("x\n")
        with open_vcf(str(vcf)) as f:
            assert not f.closed
        assert f.closed


class TestNormalizeChrom:

    def test_strips_lowercase_chr(self):
        assert normalize_chrom("chr1") == "1"

    def test_strips_mixed_case_Chr(self):
        assert normalize_chrom("Chr1") == "1"

    def test_strips_uppercase_CHR(self):
        assert normalize_chrom("CHR1") == "1"

    def test_leaves_plain_name_unchanged(self):
        assert normalize_chrom("1") == "1"

    def test_leaves_X_unchanged(self):
        assert normalize_chrom("X") == "X"

    def test_handles_chrX(self):
        assert normalize_chrom("chrX") == "X"


class TestParseInfoField:

    def test_single_key_value(self):
        assert parse_info_field("SVTYPE=DEL") == {"SVTYPE": "DEL"}

    def test_multiple_key_values(self):
        result = parse_info_field("END=84000;SVTYPE=DEL;SVLEN=-84000")
        assert result == {"END": "84000", "SVTYPE": "DEL", "SVLEN": "-84000"}

    def test_flag_entries_are_skipped(self):
        result = parse_info_field("IMPRECISE;SVTYPE=DEL")
        assert result == {"SVTYPE": "DEL"}
        assert "IMPRECISE" not in result

    def test_value_containing_equals_is_preserved(self):
        # e.g. CSQ values or base64 can contain '='
        result = parse_info_field("FOO=bar=baz")
        assert result == {"FOO": "bar=baz"}

    def test_empty_string_returns_empty_dict(self):
        assert parse_info_field("") == {}


class TestParseCi:

    def test_two_values(self):
        assert parse_ci("100,200") == (100, 200)

    def test_single_value_duplicated(self):
        assert parse_ci("300") == (300, 300)

    def test_strips_parentheses(self):
        assert parse_ci("(-50,50)") == (50, 50)

    def test_negative_values_become_absolute(self):
        assert parse_ci("-100,200") == (100, 200)

    def test_single_negative_value(self):
        assert parse_ci("-300") == (300, 300)
