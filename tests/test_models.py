import pytest

from svdb.models import MergeVariant, VCFVariant


class TestVCFVariant:

    def _make(self, chrA="1", posA=100, chrB="1", posB=200, event_type="DEL"):
        return VCFVariant(chrA=chrA, posA=posA, chrB=chrB, posB=posB,
                          event_type=event_type, info={}, fmt={})

    # --- is_interchromosomal ---

    def test_interchromosomal_true_for_different_chromosomes(self):
        v = self._make(chrA="1", chrB="2", event_type="BND")
        assert v.is_interchromosomal()

    def test_interchromosomal_false_for_same_chromosome(self):
        v = self._make(chrA="1", chrB="1", event_type="DEL")
        assert not v.is_interchromosomal()

    # --- is_insertion ---

    def test_insertion_true_for_ins(self):
        v = self._make(event_type="INS")
        assert v.is_insertion()

    def test_insertion_true_for_mobile_element(self):
        # MEINS and similar subtypes must also match
        v = self._make(event_type="MEINS")
        assert v.is_insertion()

    def test_insertion_false_for_deletion(self):
        v = self._make(event_type="DEL")
        assert not v.is_insertion()

    def test_insertion_false_for_bnd(self):
        v = self._make(chrA="1", chrB="2", event_type="BND")
        assert not v.is_insertion()

    # --- is_precise ---

    def test_precise_true_for_bnd(self):
        v = self._make(chrA="1", chrB="2", event_type="BND")
        assert v.is_precise()

    def test_precise_false_for_del(self):
        v = self._make(event_type="DEL")
        assert not v.is_precise()

    def test_precise_false_for_ins(self):
        v = self._make(event_type="INS")
        assert not v.is_precise()

    # --- immutability ---

    def test_frozen_raises_on_field_assignment(self):
        v = self._make()
        with pytest.raises(AttributeError):
            v.chrA = "2"

    def test_frozen_raises_on_new_attribute(self):
        v = self._make()
        with pytest.raises(AttributeError):
            v.new_field = "x"


class TestMergeVariant:

    def _make(self, event_type="DEL"):
        return MergeVariant(chrB="1", event_type=event_type,
                            posA=100, posB=200, source="test.vcf",
                            index=0, raw_line="raw")

    def test_insertion_true_for_ins(self):
        v = self._make(event_type="INS")
        assert v.is_insertion()

    def test_insertion_true_for_mobile_element(self):
        v = self._make(event_type="MEINS")
        assert v.is_insertion()

    def test_insertion_false_for_deletion(self):
        v = self._make(event_type="DEL")
        assert not v.is_insertion()

    def test_insertion_false_for_dup(self):
        v = self._make(event_type="DUP")
        assert not v.is_insertion()
