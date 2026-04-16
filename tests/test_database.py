import pytest

from svdb.database import DB, CREATE_TABLE_SQL


@pytest.fixture
def db(tmp_path):
    """In-memory-equivalent DB backed by a temp file, pre-populated with two rows."""
    path = str(tmp_path / "test")
    d = DB(path)
    d.create(CREATE_TABLE_SQL)
    d.insert_many([
        ("DEL", "1", "1", 100, 0, 0, 200, 0, 0, "sample_A", 0),
        ("INS", "2", "2", 500, 0, 0, 500, 0, 0, "sample_B", 1),
    ])
    return d


class TestDBQueryColumn:

    def test_returns_flat_list(self, db):
        result = db.query_column("SELECT DISTINCT var FROM SVDB")
        assert set(result) == {"DEL", "INS"}

    def test_single_row(self, db):
        result = db.query_column("SELECT sample FROM SVDB WHERE var = 'DEL'")
        assert result == ["sample_A"]

    def test_empty_result(self, db):
        result = db.query_column("SELECT sample FROM SVDB WHERE var = 'BND'")
        assert result == []


class TestDBProperties:

    def test_tables_contains_svdb(self, db):
        assert "SVDB" in db.tables

    def test_sample_ids_returns_flat_list(self, db):
        assert set(db.sample_ids) == {"sample_A", "sample_B"}

    def test_len_counts_distinct_samples(self, db):
        assert len(db) == 2
