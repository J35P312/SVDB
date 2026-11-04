from types import SimpleNamespace

import pytest

from svdb.database import DB, CREATE_TABLE_SQL
from svdb.query_module import SQDB


@pytest.fixture
def db(tmp_path):
    path = str(tmp_path / "test")
    with DB(path) as d:
        d.create(CREATE_TABLE_SQL)
        d.insert_many([
            ("DEL", "1", "1", 100, 0, 0, 200, 0, 0, "sample_A", 0),
        ])
        yield d


def _args():
    return SimpleNamespace(
        ins_distance=1000, bnd_distance=1000, overlap=0.8,
        no_ins_seq=True, ins_svlen_ratio=0.8, ins_seq_similarity=0.8,
        max_ins_seq_len=1000,
    )


def _query(chrA="1", posA=100, chrB="1", posB=200, event="DEL"):
    return [chrA, posA, chrB, posB, event, None, "line", "", None]


def test_sqdb_matches_plain_query(db):
    assert SQDB(_query(), _args(), db) == 1


def test_sqdb_quote_in_chromosome_is_literal(db):
    assert SQDB(_query(chrA="1' OR '1'='1"), _args(), db) == 0


def test_sqdb_quote_in_type_is_literal(db):
    assert SQDB(_query(event="DEL' OR '1'='1"), _args(), db) == 0


@pytest.mark.parametrize("event, db_svtypes, expected", [
    ("CNV", None, 0),        # exact: CNV != DEL
    ("CNV", ["DEL"], 1),     # hierarchical: parent matches child
    ("DUP", ["DEL"], 0),     # sibling: no db type matches, empty IN ()
])
def test_sqdb_svtype_match(db, event, db_svtypes, expected):
    assert SQDB(_query(event=event), _args(), db, db_svtypes=db_svtypes) == expected


def test_db_query_accepts_params(db):
    rows = db.query(
        "SELECT sample FROM SVDB WHERE var = ? AND chrA = ?", ("DEL", "1")
    )
    assert rows == [("sample_A",)]
    rows = db.query(
        "SELECT sample FROM SVDB WHERE chrA = ?", ("1' OR '1'='1",)
    )
    assert rows == []
