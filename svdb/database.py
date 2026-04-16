import sqlite3
from typing import Any, List, Tuple

SCHEMA_COLUMNS = (
    "var TEXT",
    "chrA TEXT",
    "chrB TEXT",
    "posA INT",
    "ci_A_lower INT",
    "ci_A_upper INT",
    "posB INT",
    "ci_B_lower INT",
    "ci_B_upper INT",
    "sample TEXT",
    "idx INT",
)

CREATE_TABLE_SQL = "CREATE TABLE SVDB ({})".format(", ".join(SCHEMA_COLUMNS))
INSERT_PLACEHOLDERS = "({})".format(", ".join("?" for _ in SCHEMA_COLUMNS))


class DB:
    def __init__(self, db: str, memory: bool = False) -> None:
        if not db.endswith('.db'):
            db += '.db'

        self.conn = sqlite3.connect(db)

        if memory:
            memory_db = sqlite3.connect(':memory:')
            db_dump = "".join(line for line in self.conn.iterdump())
            memory_db.executescript(db_dump)
            self.conn.close()
            self.cursor = memory_db.cursor()
        else:
            self.cursor = self.conn.cursor()

    def __len__(self) -> int:
        return len(self.query('SELECT DISTINCT sample FROM SVDB'))

    def query(self, query: str) -> List[Tuple[Any, ...]]:
        self.cursor.execute(query)
        res = self.cursor.fetchall()
        return res

    def query_column(self, query: str) -> List[Any]:
        """Run a single-column SELECT and return a flat list of values."""
        return [row[0] for row in self.query(query)]

    def drop(self, query: str) -> None:
        self.cursor.execute(query)

    def create(self, query: str) -> None:
        self.cursor.execute(query)
        self.conn.commit()

    def insert_many(self, data: List[Tuple[Any, ...]]) -> None:
        self.cursor.executemany(f'INSERT INTO SVDB VALUES {INSERT_PLACEHOLDERS}', data)
        self.conn.commit()

    def create_index(self, name: str, columns: str) -> None:
        query = f"CREATE INDEX {name} ON SVDB {columns}"
        self.cursor.execute(query)
        self.conn.commit()

    @property
    def tables(self) -> List[str]:
        return self.query_column("SELECT name FROM sqlite_master WHERE type='table'")

    @property
    def sample_ids(self) -> List[str]:
        return self.query_column('SELECT DISTINCT sample FROM SVDB')
