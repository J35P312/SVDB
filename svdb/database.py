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

    def drop(self, query: str) -> None:
        self.cursor.execute(query)

    def create(self, query: str) -> None:
        self.cursor.execute(query)
        self.conn.commit()

    def insert_many(self, data: List[Tuple[Any, ...]]) -> None:
        self.cursor.executemany(f'INSERT INTO SVDB VALUES {INSERT_PLACEHOLDERS}', data)
        self.conn.commit()

    def create_index(self, name: str, columns: str) -> None:
        query = "CREATE INDEX {} ON SVDB {}".format(name, columns)
        self.cursor.execute(query)
        self.conn.commit()

    @property
    def tables(self) -> List[str]:
        res = self.query("SELECT name FROM sqlite_master WHERE type=\'table\'")
        return [table[0] for table in res]

    @property
    def sample_ids(self) -> List[Tuple[str, ...]]:
        return [sample for sample in self.query('SELECT DISTINCT sample FROM SVDB')]
