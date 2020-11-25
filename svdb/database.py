import sqlite3


class DB:
    def __init__(self, db, memory=False):
        if not db.endswith('.db'):
            db += '.db'

        conn = sqlite3.connect(db)

        if memory:
            memory_db = sqlite3.connect(':memory:')
            db_dump = "".join(line for line in conn.iterdump())
            memory_db.executescript(db_dump)
            conn.close()
            self.cursor = memory_db.cursor()
        else:
            self.cursor = conn.cursor()

    def __len__(self):
        return len(self.query('SELECT DISTINCT sample FROM SVDB'))

    def query(self, query):
        self.cursor.execute(query)
        res = self.cursor.fetchall()
        return res

    def drop(self, query):
        self.cursor.execute(query)
        try:
            self.cursor.execute(query)
        except Exception:
            pass

    def create(self, query):
        self.cursor.execute(query)
        self.conn.commit()

    def insert_many(self, data):
        self.cursor.executemany('INSERT INTO SVDB VALUES (?,?,?,?,?,?,?,?,?,?,?)', data)
        self.conn.commit()

    def create_index(self, name, columns):
        self.cursor.execute("CREATE INDEX {} ON SVDB {}".format(name, columns))
        self.conn.commit()

    @property
    def tables(self):
        res = self.query("SELECT name FROM sqlite_master WHERE type=\'table\'")
        return [table[0] for table in res]

    @property
    def sample_ids(self):
        return [sample for sample in self.query('SELECT DISTINCT sample FROM SVDB')]
