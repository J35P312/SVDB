import sqlite3


class DB:
    def __init__(self, db):
        if not db.endswith('.db'):
            db += '.db'
        self.conn = sqlite3.connect(db)

    def query(self, query):
        with self.conn.cursor() as cursor:
            cursor.execute(query)
            res = cursor.fetchall()
        return res

    def drop(self, query):
        with self.conn.cursor() as cursor:
            cursor.execute(query)
            try:
                cursor.execute(query)
            except Exception:
                pass

    def create(self, query):
        with self.conn.cursor() as cursor:
            cursor.execute(query)
            self.conn.commit()

    def insert_many(self, data):
        with self.conn.cursor() as cursor:
            cursor.executemany('INSERT INTO SVDB VALUES (?,?,?,?,?,?,?,?,?,?,?)', data)
            self.conn.commit()

    def create_index(self, name, columns):
        with self.conn.cursor() as cursor:
            cursor.execute("CREATE INDEX {} ON SVDB {}".format(name, columns))
            self.conn.commit()

    @property
    def tables(self):
        res = self.query("SELECT name FROM sqlite_master WHERE type=\'table\'")
        return [table[0] for table in res]

    @property
    def sample_ids(self):
        return [sample for sample in self.query('SELECT DISTINCT sample FROM SVDB')]
