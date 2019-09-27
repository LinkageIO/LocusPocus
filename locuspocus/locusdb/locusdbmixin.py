
class LocusDBMixin(object):

    __slots__ = ('_db')
    
    def _initialize_tables(self):
        cur = self._db.cursor()
        cur.execute('''
            PRAGMA foreign_keys = ON;
        ''')
        cur.execute('''
            CREATE TABLE IF NOT EXISTS loci (
                /* Store the UUID as text  */
                /* https://stackoverflow.com/questions/11337324/how-to-efficient-insert-and-fetch-uuid-in-core-data/11337522#11337522 */
                UUID TEXT PRIMARY KEY,
                
                /* Store the locus values  */
                chromosome TEXT NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER NOT NULL,

                source TEXT,
                feature_type TEXT,
                strand TEXT,
                frame INT,

                /* Store things to make my life easier */
                name TEXT
                
            );
        ''')
        cur.execute('''
            CREATE INDEX IF NOT EXISTS locus_id ON loci (name);
            CREATE INDEX IF NOT EXISTS locus_chromosome ON loci (chromosome);
            CREATE INDEX IF NOT EXISTS locus_start ON loci (start);
            CREATE INDEX IF NOT EXISTS locus_feature_type ON loci (feature_type);
            CREATE INDEX IF NOT EXISTS locus_end ON loci (end);
            '''
        )
        cur.execute('''
            -- Create a table that contains loci attribute mapping
            CREATE TABLE IF NOT EXISTS loci_attrs (
                UUID TEXT REFERENCES loci(UUID) ON DELETE CASCADE,
                key TEXT,
                val TEXT,
                type TEXT,
                UNIQUE(UUID,key)
            );
        ''')
        cur.execute('''
            -- Create a table with parent-child relationships        
            CREATE TABLE IF NOT EXISTS relationships (
                parent INT REFERENCES loci(UUID) ON DELETE CASCADE,
                child INT REFERENCES loci(UUID) ON DELETE CASCADE
            );
            CREATE INDEX IF NOT EXISTS relationships_parent ON relationships (parent);
            CREATE INDEX IF NOT EXISTS relationships_child ON relationships (child);
        '''
        )
        cur.execute('''
            -- Create a R*Tree table so we can efficiently query by ranges
            CREATE VIRTUAL TABLE IF NOT EXISTS positions USING rtree_i32( 
                UUID, 
                start INT,
                end INT,
                +chromosome TEXT
            );
        ''')
