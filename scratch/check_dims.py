import psycopg2

def main():
    conn = psycopg2.connect(dbname='hypercluster', user='stephane', password='tallis', host='localhost')
    cursor = conn.cursor()
    
    print("Checking indexes on star_snapshots...")
    cursor.execute("""
        SELECT indexname, indexdef 
        FROM pg_indexes 
        WHERE tablename = 'star_snapshots';
    """)
    for row in cursor.fetchall():
        print(f"Index: {row[0]} -> {row[1]}")
        
    print("\nEstimating row count...")
    cursor.execute("SELECT reltuples AS estimate FROM pg_class WHERE relname = 'star_snapshots';")
    print(f"Estimated row count: {cursor.fetchone()[0]}")
    
    candidate_dims = [1, 2, 3, 4, 5, 6, 25, 50, 100]
    existing_dims = []
    
    print("\nChecking dimensions...")
    for dim in candidate_dims:
        cursor.execute("SELECT 1 FROM star_snapshots WHERE dim_space = %s LIMIT 1;", (dim,))
        row = cursor.fetchone()
        if row:
            existing_dims.append(dim)
            print(f"Dimension {dim}D exists!")
            
            # Also get the max snapshot_id and count of stars for this dim
            cursor.execute("SELECT MAX(snapshot_id) FROM star_snapshots WHERE dim_space = %s;", (dim,))
            max_snap = cursor.fetchone()[0]
            
            cursor.execute("SELECT COUNT(*) FROM star_snapshots WHERE dim_space = %s AND snapshot_id = %s;", (dim, max_snap))
            count = cursor.fetchone()[0]
            print(f"  Max snapshot: {max_snap}, Stars count: {count}")
            
    print(f"All existing dimensions: {existing_dims}")
    cursor.close()
    conn.close()

if __name__ == '__main__':
    main()
