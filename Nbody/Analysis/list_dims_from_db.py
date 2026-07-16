#!/usr/bin/env python3
import argparse
import psycopg2

def main():
    parser = argparse.ArgumentParser(description="List available dim_space simulations in the hypercluster database.")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    args = parser.parse_args()

    print(f"Connecting to database '{args.dbname}' on '{args.host}'...")
    try:
        conn = psycopg2.connect(
            dbname=args.dbname,
            user=args.user,
            password=args.password,
            host=args.host
        )
        cursor = conn.cursor()
    except Exception as e:
        print(f"Failed to connect to database: {e}")
        return

    print("Connected successfully. Retrieving available spatial dimensions...")
    
    try:
        # Get distinct dimensions
        cursor.execute("""
            SELECT DISTINCT dim_space 
            FROM star_snapshots 
            ORDER BY dim_space ASC;
        """)
        dims = [row[0] for row in cursor.fetchall()]
        
        if not dims:
            print("No simulation data found in the 'star_snapshots' table.")
            return

        print(f"\nFound {len(dims)} spatial dimension(s) in the database:")
        print(f"{'Dimension (dim_space)':<25} | {'Max Snapshot ID':<16} | {'Number of Stars':<16}")
        print("-" * 65)
        
        for dim in dims:
            # Get max snapshot_id
            cursor.execute("""
                SELECT MAX(snapshot_id) 
                FROM star_snapshots 
                WHERE dim_space = %s;
            """, (dim,))
            max_snapshot = cursor.fetchone()[0]
            
            # Get number of stars in the final snapshot
            if max_snapshot is not None:
                cursor.execute("""
                    SELECT COUNT(*) 
                    FROM star_snapshots 
                    WHERE dim_space = %s AND snapshot_id = %s;
                """, (dim, max_snapshot))
                    
                num_stars = cursor.fetchone()[0]
            else:
                num_stars = 0
                
            print(f"{dim:<25} | {max_snapshot if max_snapshot is not None else 'N/A':<16} | {num_stars:<16}")
            
        print()

    except Exception as e:
        print(f"Error querying database: {e}")
    finally:
        cursor.close()
        conn.close()

if __name__ == "__main__":
    main()
