#!/usr/bin/env python3
import argparse
import psycopg2

def main():
    parser = argparse.ArgumentParser(description="List the initial total mass and mass statistics of the stellar cluster for each dimension N.")
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

    try:
        # Get distinct dimensions
        cursor.execute("SELECT DISTINCT dim_space FROM star_snapshots ORDER BY dim_space ASC;")
        dims = [row[0] for row in cursor.fetchall()]

        if not dims:
            print("No simulation data found in table 'star_snapshots'.")
            return

        print("\n" + "=" * 95)
        print(f"{'Dimension (N)':<15} | {'Min Snap ID':<12} | {'Star Count':<12} | {'Total Mass (M_sun)':<20} | {'Mean Mass (M_sun)':<18}")
        print("=" * 95)

        for dim in dims:
            # Find minimum snapshot_id for initial state
            cursor.execute("SELECT MIN(snapshot_id) FROM star_snapshots WHERE dim_space = %s;", (dim,))
            min_snap = cursor.fetchone()[0]

            if min_snap is None:
                continue

            # Query mass statistics for initial snapshot (using DISTINCT ON star_id for single run)
            cursor.execute("""
                SELECT 
                    COUNT(mass), 
                    SUM(mass), 
                    AVG(mass), 
                    MIN(mass), 
                    MAX(mass) 
                FROM (
                    SELECT DISTINCT ON (star_id) mass
                    FROM star_snapshots
                    WHERE dim_space = %s AND snapshot_id = %s
                    ORDER BY star_id, id DESC
                ) sub;
            """, (dim, min_snap))

            count, total_mass, avg_mass, min_m, max_m = cursor.fetchone()

            if count > 0:
                print(f"{dim:<15} | {min_snap:<12} | {count:<12} | {total_mass:<20.2f} | {avg_mass:<18.4f}")

        print("=" * 95 + "\n")

    except Exception as e:
        print(f"Error querying database: {e}")
    finally:
        cursor.close()
        conn.close()

if __name__ == "__main__":
    main()
