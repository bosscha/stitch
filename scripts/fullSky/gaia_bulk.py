import os
import re
import requests
import psycopg2
import gzip
import csv
import io
import sys

# --- CONFIGURATION ---
BASE_URL = "https://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/"
MD5_URL = BASE_URL + "_MD5SUM.txt"

DB_DBNAME = os.getenv("DB_NAME", "gaiadb")
DB_USER = os.getenv("DB_USER", "stephane")
DB_PASS = os.getenv("DB_PASS", "tallis")
DB_HOST = os.getenv("DB_HOST", "127.0.0.1")
DB_CONFIG = f"dbname={DB_DBNAME} user={DB_USER} password={DB_PASS} host={DB_HOST}"

# Scientific Threshold
RUWE_THRESHOLD = 1.4

# Column mapping matching your PostgreSQL schema exactly
TARGET_COLUMNS = [
    'source_id', 'ra', 'dec', 'parallax', 'parallax_error',
    'pmra', 'pmdec', 'l', 'b', 'phot_g_mean_mag', 'phot_bp_mean_mag',
    'phot_rp_mean_mag', 'phot_g_mean_flux', 'phot_g_mean_flux_error',
    'phot_bp_mean_flux', 'phot_bp_mean_flux_error',
    'phot_rp_mean_flux', 'phot_rp_mean_flux_error',
    'ruwe', 'astrometric_excess_noise', 'astrometric_params_solved',
    'radial_velocity', 'radial_velocity_error',
    'ag_gspphot', 'azero_gspphot', 'ebpminrp_gspphot', 'mh_gspphot',
    'nu_eff_used_in_astrometry', 'pseudocolour', 'ecl_lat'
]

def get_files_from_md5():
    """Reads the static MD5 sum file list to bypass dynamic JS loading."""
    print(f"Fetching DR3 file list from: {MD5_URL}")
    headers = {'User-Agent': 'Mozilla/5.0'}
    try:
        r = requests.get(MD5_URL, headers=headers, timeout=30)
        r.raise_for_status()
        files = re.findall(r'GaiaSource_\d{6}-\d{6}\.csv\.gz', r.text)
        return sorted(list(set(files)))
    except Exception as e:
        print(f"Error reading MD5 file list: {e}")
        return []

def process_file(conn, filename, session):
    cur = conn.cursor()
    local_path = filename

    try:
        # 1. DOWNLOAD
        print(f"Downloading {filename}...")
        with session.get(BASE_URL + filename, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(local_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=1024*1024):
                    f.write(chunk)

        # 2. FILTER & TRANSFORM
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=TARGET_COLUMNS, extrasaction='ignore', delimiter='\t')

        in_count = 0
        out_count = 0

        with gzip.open(local_path, 'rt') as f:
            # Skip comments and empty lines
            clean_lines = (line for line in f if line.strip() and not line.startswith('#'))

            reader = csv.DictReader(clean_lines)

            # --- DYNAMIC HEADER MAPPING ---
            actual_keys = {k.strip(): k for k in reader.fieldnames if k}
            sid_key = actual_keys.get('source_id')
            ruwe_key = actual_keys.get('ruwe')

            if not sid_key:
                print(f"CRITICAL: 'source_id' column not found in {filename}!")
                return False

            for row in reader:
                in_count += 1

                # VALIDATION: Skip empty lines
                sid = row.get(sid_key, '').strip()
                if not sid:
                    continue

                # --- PERMISSIVE RUWE LOGIC ---
                try:
                    raw_ruwe = row.get(ruwe_key, '').strip() if ruwe_key else ''
                    # Catch literal 'null' strings just in case they appear in the RUWE column
                    if raw_ruwe.lower() in ['', 'null', 'nan']:
                        ruwe = 0.0
                    else:
                        ruwe = float(raw_ruwe)
                except:
                    ruwe = 99.0

                if ruwe <= RUWE_THRESHOLD:
                    # Clean strings and catch literal "null" text
                    clean_row = {}
                    for col_name in TARGET_COLUMNS:
                        csv_key = actual_keys.get(col_name)
                        val = row.get(csv_key, '') if csv_key else ''

                        if val is not None:
                            val = str(val).strip()
                            # CRITICAL FIX: Catch literal "null" and "NaN" text
                            if val.lower() in ['null', 'nan', '']:
                                val = None

                        clean_row[col_name] = val

                    writer.writerow(clean_row)
                    out_count += 1

        # 3. BULK COPY
        if out_count > 0:
            buffer.seek(0)
            columns_str = ','.join(TARGET_COLUMNS)
            cur.copy_expert(f"COPY gaia_source ({columns_str}) FROM STDIN WITH (DELIMITER '\t', NULL '')", buffer)
        else:
            print(f"WARNING: No stars passed filter for {filename}")

        # 4. LOG & COMMIT
        cur.execute("INSERT INTO import_log (file_name) VALUES (%s) ON CONFLICT DO NOTHING", (filename,))
        conn.commit()
        print(f"SUCCESS: {filename} (Kept {out_count}/{in_count})")
        return True

    except Exception as e:
        conn.rollback()
        print(f"FAILED {filename}: {e}")
        return False
    finally:
        cur.close()
        if os.path.exists(local_path):
            os.remove(local_path)

def main():
    file_list = get_files_from_md5()
    if not file_list:
        print("No files found. Check your internet connection.")
        return

    print(f"Found {len(file_list)} files. Starting Full DR3 ingestion...")

    try:
        conn = psycopg2.connect(DB_CONFIG)
    except Exception as e:
        print(f"Database connection failed: {e}")
        return

    session = requests.Session()
    session.headers.update({'User-Agent': 'Mozilla/5.0'})

    for filename in file_list:
        check_cur = conn.cursor()
        check_cur.execute("SELECT 1 FROM import_log WHERE file_name = %s", (filename,))
        done = check_cur.fetchone()
        check_cur.close()

        if done:
            continue

        if not process_file(conn, filename, session):
            print("Stopping due to error. Restart script to resume.")
            break

    conn.close()
    print("Process Finished.")

if __name__ == "__main__":
    main()
