import os
import json
import logging
from starlette.applications import Starlette
from starlette.responses import JSONResponse, HTMLResponse, FileResponse
from starlette.routing import Route, Mount
from starlette.staticfiles import StaticFiles
from starlette.middleware import Middleware
from starlette.middleware.cors import CORSMiddleware
import psycopg2
from psycopg2.extras import RealDictCursor

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("gaia_dashboard")

# Database Connection details with defaults matching the stitch pipeline
DB_HOST = os.getenv("DB_HOST", "127.0.0.1")
DB_USER = os.getenv("DB_USER", "stephane")
DB_PASS = os.getenv("DB_PASS", "tallis")
DB_NAME = os.getenv("DB_NAME", "gaiadb")

def get_db_connection():
    return psycopg2.connect(
        host=DB_HOST,
        user=DB_USER,
        password=DB_PASS,
        dbname=DB_NAME,
        cursor_factory=RealDictCursor
    )

async def api_stats(request):
    """Fetches high-level stats of processed pixels, cluster counts, and tracked stars."""
    try:
        conn = get_db_connection()
        with conn.cursor() as cur:
            # Check if tables exist first
            cur.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_schema = 'public' 
                    AND table_name = 'clusters_processed_pixels'
                );
            """)
            has_pixels = cur.fetchone()["exists"]

            cur.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_schema = 'public' 
                    AND table_name = 'clusters_metadata'
                );
            """)
            has_metadata = cur.fetchone()["exists"]

            cur.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_schema = 'public' 
                    AND table_name = 'clusters'
                );
            """)
            has_clusters = cur.fetchone()["exists"]

            processed_pixels = 0
            total_clusters = 0
            total_stars = 0
            avg_cluster_size = 0.0

            if has_pixels:
                cur.execute("SELECT COUNT(*) as count FROM clusters_processed_pixels;")
                processed_pixels = cur.fetchone()["count"]

            if has_metadata:
                cur.execute("SELECT COUNT(*) as count, COALESCE(AVG(nstars), 0.0) as avg_size FROM clusters_metadata;")
                row = cur.fetchone()
                total_clusters = row["count"]
                avg_cluster_size = float(row["avg_size"])

            if has_clusters:
                cur.execute("SELECT COUNT(*) as count FROM clusters;")
                total_stars = cur.fetchone()["count"]

        conn.close()
        return JSONResponse({
            "status": "success",
            "db_connected": True,
            "data": {
                "processed_pixels": processed_pixels,
                "total_pixels_grid": 12288, # Level 5 Healpix total
                "total_clusters": total_clusters,
                "total_stars": total_stars,
                "avg_cluster_size": round(avg_cluster_size, 1)
            }
        })
    except Exception as e:
        logger.error(f"Error fetching stats: {e}")
        # Return fallback mock/placeholder data so UI doesn't break
        return JSONResponse({
            "status": "partial_success",
            "db_connected": False,
            "error": str(e),
            "data": {
                "processed_pixels": 1,
                "total_pixels_grid": 12288,
                "total_clusters": 11,
                "total_stars": 16428,
                "avg_cluster_size": 1493.5
            }
        })

async def api_clusters(request):
    """Returns a list of all processed clusters from clusters_metadata."""
    try:
        conn = get_db_connection()
        with conn.cursor() as cur:
            cur.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_schema = 'public' 
                    AND table_name = 'clusters_metadata'
                );
            """)
            if not cur.fetchone()["exists"]:
                return JSONResponse({"status": "success", "data": []})

            cur.execute("""
                SELECT cluster_id, votname, datetime, nstars, ntail, ra, dec, distance, age, feh, qc,
                       l, b, vl, vb, vldisp, vbdisp, xdisp, ydisp, zdisp, entropy_core
                FROM clusters_metadata
                ORDER BY datetime DESC;
            """)
            rows = cur.fetchall()
            # Convert float fields to float
            for row in rows:
                for key in ["ra", "dec", "distance", "age", "feh", "qc", "l", "b", "vl", "vb", "vldisp", "vbdisp", "xdisp", "ydisp", "zdisp", "entropy_core"]:
                    if row.get(key) is not None:
                        row[key] = float(row[key])
        conn.close()
        return JSONResponse({"status": "success", "data": rows})
    except Exception as e:
        logger.error(f"Error fetching clusters: {e}")
        return JSONResponse({"status": "error", "error": str(e), "data": []})

async def api_pixels(request):
    """Returns details of all processed pixel IDs."""
    try:
        conn = get_db_connection()
        with conn.cursor() as cur:
            cur.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_schema = 'public' 
                    AND table_name = 'clusters_processed_pixels'
                );
            """)
            if not cur.fetchone()["exists"]:
                return JSONResponse({"status": "success", "data": []})

            cur.execute("SELECT pix, datetime FROM clusters_processed_pixels ORDER BY pix ASC;")
            rows = cur.fetchall()
        conn.close()
        return JSONResponse({"status": "success", "data": rows})
    except Exception as e:
        logger.error(f"Error fetching pixels: {e}")
        return JSONResponse({"status": "error", "error": str(e), "data": []})

async def serve_home(request):
    """Serves the main HTML dashboard template."""
    return FileResponse(os.path.join(os.path.dirname(__file__), "index.html"))

# App Routing and Middleware setup
routes = [
    Route("/", serve_home),
    Route("/api/stats", api_stats),
    Route("/api/clusters", api_clusters),
    Route("/api/pixels", api_pixels)
]

middleware = [
    Middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])
]

app = Starlette(debug=True, routes=routes, middleware=middleware)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app:app", host="0.0.0.0", port=8050, reload=True)
