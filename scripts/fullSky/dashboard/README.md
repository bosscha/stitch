# Gaia DR3 Clustering - All-Sky Observability Dashboard

A lightweight, self-contained, real-time observability dashboard for tracking processed HEALPix pixels and newly detected stellar clusters in your PostgreSQL database (`gaiadb`).

---

## Features
- **Real-Time KPI Tracking**: Total processed pixels (out of 12,288 Level 5 pixels), detected clusters, tracked stars, and average cluster sizes.
- **Interactive Sky Plot**: An RA vs. Dec scatter plot of detected open cluster centroids.
- **Scale Distribution**: A distribution histogram showing cluster populations.
- **Interactive Catalog**: A searchable and sortable table of all detected cluster candidates including their ages, metallicity (`feh`), distance, and quality contrast metrics (`qc`).
- **Resilient Fallback**: Automatically displays fallback/placeholder metrics if PostgreSQL tables have not been created yet, so the UI is always accessible.

---

## Setup & Running the Dashboard

### 1. Launch the Server
To start the dashboard backend API and host the frontend, run the following command from the `dashboard/` directory:

```bash
../../.venv/bin/uvicorn app:app --host 0.0.0.0 --port 8050 --reload
```

Once started, open your web browser and navigate to:
**[http://localhost:8050](http://localhost:8050)**

### 2. Configure Database Parameters (Optional)
By default, the dashboard connects using the credentials matching your `stitch` configuration files:
- **Host**: `127.0.0.1`
- **Database**: `gaiadb`
- **User**: `postgres`
- **Password**: `tallis`

To override these parameters, simply define the environment variables when starting the server:

```bash
DB_HOST="localhost" DB_NAME="gaiadb" DB_USER="postgres" DB_PASS="tallis" ../../.venv/bin/uvicorn app:app --host 0.0.0.0 --port 8050 --reload
```
