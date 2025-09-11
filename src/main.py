# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 17:20:26 2025
@author: Ali.MorshediShahreba
"""

import os, json, warnings, re
from pathlib import Path
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

# ── Helpers ───────────────────────────────────────────────────────────────────
from Esri_Processing import esri_envelope, query_to_gdf
from drought_indices import (list_rasters, build_spei_manifest,
    check_monthly_completeness, zonal_means_for_stack, filter_manifest)
from saving_layers import init_paths, get_or_download, load_cached_layer, save_layer, list_layers_in_gpkg

# -------------------- Layer URLs --------------------
CWS_URL  = "https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/Water_System_Boundaries/FeatureServer/0"
SSA_URL  = "https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/Sole_Source_Aquifers_August_2019/FeatureServer/0"

SUBWATERSHEDS_SERVICE = "https://watersgeo.epa.gov/arcgis/rest/services/OWOTHER/WATERS_GeoViewer_SurfaceWater/MapServer"
HUC12_URL = f"{SUBWATERSHEDS_SERVICE}/13"  # HUC12_Subwatersheds

WBD12_URL = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/6"
UST_URL   = "https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/US_Underground_Storage_Tank_2019/FeatureServer/0"

#  a percent-by-HUC12 layer for protected areas.
SWPA_URL = "https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/Percent_HUC12_SWPA/FeatureServer"
SWPA_PCT_LAYER = f"{SWPA_URL.rstrip('/')}/0"  # adjust if your percent layer index differs

# -------------------- Optional AOI filter --------------------
BBOX = None
# Example: BBOX = (-88.7, 41.3, -87.3, 42.5)

# -------------------- cache --------------------
CACHE_DIR = Path(r"C:\Users\Ali.MorshediShahreba\OneDrive - Arup\Projects\IIa\Drought IiA\Data\ModelLayers")
init_paths(CACHE_DIR)  # sets water_layers.gpkg + /parquet

# -------------------- SPEI config --------------------
TIMESCALES  = ("06", "36")
SCENARIOS   = ("rcp45", "rcp85")
OUTDIR      = r"C:\Users\Ali.MorshediShahreba\OneDrive - Arup\Projects\IIa\Drought IiA\Data\drought_outputs"
SPEI_ROOT   = r"C:\Users\Ali.MorshediShahreba\OneDrive - Arup\Projects\IIa\Drought IiA\Data\SPEI - Climate Projections\RDS-2022-0075_Data_spei_06_cnrm_cm5"

# ==================== Small utilities kept local ====================
def detect_huc12_col(gdf: gpd.GeoDataFrame) -> str:
    candidates = ["HUC12", "HUC_12", "huc12", "HUC12CODE", "HUC12_Code"]
    for c in candidates:
        if c in gdf.columns:
            return c
    for c in gdf.columns:
        try:
            if gdf[c].astype(str).str.len().median() == 12:
                warnings.warn(f"Using '{c}' as HUC12 (heuristic).")
                return c
        except Exception:
            pass
    raise ValueError("Could not find a HUC12 column. Please add/rename one.")

def pick_percent_columns(df: pd.DataFrame):
    """Return columns likely to be percent fields."""
    cols = []
    for c in df.columns:
        cl = c.lower()
        if ("pct" in cl or "percent" in cl or cl.endswith("_pc")) and pd.api.types.is_numeric_dtype(df[c]):
            cols.append(c)
    # Fallback: include any numeric column explicitly named like '*swpa*'
    if not cols:
        for c in df.columns:
            cl = c.lower()
            if "swpa" in cl and pd.api.types.is_numeric_dtype(df[c]):
                cols.append(c)
    return cols

# ==================== Download / Load (cached) ====================
cws   = get_or_download(CWS_URL,  layer_name="cws_service_areas", bbox=BBOX)
ssa   = get_or_download(SSA_URL,  layer_name="ssa",                bbox=BBOX)
huc12 = get_or_download(HUC12_URL, layer_name="huc12_epa",         bbox=BBOX)

# USGS WBD HUC12 (optional)
wbd12 = load_cached_layer("wbd12_usgs")
if wbd12 is None:
    try:
        print("Attempting USGS WBD HUC12…")
        wbd12 = query_to_gdf(WBD12_URL, bbox=BBOX)
        print(f"Downloaded WBD12: {len(wbd12)} rows")
        save_layer(wbd12, "wbd12_usgs")
    except Exception as e:
        print("USGS WBD HUC12 failed; skipping. Reason:", e)
        wbd12 = None

# UST (optional)
ust = load_cached_layer("ust")
if ust is None and UST_URL and UST_URL.rstrip("/").endswith(("/0","/1","/2","/3","/4","/5")):
    try:
        print("Downloading UST…")
        ust = query_to_gdf(UST_URL, bbox=BBOX)
        save_layer(ust, "ust")
    except Exception as e:
        print("UST request failed (continuing):", e)
        ust = None

# ----- Your percent-by-HUC12 protected areas layer -----
swpa_pct = None
if SWPA_URL:
    try:
        swpa_pct = get_or_download(SWPA_PCT_LAYER, layer_name="huc12_swpa_percent", bbox=None)
    except Exception as e:
        print("Percent-HUC12 SWPA request failed (continuing):", e)

print(f"CWS:        {len(cws)} rows")
print(f"SSA:        {len(ssa)} rows")
print(f"HUC12 EPA:  {len(huc12)} rows")
print(f"WBD12 USGS: {0 if wbd12 is None else len(wbd12)} rows")
if ust is not None:    print(f"UST:        {len(ust)} rows")
if swpa_pct is not None:
    print(f"SWPA %:     {len(swpa_pct)} rows (percent-by-HUC12)")


# === Quick plot of HUC12 ===
if not huc12.empty:
    ax = huc12.plot(figsize=(10, 10), edgecolor="black", linewidth=0.2, alpha=0.6)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()
else:
    print("HUC12 returned no features to plot.")

# ==================== Join the percent layer to HUC12 ====================
huc_col = detect_huc12_col(huc12)

# If your percent layer is a *table* (no geometry), this still works.
if swpa_pct is not None and not swpa_pct.empty:
    # Find HUC12 column on the percent layer (same heuristic)
    pct_huc_col = detect_huc12_col(swpa_pct)
    pct_cols = pick_percent_columns(swpa_pct)
    if not pct_cols:
        print("Warning: could not find obvious percent columns on SWPA percent layer. Will join only HUC12.")
        pct_cols = []
    # Keep only HUC + percent columns (avoid duplicating geometry if any)
    pct_df = swpa_pct.drop(columns=[c for c in swpa_pct.columns if c not in [pct_huc_col] + pct_cols])
    # Merge onto HUC12 polygons
    huc12_percents = huc12.merge(pct_df, left_on=huc_col, right_on=pct_huc_col, how="left")
else:
    # If percent layer is absent, just pass HUC12 through (no computed overlay)
    huc12_percents = huc12.copy()

# Cache the percent-enriched HUC12
save_layer(huc12_percents, "huc12_percents_joined")

# Optional CSV export (no geometry)
csv_out = CACHE_DIR / "huc12_percents_joined.csv"
huc12_percents.drop(columns="geometry").to_csv(csv_out, index=False)
print(f"Saved HUC12 percents (joined) → {csv_out}")

# ==================== SPEI: list rasters + build manifest ====================
print("\nSTEP 2 — Listing SPEI rasters…")
files = list_rasters(SPEI_ROOT)
print(f"Found {len(files)} raster files under SPEI_ROOT")

print("Building SPEI manifest (timescales=06/36, scenarios=RCP4.5/8.5)…")
manifest = build_spei_manifest(SPEI_ROOT, timescales=TIMESCALES, scenarios=SCENARIOS)

if manifest.empty:
    print("No SPEI files matched the filters. Check SPEI_ROOT or loosen filters.")
else:
    print("\nManifest preview (top 5):")
    try:
        print(manifest.head(5).to_string(index=False))
    except Exception:
        print(manifest.head(5))
    print("\nCounts by timescale/scenario:")
    try:
        counts = manifest.groupby(["timescale", "scenario"]).size().rename("count")
        print(counts.to_string())
    except Exception:
        print(manifest.groupby(["timescale", "scenario"]).size())

gaps = check_monthly_completeness(manifest)
print(gaps[["timescale","scenario","model","n_months","start","end","missing_count"]].to_string(index=False))

# Small test run (no upstream HUC12s anywhere)
huc12_small = huc12.sample(min(10, len(huc12)), random_state=0) if len(huc12) > 10 else huc12
mf_small = filter_manifest(manifest, timescales=("6",), scenarios=("rcp45",), start="2010-01", end="2011-12").head(3)
ts_small = zonal_means_for_stack(huc12_small, mf_small, id_field=None, all_touched=False, max_rasters=None)

# Save manifest
os.makedirs(OUTDIR, exist_ok=True)
mpath = os.path.join(OUTDIR, "spei_manifest.csv")
manifest.to_csv(mpath, index=False)
print(f"\nSaved manifest → {mpath}")
