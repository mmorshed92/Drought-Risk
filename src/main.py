# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 17:20:26 2025

@author: Ali.MorshediShahreba
"""

import os
import json
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

# ── Import ArcGIS helpers from your local module ──────────────────────────────
from Esri_Processing import esri_envelope, query_to_gdf
from drought_indices import list_rasters, build_spei_manifest

# -------------------- Layer URLs --------------------
# 1) CWS boundaries (EPA)
CWS_URL = "https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/Water_System_Boundaries/FeatureServer/0"

# 2) Sole Source Aquifers (EPA)
SSA_URL = "https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/Sole_Source_Aquifers_August_2019/FeatureServer/0"

# 3) EPA Surface Water Features (Hydrologic Units sublayers)
SUBWATERSHEDS_SERVICE = "https://watersgeo.epa.gov/arcgis/rest/services/OWOTHER/WATERS_GeoViewer_SurfaceWater/MapServer"
HYDRO_LAYERS = {
    "HUC10_Watersheds": 12,
    "HUC12_Subwatersheds": 13,
}
HUC10_URL = f"{SUBWATERSHEDS_SERVICE}/{HYDRO_LAYERS['HUC10_Watersheds']}"
HUC12_URL = f"{SUBWATERSHEDS_SERVICE}/{HYDRO_LAYERS['HUC12_Subwatersheds']}"

# 4) USGS WBD HUC12 (GeoJSON endpoint can be flaky; we handle gracefully)
WBD12_URL = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/6"

# 5) Underground Storage Tanks (example FeatureServer root; pick a sublayer id)
UST_URL = "https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/US_Underground_Storage_Tank_2019/FeatureServer/0"

# -------------------- Optional AOI filter --------------------
# set to None for nationwide pulls
BBOX = None
# Example (Chicagoland): BBOX = (-88.7, 41.3, -87.3, 42.5)

# -------------------- SPEI config (Step 2) --------------------
TIMESCALES  = ("06", "36")
SCENARIOS   = ("rcp45", "rcp85")
OUTDIR      = r"C:\Users\Ali.MorshediShahreba\OneDrive - Arup\Projects\IIa\Drought IiA\Data\drought_outputs"  # where to save the manifest


# === Download polygons ===
print("Downloading CWS…")
cws   = query_to_gdf(CWS_URL, bbox=BBOX)

print("Downloading SSA…")
ssa   = query_to_gdf(SSA_URL, bbox=BBOX)

print("Downloading EPA HUC12…")
huc12 = query_to_gdf(HUC12_URL, bbox=BBOX)

# print("Downloading EPA HUC10…")
# huc10 = query_to_gdf(HUC10_URL, bbox=BBOX)

# USGS WBD HUC12 is flaky for f=geojson; try once, then skip on error
wbd12 = None
try:
    print("Downloading USGS WBD HUC12…")
    wbd12 = query_to_gdf(WBD12_URL, bbox=BBOX)  # may raise HTTPError 500
except Exception as e:
    print("USGS WBD HUC12 failed; skipping. Reason:", e)
    print("Tip: rely on EPA HUC10/HUC12 above, or download WBD as a FileGDB and read offline.")

# UST (only if a concrete sublayer is provided)
ust = None
try:
    if UST_URL and UST_URL.rstrip("/").endswith(("/0","/1","/2","/3","/4","/5")):
        print("Downloading UST…")
        ust = query_to_gdf(UST_URL, bbox=BBOX)
except Exception as e:
    print("UST request failed (continuing):", e)

print(f"CWS:   {len(cws)} rows")
print(f"SSA:   {len(ssa)} rows")
print(f"HUC12: {len(huc12)} rows (EPA)")
# print(f"HUC10: {len(huc10)} rows (EPA)")
print(f"WBD12: {0 if wbd12 is None else len(wbd12)} rows (USGS)")
if ust is not None: print(f"UST:   {len(ust)} rows")



# === Quick plot of HUC12 (optional) ===
if not huc12.empty:
    ax = huc12.plot(figsize=(10, 10), edgecolor="black", linewidth=0.2, alpha=0.6)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()
else:
    print("HUC12 returned no features to plot.")
    
    

# === SPEI Step 2: list rasters + build manifest ===
SPEI_ROOT   = r"C:\Users\Ali.MorshediShahreba\OneDrive - Arup\Projects\IIa\Drought IiA\Data\SPEI - Climate Projections\RDS-2022-0075_Data_spei_06_cnrm_cm5"

print("\nSTEP 2 — Listing SPEI rasters…")
files = list_rasters(SPEI_ROOT)
print(f"Found {len(files)} raster files under SPEI_ROOT")

print("Building SPEI manifest (timescales=06/36, scenarios=RCP4.5/8.5)…")
manifest = build_spei_manifest(SPEI_ROOT, timescales=TIMESCALES, scenarios=SCENARIOS)

if manifest.empty:
    print("No SPEI files matched the filters. Check SPEI_ROOT or loosen filters.")
else:
    # Preview a few rows and counts
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


from drought_indices import check_monthly_completeness
gaps = check_monthly_completeness(manifest)
print(gaps[["timescale","scenario","model","n_months","start","end","missing_count"]].to_string(index=False))


from drought_indices import zonal_means_for_stack
from drought_indices import filter_manifest

huc12_small = huc12.sample(min(10, len(huc12)), random_state=0) if len(huc12) > 10 else huc12
mf_small = filter_manifest(manifest, timescales=("6",), scenarios=("rcp45",), start="2010-01", end="2011-12").head(3)
ts_small = zonal_means_for_stack(huc12_small, mf_small, id_field=None, all_touched=False, max_rasters=None)


    # Save manifest
os.makedirs(OUTDIR, exist_ok=True)
mpath = os.path.join(OUTDIR, "spei_manifest.csv")
manifest.to_csv(mpath, index=False)
print(f"\nSaved manifest → {mpath}")

