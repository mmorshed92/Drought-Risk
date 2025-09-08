# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 17:20:26 2025

@author: Ali.MorshediShahreba
"""

import os, json, pandas as pd, geopandas as gpd, requests
from urllib.parse import urlencode

def esri_envelope(bbox4326):
    xmin, ymin, xmax, ymax = bbox4326
    return json.dumps({"xmin": xmin, "ymin": ymin, "xmax": xmax, "ymax": ymax,
                       "spatialReference": {"wkid": 4326}})

def query_to_gdf(layer_url, where="1=1", bbox=None, fields="*",
                      page_size=2000, timeout=120,
                      simplify_tol_deg=None, geometry_precision=None,
                      session=None):
    """
    Faster ArcGIS query with field whitelisting, geometry simplification, and keep-alive.
    - layer_url: .../FeatureServer/<id> or .../MapServer/<id>
    - where: SQL filter (e.g., "states LIKE '%IL%'")
    - bbox: (xmin, ymin, xmax, ymax) in EPSG:4326
    - fields: e.g., "huc12,name,states"
    - simplify_tol_deg: maxAllowableOffset in degrees (e.g., 0.0002)
    - geometry_precision: int digits to keep (e.g., 6)
    """
    s = session or requests.Session()
    base = layer_url.rstrip("/") + "/query"
    common = {
        "where": where,
        "outFields": fields,
        "outSR": 4326,
        "returnGeometry": "true",
        "f": "geojson"
    }
    if bbox:
        common.update({
            "geometry": esri_envelope(bbox),
            "geometryType": "esriGeometryEnvelope",
            "spatialRel": "esriSpatialRelIntersects",
            "inSR": 4326
        })
    if simplify_tol_deg is not None:
        common["maxAllowableOffset"] = simplify_tol_deg
    if geometry_precision is not None:
        common["geometryPrecision"] = geometry_precision

    # Try count (optional)
    total = None
    try:
        rr = s.get(base, params={**common, "f": "json", "returnCountOnly": "true"},
                   timeout=timeout, headers={"Accept-Encoding":"gzip"})
        rr.raise_for_status()
        total = int(rr.json().get("count", 0))
    except Exception:
        pass

    gdfs, offset = [], 0
    while True:
        params = {**common, "resultOffset": offset, "resultRecordCount": page_size}
        url = f"{base}?{urlencode(params)}"
        r = s.get(url, timeout=timeout, headers={"Accept-Encoding":"gzip"})
        r.raise_for_status()
        txt = r.text.strip()
        if not txt:
            break
        fc = r.json()
        feats = fc.get("features", [])
        if not feats:
            break
        gdf = gpd.GeoDataFrame.from_features(fc, crs="EPSG:4326")
        gdfs.append(gdf)
        offset += page_size
        if total is not None and offset >= total:
            break

    if not gdfs:
        return gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    return gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), crs="EPSG:4326")

# -------------------- Layer URLs --------------------

# 1) CWS boundaries (EPA)
CWS_URL = "https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/Water_System_Boundaries/FeatureServer/0"

# 2) Sole Source Aquifers (EPA)
SSA_URL = "https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/Sole_Source_Aquifers_August_2019/FeatureServer/0"

# 3) EPA Surface Water Features (Hydrologic Units sublayers)
SUBWATERSHEDS_SERVICE = "https://watersgeo.epa.gov/arcgis/rest/services/OWOTHER/WATERS_GeoViewer_SurfaceWater/MapServer"
HYDRO_LAYERS = {
    "HUC2_Regions": 8,
    "HUC4_Subregions": 9,
    "HUC6_Basins": 10,
    "HUC8_Subbasins": 11,
    "HUC10_Watersheds": 12,
    "HUC12_Subwatersheds": 13
}
HUC10_URL = f"{SUBWATERSHEDS_SERVICE}/{HYDRO_LAYERS['HUC10_Watersheds']}"
HUC12_URL = f"{SUBWATERSHEDS_SERVICE}/{HYDRO_LAYERS['HUC12_Subwatersheds']}"

# 4) USGS WBD HUC12 (this one often 500s on GeoJSON; we handle gracefully)
WBD12_URL = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/6"

# 5) UST / PSC placeholders (replace with real FeatureServer layer URLs when ready)
UST_URL = "<PASTE_UST_FEATURESERVER_LAYER_URL_HERE>"            # e.g., .../FeatureServer/0
PSC_URL = "<PASTE_PSC_FEATURESERVER_OR_POINTS_LAYER_URL_HERE>"  # e.g., .../FeatureServer/0

# -------------------- Optional AOI filter --------------------
# Keep downloads small while testing; set to None for nationwide pulls
BBOX = None
# Example (Chicagoland):
# BBOX = (-88.7, 41.3, -87.3, 42.5)

# -------------------- Download --------------------

print("Downloading CWS…")
cws   = query_to_gdf(CWS_URL, bbox=BBOX)
print("Downloading SSA…")
ssa   = query_to_gdf(SSA_URL, bbox=BBOX)
print("Downloading EPA HUC12…")
huc12 = query_to_gdf(HUC12_URL, bbox=BBOX)
print("Downloading EPA HUC10…")
huc10 = query_to_gdf(HUC10_URL, bbox=BBOX)

# USGS WBD HUC12 is flaky for f=geojson; try once, then skip on error
wbd12 = None
try:
    print("Downloading USGS WBD HUC12…")
    wbd12 = query_to_gdf(WBD12_URL, bbox=BBOX)   # may raise HTTPError 500
except Exception as e:
    print("USGS WBD HUC12 failed; skipping. Reason:", e)
    print("Tip: rely on EPA HUC10/HUC12 above, or download WBD as a FileGDB and read offline.")

# Only try UST / PSC if placeholders replaced
ust = query_to_gdf(UST_URL, bbox=BBOX) if "PASTE_UST" not in UST_URL else None
psc = query_to_gdf(PSC_URL, bbox=BBOX) if "PASTE_PSC" not in PSC_URL else None

print(f"CWS:   {len(cws)} rows")
print(f"SSA:   {len(ssa)} rows")
print(f"HUC12: {len(huc12)} rows (EPA)")
print(f"HUC10: {len(huc10)} rows (EPA)")
print(f"WBD12: {0 if wbd12 is None else len(wbd12)} rows (USGS)")
if ust is not None: print(f"UST:   {len(ust)} rows")
if psc is not None: print(f"PSC:   {len(psc)} rows")

# -------------------- Save --------------------

out_gpkg = r"C:\Users\Ali.MorshediShahreba\OneDrive - Arup\Projects\IIa\Drought IiA\Code\Drought-Risk\data\dwmap_extract.gpkg"
ensure_dir_for_file(out_gpkg)

cws.to_file(out_gpkg,  layer="EPA_CWS_Boundaries",       driver="GPKG")
ssa.to_file(out_gpkg,  layer="EPA_Sole_Source_Aquifers", driver="GPKG")
huc12.to_file(out_gpkg, layer="EPA_HUC12_Subwatersheds", driver="GPKG")
huc10.to_file(out_gpkg, layer="EPA_HUC10_Watersheds",    driver="GPKG")
if wbd12 is not None:
    wbd12.to_file(out_gpkg, layer="USGS_WBD_HUC12",      driver="GPKG")
if ust is not None:
    ust.to_file(out_gpkg,  layer="UST",                  driver="GPKG")
if psc is not None:
    psc.to_file(out_gpkg,  layer="PSC",                  driver="GPKG")

# Optional: fast columnar copy for tabular workflows
cws.to_parquet(out_gpkg.replace(".gpkg", "_cws.parquet"))
