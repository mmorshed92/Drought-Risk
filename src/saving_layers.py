# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 11:18:59 2025

@author: Ali.MorshediShahreba
"""

# saving_layers.py
from pathlib import Path
import fiona
import geopandas as gpd

# Globals configured via init_paths(...)
CACHE_DIR = None
GPKG_PATH = None

def init_paths(base_dir):
    """Initialize cache directories (GeoPackage + Parquet)."""
    global CACHE_DIR, GPKG_PATH, PARQ_DIR
    CACHE_DIR = Path(base_dir)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    GPKG_PATH = CACHE_DIR / "water_layers.gpkg"

def list_layers_in_gpkg(gpkg_path=None):
    gpkg_path = Path(gpkg_path) if gpkg_path else GPKG_PATH
    if not gpkg_path or not gpkg_path.exists():
        return set()
    try:
        return set(fiona.listlayers(str(gpkg_path)))
    except Exception:
        return set()

def save_layer(gdf: gpd.GeoDataFrame, layer: str, gpkg_path=None, parq_dir=None):
    if gdf is None or gdf.empty:
        return
    gpkg_path = Path(gpkg_path) if gpkg_path else GPKG_PATH
    # GeoPackage (overwrites layer)
    gdf.to_file(gpkg_path, layer=layer, driver="GPKG")

def load_cached_layer(layer: str, gpkg_path=None):
    gpkg_path = Path(gpkg_path) if gpkg_path else GPKG_PATH
    if gpkg_path and gpkg_path.exists() and layer in list_layers_in_gpkg(gpkg_path):
        return gpd.read_file(gpkg_path, layer=layer)
    return None

def get_or_download(url: str, layer_name: str, bbox=None, gpkg_path=None):
    """Load layer from cache if present; otherwise query ArcGIS and cache it."""
    from Esri_Processing import query_to_gdf  # local import to avoid hard dependency at import time
    gdf = load_cached_layer(layer_name, gpkg_path=gpkg_path)
    if gdf is not None and not gdf.empty:
        print(f"Loaded cached {layer_name}: {len(gdf)} rows")
        return gdf
    print(f"Downloading {layer_name}…")
    gdf = query_to_gdf(url, bbox=bbox)
    print(f"Downloaded {layer_name}: {len(gdf)} rows")
    save_layer(gdf, layer_name, gpkg_path=gpkg_path)
    return gdf
