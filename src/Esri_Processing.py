# esri_utils.py
# Lightweight helpers for querying ArcGIS FeatureServer/MapServer layers into GeoPandas.

from __future__ import annotations

from typing import Optional, Tuple, Union
import time
import requests
import pandas as pd
import geopandas as gpd

# Accept either a numeric bbox tuple or any GeoPandas object with .total_bounds
BBoxLike = Union[Tuple[float, float, float, float], gpd.GeoSeries, gpd.GeoDataFrame]

HEADERS = {
    "Accept-Encoding": "gzip",
    "User-Agent": "esri-utils/1.0 (+https://arup.com)"
}


def esri_envelope(bbox4326: BBoxLike) -> str:
    """
    Convert (xmin, ymin, xmax, ymax) in EPSG:4326 (or any GeoPandas object with .total_bounds)
    to an Esri JSON envelope string, suitable for ArcGIS 'geometry' param with geometryType=envelope.
    """
    if hasattr(bbox4326, "total_bounds"):
        xmin, ymin, xmax, ymax = bbox4326.total_bounds  # type: ignore[attr-defined]
    else:
        xmin, ymin, xmax, ymax = bbox4326  # type: ignore[misc]

    if xmin >= xmax or ymin >= ymax:
        raise ValueError("Invalid bbox: expected xmin < xmax and ymin < ymax.")

    return (
        f'{{"xmin": {xmin}, "ymin": {ymin}, "xmax": {xmax}, '
        f'"ymax": {ymax}, "spatialReference": {{"wkid": 4326}}}}'
    )


def query_to_gdf(
    layer_url: str,
    where: str = "1=1",
    bbox: Optional[BBoxLike] = None,
    fields: str = "*",
    page_size: int = 2000,
    timeout: int = 120,
    simplify_tol_deg: Optional[float] = None,
    geometry_precision: Optional[int] = None,
    session: Optional[requests.Session] = None,
    max_retries: int = 3,
    backoff: float = 1.5,
) -> gpd.GeoDataFrame:
    """
    Query an ArcGIS FeatureServer/MapServer layer and return a GeoDataFrame (EPSG:4326).

    Parameters
    ----------
    layer_url : str
        .../FeatureServer/<id> or .../MapServer/<id>
    where : str
        SQL filter, e.g., "states LIKE '%IL%'"
    bbox : tuple or GeoPandas object with .total_bounds
        (xmin, ymin, xmax, ymax) in EPSG:4326
    fields : str
        Comma-separated field whitelist or "*"
    page_size : int
        Records per page (respects layer maxRecordCount if smaller)
    timeout : int
        Request timeout (seconds)
    simplify_tol_deg : float | None
        maxAllowableOffset in degrees (geometry generalization)
    geometry_precision : int | None
        geometryPrecision (digits to keep)
    session : requests.Session | None
        Provide your own keep-alive Session
    max_retries : int
        Retries per page on transient failures
    backoff : float
        Exponential backoff base between retries

    Returns
    -------
    GeoDataFrame in EPSG:4326 (empty if no features).
    """
    s = session or requests.Session()

    # Best-effort layer info to respect server paging limits
    supports_pagination = True
    try:
        info = s.get(layer_url.rstrip("/") + "?f=json", timeout=timeout, headers=HEADERS).json()
        max_rc = int(info.get("maxRecordCount", page_size))
        page_size = min(page_size, max_rc if max_rc > 0 else page_size)
        supports_pagination = bool(info.get("supportsPagination", True))
    except Exception:
        pass  # continue optimistically

    base = layer_url.rstrip("/") + "/query"

    common = {
        "where": where,
        "outFields": fields,
        "outSR": 4326,
        "returnGeometry": "true",
        "f": "geojson",
    }
    if bbox is not None:
        common.update(
            {
                "geometry": esri_envelope(bbox),
                "geometryType": "esriGeometryEnvelope",
                "spatialRel": "esriSpatialRelIntersects",
                "inSR": 4326,
            }
        )
    if simplify_tol_deg is not None:
        common["maxAllowableOffset"] = simplify_tol_deg
    if geometry_precision is not None:
        common["geometryPrecision"] = geometry_precision

    # Optional fast count
    total = None
    try:
        rr = s.get(
            base,
            params={**common, "f": "json", "returnCountOnly": "true"},
            timeout=timeout,
            headers=HEADERS,
        )
        rr.raise_for_status()
        total = int(rr.json().get("count", 0))
    except Exception:
        pass

    gdfs = []
    offset = 0
    while True:
        params = {**common, "resultOffset": offset, "resultRecordCount": page_size}

        # Robust fetch with retries
        last_err = None
        for attempt in range(max_retries):
            try:
                r = s.get(base, params=params, timeout=timeout, headers=HEADERS)
                r.raise_for_status()
                fc = r.json()
                break
            except Exception as e:
                last_err = e
                if attempt == max_retries - 1:
                    raise
                time.sleep(backoff ** attempt)

        feats = fc.get("features") or []
        if not feats:
            # If server doesn't support pagination and didn't exceed transfer limit, stop.
            if not supports_pagination and not fc.get("exceededTransferLimit", False):
                break
            break

        gdf = gpd.GeoDataFrame.from_features(fc, crs="EPSG:4326")
        if not gdf.empty:
            gdfs.append(gdf)

        offset += page_size

        # Stop conditions
        if total is not None and offset >= total:
            break
        if not fc.get("exceededTransferLimit", False) and not supports_pagination:
            break

    if not gdfs:
        return gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")

    out = pd.concat(gdfs, ignore_index=True)
    return gpd.GeoDataFrame(out, crs="EPSG:4326")
