# -*- coding: utf-8 -*-
"""
Tools to build a raster manifest (SPEI), compute HUC12 zonal means, and aggregate drought stats.
Step 1–2 only: filename parsing + recursive file listing + manifest builder.
"""

from __future__ import annotations
import os, re
from pathlib import Path
import pandas as pd


# ---------------- Step 1: filename parser ----------------

def _parse_spei_filename(path: str) -> dict:
    """
    Parse common SPEI filenames to extract:
      timescale ('06','36'), scenario ('rcp45','rcp85' or 'unknown'),
      model (best-effort), year (int), month (int).
    Examples:
      SPEI06_RCP45_GFDL_195001.tif
      spei_36_rcp8.5_CanESM2_2003-07.tif
    """
    name = os.path.basename(path).lower()

    # timescale (e.g., 06 or 36)
    m_ts = re.search(r'spei[^0-9]*0?([0-9]{1,2})', name)
    timescale = m_ts.group(1).zfill(2) if m_ts else None

    # scenario
    if re.search(r'rcp4\.?5|rcp45', name):
        scenario = 'rcp45'
    elif re.search(r'rcp8\.?5|rcp85', name):
        scenario = 'rcp85'
    else:
        scenario = 'unknown'

    # year-month (YYYYMM or YYYY-MM)
    m_ym = re.search(r'((?:19|20)\d{2})[^\d]?((?:0[1-9]|1[0-2]))', name)
    year  = int(m_ym.group(1)) if m_ym else None
    month = int(m_ym.group(2)) if m_ym else None

    # crude model guess = text between scenario and first year token
    model = 'unknown'
    try:
        parts = re.split(r'rcp4\.?5|rcp45|rcp8\.?5|rcp85', name)
        if len(parts) > 1:
            tail = parts[1]
            tail = re.split(r'(?:19|20)\d{2}', tail)[0]
            model = re.sub(r'[^a-z0-9]+', '', tail) or 'unknown'
    except Exception:
        pass

    return dict(timescale=timescale, scenario=scenario, model=model, year=year, month=month)

# ---------------- Step 2: list files + build manifest ----------------

def list_rasters(root: str, extensions=(".tif", ".tiff")) -> list[str]:
    """
    Recursively list raster files under `root`. Uses pathlib instead of glob.
    """
    exts = {e.lower() for e in extensions}
    return [
        str(p)
        for p in Path(root).rglob("*")
        if p.is_file() and p.suffix.lower() in exts
    ]


def build_spei_manifest(
    spei_root: str,
    timescales=("06","36"),        # set to None to accept any timescale
    scenarios=("rcp45","rcp85"),   # set to None to accept any scenario
) -> pd.DataFrame:
    """
    Scan `spei_root` for GeoTIFFs and parse filename metadata.
    Returns columns: path, date, timescale, scenario, model, year, month
    """
    tif_paths = list_rasters(spei_root)

    rows = []
    for p in tif_paths:
        meta = _parse_spei_filename(p)
        # must have a parsable year-month and timescale
        if not (meta["timescale"] and meta["year"] and meta["month"]):
            continue

        if timescales:
            keep_ts = {str(ts).zfill(2) for ts in timescales}
            if meta["timescale"] not in keep_ts:
                continue

        if scenarios:
            if meta["scenario"] not in set(scenarios):
                continue

        rows.append({
            **meta,
            "path": p,
            "date": pd.Timestamp(meta["year"], meta["month"], 1),
        })

    if not rows:
        return pd.DataFrame(columns=["path","date","timescale","scenario","model","year","month"])

    return pd.DataFrame(rows).sort_values("date").reset_index(drop=True)



def check_monthly_completeness(mf: pd.DataFrame,
                               by=("timescale","scenario","model"),
                               freq="MS") -> pd.DataFrame:
    """
    For each group, report #months, span, and missing months within the span.
    """
    rows = []
    if mf.empty:
        return pd.DataFrame(rows)

    for keys, g in mf.groupby(list(by), dropna=False):
        s = g.set_index("date").sort_index()
        full_idx = pd.date_range(s.index.min(), s.index.max(), freq=freq)
        present = pd.Series(True, index=s.index)
        mask = pd.Series(False, index=full_idx)
        mask.loc[present.index] = True
        missing = mask[~mask].index
        rows.append({
            **{k:v for k,v in zip(by, keys)},
            "n_months": len(s),
            "start": s.index.min(),
            "end": s.index.max(),
            "missing_count": int(len(missing)),
            # cap display to first 24 to keep tables readable
            "missing_months": [d.strftime("%Y-%m") for d in missing[:24]],
        })
    return pd.DataFrame(rows).sort_values(list(by))



def filter_manifest(mf: pd.DataFrame,
                    timescales=None,        # e.g., ("06","36")
                    scenarios=None,         # e.g., ("rcp45","rcp85")
                    models=None,            # e.g., ("GFDL","CanESM2")
                    start=None,             # e.g., "1980-01"
                    end=None,               # e.g., "2060-12"
                    months=None             # subset months, e.g., (6,7,8)
                   ) -> pd.DataFrame:
    """
    Flexible filtering for downstream work.
    """
    df = mf.copy()
    if timescales:
        keep_ts = {str(ts).zfill(2) for ts in timescales}
        df = df[df["timescale"].isin(keep_ts)]
    if scenarios:
        df = df[df["scenario"].isin(scenarios)]
    if models:
        df = df[df["model"].isin(models)]
    if start:
        df = df[df["date"] >= pd.to_datetime(start)]
    if end:
        df = df[df["date"] <= pd.to_datetime(end)]
    if months:
        df = df[df["month"].isin(months)]
    return df.reset_index(drop=True)





# --- Step 3: Zonal means per HUC12 from a raster stack -----------------------
from typing import Optional
import numpy as np

def zonal_means_for_stack(
    huc12_gdf: gpd.GeoDataFrame,
    manifest: pd.DataFrame,
    id_field: Optional[str] = None,
    all_touched: bool = False,
    use_spatial_index: bool = True,
    max_rasters: Optional[int] = None,
) -> pd.DataFrame:
    """
    Compute monthly HUC12 mean SPEI for each raster in `manifest`.

    Returns a long DataFrame:
      [huc12_id, date, timescale, scenario, model, spei]

    Parameters
    ----------
    huc12_gdf : GeoDataFrame (EPSG:4326 or any CRS)
    manifest  : DataFrame with columns ['path','date','timescale','scenario','model',...]
    id_field  : column name to use as HUC12 id (auto-detects common names)
    all_touched : pass-through to rasterstats.zonal_stats (False = conservative)
    use_spatial_index : build a per-raster spatial index to skip non-overlapping polygons
    max_rasters : if set, process only the first N rasters (for quick tests)
    """
    try:
        import rasterio
        from rasterstats import zonal_stats
        from shapely.geometry import box
    except Exception as e:
        raise ImportError("Install rasterio and rasterstats first.") from e

    if huc12_gdf.empty:
        return pd.DataFrame(columns=["huc12_id","date","timescale","scenario","model","spei"])

    # Pick an id field if none specified
    if id_field is None:
        candidates = [c for c in huc12_gdf.columns if str(c).lower() in ("huc12", "huc_12", "huc")]
        id_field = candidates[0] if candidates else huc12_gdf.columns[0]

    records = []
    # Optional limit for quick runs
    if max_rasters is not None:
        manifest = manifest.head(int(max_rasters)).copy()

    for (ts, scen, mdl), grp in manifest.groupby(["timescale", "scenario", "model"], dropna=False):
        # Iterate in chronological order for nice logs
        grp = grp.sort_values("date")
        for _, row in grp.iterrows():
            rpath = row["path"]; date = row["date"]
            try:
                with rasterio.open(rpath) as src:
                    # Reproject polygons to raster CRS
                    shapes = huc12_gdf.to_crs(src.crs)

                    # Limit to polygons that intersect raster bounds (faster)
                    idx = shapes.index
                    if use_spatial_index:
                        # Build index per raster CRS
                        sindex = shapes.sindex
                        rb = box(*src.bounds)
                        cand = list(sindex.intersection(rbinds := rb.bounds))
                        if len(cand) == 0:
                            # no overlap → all NaN for this raster
                            continue
                        shapes = shapes.iloc[cand].copy()
                        idx = shapes.index  # remember which rows we kept

                    # Zonal mean
                    zs = zonal_stats(
                        vectors=list(shapes.geometry),
                        raster=rpath,
                        stats=["mean"],
                        nodata=src.nodata,
                        all_touched=all_touched,
                        geojson_out=False,
                    )
                    means = [d["mean"] if d["mean"] is not None else np.nan for d in zs]

                # Build a frame for just the intersecting HUC12s
                rec = pd.DataFrame({
                    "huc12_id": shapes[id_field].astype(str).values,
                    "date": date,
                    "timescale": ts,
                    "scenario": scen,
                    "model": mdl,
                    "spei": means,
                })
                records.append(rec)

            except Exception as e:
                print(f"[WARN] Zonal means failed for {rpath}: {e}")
                continue

    if not records:
        return pd.DataFrame(columns=["huc12_id","date","timescale","scenario","model","spei"])

    out = pd.concat(records, ignore_index=True)
    # Some HUC12s may never intersect → absent rows, which is fine (saves space).
    return out


# --- Step 4: Decadal drought metrics -----------------------------------------
import numpy as np
import pandas as pd

def _run_stats(bool_series: pd.Series) -> tuple[int, float]:
    """Return (max_run, mean_run) of consecutive True values in a boolean series."""
    if bool_series.empty:
        return 0, 0.0
    runs, run = [], 0
    for v in bool_series.values:
        if v:
            run += 1
        else:
            if run: runs.append(run); run = 0
    if run: runs.append(run)
    return (max(runs) if runs else 0, float(np.mean(runs)) if runs else 0.0)


def compute_decadal_metrics(
    ts_df: pd.DataFrame,
    thresholds = (-1.0, -1.5)
) -> pd.DataFrame:
    """
    For each HUC12/timescale/scenario/model/decade, compute:
      - freq_{thr}: # of months with SPEI <= thr
      - maxrun_{thr}, meanrun_{thr}: lengths of consecutive drought months
      - intensity_mean_{thr}: mean SPEI during drought months (<= thr)
    """
    if ts_df.empty:
        return pd.DataFrame(columns=[
            "huc12_id","timescale","scenario","model","decade","months",
            "freq_le_1","maxrun_le_1","meanrun_le_1","intensity_mean_le_1",
            "freq_le_1.5","maxrun_le_1.5","meanrun_le_1.5","intensity_mean_le_1.5"
        ])

    df = ts_df.copy()
    # ensure datetime
    if not np.issubdtype(df["date"].dtype, np.datetime64):
        df["date"] = pd.to_datetime(df["date"])
    df["year"] = df["date"].dt.year
    df["decade"] = (df["year"] // 10) * 10

    groups = ["huc12_id","timescale","scenario","model","decade"]
    rows = []
    for keys, g in df.groupby(groups, sort=False):
        row = dict(zip(groups, keys))
        row["months"] = int(len(g))
        for thr in thresholds:
            mask = g["spei"] <= thr
            freq = int(mask.sum())
            maxrun, meanrun = _run_stats(mask)
            intensity = float(g.loc[mask, "spei"].mean()) if freq > 0 else np.nan
            lab = f"{abs(thr):.1f}".replace(".0","")
            row[f"freq_le_{lab}"] = freq
            row[f"maxrun_le_{lab}"] = maxrun
            row[f"meanrun_le_{lab}"] = meanrun
            row[f"intensity_mean_le_{lab}"] = intensity
        rows.append(row)

    out = pd.DataFrame(rows).sort_values(groups)
    return out


