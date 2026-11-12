#!/usr/bin/env python3
"""One table of in-situ surface chlorophyll-a for lakes that have Lakes_cci series (input of insitu_timing.py).

    python3 analysis/lakes_cci/insitu_chl.py <Lakes_cci lake list csv (hylak_id, central lat/lon)> <out csv> [<wqp csv> <polygons .geojson> [<more csv> ...]]

Sources (each value: surface sample, depth <= 3 m, or an epilimnion/0-5 m integrated sample):
  - EDI: Lewis et al. (2023), edi.1530.1, CC BY 4.0 (EDI_DIR, default ~/data/edi.1530.1: productivity_EDI.csv and
    lake_metadata_EDI.csv; HydroLAKES id given there); skipped when absent
  - local: project observations (LOCAL_OBS, default calibration/obs_wq/obs_wq.csv; HydroLAKES id from LOCAL_IDS,
    default global_test/local_vs_global.csv); not public, skipped when absent
  - WQP: US Water Quality Portal (insitu_wqp.py output; site inside the HydroLAKES polygon)
  - EEA Waterbase, GEMStat, SYKE (insitu_more.py outputs, <more csv>; national monitoring archives, site inside the polygon)
A lake found in several sources keeps only the source with the most sampling dates (no double counting; e.g. Finnish
sites are in both SYKE and EEA Waterbase). Lake list: Lakes_cci v2.1 per-lake set (Shahvaran et al.) or, for v3.0, a csv
with a hylak_id column (lakes without central coordinates take the mean site position).
Only lakes whose HydroLAKES id is in the Lakes_cci list are kept. Output: source, lake, hylak_id, lat, lon, date, chl
(ug/L; daily mean over sites and depths of that source), area_km2 (EDI metadata, HydroLAKES for local and WQP lakes).
"""
import json
import sys
from pathlib import Path

import pandas as pd

import os

ROOT = Path(__file__).resolve().parents[2]
EDI = Path(os.environ.get("EDI_DIR", Path.home() / "data/edi.1530.1"))           # Lewis et al. (2023) files
LOCAL_OBS = Path(os.environ.get("LOCAL_OBS", ROOT / "calibration/obs_wq/obs_wq.csv"))   # project data (not public)
LOCAL_IDS = Path(os.environ.get("LOCAL_IDS", ROOT / "global_test/local_vs_global.csv"))


def edi():
    p = pd.read_csv(EDI / "productivity_EDI.csv", usecols=["LakeID", "Date", "Depth_m", "Interval", "Chla_ugL"])
    p = p.dropna(subset=["Chla_ugL"])
    iv = p.Interval.astype(str).str.strip()
    p = p[(p.Depth_m <= 3) | iv.isin(["EPILIMNION", "0 - 5 int"])]
    m = pd.read_csv(EDI / "lake_metadata_EDI.csv").dropna(subset=["hylak_id"])
    p = p.merge(m[["LakeID", "LakeName", "hylak_id", "Latitude_DD", "Longitude_DD", "SurfaceArea_ha"]], on="LakeID")
    return pd.DataFrame({"source": "EDI", "lake": p.LakeName, "hylak_id": p.hylak_id.astype(int),
                         "lat": p.Latitude_DD, "lon": p.Longitude_DD, "date": pd.to_datetime(p.Date), "chl": p.Chla_ugL,
                         "area_km2": p.SurfaceArea_ha / 100})


def local():
    o = pd.read_csv(LOCAL_OBS)
    o = o[(o.variable == "chl") & (o.depth <= 3)]
    g = pd.read_csv(LOCAL_IDS).dropna(subset=["hylak_id_final"])
    o = o.merge(g[["lake", "hylak_id_final", "lat", "lon", "hl_area_km2"]], on="lake")
    return pd.DataFrame({"source": "local", "lake": o.lake, "hylak_id": o.hylak_id_final.astype(int), "lat": o.lat,
                         "lon": o.lon, "date": pd.to_datetime(o.datetime).dt.normalize(), "chl": o.value, "area_km2": o.hl_area_km2})


def names(geojson):
    """HydroLAKES name, area and a point inside each polygon (position where the lake list gives none)."""
    from shapely.geometry import shape
    F = json.loads(Path(geojson).read_text())["features"]
    pt = [shape(x["geometry"]).representative_point() for x in F]
    pr = pd.DataFrame([x["properties"] for x in F]).assign(plat=[q.y for q in pt], plon=[q.x for q in pt])
    return pr.drop_duplicates("Hylak_id").set_index("Hylak_id")


def wqp(f, pr, cci):
    w = pd.read_csv(f)
    c = cci.set_index("hylak_id")
    name = w.hylak_id.map(pr.Lake_name).fillna("") + " (" + w.hylak_id.astype(str) + ")"
    return pd.DataFrame({"source": "WQP", "lake": name.str.strip(), "hylak_id": w.hylak_id,
                         "area_km2": w.hylak_id.map(pr.Lake_area),
                         "lat": w.hylak_id.map(c.central_lat_wgs84_dd), "lon": w.hylak_id.map(c.central_lon_wgs84_dd),
                         "date": pd.to_datetime(w.date, errors="coerce"), "chl": w.chl})


NAMES = {16679: "Ohau", 16680: "Benmore"}   # HydroLAKES polygons without a name, identified from their position


def more(f, pr):
    w = pd.read_csv(f)
    name = w.hylak_id.map(pr.Lake_name).fillna(w.hylak_id.map(NAMES)).fillna("") + " (" + w.hylak_id.astype(str) + ")"
    return pd.DataFrame({"source": w.source, "lake": name.str.strip(), "hylak_id": w.hylak_id, "area_km2": w.hylak_id.map(pr.Lake_area),
                         "lat": w.lat, "lon": w.lon, "date": pd.to_datetime(w.date, errors="coerce"), "chl": w.chl})


def main(cci_csv, out, wqp_csv=None, geojson=None, *more_csv):
    cci = pd.read_csv(cci_csv)
    if "hylak_id" not in cci:
        cci = cci.rename(columns={"Hylak_id": "hylak_id"}).dropna(subset=["hylak_id"]).astype({"hylak_id": int})
    for k in ("central_lat_wgs84_dd", "central_lon_wgs84_dd"):
        if k not in cci:
            cci[k] = float("nan")
    pr = names(geojson) if geojson else None
    parts = ([edi()] if (EDI / "productivity_EDI.csv").exists() else []) + ([local()] if LOCAL_OBS.exists() else []) + ([wqp(wqp_csv, pr, cci)] if wqp_csv else []) + [more(f, pr) for f in more_csv]
    d = pd.concat(parts, ignore_index=True).dropna(subset=["date", "chl"])
    d = d[d.hylak_id.isin(set(cci.hylak_id)) & (d.date >= "2002-01-01")]
    # one position per lake and source (site means; WQP, EDI and local already give one)
    if pr is not None:
        d["lat"], d["lon"] = d.lat.fillna(d.hylak_id.map(pr.plat)), d.lon.fillna(d.hylak_id.map(pr.plon))
    d[["lat", "lon"]] = d.groupby(["source", "hylak_id"])[["lat", "lon"]].transform("mean")
    d = d.groupby(["source", "lake", "hylak_id", "date"], as_index=False).agg(lat=("lat", "first"), lon=("lon", "first"),
        area_km2=("area_km2", "first"), chl=("chl", "mean"))
    n = d.groupby(["hylak_id", "source"]).date.nunique().reset_index()
    best = n.sort_values(["date", "source"], ascending=[False, True], kind="stable").drop_duplicates("hylak_id")   # ties: source name
    dup = n[~n.set_index(["hylak_id", "source"]).index.isin(best.set_index(["hylak_id", "source"]).index)]
    if len(dup):
        print(f"{dup.hylak_id.nunique()} lakes in more than one source; kept the source with most dates, dropped:",
              dup.groupby("source").size().to_dict())
    d = d.merge(best[["hylak_id", "source"]], on=["hylak_id", "source"])
    d.to_csv(out, index=False, float_format="%.4g")
    s = d.groupby(["source", "hylak_id"]).agg(n=("chl", "size"), months=("date", lambda x: x.dt.month.nunique()))
    print(s.groupby("source").agg(lakes=("n", "size"), dates=("n", "sum"), with_20_dates_5_months=("n", lambda x: ((x >= 20) & (s.loc[x.index, "months"] >= 5)).sum())))


if __name__ == "__main__":
    main(*sys.argv[1:])
