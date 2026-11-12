#!/usr/bin/env python3
"""In-situ surface chlorophyll-a from the US Water Quality Portal (WQP) for the Lakes_cci lakes in the USA.

    python3 analysis/lakes_cci/insitu_wqp.py <HydroLAKES polygons .geojson> <out csv>

<polygons>: the Lakes_cci lakes to look up (property Hylak_id), cut from LakeATLAS_v10_pol_Global.gpkg with
ogr2ogr -where "Hylak_id IN (...)". Per lake: HydroLAKES polygon -> WQP stations with
chlorophyll-a inside its bounding box (siteType Lake, Reservoir, Impoundment) -> stations inside the polygon (~50 m
buffer, 0.0005 degrees) -> their results 2002-2023 (POST by site list, 200 sites per request). Kept: laboratory chlorophyll-a (CharacteristicName "Chlorophyll a", "...corrected
for pheophytin", "...uncorrected for pheophytin"; probe/fluorescence values left out), water samples, no quality-control
activities, units ug/L, mg/m3, ppb or mg/L (x1000), sample depth <= 3 m or not given (surface grabs). Output: one row per
sample: hylak_id, site, date, chl (ug/L), char, depth_m. Lakes already in <out csv> are skipped (rerun = resume).
Source: Water Quality Portal (USGS, EPA, NWQMC), https://www.waterqualitydata.us (public domain US data).
"""
import io
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from shapely import contains_xy
from shapely.geometry import shape

WQP = "https://www.waterqualitydata.us/data"
CHARS = ["Chlorophyll a", "Chlorophyll a, corrected for pheophytin", "Chlorophyll a, uncorrected for pheophytin"]
UNITS = {"ug/l": 1.0, "ug/L": 1.0, "mg/m3": 1.0, "mg/m^3": 1.0, "ppb": 1.0, "mg/l": 1000.0, "mg/L": 1000.0}


def get(kind, bbox=None, sites=None):
    """Station search by bounding box (GET), or Result search for a list of sites (POST, JSON body)."""
    body = None
    if sites is None:
        q = [("bBox", ",".join(f"{v:.5f}" for v in bbox)), ("siteType", "Lake, Reservoir, Impoundment"),
             ("startDateLo", "01-01-2002"), ("startDateHi", "12-31-2023"), ("mimeType", "csv")]
        q += [("characteristicName", c) for c in CHARS]
        url = f"{WQP}/{kind}/search?" + urllib.parse.urlencode(q)
    else:
        url = f"{WQP}/{kind}/search?mimeType=csv"
        body = json.dumps({"siteid": sorted(sites), "characteristicName": CHARS, "startDateLo": "01-01-2002",
                           "startDateHi": "12-31-2023"}).encode()
    for k in range(4):
        try:
            req = urllib.request.Request(url, data=body, headers={"Content-Type": "application/json"} if body else {})
            with urllib.request.urlopen(req, timeout=900) as r:
                return pd.read_csv(io.BytesIO(r.read()), low_memory=False)
        except Exception as e:  # noqa: BLE001  (server errors and timeouts: retry, then give up on the lake)
            err = e; time.sleep(30 * (k + 1))
    raise RuntimeError(f"{kind} {bbox}: {err}")


def depth_m(df):
    d = pd.to_numeric(df["ActivityDepthHeightMeasure/MeasureValue"], errors="coerce")
    u = df["ActivityDepthHeightMeasure/MeasureUnitCode"].astype(str)
    top = pd.to_numeric(df["ActivityTopDepthHeightMeasure/MeasureValue"], errors="coerce")
    tu = df["ActivityTopDepthHeightMeasure/MeasureUnitCode"].astype(str)
    d = d.where(~u.str.startswith("ft"), d * 0.3048)
    top = top.where(~tu.str.startswith("ft"), top * 0.3048)
    return d.fillna(top)


def lake(hy, poly):
    b = poly.bounds
    st = get("Station", b)
    if st.empty:
        return None
    inside = contains_xy(poly.buffer(0.0005), st.LongitudeMeasure.values, st.LatitudeMeasure.values)
    sites = set(st.MonitoringLocationIdentifier[inside])
    if not sites:
        return None
    s = sorted(sites)
    r = pd.concat([get("Result", sites=s[i:i + 200]) for i in range(0, len(s), 200)], ignore_index=True)
    r = r[r.MonitoringLocationIdentifier.isin(sites) & r.CharacteristicName.isin(CHARS)
          & (r.ActivityMediaName.astype(str) == "Water")
          & ~r.ActivityTypeCode.astype(str).str.contains("Quality Control", case=False)]
    f = r["ResultMeasure/MeasureUnitCode"].astype(str).map(UNITS)
    v = pd.to_numeric(r.ResultMeasureValue, errors="coerce") * f
    dep = depth_m(r)
    k = v.notna() & (v >= 0) & (dep.isna() | (dep <= 3))
    return pd.DataFrame({"hylak_id": hy, "site": r.MonitoringLocationIdentifier[k], "date": r.ActivityStartDate[k],
                         "chl": v[k].round(3), "char": r.CharacteristicName[k], "depth_m": dep[k].round(2)})


def main(geojson, out):
    out = Path(out)
    P = {int(f["properties"]["Hylak_id"]): shape(f["geometry"]) for f in json.loads(Path(geojson).read_text())["features"]}
    ids = sorted(P)
    done = set(pd.read_csv(out).hylak_id) if out.exists() else set()
    tried = out.with_suffix(".tried")
    done |= set(int(x) for x in tried.read_text().split()) if tried.exists() else set()
    todo = [h for h in ids if h not in done]
    print(f"{len(ids)} lakes, {len(todo)} to do", flush=True)
    for hy in todo:
        t0 = time.time()
        try:
            d = lake(hy, P[hy])
        except RuntimeError as e:
            print(f"{hy}: {e}", flush=True); continue
        if d is not None and len(d):
            d.to_csv(out, mode="a", header=not out.exists(), index=False)
        with tried.open("a") as f:
            f.write(f"{hy}\n")
        print(f"{hy}: {0 if d is None else len(d)} samples, {0 if d is None else d.site.nunique()} sites, "
              f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
