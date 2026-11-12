#!/usr/bin/env python3
"""In-situ surface chlorophyll-a from national monitoring archives, matched to the Lakes_cci lakes (input of insitu_chl.py).

    python3 analysis/lakes_cci/insitu_more.py <source> <CCI lake polygons .geojson> <out csv> [<source dir>]

Every source also writes <out>_sites.csv: all sites with chl-a it holds (inside a Lakes_cci lake or not; for gemstat,
cogerh, epavic, ecan), the inventory shown on the report map.
<source>: eea (EEA Waterbase Water Quality ICM 2026, <dir> with WISE6_DisaggregatedData-csv.zip and
WISE6_SpatialObjects_DerivedData-csv.zip), gemstat (UNEP GEMStat open subset GFQA v3, Zenodo 10.5281/zenodo.18459694,
<dir> with Pigment.csv and GEMStat_station_metadata.csv), syke (Finnish Environment Institute VESLA open API, fetched
here), cogerh (Ceara state reservoir monitoring, Brazil, fetched here; <dir> caches the campaign pages; no Ceara reservoir
is a Lakes_cci lake, kept for model checks), epavic (EPA Victoria Gippsland Lakes, <dir> with the DataVic xlsx), ecan
(Environment Canterbury Hilltop server, fetched here). EDI/NTL-LTER, Sweden's MVM API and DataStream (Canada) need a log-in or key and are not fetched here; Ontario's
open inland-lake chemistry has no Lakes_cci lake sampled more than once a year. <polygons>: HydroLAKES polygons of the Lakes_cci lakes (property Hylak_id), cut from LakeATLAS_v10_pol_Global.gpkg.

Independence from Lakes_cci: Lakes_cci chl-a algorithms were tuned and validated against LIMNADES, a research archive of
bio-optical campaign data, and "very few data sets from national monitoring programmes" were used (Lakes_cci PVP
CCI-LAKES-0030 v1.3, section 6.3; ATBD D2.2). All three sources here are statutory monitoring archives of water-sample
chlorophyll-a; no value comes from satellite or model data. SYKE values are laboratory extractions (analyte CP; in-field fluorescence left out); EEA and
most GEMStat values do not report the analysis method. Per sample: site inside the lake polygon (~50 m buffer, 0.0005 degrees), 2002-2023, sample depth <= 3 m or not
given, chl-a in ug/L; values below the limit of quantification are set to half of it (EEA). Output, one row per sample:
source, site, hylak_id, lat, lon, date, chl, depth_m, method.
"""
import io
import json
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from shapely import STRtree, points
from shapely.geometry import shape

Y0, Y1 = "2002-01-01", "2023-12-31"


def polygons(geojson):
    F = json.loads(Path(geojson).read_text())["features"]
    ids = np.array([int(f["properties"]["Hylak_id"]) for f in F])
    return ids, [shape(f["geometry"]) for f in F]


def match(sites, geojson):
    """sites: DataFrame with site, lat, lon -> adds hylak_id (sites outside every CCI lake are dropped)."""
    ids, P = polygons(geojson)
    s = sites.dropna(subset=["lat", "lon"]).drop_duplicates("site").copy()
    tree = STRtree(P)
    pi, gi = tree.query(points(s.lon.values, s.lat.values), predicate="dwithin", distance=0.0005)
    s = s.iloc[pi].assign(hylak_id=ids[gi])
    return s.drop_duplicates("site")[["site", "hylak_id"]]


def eea(d, geojson):
    d = Path(d)
    with zipfile.ZipFile(d / "WISE6_SpatialObjects_DerivedData-csv.zip") as z:
        sp = pd.read_csv(z.open([n for n in z.namelist() if n.endswith(".csv")][0]), low_memory=False)
    sp = sp.dropna(subset=["monitoringSiteIdentifier"]).rename(columns={"monitoringSiteIdentifier": "site"})
    hy = match(sp[["site", "lat", "lon"]], geojson)
    print(f"EEA: {len(sp)} sites, {len(hy)} inside CCI lakes ({hy.hylak_id.nunique()} lakes)", flush=True)
    keep = set(hy.site)
    cols = ["monitoringSiteIdentifier", "parameterWaterBodyCategory", "observedPropertyDeterminandCode", "procedureAnalysedMatrix",
            "resultUom", "phenomenonTimeSamplingDate", "resultObservedValue", "resultQualityObservedValueBelowLOQ", "procedureLOQValue",
            "parameterSampleDepth", "metadata_observationStatus", "metadata_statusCode"]
    raw = d / "eea_chl_raw.csv"   # cache of the chl-a rows at CCI sites (one pass over the 28 GB table)
    if raw.exists():
        r = pd.read_csv(raw, low_memory=False)
        return eea_clean(r, hy, sp)
    out = []
    # the 28 GB table is Deflate64-zipped (not readable by zipfile): stream it through unzip, keeping only chl-a lines
    cmd = f"unzip -p '{d / 'WISE6_DisaggregatedData-csv.zip'}' | awk 'NR == 1 || /EEA_3164-01-0/'"
    with subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE) as pr:
        for ch in pd.read_csv(pr.stdout, usecols=cols, chunksize=1_000_000, low_memory=False):
            ch = ch[(ch.observedPropertyDeterminandCode == "EEA_3164-01-0") & ch.monitoringSiteIdentifier.isin(keep)]
            if len(ch):
                out.append(ch)
    r = pd.concat(out, ignore_index=True)
    r.to_csv(raw, index=False)
    return eea_clean(r, hy, sp)


def eea_clean(r, hy, sp):
    n0 = len(r)
    f = r.resultUom.astype(str).map({"ug/L": 1.0, "ug/l": 1.0, "mg/m3": 1.0, "mg/L": 1000.0, "mg/l": 1000.0})
    v = pd.to_numeric(r.resultObservedValue, errors="coerce")
    loq = r.resultQualityObservedValueBelowLOQ.astype(str).str.lower().isin(["true", "1"])
    v = v.where(~loq, pd.to_numeric(r.procedureLOQValue, errors="coerce").fillna(v) / 2) * f
    ok = (r.procedureAnalysedMatrix.astype(str) == "W") & (r.metadata_observationStatus.astype(str).isin(["A", "nan"])) \
        & (r.metadata_statusCode.astype(str) != "invalid") & v.notna() & (v >= 0) \
        & (r.parameterSampleDepth.isna() | (r.parameterSampleDepth <= 3))
    r = r[ok].assign(chl=v[ok])
    print(f"EEA: {n0} chl-a rows at those sites, {len(r)} kept", flush=True)
    r = r.rename(columns={"monitoringSiteIdentifier": "site", "phenomenonTimeSamplingDate": "date", "parameterSampleDepth": "depth_m"})
    r["date"] = pd.to_datetime(r.date.astype("Int64").astype(str), format="%Y%m%d", errors="coerce")   # stored as yyyymmdd
    r = r.merge(hy, on="site").merge(sp[["site", "lat", "lon"]].drop_duplicates("site"), on="site")
    return r.assign(source="EEA", method="lab (WISE SoE determinand EEA_3164-01-0)")


def gemstat(d, geojson):
    d = Path(d)
    st = pd.read_csv(d / "GEMStat_station_metadata.csv", encoding="latin1")
    st = st.rename(columns={"GEMS Station Number": "site", "Latitude": "lat", "Longitude": "lon"}).drop_duplicates("site")
    hy = match(st[["site", "lat", "lon"]], geojson)
    print(f"GEMStat: {len(st)} stations, {len(hy)} inside CCI lakes ({hy.hylak_id.nunique()} lakes)", flush=True)
    p = pd.read_csv(d / "Pigment.csv", encoding="latin1", low_memory=False)
    p = p[p["Parameter Code"] == "Chl-a"]
    lk = st[st["Water Type"].isin(["Lake station", "Reservoir station"])]
    found(p.rename(columns={"GEMS Station Number": "site", "Sample Date": "date"}).merge(lk[["site", "lat", "lon"]], on="site"),
          "GEMStat", "lake or reservoir station", geojson)
    p = p[p["GEMS Station Number"].isin(set(hy.site))]
    f = p.Unit.astype(str).str.lower().map({"mg/l": 1000.0, "µg/l": 1.0, "ug/l": 1.0, "mg/m3": 1.0})
    v = pd.to_numeric(p.Value, errors="coerce") * f
    dep = pd.to_numeric(p.Depth, errors="coerce")
    ok = v.notna() & (v >= 0) & (dep.isna() | (dep <= 3)) & (p["Data Quality"].astype(str) != "Suspect") \
        & (p["Value Flags"].astype(str).str.strip() != "<")   # "<": below detection, value is the limit
    p = p[ok].assign(chl=v[ok], depth_m=dep[ok])
    p = p.rename(columns={"GEMS Station Number": "site", "Sample Date": "date", "Analysis Method Code": "method"})
    p = p.merge(hy, on="site").merge(st[["site", "lat", "lon"]], on="site")
    return p.assign(source="GEMStat", method="lab (GEMStat " + p.method.astype(str) + ")")


SYKE = "https://rajapinnat.ymparisto.fi/api/vesla/2.0/odata"


def odata(path):
    url = f"{SYKE}/{path}"
    for k in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=300) as r:
                return json.load(r)["value"]
        except Exception as e:  # noqa: BLE001  (server errors and timeouts: retry, then stop)
            err = e; time.sleep(20 * (k + 1))
    raise RuntimeError(f"{url}: {err}")


def syke(_, geojson):
    rows, k = [], 0
    while True:   # all lake sites (paged, 1000 per page)
        v = odata("Site_Wide?$filter=EnvironmentType%20eq%20'Lake'&$select=Site_Id,Name,Lake,CoordEUREFFIN_WGS84_Lat,"
                  f"CoordEUREFFIN_WGS84_Long&$top=1000&$skip={k}")
        rows += v; k += 1000
        if len(v) < 1000:
            break
    st = pd.DataFrame(rows).rename(columns={"Site_Id": "site", "CoordEUREFFIN_WGS84_Lat": "lat", "CoordEUREFFIN_WGS84_Long": "lon"})
    st["lat"], st["lon"] = pd.to_numeric(st.lat, errors="coerce"), pd.to_numeric(st.lon, errors="coerce")
    hy = match(st[["site", "lat", "lon"]], geojson)
    print(f"SYKE: {len(st)} lake sites, {len(hy)} inside CCI lakes ({hy.hylak_id.nunique()} lakes)", flush=True)
    out = []
    for s in hy.site:   # laboratory chlorophyll-a (analyte CP: acetone/ethanol extraction), 2002-2023
        k = 0
        while True:
            v = odata(f"Result_Wide?$filter=Site_Id%20eq%20{int(s)}%20and%20AnalyteCode%20eq%20'CP'%20and%20Time%20ge%20"
                      f"datetimeoffset'{Y0}T00:00:00Z'&$select=Site_Id,Time,SampleDepth_m,SampleDepthUpper_m,SampleDepthLower_m,Flag,"
                      f"Value,Unit,DeterminationName&$top=1000&$skip={k}")
            out += v; k += 1000
            if len(v) < 1000:
                break
    r = pd.DataFrame(out)
    n0 = len(r)
    dep = pd.to_numeric(r.SampleDepthLower_m.fillna(r.SampleDepth_m), errors="coerce")   # composite samples: their lower depth
    v = pd.to_numeric(r.Value, errors="coerce") * r.Unit.astype(str).str.lower().map({"µg/l": 1.0, "ug/l": 1.0, "mg/m3": 1.0, "mg/l": 1000.0})
    field = r.DeterminationName.astype(str).str.contains("in field", case=False)   # in-situ fluorescence probes: left out (lab only)
    ok = v.notna() & (v >= 0) & (dep.isna() | (dep <= 3)) & ~r.Flag.astype(str).str.contains("<") & ~field
    r = r[ok].assign(chl=v[ok], depth_m=dep[ok], date=r.Time.str[:10])
    print(f"SYKE: {n0} chl-a results, {len(r)} kept", flush=True)
    r = r.rename(columns={"Site_Id": "site", "DeterminationName": "method"}).merge(hy, on="site").merge(st[["site", "lat", "lon"]], on="site")
    return r.assign(source="SYKE", method="lab (" + r.method.astype(str).str.strip() + ")")


COGERH = "http://www.hidro.ce.gov.br/acude/eutrofizacao"


def cogerh(d, geojson):
    """Ceara (Brazil) reservoirs: COGERH quarterly trophic-state campaigns (Feb, May, Aug, Nov; 2008-), shown on the state's
    Portal Hidrologico; one page per campaign with the sample points as embedded JSON (cached as <dir>/<campaign>.json)."""
    d = Path(d); d.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(COGERH, timeout=120) as r:
        camps = sorted(set(re.findall(r'option value="([A-Z]{3}/\d{4})"', r.read().decode("utf-8", "replace"))))
    rows = []
    for c in camps:
        f = d / (c.replace("/", "_") + ".json")
        if not f.exists():
            with urllib.request.urlopen(f"{COGERH}?campanha={urllib.parse.quote(c)}", timeout=120) as r:
                m = re.search(r"var pontos_eutrofizados = (\[.*?\]);", r.read().decode("utf-8", "replace"), re.S)
            f.write_text(m.group(1) if m else "[]"); time.sleep(2)
        rows += [x["reservatorio"] | {"campanha": c} for x in json.loads(f.read_text())]
    r = pd.DataFrame(rows)
    r["site"] = r.nome.astype(str) + " " + r.ponto.astype(str)
    hy = match(r.rename(columns={"latitude": "lat", "longitude": "lon"})[["site", "lat", "lon"]], geojson)
    print(f"COGERH: {len(camps)} campaigns, {r.site.nunique()} points, {len(hy)} inside CCI lakes ({hy.hylak_id.nunique()} lakes)", flush=True)
    v = pd.to_numeric(r.clorofila, errors="coerce"); dep = pd.to_numeric(r.profundidade, errors="coerce")
    ok = v.notna() & (v >= 0) & (dep.isna() | (dep <= 3))
    r = r[ok].assign(chl=v[ok], depth_m=dep[ok]).rename(columns={"data": "date", "latitude": "lat", "longitude": "lon"})
    found(r, "COGERH", "reservoir point", geojson)
    return r.merge(hy, on="site").assign(source="COGERH", method="lab (COGERH trophic-state monitoring)")


def epavic(d, geojson):
    """EPA Victoria, Gippsland Lakes water quality 1990-2026 (DataVic, CC BY 4.0; monthly lab chl-a by spectrophotometry).
    <dir>: folder with 1990_01-2026_05_Gippsland_Lakes_Water_quality_data.xlsx (from apps.epa.vic.gov.au/datavic)."""
    f = sorted(Path(d).glob("*Gippsland_Lakes_Water_quality_data.xlsx"))[-1]
    x = pd.read_excel(f, sheet_name=None)
    st = x["Site Metadata"].rename(columns={"site_id": "site", "latitude": "lat", "longitude": "lon"})
    hy = match(st[["site", "lat", "lon"]], geojson)
    v = x["Data"]
    v = v[v.Type.astype(str).str.lower().str.startswith("surf") & v.CHL_A.notna() & (v.CHL_A >= 0)]
    found(v.assign(site=v.site_id).merge(st, on="site"), "EPA Victoria", "lake", geojson)
    print(f"EPA Victoria: {len(st)} sites, {len(hy)} inside CCI lakes ({hy.hylak_id.nunique()} lakes)", flush=True)
    r = v.rename(columns={"site_id": "site", "CHL_A": "chl"}).merge(hy, on="site").merge(st[["site", "lat", "lon"]], on="site")
    return r.assign(depth_m=np.nan, source="EPA Victoria", method="lab (spectrophotometry)")


HILLTOP = {"ECan": ("http://wateruse.ecan.govt.nz/WQAll.hts", "Chlorophyll a (planktonic)")}


def hilltop(agency, geojson, cache):
    """New Zealand regional council Hilltop servers (public web service behind LAWA): sites with planktonic chl-a and
    coordinates, their samples 2002-2023 (only time, value and method are kept; samples flagged '<' are left out)."""
    import xml.etree.ElementTree as ET
    url, meas = HILLTOP[agency]
    q = lambda **k: url + "?service=Hilltop&" + urllib.parse.urlencode(k, quote_via=urllib.parse.quote)   # %20, not +
    with urllib.request.urlopen(q(request="SiteList", Location="LatLong", Measurement=meas), timeout=300) as r:
        root = ET.fromstring(r.read())
    st = pd.DataFrame([{"site": s.get("Name"), "lat": float(s.findtext("Latitude") or "nan"), "lon": float(s.findtext("Longitude") or "nan")}
                       for s in root.iter("Site")]).dropna()
    hy = match(st, geojson)
    print(f"{agency}: {len(st)} sites with {meas}, {len(hy)} inside CCI lakes ({hy.hylak_id.nunique()} lakes)", flush=True)
    cache = Path(cache); rows = []
    for s in ([] if cache.exists() else st.site):   # all sites with chl-a (also outside CCI lakes: inventory of what exists)
        with urllib.request.urlopen(q(request="GetData", Site=s, Measurement=meas, From="2002-01-01", To="2023-12-31"), timeout=300) as r:
            for e in ET.fromstring(r.read()).iter("E"):
                m = {p.get("Name"): p.get("Value") for p in e.iter("Parameter")}
                rows.append({"site": s, "date": e.findtext("T"), "val": e.findtext("Value"), "method": m.get("MethodText", "")})
        time.sleep(0.5)
    if rows:
        pd.DataFrame(rows).to_csv(cache, index=False)
    r = pd.read_csv(cache)
    r = r[~r.val.astype(str).str.contains("<")]
    r["chl"] = pd.to_numeric(r.val, errors="coerce"); r = r[r.chl.notna() & (r.chl >= 0)]
    found(r.merge(st, on="site"), agency, "site with planktonic chl-a", geojson)
    r = r.merge(hy, on="site").merge(st, on="site")
    return r.assign(depth_m=np.nan, source=agency, method="lab (" + r.method.astype(str) + ")")


FOUND = []


def found(r, source, kind, geojson):
    """Inventory of every site with in-situ chl-a seen in a source, inside a Lakes_cci lake or not (for the map)."""
    r = r.assign(date=pd.to_datetime(r.date, errors="coerce")).dropna(subset=["date", "lat", "lon"])
    r = r[(r.date >= Y0) & (r.date <= Y1)]
    g = r.groupby("site").agg(lat=("lat", "first"), lon=("lon", "first"), n_dates=("date", "nunique"),
                              months=("date", lambda x: x.dt.month.nunique())).reset_index()
    hy = match(g[["site", "lat", "lon"]], geojson)
    FOUND.append(g.merge(hy, on="site", how="left").assign(source=source, kind=kind))


def main(source, geojson, out, d=None):
    F = {"eea": eea, "gemstat": gemstat, "syke": syke, "cogerh": cogerh, "epavic": epavic,
         "ecan": lambda _, g: hilltop("ECan", g, Path(out).with_name("ecan_raw.csv"))}
    r = F[source](d, geojson)
    if FOUND:   # all sites with chl-a in this source (also outside Lakes_cci lakes)
        pd.concat(FOUND).to_csv(Path(out).with_name(Path(out).stem + "_sites.csv"), index=False, float_format="%.5g")
    if r.empty:
        print(f"{source}: no samples inside Lakes_cci lakes"); return
    r["date"] = pd.to_datetime(r.date, errors="coerce").dt.normalize()
    r = r[(r.date >= Y0) & (r.date <= Y1)]
    cols = ["source", "site", "hylak_id", "lat", "lon", "date", "chl", "depth_m", "method"]
    r[cols].to_csv(out, index=False, float_format="%.4g")
    g = r.groupby("hylak_id").agg(n=("date", "nunique"), months=("date", lambda x: x.dt.month.nunique()))
    print(f"{source}: {len(r)} samples, {len(g)} lakes, {(g.n >= 20).sum()} with >= 20 dates, {((g.n >= 20) & (g.months >= 5)).sum()} "
          f"of them over >= 5 months -> {out}")


if __name__ == "__main__":
    main(*sys.argv[1:5])
