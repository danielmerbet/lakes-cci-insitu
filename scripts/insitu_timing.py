#!/usr/bin/env python3
"""Does Lakes_cci chlorophyll-a get the seasonal timing right? Comparison with in-situ surface chlorophyll-a.

    python3 analysis/lakes_cci/insitu_timing.py <insitu_chl.csv> <cci source> <out dir> [<*_sites.csv> ...]

<insitu_chl.csv>: insitu_chl.py output (source, lake, hylak_id, lat, lon, date, chl).
<cci source>: either a folder of per-lake daily CSVs (Shahvaran et al., ESSD preprint essd-2026-528, Lakes_cci v2.1,
  <hylak_id>.csv) -> label "v2.1"; or "v3:<kd dir>:<cci_registry_match.csv>" (extract_kd.py tables, Lakes_cci v3.0,
  chla_median per lake and day, mapped to HydroLAKES with the match table) -> label "v3.0".
CCI days: coverage >= 10 % of the lake (v2.1: chla_mean_coverage_ratio; v3.0: n_px / n_lake_px), 2002-2022.
Lakes: >= 20 in-situ dates in 2002-2022 spread over >= 5 calendar months, and >= 100 CCI days (tier "standard").
Tropical lakes (|lat| < 23.5) with fewer data: >= 10 dates over >= 4 months, in-situ months kept with data in one year,
cycles compared on >= 4 common months (tier "tropical low-data"; few open tropical series exist).

Per lake (lakes.csv):
  matchups (same day): n_match, bias (median log10 CCI/in-situ), r_match (log10, n >= 10);
  seasonal cycles: monthly means per year, then the mean over years per calendar month (month kept if it has data
  in >= 2 years for in-situ), on the months both have (>= 5): peak month of each, dpeak (CCI - in-situ, circular,
  -6..6), r_season (log10 climatologies), and the same from matchup days only (dpeak_m, r_season_m), which removes
  the effect of when each source samples and leaves the retrieval.
Southern-hemisphere months are shifted by 6 in the pooled monthly figures (month 1 = January in the north, July in the south).
Writes lakes.csv, matchups.csv, monthly_bias.csv, cycles.json (per-lake monthly cycles, for the report) and
fig1..fig7 PNGs to <out dir> (fig7: tropical lakes). <*_sites.csv> (insitu_more.py site inventories) add all other
sites with in-situ chl-a to the map (fig1) and to sites_found.csv.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plotstyle as PL  # noqa: E402
plt = PL.plt

FIRST, LAST, MIN_COVER = "2002-01-01", "2022-12-31", 10.0
SEASON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
          9: "SON", 10: "SON", 11: "SON"}
SCOL = {"DJF": PL.CAT[0], "MAM": PL.CAT[2], "JJA": PL.CAT[1], "SON": PL.CAT[3]}
lg = lambda x: np.log10(np.maximum(np.asarray(x, float), 0.05))
MON = "JFMAMJJASOND"


# ---------------------------------------------------------------- CCI series
class CCI:
    def __init__(self, spec):
        if spec.startswith("v3:"):
            _, kd, match = spec.split(":", 2)
            m = pd.read_csv(match).dropna(subset=["Hylak_id"])
            cmap = dict(zip(m.id.astype(int), m.Hylak_id.astype(int)))
            K = pd.concat([pd.read_csv(f, parse_dates=["date"]) for f in sorted(Path(kd).glob("kd_*.csv"))])
            K = K[K.chla_median.notna() & (100 * K.n_px / K.n_lake_px >= MIN_COVER)]
            K["hylak_id"] = K.cci_id.map(cmap)
            self.g = {h: d.set_index("date").chla_median for h, d in K.dropna(subset=["hylak_id"]).groupby("hylak_id")}
            self.label = "v3.0"
        else:
            self.dir, self.label, self.g = Path(spec), "v2.1", None

    def get(self, hy):
        if self.g is not None:
            s = self.g.get(hy)
        else:
            f = self.dir / f"{int(hy)}.csv"
            if not f.exists():
                return None
            o = pd.read_csv(f, parse_dates=["date"]).set_index("date")
            s = o.chla_mean_mean.where(o.chla_mean_coverage_ratio >= MIN_COVER).dropna()
        if s is None:
            return None
        s = s[(s.index >= FIRST) & (s.index <= LAST)]
        return s.groupby(s.index).mean()


# ---------------------------------------------------------------- metrics
def clim(s, min_years=1):
    g = s.groupby([s.index.year, s.index.month]).mean()
    c = g.groupby(level=1).agg(["mean", "size"])
    return c["mean"][c["size"] >= min_years]


def cdiff(a, b):
    return int((a - b + 6) % 12 - 6)


def season_stats(o, c, tag="", min_years=2, min_months=5):
    co, cc = clim(o, min_years if not tag else 1), clim(c)
    m = co.index.intersection(cc.index)
    if len(m) < min_months:
        return {}
    co, cc = co[m], cc[m]
    return {f"months{tag}": len(m), f"peak_insitu{tag}": int(co.idxmax()), f"peak_cci{tag}": int(cc.idxmax()),
            f"dpeak{tag}": cdiff(cc.idxmax(), co.idxmax()), f"r_season{tag}": float(np.corrcoef(lg(co), lg(cc))[0, 1])}


def analyse(ins, cci, out):
    ins = ins[(ins.date >= FIRST) & (ins.date <= LAST)]
    rows, M, cyc = [], [], {}
    for (src, lake, hy), d in ins.groupby(["source", "lake", "hylak_id"]):
        o = d.set_index("date").chl.groupby(level=0).mean()
        nm = o.index.month.nunique()
        tier = "standard" if len(o) >= 20 and nm >= 5 else \
            "tropical low-data" if abs(d.lat.iloc[0]) < 23.5 and len(o) >= 10 and nm >= 4 else None
        if tier is None:
            continue
        my, mmin = (2, 5) if tier == "standard" else (1, 4)
        c = cci.get(hy)
        if c is None or len(c) < 100:
            continue
        mm = pd.DataFrame({"insitu": o, "cci": c}).dropna()
        r = {"source": src, "lake": lake, "hylak_id": hy, "tier": tier, "lat": d.lat.iloc[0], "lon": d.lon.iloc[0], "area_km2": d.area_km2.iloc[0],
             "n_insitu": len(o), "n_cci": len(c), "n_match": len(mm), "insitu_mean": o.mean(), "cci_mean": c.mean()}
        if len(mm):
            e = lg(mm.cci) - lg(mm.insitu)
            r["bias"] = float(np.median(e))
            r["r_match"] = float(np.corrcoef(lg(mm.insitu), lg(mm.cci))[0, 1]) if len(mm) >= 10 else np.nan
            M.append(mm.assign(hylak_id=hy, lake=lake, lat=r["lat"], logratio=e).rename_axis("date").reset_index())
        r.update(season_stats(o, c, min_years=my, min_months=mmin))
        if len(mm) >= 15:
            r.update(season_stats(mm.insitu, mm.cci, "_m"))
        rows.append(r)
        cyc[hy] = (clim(o, my), clim(c), clim(mm.insitu) if len(mm) >= 15 else None, clim(mm.cci) if len(mm) >= 15 else None)
    L = pd.DataFrame(rows)
    M = pd.concat(M, ignore_index=True) if M else pd.DataFrame()
    if len(M):
        M["month"] = M.date.dt.month
        M["month_nh"] = np.where(M.lat < 0, (M.month + 5) % 12 + 1, M.month)
        M["season"] = M.month_nh.map(SEASON)
    L.to_csv(out / "lakes.csv", index=False, float_format="%.4g")
    arr = lambda s: None if s is None else [None if m not in s.index else round(float(s[m]), 3) for m in range(1, 13)]
    J = []
    for _, r in L.iterrows():
        co, cc, mo, mc = cyc[r.hylak_id]
        J.append({k: (None if pd.isna(v) else (round(v, 3) if isinstance(v, float) else v)) for k, v in r.items()} |
                 {"insitu": arr(co), "cci": arr(cc), "cci_m": arr(mc), "insitu_m": arr(mo)})
    (out / "cycles.json").write_text(json.dumps(J, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    M.to_csv(out / "matchups.csv", index=False, float_format="%.4g")
    return L, M, cyc


# ---------------------------------------------------------------- figures
def fig_map(L, out, label, ins=None, found=None):
    """All lakes: used (colour = source, ring = tropical low-data tier), Lakes_cci lakes with too few in-situ dates (open),
    other sites with in-situ chl-a outside Lakes_cci lakes (grey; sources with a site inventory)."""
    import cartopy.crs as ccrs
    import cartopy.feature as cf
    fig = plt.figure(figsize=(12, 6.2))
    ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=10))
    ax.set_global(); ax.add_feature(cf.LAND, facecolor="#ecebe4", edgecolor="none")
    ax.add_feature(cf.COASTLINE, linewidth=0.3, edgecolor=PL.AXIS)
    pc = ccrs.PlateCarree()
    if found is not None and len(found):
        o = found[found.hylak_id.isna()]
        ax.scatter(o.lon, o.lat, s=5, color="#8a8f8c", alpha=0.55, linewidth=0, transform=pc, zorder=2,
                   label=f"in-situ chl-a, not a Lakes_cci lake ({len(o):,} sites)")
    if ins is not None:
        few = ins[~ins.hylak_id.isin(L.hylak_id)].groupby("hylak_id").agg(lat=("lat", "mean"), lon=("lon", "mean"))
        ax.scatter(few.lon, few.lat, s=14, facecolor="none", edgecolor=PL.INK2, linewidth=0.8, transform=pc, zorder=3,
                   label=f"Lakes_cci lake, too few in-situ dates ({len(few)})")
    for i, (s, d) in enumerate(L.groupby("source")):
        ax.scatter(d.lon, d.lat, s=8 + 4 * np.sqrt(d.n_insitu), color=PL.CAT[i % len(PL.CAT)], alpha=0.8, edgecolor=PL.SURFACE,
                   linewidth=0.6, transform=pc, label=f"{s} ({len(d)} lakes)", zorder=4)
    if "tier" in L:
        lo = L[L.tier != "standard"]
        ax.scatter(lo.lon, lo.lat, s=60, facecolor="none", edgecolor=PL.INK, linewidth=1.2, transform=pc, zorder=5,
                   label=f"tropical low-data tier ({len(lo)})")
    ax.legend(loc="lower left", fontsize=7.5, ncol=2, framealpha=0.9)
    ax.set_title(f"Lakes with in-situ chlorophyll-a; coloured: compared with Lakes_cci {label} (size: in-situ dates)")
    fig.savefig(out / "fig1_map.png", dpi=130, bbox_inches="tight"); plt.close(fig)


def fig_matchups(L, M, out, label):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.8))
    for s in ("DJF", "MAM", "JJA", "SON"):
        d = M[M.season == s]
        ax[0].scatter(d.insitu, d.cci, s=6, alpha=0.35, color=SCOL[s], edgecolor="none", label=f"{s} ({len(d)})")
    lim = [0.1, 300]
    for a in ax:
        a.plot(lim, lim, color=PL.INK2, lw=1); a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(lim); a.set_ylim(lim)
        a.set_xlabel("in-situ chl-a (µg/L)"); a.set_ylabel(f"Lakes_cci {label} chl-a (µg/L)")
    e = M.logratio
    ax[0].set_title(f"Same-day matchups: {len(M)} pairs in {M.hylak_id.nunique()} lakes")
    ax[0].text(0.02, 0.97, f"median bias ×{10 ** e.median():.2f}\nr (log) {np.corrcoef(lg(M.insitu), lg(M.cci))[0, 1]:.2f}",
               transform=ax[0].transAxes, va="top", color=PL.INK)
    ax[0].legend(loc="lower right", markerscale=2, title="season (NH)")
    d = L.dropna(subset=["bias"])
    g = M.groupby("hylak_id")[["insitu", "cci"]].median()
    ax[1].scatter(g.insitu, g.cci, s=22, color=PL.CAT[0], edgecolor=PL.SURFACE, linewidth=0.6)
    ax[1].set_title(f"Per-lake median of the matchups ({len(g)} lakes)")
    rr = np.corrcoef(lg(g.insitu), lg(g.cci))[0, 1]
    ax[1].text(0.02, 0.97, f"between-lake r (log) {rr:.2f}\nmedian lake bias ×{10 ** d.bias.median():.2f}",
               transform=ax[1].transAxes, va="top", color=PL.INK)
    fig.tight_layout(); fig.savefig(out / "fig2_matchups.png", dpi=130); plt.close(fig)


def fig_month_bias(M, out, label):
    # per lake and month first (so lakes with many samples do not dominate), then across lakes
    g = M.groupby(["hylak_id", "month_nh"]).logratio.median().reset_index()
    g["dev"] = g.logratio - g.groupby("hylak_id").logratio.transform("median")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
    for a, k, t in ((ax[0], "logratio", "CCI / in-situ, same day (log10)"),
                    (ax[1], "dev", "same, minus each lake's own median bias")):
        q = g.groupby("month_nh")[k].quantile([0.25, 0.5, 0.75]).unstack()
        n = g.groupby("month_nh").size()
        a.fill_between(q.index, q[0.25], q[0.75], color=PL.CAT[0], alpha=0.18, lw=0)
        a.plot(q.index, q[0.5], "o-", color=PL.CAT[0], ms=5)
        a.axhline(0, color=PL.INK2, lw=1)
        a.set_title(t); a.set_xticks(range(1, 13))
        a.set_xticklabels([f"{MON[m - 1]}\n{n.get(m, 0)}" for m in range(1, 13)])
    ax[0].set_ylabel("log10 ratio (0 = no bias; 0.3 = ×2)")
    fig.suptitle(f"Does the Lakes_cci {label} bias change through the year? Median and IQR over lakes "
                 "(numbers: lakes per month; southern lakes shifted 6 months)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig3_bias_by_month.png", dpi=130); plt.close(fig)
    g.to_csv(out / "monthly_bias.csv", index=False, float_format="%.4g")


def fig_cycles(L, cyc, out, label, n=24, name="fig4_cycles", title="lakes with most in-situ dates"):
    top = L.dropna(subset=["r_season"]).sort_values("n_insitu", ascending=False).head(n)
    if top.empty:
        return
    nc = 6; nr = int(np.ceil(len(top) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(2.3 * nc, 2.0 * nr), sharex=True)
    for a, (_, r) in zip(axs.flat, top.iterrows()):
        co, cc, mo, mc = cyc[r.hylak_id]
        nm = lambda s: s / s.mean()
        a.plot(co.index, nm(co), "o-", color=PL.INK, ms=3, lw=1.6, label="in-situ")
        a.plot(cc.index, nm(cc), "-", color=PL.CAT[0], lw=1.8, label=f"CCI {label}, all days")
        if mc is not None:
            a.plot(mc.index, nm(mc), "--", color=PL.CAT[1], lw=1.2, label="CCI, matchup days")
        nmm = str(r.lake)[:22]
        a.set_title(f"{nmm} {r.lat:.0f}°  r={r.r_season:.2f}" + ("  (low data)" if r.get("tier") == "tropical low-data" else ""), fontsize=7.5)
        a.set_xticks([1, 4, 7, 10]); a.set_xticklabels(["J", "A", "J", "O"]); a.tick_params(labelsize=7)
    for a in list(axs.flat)[len(top):]:
        a.axis("off")
    axs.flat[0].legend(fontsize=6.5, loc="upper left")
    fig.suptitle(f"Monthly cycles relative to each series' own mean ({title})", x=0.01, ha="left")
    fig.tight_layout(); fig.savefig(out / f"{name}.png", dpi=130); plt.close(fig)


def fig_peaks(L, out, label):
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    bins = np.arange(-6.5, 6.6, 1)
    for i, (k, t) in enumerate((("dpeak", "all days"), ("dpeak_m", "matchup days only"))):
        v = L[k].dropna()
        ax[0].hist(v, bins=bins + (i - 0.5) * 0.35, rwidth=0.33, color=PL.CAT[i], label=f"{t} (n={len(v)}, "
                   f"|Δ|≤1: {(v.abs() <= 1).mean():.0%})")
    ax[0].set_title("Peak month, CCI minus in-situ"); ax[0].set_xlabel("months (+: CCI peaks later)"); ax[0].legend(fontsize=7)
    for i, (k, t) in enumerate((("r_season", "all days"), ("r_season_m", "matchup days only"))):
        v = L[k].dropna()
        ax[1].hist(v, bins=np.linspace(-1, 1, 11) + (i - 0.5) * 0.06, rwidth=0.3, color=PL.CAT[i],
                   label=f"{t}: median {v.median():.2f}")
    ax[1].set_title("Seasonal correlation (log10 monthly cycles)"); ax[1].set_xlabel("r"); ax[1].legend(fontsize=7)
    d = L.dropna(subset=["peak_insitu", "peak_cci"])
    j = np.random.default_rng(1).uniform(-0.25, 0.25, (2, len(d)))
    pi = np.where(d.lat < 0, (d.peak_insitu + 5) % 12 + 1, d.peak_insitu)
    pc = np.where(d.lat < 0, (d.peak_cci + 5) % 12 + 1, d.peak_cci)
    ax[2].scatter(pi + j[0], pc + j[1], s=18, color=PL.CAT[0], alpha=0.7, edgecolor=PL.SURFACE, linewidth=0.5)
    ax[2].plot([1, 12], [1, 12], color=PL.INK2, lw=1)
    ax[2].set_xticks(range(1, 13)); ax[2].set_xticklabels(list(MON)); ax[2].set_yticks(range(1, 13))
    ax[2].set_yticklabels(list(MON)); ax[2].set_xlabel("in-situ peak month (NH-aligned)"); ax[2].set_ylabel("CCI peak month")
    ax[2].set_title("Peak months, all days")
    fig.tight_layout(); fig.savefig(out / "fig5_peaks.png", dpi=130); plt.close(fig)


def fig_level(L, out, label):
    d = L.dropna(subset=["bias"]).copy()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].scatter(d.insitu_mean, 10 ** d.bias, s=22, color=PL.CAT[0], edgecolor=PL.SURFACE, linewidth=0.6)
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].axhline(1, color=PL.INK2, lw=1)
    ax[0].set_xlabel("in-situ mean chl-a (µg/L)"); ax[0].set_ylabel("CCI / in-situ (median of matchups)")
    ax[0].set_title("Bias against lake chlorophyll level")
    ax[1].scatter(d.area_km2, 10 ** d.bias, s=22, color=PL.CAT[2], edgecolor=PL.SURFACE, linewidth=0.6)
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].axhline(1, color=PL.INK2, lw=1)
    ax[1].set_xlabel("lake area (km²)"); ax[1].set_title("Bias against lake size")
    fig.tight_layout(); fig.savefig(out / "fig6_level.png", dpi=130); plt.close(fig)


def main(ins_csv, cci_spec, out, *sites_csv):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    ins = pd.read_csv(ins_csv, parse_dates=["date"])
    cci = CCI(cci_spec)
    L, M, cyc = analyse(ins, cci, out)
    print(f"{len(L)} lakes ({L.groupby('source').size().to_dict()}), {len(M)} matchups in {M.hylak_id.nunique()} lakes")
    cols = ["n_match", "bias", "r_match", "dpeak", "r_season", "dpeak_m", "r_season_m"]
    print(L[[c for c in cols if c in L]].describe().round(2).to_string())
    found = pd.concat([pd.read_csv(f) for f in sites_csv], ignore_index=True) if sites_csv else None
    if found is not None:
        found.to_csv(out / "sites_found.csv", index=False)
    ins.groupby(["source", "hylak_id"]).date.nunique().rename("n_dates").reset_index().to_csv(out / "insitu_lakes_all.csv", index=False)
    fig_map(L, out, cci.label, ins, found); fig_matchups(L, M, out, cci.label); fig_month_bias(M, out, cci.label)
    fig_cycles(L, cyc, out, cci.label)
    fig_cycles(L[L.lat.abs() < 23.5], cyc, out, cci.label, name="fig7_tropics", title="tropical lakes, calendar months")
    fig_peaks(L, out, cci.label); fig_level(L, out, cci.label)


if __name__ == "__main__":
    main(*sys.argv[1:])
