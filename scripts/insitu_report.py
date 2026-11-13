#!/usr/bin/env python3
"""HTML report: Lakes_cci chlorophyll-a seasonality against in-situ data (insitu_timing.py outputs).

    python3 analysis/lakes_cci/insitu_report.py <insitu_timing out dir> <out.html> [--standalone]

--standalone writes a full HTML document (GitHub Pages); without it the page is a fragment for the Artifact viewer.
INSITU_SOURCES: the sources table (default calibration/obs_lakes_cci/insitu_sources.csv, or data/sources.csv next to
the scripts' folder); REPORT_REPO: the code repository linked from the page.

Figures are embedded as PNG data URIs; the lake explorer draws the per-lake monthly cycles (cycles.json) in SVG.
"""
import base64
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def img(p):
    return "data:image/png;base64," + base64.b64encode(Path(p).read_bytes()).decode()


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
# table of every in-situ database used or checked: GLM4-AED3 layout, or data/sources.csv in the standalone repository
SOURCES = Path(os.environ.get("INSITU_SOURCES") or next(
    (f for f in (ROOT / "calibration/obs_lakes_cci/insitu_sources.csv", HERE.parent / "data/sources.csv") if f.exists()),
    HERE.parent / "data/sources.csv"))
REPO = os.environ.get("REPORT_REPO", "https://github.com/danielmerbet/lakes-cci-insitu")


def esc(s):
    s = "" if pd.isna(s) else str(s)
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def sources_table(d, L):
    """Every in-situ database used or checked (calibration/obs_lakes_cci/insitu_sources.csv), with lakes compared,
    Lakes_cci lakes with any data, and sites found outside Lakes_cci lakes."""
    S = pd.read_csv(SOURCES)
    A = pd.read_csv(d / "insitu_lakes_all.csv") if (d / "insitu_lakes_all.csv").exists() else pd.DataFrame(columns=["source", "hylak_id"])
    F = pd.read_csv(d / "sites_found.csv") if (d / "sites_found.csv").exists() else pd.DataFrame(columns=["source", "hylak_id"])
    used = L.source.value_counts(); anyl = A.groupby("source").hylak_id.nunique()
    out = F[F.hylak_id.isna()].groupby("source").size()
    rows = []
    for st, lab in (("used", "Used"), ("inventory only", "On the map only"), ("not yet", "Not yet available"),
                    ("checked", "Checked, nothing usable"), ("excluded", "Excluded")):
        x = S[S.status == st]
        if x.empty:
            continue
        rows.append(f'<tr class="grp"><td colspan="6">{lab}</td></tr>')
        for _, r in x.iterrows():
            nm = esc(r.source)
            link = f'<a href="{esc(r.url)}">{nm}</a>' if isinstance(r.url, str) and r.url.startswith("http") else nm
            n = lambda s, k: "" if st not in ("used", "inventory only") else f"{int(s.get(k, 0)):,}"
            note = "; ".join(esc(v) for v in (r.access, r.independence if st != "used" else "", r.notes) if isinstance(v, str) and v)
            rows.append(f"<tr><td>{link}<br><span class='meta'>{esc(r.provider)} · {esc(r.region)}</span></td>"
                        f"<td>{esc(r.method)}</td><td>{esc(r.licence)}</td><td>{n(used, r.source)}</td><td>{n(anyl, r.source)}</td>"
                        f"<td>{n(out, r.source) if r.source in set(F.source) else ('searched inside Lakes_cci lakes only' if st == 'used' else '')}</td></tr>" + (f"<tr class='note'><td colspan='6'>{note}</td></tr>" if note else ""))
    return "".join(rows)


def numbers(d):
    L = pd.read_csv(d / "lakes.csv"); M = pd.read_csv(d / "matchups.csv"); B = pd.read_csv(d / "monthly_bias.csv")
    lg = lambda x: np.log10(np.maximum(np.asarray(x, float), 0.05))
    g = M.groupby("hylak_id")[["insitu", "cci"]].median()
    S = L.dropna(subset=["dpeak"]); b = L.dropna(subset=["bias"])
    n = {"lakes": len(L), "lakes_m": int((L.n_match > 0).sum()), "pairs": len(M), "sources": L.source.value_counts().to_dict(),
         "r_pairs": np.corrcoef(lg(M.insitu), lg(M.cci))[0, 1], "bias_pairs": 10 ** np.median(lg(M.cci) - lg(M.insitu)),
         "r_lakes": np.corrcoef(lg(g.insitu), lg(g.cci))[0, 1], "bias_lakes": 10 ** b.bias.median(),
         "n_season": len(S), "peak1": (S.dpeak.abs() <= 1).mean(), "later": (S.dpeak > 1).mean(), "earlier": (S.dpeak < -1).mean(),
         "r_season": S.r_season.median(), "r_season_m": L.r_season_m.median(), "n_season_m": int(L.r_season_m.notna().sum()),
         "lat0": L.lat.min(), "lat1": L.lat.max()}
    w = B[B.month_nh.isin([12, 1, 2, 3, 4, 5])].dev.median(); s = B[B.month_nh.isin([6, 7, 8, 9])].dev.median()
    n["ws_vs_summer"] = 10 ** (w - s)
    n["troph"] = [(lab, len(x), 10 ** x.bias.median()) for lab, x in
                  (("&lt; 5", b[b.insitu_mean < 5]), ("5–10", b[(b.insitu_mean >= 5) & (b.insitu_mean < 10)]),
                   ("10–30", b[(b.insitu_mean >= 10) & (b.insitu_mean < 30)]), ("≥ 30", b[b.insitu_mean >= 30]))]
    n["size"] = [(lab, len(x), 10 ** x.bias.median(), 10 ** x.bias.abs().median()) for lab, x in
                 (("&lt; 20", b[b.area_km2 < 20]), ("20–100", b[(b.area_km2 >= 20) & (b.area_km2 < 100)]),
                  ("100–1,000", b[(b.area_km2 >= 100) & (b.area_km2 < 1000)]), ("≥ 1,000", b[b.area_km2 >= 1000]))]
    P = L.dropna(subset=["peak_insitu"]); a = P.lat.abs()
    zones = (("≥ 55°", P[a >= 55]), ("40–55°", P[(a >= 40) & (a < 55)]), ("23.5–40°", P[(a >= 23.5) & (a < 40)]),
             ("&lt; 23.5° (tropics)", P[a < 23.5]))
    n["lat"] = [(lab, len(x), x.peak_insitu.between(7, 10).mean(), x.peak_cci.between(7, 10).mean(),
                 x.peak_insitu.between(4, 6).mean(), x.peak_cci.between(4, 6).mean(), x.dpeak.median(), x.dpeak_m.median(),
                 x.r_season.median(), x.r_season_m.median()) for lab, x in zones]
    n["north"] = [(x.peak_insitu.between(7, 10).mean(), x.peak_cci.between(7, 10).mean()) for x in (P[P.lat >= 40],)][0]
    B2 = P[a >= 55]
    n["boreal"] = (len(B2), B2.dpeak.median(), B2.dpeak_m.median(), B2.r_season_m.median(), (B2.dpeak <= -2).mean())
    T = L[L.lat.abs() < 23.5]; Ts = T.dropna(subset=["dpeak"])
    n["trop"] = {"n": len(T), "low": int((T.tier == "tropical low-data").sum()), "n_season": len(Ts),
                 "peak1": (Ts.dpeak.abs() <= 1).mean(), "r": Ts.r_season.median(), "bias": 10 ** T.bias.median(),
                 "dates": T.n_insitu.median(), "src": T.source.value_counts().to_dict()}
    m8 = L[L.months >= 8]
    n["m8"] = (len(m8), m8.r_season.median(), (m8.dpeak.abs() <= 1).mean())
    return n


CSS = """
:root{
  /* Layout: one reading column (68ch) with figure plates that break out to 1120px; findings strip under the title. */
  --bg:#f5f7f6; --surface:#ffffff; --plate:#fcfcfb; --ink:#10282c; --ink2:#3e5459; --muted:#6c7f83; --rule:#d5dedc;
  --accent:#2a78d6; --insitu:#10282c; --cci:#2a78d6; --ccim:#eb6834; --chl:#1b8f63; --chip:#e6efec;
  --display:"Bricolage Grotesque",ui-sans-serif,system-ui,sans-serif; --body:"Source Serif 4",Georgia,serif;
  --mono:"JetBrains Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --bg:#0e1719; --surface:#142124; --plate:#fcfcfb; --ink:#e4eeec; --ink2:#b4c6c3; --muted:#8aa09d; --rule:#2a3b3e;
  --accent:#6da7ec; --insitu:#e4eeec; --cci:#6da7ec; --ccim:#f4935f; --chl:#4cc79a; --chip:#1d2f32; color-scheme:dark}}
:root[data-theme="dark"]{
  --bg:#0e1719; --surface:#142124; --plate:#fcfcfb; --ink:#e4eeec; --ink2:#b4c6c3; --muted:#8aa09d; --rule:#2a3b3e;
  --accent:#6da7ec; --insitu:#e4eeec; --cci:#6da7ec; --ccim:#f4935f; --chl:#4cc79a; --chip:#1d2f32; color-scheme:dark}
*{box-sizing:border-box}
body{background:var(--bg);color:var(--ink);font:17px/1.62 var(--body);margin:0}
.wrap{max-width:1120px;margin:0 auto;padding-block:48px 80px;padding-left:max(16px,env(safe-area-inset-left,0px));
  padding-right:max(16px,env(safe-area-inset-right,0px));display:flex;flex-direction:column;gap:28px}
.col{max-width:68ch}
h1,h2,h3{font-family:var(--display);text-wrap:balance;letter-spacing:-0.01em;margin:0}
h1{font-size:clamp(2rem,4.6vw,3.1rem);line-height:1.08;font-weight:700}
h2{font-size:1.55rem;line-height:1.2;font-weight:650;padding-top:18px;border-top:1px solid var(--rule)}
h3{font-size:1.08rem;font-weight:650}
p{margin:0}
.stack{display:flex;flex-direction:column;gap:14px}
.eyebrow{font:600 .74rem/1 var(--mono);letter-spacing:.12em;text-transform:uppercase;color:var(--chl)}
.lede{font-size:1.18rem;color:var(--ink2)}
.meta{font:.82rem/1.5 var(--mono);color:var(--muted)}
.findings{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:12px}
.find{background:var(--surface);border:1px solid var(--rule);border-radius:10px;padding:16px;display:flex;flex-direction:column;gap:6px;min-width:0}
.find b{font:700 1.9rem/1 var(--display);color:var(--ink);font-variant-numeric:tabular-nums}
.find span{font-size:.92rem;color:var(--ink2);line-height:1.45}
.plate{background:var(--plate);border:1px solid var(--rule);border-radius:10px;padding:10px;overflow-x:auto}
.plate img{display:block;width:100%;min-width:640px;height:auto}
figure{margin:0;display:flex;flex-direction:column;gap:10px}
figcaption{font-size:.93rem;color:var(--ink2);max-width:80ch}
figcaption b{color:var(--ink)}
table{border-collapse:collapse;font:.88rem/1.4 var(--mono);font-variant-numeric:tabular-nums}
th,td{padding:6px 12px;text-align:right;border-bottom:1px solid var(--rule)}
th:first-child,td:first-child{text-align:left}
th{color:var(--muted);font-weight:600}
.srcs td{text-align:left;vertical-align:top} .srcs td:nth-child(n+4){text-align:right}
.srcs .grp td{font:600 .78rem var(--mono);letter-spacing:.08em;text-transform:uppercase;color:var(--chl);padding-top:16px}
.srcs .note td{font-size:.8rem;color:var(--muted);border-bottom:1px solid var(--rule);padding-top:0}
.srcs tr:not(.note):not(.grp) td{border-bottom:none}
.tables{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:20px}
.tables>div{overflow-x:auto;min-width:0}
ul{margin:0;padding-left:1.2em;display:flex;flex-direction:column;gap:8px}
a{color:var(--accent)}
code{font:.86em var(--mono);background:var(--chip);padding:1px 5px;border-radius:4px}
.ask{background:var(--surface);border:1px solid var(--rule);border-left:4px solid var(--chl);border-radius:6px;padding:16px 18px}
.explorer{background:var(--surface);border:1px solid var(--rule);border-radius:10px;padding:16px;display:grid;
  grid-template-columns:minmax(0,1fr) 260px;gap:18px}
@media (max-width:760px){.explorer{grid-template-columns:1fr}}
.explorer select,.explorer input{font:.9rem var(--mono);padding:7px 9px;border:1px solid var(--rule);border-radius:6px;
  background:var(--bg);color:var(--ink);width:100%}
.controls{display:flex;gap:8px;flex-wrap:wrap}
.controls>*{flex:1 1 200px}
select:focus-visible,input:focus-visible{outline:2px solid var(--accent);outline-offset:1px}
.stats{display:grid;grid-template-columns:auto auto;gap:4px 12px;font:.85rem/1.45 var(--mono);align-content:start}
.stats dt{color:var(--muted)} .stats dd{margin:0;text-align:right;font-variant-numeric:tabular-nums}
.legend{display:flex;gap:16px;flex-wrap:wrap;font:.8rem var(--mono);color:var(--ink2)}
.legend i{display:inline-block;width:18px;height:0;border-top:2px solid;vertical-align:middle;margin-right:6px}
svg text{font-family:var(--mono)}
"""

JS = """
const D = __DATA__;
const MON = "JFMAMJJASOND";
const sel = document.getElementById('lake'), q = document.getElementById('q'), svg = document.getElementById('chart');
function fill(){
  const f = q.value.toLowerCase();
  const cur = sel.value;
  sel.innerHTML = '';
  D.filter(d => !f || d.lake.toLowerCase().includes(f)).forEach(d => {
    const o = document.createElement('option'); o.value = d.hylak_id;
    o.textContent = `${d.lake} · ${d.lat.toFixed(1)}° · ${d.n_insitu} in-situ dates`; sel.appendChild(o);
  });
  if ([...sel.options].some(o => o.value == cur)) sel.value = cur;
  draw();
}
function norm(a){ if(!a) return null; const v = a.filter(x => x != null); if(!v.length) return null;
  const m = v.reduce((s,x)=>s+x,0)/v.length; return a.map(x => x == null ? null : x/m); }
function draw(){
  const d = D.find(x => x.hylak_id == sel.value); if(!d){ svg.innerHTML=''; return; }
  const W=640,H=300,l=44,r=14,t=14,b=30;
  const S = [[norm(d.insitu),'var(--insitu)','',3.5],[norm(d.cci),'var(--cci)','',0],[norm(d.cci_m),'var(--ccim)','6 4',0]];
  let ymax = 0; S.forEach(s => s[0] && s[0].forEach(v => { if(v!=null) ymax=Math.max(ymax,v); }));
  ymax = Math.max(1.5, Math.ceil(ymax*2)/2);
  const x = m => l + (m-0.5)*(W-l-r)/12, y = v => t + (1 - v/ymax)*(H-t-b);
  let h = '';
  for(let v=0; v<=ymax+1e-9; v+=0.5){ h += `<line x1="${l}" x2="${W-r}" y1="${y(v)}" y2="${y(v)}" stroke="var(--rule)" stroke-width="${v==1?1.2:0.6}"/>`;
    h += `<text x="${l-8}" y="${y(v)+4}" text-anchor="end" font-size="11" fill="var(--muted)">${v.toFixed(1)}</text>`; }
  for(let m=1;m<=12;m++) h += `<text x="${x(m)}" y="${H-10}" text-anchor="middle" font-size="11" fill="var(--muted)">${MON[m-1]}</text>`;
  S.forEach(([a,c,dash,dot]) => { if(!a) return; let p='', seg=false;
    a.forEach((v,i)=>{ if(v==null){seg=false;return;} p += (seg?'L':'M')+x(i+1)+' '+y(v)+' '; seg=true; });
    h += `<path d="${p}" fill="none" stroke="${c}" stroke-width="2.2" stroke-dasharray="${dash}" stroke-linejoin="round"/>`;
    if(dot) a.forEach((v,i)=>{ if(v!=null) h += `<circle cx="${x(i+1)}" cy="${y(v)}" r="${dot}" fill="${c}" stroke="var(--surface)" stroke-width="1.5"><title>${MON[i]}: ${(v).toFixed(2)} × mean</title></circle>`; });
  });
  svg.innerHTML = h;
  const f = (v,k=2) => v==null ? '–' : (+v).toFixed(k);
  const mn = v => v==null ? '–' : MON[v-1] + ' (' + ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'][v-1] + ')';
  document.getElementById('stats').innerHTML = [
    ['source', d.source], ['HydroLAKES id', d.hylak_id], ['latitude', f(d.lat,2)+'°'], ['area', d.area_km2==null?'–':Math.round(d.area_km2).toLocaleString()+' km²'],
    ['in-situ dates', d.n_insitu], ['CCI days', d.n_cci], ['same-day pairs', d.n_match],
    ['in-situ mean', f(d.insitu_mean,1)+' µg/L'], ['CCI mean', f(d.cci_mean,1)+' µg/L'],
    ['CCI / in-situ (pairs)', d.bias==null?'–':'×'+Math.pow(10,d.bias).toFixed(2)],
    ['peak in-situ', mn(d.peak_insitu)], ['peak CCI', mn(d.peak_cci)], ['seasonal r', f(d.r_season)],
  ].map(([k,v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
}
q.addEventListener('input', fill); sel.addEventListener('change', draw);
fill();
"""


def html(n, d):
    J = json.loads((d / "cycles.json").read_text())
    J = sorted(J, key=lambda r: -r["n_insitu"])
    src = ", ".join(f"{v} {k}" for k, v in n["sources"].items())
    tro = "".join(f"<tr><td>{a}</td><td>{k}</td><td>×{v:.2f}</td></tr>" for a, k, v in n["troph"])
    siz = "".join(f"<tr><td>{a}</td><td>{k}</td><td>×{v:.2f}</td><td>×{e:.2f}</td></tr>" for a, k, v, e in n["size"])
    f = lambda v, k=0: "–" if v != v else f"{v:+.{k}f}" if k == 0 else f"{v:.{k}f}"
    lat = "".join(f"<tr><td>{a}</td><td>{k}</td><td>{pi:.0%}</td><td>{pc:.0%}</td><td>{si:.0%}</td><td>{sc:.0%}</td>"
                  f"<td>{f(dp)}</td><td>{f(dm)}</td><td>{f(rs, 2)}</td><td>{f(rm, 2)}</td></tr>"
                  for a, k, pi, pc, si, sc, dp, dm, rs, rm in n["lat"])
    bn, bdp, bdm, brm, bear = n["boreal"]
    t5 = n["troph"][0][2]
    F = {k: img(d / f"{k}.png") for k in ("fig1_map", "fig2_matchups", "fig3_bias_by_month", "fig4_cycles", "fig5_peaks", "fig6_level", "fig7_tropics")}
    m8n, m8r, m8p = n["m8"]
    LL = pd.read_csv(d / "lakes.csv")
    srct = sources_table(d, LL)
    oz = LL[(LL.lat < -30) & LL.lon.between(110, 180)]
    mo = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    ozt = "; ".join(f"{str(r.lake).split(' (')[0] if not str(r.lake).startswith('(') else 'lake ' + str(r.hylak_id)} "
                    f"({r.source}, {int(r.n_insitu)} dates): field peak {mo[int(r.peak_insitu) - 1] if r.peak_insitu == r.peak_insitu else '–'}, "
                    f"CCI {mo[int(r.peak_cci) - 1] if r.peak_cci == r.peak_cci else '–'}, r {r.r_season:.2f}"
                    for _, r in oz.sort_values("lat").iterrows())
    return f"""<title>Lakes_cci Chlorophyll Seasonality Check</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,500..700&family=JetBrains+Mono:wght@400;600&family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&display=swap">
<style>{CSS}</style>
<div class="wrap">
<header class="stack col">
  <span class="eyebrow">ESA Lakes_cci · chlorophyll-a · in-situ check</span>
  <h1>Does Lakes_cci get the chlorophyll season right?</h1>
  <p class="lede">We compared Lakes_cci chlorophyll-a with {n['pairs']:,} same-day field samples and with the monthly
  cycles of {n['lakes']} lakes that have both, from North America, Europe and Mexico. In temperate lakes the seasonal
  timing holds up. In boreal lakes Lakes_cci peaks too early, just after ice-off, while the field data peak in July–August.
  The level is compressed: clear lakes read too high, very green lakes about right.</p>
  <p class="meta">Lakes_cci v2.1 via the per-lake tables of Shahvaran et al. (ESSD preprint essd-2026-528) · in-situ:
  national and regional monitoring (US Water Quality Portal, EEA Waterbase, Finnish SYKE, UNEP GEMStat, EPA Victoria,
  Environment Canterbury), Lewis et al. (2023)
  EDI synthesis, local monitoring · none used to build or validate Lakes_cci · 2002–2022 · prepared November 2026 in the
  GLM4-AED3 global lake modelling work · v3.0 comparison to follow</p>
</header>

<section class="findings" aria-label="Main numbers">
  <div class="find"><b>{n['r_lakes']:.2f}</b><span>correlation between lakes of the per-lake median chl-a (log), {n['lakes_m']} lakes</span></div>
  <div class="find"><b>×{n['bias_lakes']:.2f}</b><span>median CCI / in-situ ratio per lake, same-day pairs</span></div>
  <div class="find"><b>{n['peak1']:.0%}</b><span>of {n['n_season']} lakes where CCI peaks within one month of the field data</span></div>
  <div class="find"><b>{n['north'][0]:.0%}</b><span>of field series north of 40° N peak in Jul–Oct (CCI: {n['north'][1]:.0%})</span></div>
  <div class="find"><b>{bdm:+.0f} mo</b><span>CCI peak minus field peak on the same days, {bn} lakes north of 55°, median</span></div>
</section>

<section class="stack col">
  <h2>Why we looked</h2>
  <p>We run the GLM-AED lake model for about 40,000 lakes worldwide (ISIMIP) and validate surface chlorophyll-a against
  Lakes_cci. In lakes north of 40° N the model peaks in April–June, while Lakes_cci peaks in July–October. Tests with the
  model's growth, nutrient, grazing, nitrogen-fixation and sediment-phosphorus settings did not move the model peak.
  Before changing the model structure we needed to know whether the satellite timing is right. This page is that check,
  shared with the Lakes_cci team because the same comparison says something about the product.</p>
</section>

<section class="stack">
  <h2 class="col">Where the data are</h2>
  <figure>
    <div class="plate"><img src="{F['fig1_map']}" alt="World map of every lake with in-situ chlorophyll-a found: lakes compared with Lakes_cci coloured by source, Lakes_cci lakes with too few samples as open circles, other sites in grey."></div>
    <figcaption>Coloured: <b>{n['lakes']} lakes</b> compared with Lakes_cci ({src}), {n['lat0']:.0f}–{n['lat1']:.0f}°.
    Open circles: Lakes_cci lakes with field data, but too few dates. Grey: other sites with field chlorophyll-a that are
    not Lakes_cci lakes, from the sources that list all their sites (GEMStat, Ceará, EPA Victoria, Canterbury); EEA,
    Finland and the US portal were searched only inside Lakes_cci lakes. A lake is used when it has at
    least 20 in-situ sampling dates spread over at least 5 calendar months, and at least 100 Lakes_cci days with 10 % or
    more of the lake covered; tropical lakes with fewer data (at least 10 dates over 4 months) are kept as a separate
    low-data tier ({n['trop']['low']} lakes). A lake found in several archives is taken from the one with most dates. Canada, Russia and
    China, where most Lakes_cci lakes are, are still missing.</figcaption>
  </figure>
</section>

<section class="stack">
  <h2 class="col">Same-day pairs</h2>
  <figure>
    <div class="plate"><img src="{F['fig2_matchups']}" alt="Scatter plots of Lakes_cci against in-situ chlorophyll-a, all pairs on the left and per-lake medians on the right."></div>
    <figcaption><b>Left:</b> {n['pairs']:,} pairs in {n['lakes_m']} lakes (r = {n['r_pairs']:.2f} in log space,
    median ratio ×{n['bias_pairs']:.2f}). <b>Right:</b> the per-lake medians line up well (r = {n['r_lakes']:.2f}),
    so CCI ranks lakes correctly. Single pairs scatter widely: a lake-mean satellite value against one point sample, often
    taken near shore.</figcaption>
  </figure>
  <div class="tables col">
    <div><h3>Bias by lake chlorophyll level</h3><table><thead><tr><th>in-situ mean chl-a (µg/L)</th><th>lakes</th><th>CCI / in-situ</th></tr></thead><tbody>{tro}</tbody></table></div>
    <div><h3>Bias by lake size</h3><table><thead><tr><th>area (km²)</th><th>lakes</th><th>median</th><th>typical error</th></tr></thead><tbody>{siz}</tbody></table></div>
  </div>
  <figure>
    <div class="plate"><img src="{F['fig6_level']}" alt="Per-lake CCI over in-situ ratio against in-situ mean chlorophyll-a and against lake area."></div>
    <figcaption>The ratio falls as the lake gets greener (higher mean chl-a): lakes under 5 µg/L read ×{t5:.1f} too high, lakes above 30 µg/L
    about right. Lakes under 20 km² have the largest errors in both directions (typical factor in the table: median of
    the absolute log ratio).</figcaption>
  </figure>
</section>

<section class="stack">
  <h2 class="col">Does the bias change through the year?</h2>
  <figure>
    <div class="plate"><img src="{F['fig3_bias_by_month']}" alt="Median CCI over in-situ ratio by calendar month, raw and relative to each lake's own bias."></div>
    <figcaption>Each lake's pairs are first summarised per month, so lakes with many samples do not dominate.
    <b>Right panel:</b> after removing each lake's own bias, CCI reads about ×{n['ws_vs_summer']:.2f} higher in
    December–May than in June–September. A seasonal bias of this kind flattens a summer peak; it does not create one.
    Winter months rest on fewer lakes (numbers under the axis).</figcaption>
  </figure>
</section>

<section class="stack">
  <h2 class="col">Seasonal cycles</h2>
  <figure>
    <div class="plate"><img src="{F['fig4_cycles']}" alt="Small multiples of monthly chlorophyll-a cycles for the 24 lakes with most in-situ dates."></div>
    <figcaption>Monthly cycles, each relative to its own mean. <b>Black:</b> in-situ. <b>Blue:</b> Lakes_cci, all
    valid days. <b>Orange dashed:</b> Lakes_cci on the same days as the field samples. Where blue and orange agree, the
    difference from black comes from the retrieval, not from when the satellite can see the lake. For the largest
    lakes (Michigan, Superior) the field series pool a few sites that need not represent the whole lake, while CCI
    averages the whole surface.</figcaption>
  </figure>
  <figure>
    <div class="plate"><img src="{F['fig5_peaks']}" alt="Histograms of the peak-month difference and of the seasonal correlation, and a scatter of peak months."></div>
    <figcaption>Peak months agree within one month in {n['peak1']:.0%} of lakes; CCI peaks later in {n['later']:.0%} and
    earlier in {n['earlier']:.0%}. Median seasonal correlation {n['r_season']:.2f} (all days) and
    {n['r_season_m']:.2f} (matchup days, {n['n_season_m']} lakes). With 8 or more months in common ({m8n} lakes) the
    agreement is weaker (r {m8r:.2f}, peak within a month in {m8p:.0%}): winter and spring are where the two
    disagree most.</figcaption>
  </figure>
  <div class="col" style="overflow-x:auto"><h3>When is the peak?</h3><table><thead><tr><th>latitude</th><th>n</th><th>in-situ Jul–Oct</th><th>CCI Jul–Oct</th><th>in-situ Apr–Jun</th><th>CCI Apr–Jun</th><th>Δpeak (mo)</th><th>Δpeak, same days</th><th>r</th><th>r, same days</th></tr></thead><tbody>{lat}</tbody></table></div>
  <p class="col">North of 40° N, {n['north'][0]:.0%} of the field series peak in July–October. For our model this settles
  the question: the late-summer maximum is real, and the model's spring peak is a model problem. Lakes_cci follows the
  field data in temperate lakes, but north of 55° it peaks a median {-bdp:.1f} months early ({bear:.0%} of {bn} lakes
  two or more months early). The shift stays when only the field-sample days are used ({bdm:+.0f} months, r {brm:.2f}),
  so it comes from the retrieval, not from when the satellite sees the lake: CCI reads high in May, after ice-off, in
  humic Finnish lakes such as Päijänne, Pielinen and Keitele. In boreal lakes Lakes_cci should not be used for timing.</p>
  <p class="col"><b>Tropics.</b> {n['trop']['n']} lakes, nearly all Mexican reservoirs (UNEP GEMStat) plus Bosumtwi in
  Ghana, with a median of {n['trop']['dates']:.0f} field dates each (about three a year). CCI peaks within a month of the
  field data in {n['trop']['peak1']:.0%} of {n['trop']['n_season']} lakes, median seasonal r {n['trop']['r']:.2f}, and reads
  ×{n['trop']['bias']:.1f} the field values. The field cycles look ten times stronger than CCI's, but with one sample per
  month and year most of that is sampling noise; tropical timing stays an open question. No open in-situ series exists
  for the large African lakes (Victoria, Tanganyika, Malawi, Kivu) or Titicaca. Brazil holds 361 tropical Lakes_cci lakes;
  São Paulo (CETESB) and Minas Gerais (IGAM) monitoring would add reservoirs but cannot be fetched automatically.</p>
  <p class="col"><b>Australia and New Zealand.</b> {len(oz)} lakes so far: {ozt}. Lake Wellington follows the field
  cycle; the clear glacial lakes of Canterbury (0.5–1 µg/L) do not, as for other clear lakes. Lake Alexandrina, Taupō,
  Rotorua, Te Anau and Waikaremoana have monitoring that could not be fetched automatically (see Data sources).</p>
  <figure>
    <div class="plate"><img src="{F['fig7_tropics']}" alt="Small multiples of monthly chlorophyll-a cycles for the tropical lakes, in-situ against Lakes_cci."></div>
    <figcaption>Tropical lakes, same layout as above. In-situ points are often a single sample per month over the whole
    record, so their cycles are noisy.</figcaption>
  </figure>
</section>

<section class="stack">
  <h2 class="col">Explore a lake</h2>
  <div class="explorer">
    <div class="stack" style="min-width:0">
      <div class="controls">
        <label class="stack" style="gap:4px"><span class="meta">filter by name</span><input id="q" type="search" placeholder="e.g. Champlain" autocomplete="off"></label>
        <label class="stack" style="gap:4px"><span class="meta">lake</span><select id="lake"></select></label>
      </div>
      <svg id="chart" viewBox="0 0 640 300" role="img" aria-label="Monthly chlorophyll-a cycle of the selected lake, relative to each series' mean" style="width:100%;height:auto"></svg>
      <div class="legend"><span><i style="border-color:var(--insitu)"></i>in-situ</span><span><i style="border-color:var(--cci)"></i>CCI, all days</span><span><i style="border-color:var(--ccim);border-top-style:dashed"></i>CCI, matchup days</span></div>
    </div>
    <dl class="stats" id="stats"></dl>
  </div>
</section>

<section class="stack">
  <h2 class="col">Data sources</h2>
  <p class="col">Every in-situ database used or checked for this comparison. Lakes_cci chlorophyll-a was tuned and
  validated on LIMNADES, so it is excluded, and only statutory monitoring or long-term research data are used. The same list is a
  table in the <a href="{REPO}">code repository</a> (<code>{SOURCES.name}</code>).</p>
  <div style="overflow-x:auto"><table class="srcs"><thead><tr><th>source</th><th>chlorophyll-a method</th><th>licence</th>
  <th>lakes compared</th><th>Lakes_cci lakes with data</th><th>sites outside Lakes_cci</th></tr></thead><tbody>{srct}</tbody></table></div>
</section>

<section class="stack col">
  <h2>Method</h2>
  <ul>
    <li><b>Lakes_cci:</b> v2.1 merged L3S daily chlorophyll-a, summarised per HydroLAKES lake by Shahvaran et al.
    (Zenodo <a href="https://doi.org/10.5281/zenodo.20347115">10.5281/zenodo.20347115</a>, CC BY 4.0): daily lake mean,
    days with ≥ 10 % of the lake covered.</li>
    <li><b>In-situ:</b> laboratory chlorophyll-a (corrected and uncorrected for pheophytin; probe fluorescence left out),
    surface samples ≤ 3 m or epilimnion composites, daily mean over sites. US Water Quality Portal
    (<a href="https://www.waterqualitydata.us">waterqualitydata.us</a>): only stations inside the HydroLAKES polygon of
    a Lakes_cci lake. EEA Waterbase Water Quality ICM 2026 (WISE SoE, determinand EEA_3164-01-0; values below the
    quantification limit set to half of it). Finnish Environment Institute VESLA open API (analyte CP, ethanol or acetone
    extraction). UNEP GEMStat open subset GFQA v3 (Zenodo <a href="https://doi.org/10.5281/zenodo.18459694">10.5281/zenodo.18459694</a>,
    CC BY 4.0). Lewis et al. (2023), EDI <a href="https://doi.org/10.6073/pasta/2cd6628a942de2a8b12d2b19962712a0">edi.1530.1</a>, CC BY 4.0. Local monitoring (Erken).</li>
    <li><b>Independence:</b> Lakes_cci chlorophyll-a algorithms were tuned and validated against LIMNADES, a research
    archive of bio-optical campaign data (ATBD D2.2; PVP CCI-LAKES-0030 v1.3, section 6.3, which notes that "very few data
    sets from national monitoring programmes" were used). LIMNADES is not used here; all values are statutory monitoring
    or long-term research samples, none satellite- or model-derived.</li>
    <li><b>Pairs:</b> same calendar day. Bias = median of log10(CCI / in-situ) per lake.</li>
    <li><b>Seasonal cycles:</b> monthly means per year, then the mean over years per calendar month; in-situ months need
    data in at least 2 years. Compared on the months both have (at least 5). Southern lakes are shifted 6 months in the
    pooled monthly figures.</li>
    <li><b>Code and data:</b> <a href="{REPO}">{REPO.split("github.com/")[-1]}</a>: <code>insitu_wqp.py</code>,
    <code>insitu_more.py</code>, <code>insitu_chl.py</code>, <code>insitu_timing.py</code>, <code>insitu_report.py</code>
    (code GPL-3; data tables CC BY 4.0 with the source attributions).</li>
  </ul>
</section>

<section class="stack col">
  <h2>Limits of this check</h2>
  <ul>
    <li>The tropics rest on {n['trop']['n']} lakes sampled about three times a year; Africa has one lake, South America none.</li>
    <li>Canada, Russia and China hold most Lakes_cci lakes but no open, independent series here yet. Sweden (MVM) and
    Canada (DataStream) need an access key; NTL-LTER needs a log-in.</li>
    <li>EEA and most GEMStat records do not state the chlorophyll-a method; Finnish values are laboratory extractions.</li>
    <li>Field programmes sample mostly in summer and often near shore, while CCI reports a lake mean.</li>
    <li>Pheophytin-corrected and uncorrected values are pooled; the difference is usually 10–30 %.</li>
    <li>This uses v2.1. The same comparison with v3.0 is being extracted now and will be added here.</li>
  </ul>
</section>

<section class="ask col stack">
  <h3>Questions for the Lakes_cci team</h3>
  <ul>
    <li>In boreal humic lakes CCI peaks in May, two to four months before the field data, even on the same days. Is a
    post-ice-off bias known (CDOM, resuspended sediment, ice edges or the optical water type assignment)?</li>
    <li>Is the higher winter–spring reading (×{n['ws_vs_summer']:.2f}) known, for example from low sun angles, ice
    edges or the optical water type mix in those months?</li>
    <li>Is the compression of the range (clear lakes ×1.6–1.7 too high) expected from the OWT-blended algorithm, and
    does v3.0 change it?</li>
    <li>Are per-pixel uncertainties for chlorophyll-a and Kd planned? They would let users weight days instead of
    using a fixed coverage threshold.</li>
  </ul>
</section>
</div>
<script>{JS.replace("__DATA__", json.dumps(J))}</script>
"""


def standalone(page):
    """Full HTML document (GitHub Pages); the Artifact viewer adds this skeleton itself."""
    return ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
            + page.replace('\n<div class="wrap">', '\n</head>\n<body>\n<div class="wrap">', 1)
            + "</body>\n</html>\n")


def main(d, out, mode=""):
    d = Path(d)
    n = numbers(d)
    page = html(n, d)
    Path(out).write_text(standalone(page) if mode == "--standalone" else page)
    print(f"{out}: {Path(out).stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main(*sys.argv[1:4])
