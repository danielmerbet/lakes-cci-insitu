# Lakes_cci chlorophyll-a against independent in-situ data

Does the ESA Lakes_cci satellite chlorophyll-a get the level and the seasonal cycle right? This repository compares
Lakes_cci v2.1 with field chlorophyll-a from national and regional monitoring programmes that were **not** used to
build or validate Lakes_cci.

**Interactive report:** https://danielmerbet.github.io/lakes-cci-insitu/ (`index.html`, same file as `docs/index.html`; Pages can serve either the root or `/docs`)

## Main results (Lakes_cci v2.1, 2002–2022)

- **149 lakes**, 5,463 same-day pairs, from the US, Europe, Finland, Mexico, Australia, New Zealand and Ghana.
- **Level:** per lake, Lakes_cci reads ×1.49 the field value (median). The range is compressed: clear lakes
  (< 5 µg/L) read about ×2 too high, very green lakes (> 30 µg/L) about right.
- **Timing, temperate lakes:** the peak month agrees within one month in about half of the lakes. North of 40° N,
  80 % of the field series peak in July–October (Lakes_cci: 57 %).
- **Timing, boreal lakes (≥ 55°):** Lakes_cci peaks a median 2 months early (high in May, after ice-off, mostly
  in humic Finnish lakes), also on same-day pairs, so the shift comes from the retrieval. Do not use Lakes_cci
  for seasonal timing there.
- **Tropics:** 18 lakes, about three field samples a year; timing unresolved, Lakes_cci ×2.2 the field values.

## Independence from Lakes_cci

Lakes_cci chlorophyll-a algorithms were tuned and validated on LIMNADES, a research archive of bio-optical campaign
data (ATBD D2.2; PVP CCI-LAKES-0030 v1.3, section 6.3, which notes that "very few data sets from national monitoring
programmes" were used). LIMNADES is not used here. All values are statutory monitoring or long-term research
samples; none are satellite- or model-derived. Field probes (in-situ fluorescence) are left out where the source
reports the method.

## Data sources

Every database used or checked, with access route, licence, method and why it was or was not used:
[`data/sources.csv`](data/sources.csv) (also the "Data sources" table of the report).

| Source | Region | Licence |
|---|---|---|
| [US Water Quality Portal](https://www.waterqualitydata.us) (USGS, US EPA, NWQMC) | USA | public domain |
| [EEA Waterbase Water Quality ICM 2026](https://sdi.eea.europa.eu/datashare/s/ZcJ37WK6oaWonFJ) | Europe | CC BY 4.0 |
| [SYKE VESLA open data](https://rajapinnat.ymparisto.fi/api/vesla/2.0/odata) | Finland | CC BY 4.0 |
| [UNEP GEMStat GFQA v3](https://doi.org/10.5281/zenodo.18459694) | global | CC BY 4.0 |
| [EPA Victoria, Gippsland Lakes](https://discover.data.vic.gov.au/en_AU/dataset/gippsland-lakes-water-quality-data) | Australia | CC BY 4.0 |
| [Environment Canterbury](http://wateruse.ecan.govt.nz/WQAll.hts?service=Hilltop&request=SiteList) (Hilltop service) | New Zealand | CC BY 4.0 |
| [Lewis et al. (2023), EDI edi.1530.1](https://doi.org/10.6073/pasta/2cd6628a942de2a8b12d2b19962712a0) | North America | CC BY 4.0 |
| Local lake observations (Erken, Bosumtwi; ISIMIP local lakes) | Sweden, Ghana | not public: in the results, not in `data/` |
| [COGERH Ceará](http://www.hidro.ce.gov.br/acude/eutrofizacao) | Brazil | not stated: sites only (map), re-fetch with the script |

Satellite data: Lakes_cci v2.1 per-lake daily tables of Shahvaran et al. (ESSD preprint essd-2026-528,
[Zenodo 10.5281/zenodo.20347115](https://doi.org/10.5281/zenodo.20347115), CC BY 4.0). Lake outlines: LakeATLAS v1.0
/ HydroLAKES (Lehner et al. 2022, [HydroATLAS](https://www.hydrosheds.org/hydroatlas), CC BY 4.0).

## Contents

| Path | What |
|---|---|
| `scripts/cci_polygons.py` | HydroLAKES outlines of the Lakes_cci lakes (`data/cci_hylak_ids.txt`: Lakes_cci v2.1 and v3.0 lakes) |
| `scripts/insitu_wqp.py` | US Water Quality Portal: lab chlorophyll-a at stations inside the lake outlines |
| `scripts/insitu_more.py` | EEA, SYKE, GEMStat, EPA Victoria, Environment Canterbury, COGERH: samples inside the outlines, and an inventory of every site with chlorophyll-a (`*_sites.csv`) |
| `scripts/insitu_chl.py` | one daily table per lake; a lake found in several sources keeps the one with most dates |
| `scripts/insitu_timing.py` | same-day pairs, monthly cycles, peak months, figures |
| `scripts/insitu_report.py` | the HTML report |
| `data/insitu/` | per-source samples (one row per sample) and site inventories |
| `data/insitu_chl.csv` | merged daily in-situ table (without the non-public local lakes) |
| `data/results/v21/` | per-lake results (`lakes.csv`), pairs, monthly bias, monthly cycles, figures |
| `docs/index.html` | the report (GitHub Pages) |

## Reproduce

Python ≥ 3.10 with the packages in `requirements.txt`, and GDAL's `ogr2ogr`.

```bash
# inputs (downloaded once)
export CCI_V21=/path/to/shahvaran     # hylak_id_description.csv and daily/<hylak_id>.csv (Zenodo 20347115)
export LAKEATLAS=/path/to/LakeATLAS_v10_pol_Global.gpkg
./run_all.sh                          # outlines, merge, comparison, report from the data in data/insitu
FETCH=1 EEA_DIR=... GEMSTAT_DIR=... EPAVIC_DIR=... ./run_all.sh   # also fetch the in-situ data again
```

Without the local lakes (not public) the comparison has 147 lakes instead of 149.

## Licences

Code: GNU GPL v3 (`LICENSE`). Data tables in `data/` and the report: CC BY 4.0, derived from the sources above,
which must be credited (see `DATA_LICENSE.md`).

Part of the GLM4-AED3 global lake modelling work (ISIMIP).
