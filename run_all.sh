#!/bin/bash
# Lakes_cci v2.1 chlorophyll-a vs independent in-situ data: the whole chain.
#   CCI_V21=<Shahvaran folder> LAKEATLAS=<LakeATLAS_v10_pol_Global.gpkg> ./run_all.sh
#   FETCH=1 also fetches the in-situ data again (WQP, SYKE, Environment Canterbury, COGERH over the internet;
#   EEA_DIR, GEMSTAT_DIR, EPAVIC_DIR = folders with the downloaded EEA Waterbase zips, GEMStat files and EPA xlsx).
# Outputs: data/insitu/ (if FETCH=1), data/insitu_chl.csv, data/results/v21/, docs/index.html. Work files in work/.
set -euo pipefail
cd "$(dirname "$0")"
: "${CCI_V21:?folder with hylak_id_description.csv and daily/ (Zenodo 10.5281/zenodo.20347115)}"
PY="python3 -I"
mkdir -p work data/insitu data/results/v21
POL=work/cci_polygons.geojson

if [ ! -s "$POL" ]; then
  : "${LAKEATLAS:?LakeATLAS_v10_pol_Global.gpkg (https://www.hydrosheds.org/hydroatlas)}"
  $PY scripts/cci_polygons.py "$LAKEATLAS" data/cci_hylak_ids.txt "$POL"
fi

if [ "${FETCH:-0}" = 1 ]; then
  $PY scripts/insitu_wqp.py "$POL" data/insitu/wqp_chl.csv
  $PY scripts/insitu_more.py syke "$POL" data/insitu/syke_chl.csv
  $PY scripts/insitu_more.py ecan "$POL" data/insitu/ecan_chl.csv
  $PY scripts/insitu_more.py cogerh "$POL" data/insitu/cogerh_chl.csv work/cogerh_pages
  [ -n "${EEA_DIR:-}" ] && $PY scripts/insitu_more.py eea "$POL" data/insitu/eea_chl.csv "$EEA_DIR"
  [ -n "${GEMSTAT_DIR:-}" ] && $PY scripts/insitu_more.py gemstat "$POL" data/insitu/gemstat_chl.csv "$GEMSTAT_DIR"
  [ -n "${EPAVIC_DIR:-}" ] && $PY scripts/insitu_more.py epavic "$POL" data/insitu/epavic_chl.csv "$EPAVIC_DIR"
  rm -f data/insitu/ecan_raw.csv data/insitu/cogerh_chl.csv   # cache; COGERH raw values are not redistributed
fi

$PY scripts/insitu_chl.py "$CCI_V21/hylak_id_description.csv" data/insitu_chl.csv data/insitu/wqp_chl.csv "$POL" \
  data/insitu/gemstat_chl.csv data/insitu/eea_chl.csv data/insitu/syke_chl.csv data/insitu/epavic_chl.csv \
  data/insitu/ecan_chl.csv
$PY scripts/insitu_timing.py data/insitu_chl.csv "$CCI_V21/daily" data/results/v21 data/insitu/*_sites.csv
$PY scripts/insitu_report.py data/results/v21 docs/index.html --standalone
cp docs/index.html index.html   # GitHub Pages from the repository root or from /docs
