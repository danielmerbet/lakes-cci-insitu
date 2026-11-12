#!/usr/bin/env python3
"""HydroLAKES polygons of the Lakes_cci lakes (input of insitu_wqp.py, insitu_more.py and insitu_chl.py).

    python3 cci_polygons.py <LakeATLAS_v10_pol_Global.gpkg> <Hylak_id list, one per line> <out .geojson>

LakeATLAS v1.0 polygons (HydroLAKES v1.0 geometry; Lehner et al. 2022, https://www.hydrosheds.org/hydroatlas, CC BY 4.0).
Keeps Hylak_id, Lake_name, Country and Lake_area. ogr2ogr's SQL parser runs out of memory with thousands of ids in
one IN (...) list, so the ids go in chunks of 800 and the pieces are merged.
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main(gpkg, ids_file, out):
    ids = [int(x) for x in Path(ids_file).read_text().split()]
    feats = []
    with tempfile.TemporaryDirectory() as tmp:
        for i in range(0, len(ids), 800):
            part = Path(tmp) / f"part{i}.geojson"
            where = "Hylak_id IN (" + ",".join(map(str, ids[i:i + 800])) + ")"
            subprocess.run(["ogr2ogr", "-f", "GeoJSON", "-select", "Hylak_id,Lake_name,Country,Lake_area",
                            "-where", where, str(part), gpkg], check=True)
            feats += json.loads(part.read_text())["features"]
    Path(out).write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    print(f"{len(feats)} of {len(ids)} lakes -> {out}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
