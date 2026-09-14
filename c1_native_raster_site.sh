#!/bin/bash
# Native 20 m trait rasters for Blois: one fresh R process per tile (OOM guard),
# 3 in parallel, resumable (skips finished tiles). Driver over c1_native_raster_tile.R.
set -u
CTG=/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm
OUT=out_files/native20/tiles
LOG=out_files/native20/driver.log
RES=20
mkdir -p "$OUT"
: > "$LOG"
echo "START $(date) — res=$RES" >> "$LOG"
ls "$CTG"/*.las | \
  xargs -P 2 -I{} bash -c '
    f="{}"; b=$(basename "$f" .las); o="'"$OUT"'/${b}.tif"
    if [ -f "$o" ] && [ $(stat -c%s "$o") -gt 2000 ]; then echo "skip $b" >> "'"$LOG"'"; exit 0; fi
    Rscript c1_native_raster_tile.R "$f" "$o" '"$RES"' >/dev/null 2>&1
    if [ -f "$o" ]; then echo "done $b" >> "'"$LOG"'"; else echo "EMPTY $b" >> "'"$LOG"'"; fi
  '
n=$(ls "$OUT"/*.tif 2>/dev/null | wc -l)
echo "ALL_TILES_DONE $(date) — $n tifs" >> "$LOG"
