#!/usr/bin/env bash
# Download MERRA-2 PBLH (planetary boundary layer height, m) hourly at the Blois
# grid cell (lat idx 275 = 47.5N, lon idx 290 = 1.25E) for the full forcing span
# 2021-01 .. 2022-07, via GES DISC OPeNDAP ascii (NASA Earthdata .netrc auth).
# Resumable (skips dates already in the CSV). Output: long CSV datetime,pblh.
set -u
cd /home/corroyez/Documents/z_Example_rmusica_31012025
OUT=out_files/Chapter3/merra2_pblh.csv
COOK=/tmp/urs_cookies
BASE="https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXFLX.5.12.4"
LAT=275; LON=290
mkdir -p out_files/Chapter3
[ -f "$OUT" ] || echo "datetime,pblh" > "$OUT"

months="2021/01 2021/02 2021/03 2021/04 2021/05 2021/06 2021/07 2021/08 2021/09 2021/10 2021/11 2021/12 2022/01 2022/02 2022/03 2022/04 2022/05 2022/06 2022/07"

fetch_day() {
  local ym="$1" fname="$2"
  local ymd=$(echo "$fname" | grep -oE "[0-9]{8}\.nc4" | cut -c1-8)
  local d="${ymd:0:4}-${ymd:4:2}-${ymd:6:2}"
  grep -q "^${d}T00:30" "$OUT" 2>/dev/null && return 0   # resume: already done
  local url="$BASE/$ym/${fname}.ascii?PBLH%5B0:1:23%5D%5B${LAT}:1:${LAT}%5D%5B${LON}:1:${LON}%5D"
  local resp h v line
  resp=$(curl -gsS -L -n -c "$COOK" -b "$COOK" "$url" 2>/dev/null)
  # parse the 24 "...PBLH.PBLH[...time=MM]..., value" lines, in order
  echo "$resp" | grep -E "PBLH\.PBLH\[" | awk -F',' '{gsub(/ /,"",$2); print $2}' > /tmp/pblh_$$_${ymd}.tmp
  local n=$(wc -l < /tmp/pblh_$$_${ymd}.tmp)
  if [ "$n" -ne 24 ]; then echo "WARN $d: got $n values (skip)"; rm -f /tmp/pblh_$$_${ymd}.tmp; return 1; fi
  h=0
  while read v; do printf "%sT%02d:30,%s\n" "$d" "$h" "$v"; h=$((h+1)); done < /tmp/pblh_$$_${ymd}.tmp >> "$OUT"
  rm -f /tmp/pblh_$$_${ymd}.tmp
}
export -f fetch_day; export OUT COOK BASE LAT LON

for ym in $months; do
  files=$(curl -gsS -L -n -c "$COOK" -b "$COOK" "$BASE/$ym/contents.html" 2>/dev/null \
          | grep -oE "MERRA2_[0-9]+\.tavg1_2d_flx_Nx\.[0-9]{8}\.nc4" | sort -u)
  [ -z "$files" ] && { echo "no files for $ym"; continue; }
  echo "$files" | xargs -I{} -P 6 bash -c 'fetch_day "$0" "$1"' "$ym" {}
  echo "done $ym : total rows $(( $(wc -l < "$OUT") - 1 ))"
done
echo "ALL DONE. rows=$(( $(wc -l < "$OUT") - 1 )) (expect ~13848 for 577 days)"
