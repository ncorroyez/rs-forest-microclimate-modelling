#!/bin/bash
# Surgically patch musica.nml to match the legacy Oct 2025 functional params,
# while keeping path-related lines intact so MuSICA still finds its inputs
# through callmusica's tempdir.
#
# Backup → musica.nml.before_nmlpatch
# Restore : run with `restore` argument.

set -e
HERE="$(cd "$(dirname "$0")/.." && pwd)"
NML="$HERE/musica.nml"
BACKUP="$HERE/musica.nml.before_nmlpatch"

case "${1:-patch}" in
  patch)
    if [ -f "$BACKUP" ]; then
      echo "ERROR: backup already exists at $BACKUP — refusing to overwrite."
      echo "Run \`$0 restore\` first, or delete the backup manually."
      exit 1
    fi
    cp "$NML" "$BACKUP"
    echo "[backup] $BACKUP  md5=$(md5sum "$NML" | cut -c1-8)"

    # Surgical edits (in-place, GNU sed)
    sed -i 's/^[[:space:]]*FORCING_TIMESTEP[[:space:]]*=[[:space:]]*3600.*$/  FORCING_TIMESTEP = 1800/'              "$NML"
    sed -i 's/^[[:space:]]*TEST_CASE_ID[[:space:]]*=[[:space:]]*-1.*$/  TEST_CASE_ID = 2/'                          "$NML"
    sed -i 's/^[[:space:]]*WATER_VOLUME_PERLAI[[:space:]]*=.*$/  WATER_VOLUME_PERLAI = 10./'                        "$NML"
    sed -i 's/^[[:space:]]*CO2_H2O_EQ_DEGREE_INLEAF[[:space:]]*=.*$/  CO2_H2O_EQ_DEGREE_INLEAF = 1.0/'              "$NML"
    sed -i 's/^[[:space:]]*PECLET_INLEAF[[:space:]]*=[[:space:]]*0\.2,[[:space:]]*0\.2.*$/  PECLET_INLEAF = 0.2/'   "$NML"
    # Activate MIXING_LENGTH_INLEAF (legacy had it set)
    sed -i 's|^[[:space:]]*!\s*MIXING_LENGTH_INLEAF[[:space:]]*=.*$|  MIXING_LENGTH_INLEAF = 0.2|'                  "$NML"
    # Add ABL_flag = 'none' just after DEBUG = .true. (if not already present)
    if ! grep -q "^[[:space:]]*ABL_flag" "$NML"; then
      sed -i "/^[[:space:]]*DEBUG[[:space:]]*=[[:space:]]*\\.true\\./a\\  ABL_flag = 'none'" "$NML"
    fi

    echo "[patch] applied to $NML  new md5=$(md5sum "$NML" | cut -c1-8)"
    echo
    echo "Diff vs backup :"
    diff "$BACKUP" "$NML" | head -30
    ;;
  restore)
    if [ ! -f "$BACKUP" ]; then
      echo "ERROR: no backup found at $BACKUP — nothing to restore."
      exit 1
    fi
    cp "$BACKUP" "$NML"
    md5_now=$(md5sum "$NML" | cut -c1-8)
    md5_bak=$(md5sum "$BACKUP" | cut -c1-8)
    if [ "$md5_now" = "$md5_bak" ]; then
      echo "[restore] OK : $NML  md5=$md5_now"
      rm "$BACKUP"
    else
      echo "ERROR: restore md5 mismatch ($md5_now vs $md5_bak), backup NOT removed"
      exit 1
    fi
    ;;
  *) echo "Usage: $0 [patch|restore]"; exit 1 ;;
esac
