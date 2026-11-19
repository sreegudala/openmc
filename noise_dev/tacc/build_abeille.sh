#!/bin/bash
# Build Abeille (https://github.com/HunterBelanger/abeille, GPL-3.0) outside
# the repo for the Stage 6c reference. Run on a login node (CMake downloads
# Abeille's dependencies, including HDF5):
#
#   bash tacc/build_abeille.sh            # from noise_dev/
#
# Destination: $ABEILLE_DIR, default <parent of $OPENMC_DEV_SRC>/abeille.
#
# Pinned to commit c6aaa5f (the one read in Stage 0) plus one local fix:
# Abeille's noise driver restores kcol after each noise batch but records the
# batch's (meaningless) complex-weight kcol into kcol_avg first. After two
# noise batches KeffFissionOperator switches to kcol_avg, and with 1e5
# particles the static power iteration then produced no fission neutrons
# ("No fission neutrons were produced", reproducible). The fix restores kcol
# before record_generation. It changes only the population normalisation of
# the static power iteration, not the noise transport (which divides by the
# input keff).
set -euo pipefail
: "${OPENMC_DEV_SRC:?source ~/openmc-dev.bashrc first}"
DEST=${ABEILLE_DIR:-$(dirname "$OPENMC_DEV_SRC")/abeille}

if [ ! -d "$DEST/.git" ]; then
  git clone https://github.com/HunterBelanger/abeille.git "$DEST"
fi
cd "$DEST"
git fetch -q origin || true
git checkout -q -B local-noise-kcol-fix c6aaa5f

python3 - <<'EOF'
from pathlib import Path
p = Path("src/noise.cpp")
s = p.read_text()
old = """  // Get new values
  Tallies::instance().calc_gen_values();

  // Keep values
  Tallies::instance().record_generation(avg_wgt_mag);"""
new = """  // Get new values
  Tallies::instance().calc_gen_values();

  // Local fix: keep kcol_avg static (see noise_dev/tacc/build_abeille.sh)
  Tallies::instance().set_kcol(original_kcol);

  // Keep values
  Tallies::instance().record_generation(avg_wgt_mag);"""
if new not in s:
    assert s.count(old) == 1, "Abeille source does not match c6aaa5f"
    p.write_text(s.replace(old, new))
EOF

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DABEILLE_USE_LTO=OFF \
  -DABEILLE_GUI_PLOT=OFF -DABEILLE_USE_MPI=OFF
cmake --build build -j "${BUILD_JOBS:-8}"
echo "Abeille built: $DEST/build/abeille"
"$DEST/build/abeille" --version | head -3
