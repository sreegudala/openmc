#!/bin/bash
# Build FEMFFUSION (https://github.com/Zonni/FEMFFUSION, GPL-3.0) outside
# the repo, for the deterministic cross-checks (Stage 9a C5G7, and the
# informational 4b/6c/8c comparisons). Run on a login node:
#
#   bash tacc/build_femffusion.sh            # from noise_dev/
#
# FEMFFUSION needs deal.II (>= 8.0) configured with PETSc and SLEPc, real
# scalars: its noise solver works on 2 x 2 real blocks. It is tested with
# deal.II 9.5.1 / PETSc 3.19 (Ubuntu 24.04 packages, the local build). LS6
# has no 9.5; 9.4.0-real is the closest older release (9.6 removed many
# functions 9.5 deprecated). TACC's PETSc modules include SLEPc.
#
# Environment:
#   FEMFFUSION_DIR       clone and build here (default: <parent of
#                        $OPENMC_DEV_SRC>/FEMFFUSION, else $SCRATCH/FEMFFUSION)
#   FEMFFUSION_PREREQ    modules to load before deal.II (compiler, MPI), as
#                        listed by `module spider <deal.II module>`; for
#                        dealii/9.4.0-real on LS6: "intel/19.1.1 impi/19.0.9"
#   FEMFFUSION_DEALII    deal.II module (default dealii/9.4.0-real)
#
# The same modules must be loaded to run femffusion.exe; the script writes
# them to $FEMFFUSION_DIR/modules.sh.
set -eo pipefail

DEALII=${FEMFFUSION_DEALII:-dealii/9.4.0-real}
PREREQ=${FEMFFUSION_PREREQ:-}
if [ -n "${OPENMC_DEV_SRC:-}" ]; then
  DEST=${FEMFFUSION_DIR:-$(dirname "$OPENMC_DEV_SRC")/FEMFFUSION}
else
  DEST=${FEMFFUSION_DIR:-$SCRATCH/FEMFFUSION}
fi

# `module` is a shell function; make sure it exists in this non-interactive
# shell (Lmod exports LMOD_PKG)
if ! type module > /dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "${LMOD_PKG:-/opt/apps/lmod/lmod}/init/bash"
fi

module reset > /dev/null 2>&1 || true
for m in $PREREQ; do
  module load "$m"
done
if ! module load "$DEALII"; then
  echo
  echo "Could not load $DEALII. Its prerequisites are:"
  module spider "$DEALII" 2>&1 | sed -n '/You will need/,/^$/p'
  echo "Set FEMFFUSION_PREREQ to one of those lines (space-separated) and rerun."
  exit 1
fi
module list 2>&1 | sed -n '1,20p'

if [ ! -d "$DEST/.git" ]; then
  git clone https://github.com/Zonni/FEMFFUSION.git "$DEST"
fi
cd "$DEST"
git fetch -q origin || true
git checkout -q 1611f90
rm -f CMakeCache.txt.tmp*

# deal.II's autopilot build, as in the FEMFFUSION README.
# CMAKE_POLICY_VERSION_MINIMUM: CMake >= 4 rejects FEMFFUSION's
# cmake_minimum_required(VERSION 2.8.8) (as for Abeille's docopt.cpp)
DEAL_II_DIR=${TACC_DEALII_DIR:-${DEAL_II_DIR:-}} cmake . \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  ${TACC_DEALII_DIR:+-DDEAL_II_DIR=$TACC_DEALII_DIR}
make release
make -j "${BUILD_JOBS:-16}"

{
  echo "# modules for $DEST/femffusion.exe (written by build_femffusion.sh)"
  echo "module reset > /dev/null 2>&1 || true"
  for m in $PREREQ; do echo "module load $m"; done
  echo "module load $DEALII"
  echo "export OMP_NUM_THREADS=1"
} > "$DEST/modules.sh"

echo "FEMFFUSION built: $DEST/femffusion.exe"
echo "Running its test suite (about a minute) ..."
OMP_NUM_THREADS=1 ./femffusion.exe -t 2>&1 | tail -4
