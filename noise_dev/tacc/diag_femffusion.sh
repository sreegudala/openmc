#!/bin/bash
# Diagnose a femffusion.exe crash on LS6 (built by build_femffusion.sh).
# Run on a login node from noise_dev/ and paste the whole output:
#
#   bash tacc/diag_femffusion.sh
#
# Each step's output is line-buffered (stdbuf), so a crash keeps it.
if [ -n "${OPENMC_DEV_SRC:-}" ]; then
  DEST=${FEMFFUSION_DIR:-$(dirname "$OPENMC_DEV_SRC")/FEMFFUSION}
else
  DEST=${FEMFFUSION_DIR:-$SCRATCH/FEMFFUSION}
fi
if ! type module > /dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "${LMOD_PKG:-/opt/apps/lmod/lmod}/init/bash"
fi
cd "$DEST" || exit 1
# shellcheck disable=SC1091
source modules.sh
export OMP_NUM_THREADS=1

step() { echo; echo "=== $*"; }
run() {  # run with a time limit; print the first lines and the status
  timeout 120 stdbuf -oL -eL "$@" > diag.out 2>&1
  local rc=$?
  head -12 diag.out
  echo "[exit status $rc]"
}

step "modules"
module list 2>&1 | sed -n '1,12p'

step "unresolved or key shared libraries"
ldd ./femffusion.exe | grep -Ei 'not found|tbb|mpi|petsc|slepc|deal|p4est'

step "1. femffusion.exe -t (singleton MPI)"
run ./femffusion.exe -t

step "2. a single static case, no test harness"
run ./femffusion.exe -f test/1D_hom_slab_2cm/1D_hom_zerocurrent.prm

step "3. -t with Intel MPI on shared memory only"
I_MPI_FABRICS=shm run ./femffusion.exe -t

step "4. -t under mpirun -n 1"
run mpirun -n 1 ./femffusion.exe -t

step "5. backtrace of the crash (catchsegv)"
if command -v catchsegv > /dev/null; then
  timeout 120 catchsegv ./femffusion.exe -t > diag.out 2>&1
  grep -A30 'Backtrace' diag.out | head -30
else
  echo "catchsegv not available"
fi
rm -f diag.out
