#!/bin/bash
# Diagnose a femffusion.exe crash on LS6 (built by build_femffusion.sh).
# Run on a login node from noise_dev/ and paste the whole output:
#
#   bash tacc/diag_femffusion.sh
#
# History: with dealii/9.4.0-real, femffusion.exe segfaults before printing
# anything (-t, a single case, I_MPI_FABRICS=shm, mpirun -n 1 alike), i.e.
# during library loading or MPI/PETSc start-up. Its deal.II pulls in all of
# Trilinos 13 (with libpytrilinos) and two Kokkos builds (Trilinos's and
# PETSc 3.16's). This script locates the crash and checks the container
# fallback.
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

step "1. last libraries initialised before the crash (ld.so trace)"
# ld.so writes this trace unbuffered to stderr, so it survives the crash.
# The last 'calling init' names the library whose constructor crashed;
# 'transferring control' means all loaded and the crash is in main
# (MPI_Init / PetscInitialize).
LD_DEBUG=libs timeout 120 ./femffusion.exe -t > diag.out 2>&1
echo "[exit status $?]"
grep -E 'calling init|transferring control' diag.out | tail -6 |
  sed 's/^ *[0-9]*:[[:space:]]*//'
if grep -q 'transferring control: ./femffusion.exe' diag.out; then
  echo "(main was reached: the crash is after start-up loading)"
else
  echo "(main was NOT reached: the crash is in a library constructor)"
fi

step "2. duplicated libraries (same name, two builds)"
ldd ./femffusion.exe | awk '{print $1}' | sed 's/\.so.*//' | sort | uniq -c |
  awk '$1 > 1'
ldd ./femffusion.exe | grep -E 'kokkoscore|libmpi(_ilp64)?\.so|pytrilinos' |
  awk '{print $1, $3}'

step "3. a deal.II-free MPI program with the same modules"
cat > diag_mpi.c <<'EOF'
#include <mpi.h>
#include <stdio.h>
int main(int argc, char **argv)
{ MPI_Init(&argc, &argv); printf("MPI ok\n"); MPI_Finalize(); return 0; }
EOF
mpicc diag_mpi.c -o diag_mpi > /dev/null 2>&1 && timeout 60 ./diag_mpi
echo "[exit status $?]"

step "4. container tools (fallback: Ubuntu 24.04 deal.II 9.5.1 image)"
command -v apptainer singularity 2> /dev/null || echo "none on PATH"
module spider apptainer 2>&1 | grep -Ei 'apptainer/|versions' | head -4
module spider tacc-apptainer 2>&1 | grep -Ei 'apptainer/|versions' | head -4

rm -f diag.out diag_mpi diag_mpi.c
