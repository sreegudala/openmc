# Running the noise_dev stage drivers on TACC

Each stage driver (`stageN/run_stageN.py`) runs 20 independent replicas per
case. `common/harness.py` runs them as separate OpenMC processes, several at
a time, so one node does a whole case in about the time of one replica.

| Variable | Default | Meaning |
|---|---|---|
| `NOISE_THREADS` | 12 locally; cores / `NOISE_JOBS` in `noise_job.slurm` | OpenMP threads per OpenMC process |
| `NOISE_JOBS` | 1 locally; 20 in `noise_job.slurm` | replicas run at once (capped at the replica count) |
| `NOISE_OPENMC_CMD` | `openmc` | launch command, e.g. `"ibrun -n 1 openmc"` if the MPI build refuses to start as a singleton |

No MPI ranks are used: the noise solver's per-generation loop is threaded,
and replica-level concurrency is the efficient way to fill a node.

## Setup (once per pull)

```bash
source ~/openmc-dev.bashrc
cd $OPENMC_DEV_SRC
git pull
odev_build                # C++ changed: rebuild and install into the venv
cd noise_dev
```

## Submitting

Submit from `noise_dev/`; the driver and its arguments follow the script.
Logs land in `noise_dev/slurm_<name>_<jobid>.out`, the stage log in
`stageN/stageN_log_<only>[_<tag>].md`, and run directories in
`stageN/work[_<tag>]/`.

```bash
# Stage 3b at 4x particles with fresh seeds (follow-up of the 1 Hz miss)
sbatch -J s3b_p20k tacc/noise_job.slurm stage3/run_stage3.py \
    --only 3b --particles 20000 --base-seed 100001 --tag p20k

# Stage 4 (all parts in one job: they share libraries and work dirs)
sbatch -J s4 tacc/noise_job.slurm stage4/run_stage4.py
```

`run_stage3.py` options for reruns: `--freqs 1,10` (subset of the F7 grid
for 3a/3b/3d), `--base-seed N` (replica r uses seed N + 7919 r), `--tag T`
(work dir `work_T/`, log suffix `_T`, so it never overwrites the default
run).

## Stage 6c: Abeille reference

Abeille (GPL-3.0) is built outside the repo, in `$ABEILLE_DIR` (default
`<parent of $OPENMC_DEV_SRC>/abeille`), on a login node because CMake
downloads its dependencies:

```bash
bash tacc/build_abeille.sh
```

The script pins commit c6aaa5f and applies one local fix to Abeille's noise
driver (see the script header). Then, from `noise_dev/`:

```bash
sbatch -J s6c_ab tacc/noise_job.slurm stage6/run_stage6c.py abeille
sbatch -J s6c_omc tacc/noise_job.slurm stage6/run_stage6c.py openmc
```

`abeille` runs 20 replicas of 600 noise batches (1e5 particles, Abeille's
cancellation) and stops each cleanly after `--max-minutes` (default 1260).
The Abeille static flux used for normalisation is committed
(`stage6c_abeille_static.npz`); `analyze` runs locally.

## FEMFFUSION (deterministic cross-checks)

FEMFFUSION (GPL-3.0) is built outside the repo, in `$FEMFFUSION_DIR`
(default `<parent of $OPENMC_DEV_SRC>/FEMFFUSION`), against TACC's deal.II
module (real scalars; TACC's PETSc modules include SLEPc). Find the
prerequisites of the deal.II module, then build on a login node:

```bash
module spider dealii/9.4.0-real          # "You will need": intel/19.1.1 impi/19.0.9
FEMFFUSION_PREREQ="intel/19.1.1 impi/19.0.9" bash tacc/build_femffusion.sh
```

The script pins commit 1611f90 (the one in the design doc), runs
FEMFFUSION's test suite, and writes the run-time modules to
`$FEMFFUSION_DIR/modules.sh`. Locally FEMFFUSION is built against the
Ubuntu 24.04 packages (`libdeal.ii-dev petsc-dev slepc-dev`), the versions
it is tested with (99/99 tests).

TACC's `dealii/9.4.0-real` needs workarounds, all in the script and none
changing results:
- it compiles in a pre-C++17 mode: `std::filesystem` (output folder
  creation) is replaced by POSIX `stat`/`mkdir`;
- its link line names `/usr/lib64/libtbb.so`, absent on LS6: pointed at the
  installed `libtbb.so.2`;
- it links PyTrilinos, whose constructor segfaults: removed from the
  executable and a private copy of `libdeal_II` (patchelf);
- it is built without LAPACK: FEMFFUSION is compiled with
  `DEAL_II_WITH_LAPACK` so its own eigenvalue calls reach MKL, and the
  three ROM unit tests, which need LAPACK inside `libdeal_II`, are skipped.
  The ROM features do not work on LS6; nothing in noise_dev uses them.

`tacc/diag_femffusion.sh` locates a start-up crash (ld.so trace). The
reference runs for Stage 9a are small and can equally be done locally.

## Bringing results back

From the local machine (adjust the host), copy logs and summaries but not the
statepoint-heavy work directories:

```bash
rsync -av --exclude 'work*/' --exclude '__pycache__/' \
    <tacc-host>:/home1/07667/sg44769/scratch_secure/sree/OpenMC/openmc-dev/noise_dev/ \
    ~/dev/openmc-dev/noise_dev/
```

If a run fails, also copy the `openmc.out` of the failing replica, e.g.
`stage3/work_p20k/3b_1.0Hz_eta1.0/rep07/openmc.out`.
