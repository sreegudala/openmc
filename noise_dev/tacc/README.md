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
