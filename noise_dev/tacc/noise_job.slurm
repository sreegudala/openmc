#!/bin/bash
# Run one noise_dev stage driver on a single node.
#
# Submit from noise_dev/ with the driver and its arguments after the script:
#   sbatch -J s3b_p20k tacc/noise_job.slurm stage3/run_stage3.py --only 3b
#
# Replicas run as independent OpenMC processes, NOISE_JOBS at a time with
# NOISE_THREADS OpenMP threads each (no MPI). Override either with
#   sbatch --export=ALL,NOISE_THREADS=8 ...
#SBATCH -p NuclearEnergy
#SBATCH -A ASC23046
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -o slurm_%x_%j.out
#SBATCH -e slurm_%x_%j.err

set -eo pipefail
source ~/openmc-dev.bashrc
set -u

cd "$SLURM_SUBMIT_DIR"

# Several OpenMC processes share the node: do not pin their threads (each
# process would bind to the same cores starting at core 0).
unset OMP_PROC_BIND OMP_PLACES

export NOISE_THREADS=${NOISE_THREADS:-4}
NCPU=$(nproc)
export NOISE_JOBS=${NOISE_JOBS:-$(( NCPU / NOISE_THREADS > 0 ? NCPU / NOISE_THREADS : 1 ))}

echo "host:     $(hostname)"
echo "cpus:     $NCPU   jobs: $NOISE_JOBS   threads/job: $NOISE_THREADS"
echo "commit:   $(git -C "$SLURM_SUBMIT_DIR" log -1 --format='%h %s')"
echo "command:  python -u $*"
echo "start:    $(date)"

status=0
python -u "$@" || status=$?
echo "end:      $(date)   exit: $status"
exit $status
