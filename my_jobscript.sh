#!/usr/bin/env bash
#SBATCH --partition=short
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=simon.moulds@ouce.ox.ac.uk
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=esmvaltool
#SBATCH --time=12:00:00

# cd $SCRATCH || exit 1

module load Anaconda3/2022.05
module load NCO/5.0.1-foss-2021a
source activate $DATA/envs/snakemake

# PYTHONPATH sometimes causes issues
export PYTHONPATH=

# Make sure R scripts can access local libraries
export R_LIBS_USER=$HOME/local/rlibs

snakemake \
    --snakefile workflow/Snakefile_1 \
    --profile $HOME/.config/snakemake/slurm.arc \
    --cores 1 \
    --config input_data_root=/data/ouce-drift/cenv0857 \
    --use-conda \
    --conda-base-path /data/ouce-drift/cenv0857 \
    --rerun-incomplete
