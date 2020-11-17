#!/bin/bash
#SBATCH -J run_souporcell
#SBATCH -p compute  # partition (queue)
#SBATCH -N 1        # number of nodes
#SBATCH -n 20        # number of cores
#SBATCH --mem 100   # memory pool for all cores
#SBATCH -t 0-24:00   # time (D-HH:MM)
#SBATCH -o %x.o%J   # STDOUT and STDERR
#SBATCH -e %x.o%J
module load singularity
module list

CWD=$(pwd)
echo "Current dir is $CWD"

echo "Running souporcell:"
singularity exec ../souporcell.sif souporcell_pipeline.py -i A.merged.bam \
    -b GSM2560245_barcodes.tsv -f refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
    -t 20 -o demux_data_test -k 4

# run this script by calling sbatch run_souporcell.slurm
