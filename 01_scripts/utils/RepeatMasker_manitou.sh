#!/bin/sh

# Run on Manitou
# srun -p medium -c 20 --mem=50G --time=2-00:00:00 -J 01_RepeatMasker -o log/01_RepeatMasker_%j.log /bin/sh 01_scripts/01_RepeatMasker_manitou.sh &


# VARIABLES
GENOME="03_genome/genome.fasta"
RM_DIR="04_RepeatMasker_manitou"
CPU=20

# LOAD REQUIRED MODULES 
module load gnu-openmpi/4.0.5
module load exonerate/2.4.0
module load RepeatMasker/4.0.8
module load ncbiblast/2.6.0
module load python/2.7
module load maker/2.31.10


# 1. Run RepeatMasker
RepeatMasker -pa $CPU -species 'salmo salar' $GENOME -dir $RM_DIR -gff
