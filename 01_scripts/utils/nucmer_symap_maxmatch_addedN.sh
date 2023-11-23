#!/bin/bash

# Prepare files required for running SyMAP to identify syntenic regions, because it cannot handle self-alignements for huge genomes. We first map the genome to itself with nucmer, then we trick symap to use the .mum file outputted by nucmer.
# Input is the masked reference genome outputted by RepeatMasker (script utils/RepeatMasker_manitou.sh), with chr names corrected to suit SyMAP (e.g. no special characters, ., _, -, ..)
# Trials resulted in nucmer failure due to a problematic region on chr OV354449, so we had to "patch" this region by masking it with Ns prior to running nucmer.

# We first "patch" the problematic region in masked reference genome, then run nucmer. The .mum file was then transfered on the machine on which SyMAP is installed. 

# srun -p large -c 8 --mem=150G --time=21-00:00:00 -J addedN_mummer_self_synteny -o log/nucmer_self_synteny_addedN_%j.log /bin/sh 01_scripts/nucmer_symap_maxmatch_addedN.sh &


# VARIABLES
NB_CPU=8
MASKED_DIR="04_RepeatMasker_manitou"
INI_REF="genome_chrs_corrected.masked.fasta"           # Original reference genome, masked using RepeatMasker, with chr names corrected for running Symap (no ., _, -, ...)

MASKED_REF="genome_chrs_corrected.masked_addedN.fasta" # Masked reference with "patched" region, output of current script's first step 

OUTPUT_DIR="05_self_align"


# LOAD REQUIRED_MODULES
#module load mummer/3.23

if [[ ! -d $OUTPUT_DIR ]]
then
  mkdir $OUTPUT_DIR
fi


# 1. Trunc ref fasta to mask 40 kb around a problematic region around position 29319494 of chr OV354449. 
## Previous trials crashed and threw the message ERROR: failed to merge alignments at position 29319494

## Make a bed with sites to mask
echo -e "ChrOV354449\t29309494\t29329494\n" > 02_infos/ChrOV354449_to_mask.bed

## Mask with N using bedtools
bedtools maskfasta -mc N -fi $MASKED_DIR/$INI_REF -bed 02_infos/ChrOV354449_to_mask.bed -fo $MASKED_DIR/$MASKED_REF

# 2. Run nucmer
nucmer -t $NB_CPU --maxmatch \
   $MASKED_DIR/$MASKED_REF $MASKED_DIR/$MASKED_REF \
   -p "$OUTPUT_DIR"/"${MASKED_REF%.fasta}"vs"${MASKED_REF%.fasta}".maxmatch

show-coords -dlTH "$OUTPUT_DIR"/"${MASKED_REF%.fasta}"vs"${MASKED_REF%.fasta}".maxmatch.delta > "$OUTPUT_DIR"/"${MASKED_REF%.fasta}"vs"${MASKED_REF%.fasta}".maxmatch.mum
