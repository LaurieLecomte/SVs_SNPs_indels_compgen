#!/bin/bash

# Compute SAF and MAF by pop

# manitou
# srun -p small -c 4 -J 01.5_SVs_saf_maf_bypop -o log/01.5_SVs_saf_maf_bypop_%j.log /bin/sh 01_scripts/01.5_SVs_saf_maf_bypop.sh &

# valeria
# srun -p ibis_small -c 4 -J 01.5_SVs_saf_maf_bypop -o log/01.5_SVs_saf_maf_bypop_%j.log /bin/sh 01_scripts/01.5_SVs_saf_maf_bypop.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"

RAW_VCF_DIR="04_vcf/SVs"
RAW_SV_VCF="$RAW_VCF_DIR/merged_SUPP2_MAF0.05_FMISS0.5.vcf.gz"

ANGSD_INPUT_DIR="05_angsd_inputs/SVs"
ANGSD_STATS_DIR="06_angsd_stats/SVs"
ANGSD_BYPOP_DIR="07_angsd_bypop/SVs"
ANGSD_FST_DIR="08_angsd_fst/SVs"

PCA_DIR="09_pca/SVs"
RDA_DIR="10_rda/SVs"
ANNOT_DIR="11_annotation/SVs"
GO_DIR="12_go/SVs"

GENOME_ANNOT="11_annotation/genome_annotation/genome_annotation_table_simplified_1.5.tsv"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

SV_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".recoded.vcf.gz"

N_IND="$(less $ID_SEX_POP | wc -l)"

CPU=4

#MIN_MAF=0.05
#MAX_MAF=0.95

POPS_FILE="02_infos/pops.txt"

#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"
#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".vcf.gz"

REALSFS_PATH="/prg/angsd/0.937/misc/realSFS"

# LOAD REQUIRED MODULES
module load python/3.7
module load angsd/0.937
module load pcangsd/1.10
module load bcftools/1.13

# VARIABLES FOR ANGSD
WINDOW=1000000
WIN_STEP=10000 

# 1. Calculate saf and maf for each population
cat $POPS_FILE | while read POP
do
	# Extract samples belonging to POP
	less $ID_SEX_POP | grep $POP | cut -f1 > 02_infos/"$POP"_IDs.txt
	N_IND=$(less 02_infos/"$POP"_IDs.txt | wc -l)
	echo "working on pop $POP, $N_IND individuals, will use the sites file provided"
	
	# Filter vcf
	bcftools view -S 02_infos/"$POP"_IDs.txt $SV_VCF_ANGSD --threads $CPU | bcftools sort -Oz > $ANGSD_BYPOP_DIR/$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)_"$POP".vcf.gz
	tabix -p vcf $ANGSD_BYPOP_DIR/$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)_"$POP".vcf.gz -f
	
	echo "Calculate the SAF, MAF for $N_IND in vcf $SV_VCF_ANGSD"
	
	# Compute saf and maf
	angsd -vcf-pl $ANGSD_BYPOP_DIR/$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)_"$POP".vcf.gz  \
	-nind $N_IND -anc $GENOME -fai "$GENOME".fai \
	-dosaf 1 -domaf 1 -doMajorMinor 5 \
	-out $ANGSD_BYPOP_DIR/"$POP"
  
  # Check idx file
  $REALSFS_PATH check $ANGSD_BYPOP_DIR/"$POP".saf.idx
done


