#!/bin/bash

# Compute SAF and MAF by pop, done after Fst script

# manitou
# srun -p small -c 1 -J SVs_get_maf_nind_bypop -o log/SVs_get_maf_nind_bypop_%j.log /bin/sh 01_scripts/utils/SVs_get_maf_nind_bypop.sh &

# valeria
# srun -p ibis_small -c 1 -J SVs_get_maf_nind_bypop -o log/SVs_get_maf_nind_bypop_%j.log /bin/sh 01_scripts/utils/SVs_get_maf_nind_bypop.sh &

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

MIN_MAF=0.05
MAX_MAF=0.95

POPS_FILE="02_infos/pops.txt"

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

POP1='RO'
POP2='PU'

# LOAD REQUIRED MODULES
module load R/4.1

# 1. Extract MAF and nb of samples with data from .mafs.gz file
cat $POPS_FILE | while read POP
do
  zless $ANGSD_BYPOP_DIR/"$POP".mafs.gz | cut -f6,7 | tail -n+2 > $ANGSD_BYPOP_DIR/"$POP".knownEM_nInd
done
  
# 2. Extract sites   
zless $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".bed | cut -f1,3 > $ANGSD_BYPOP_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".chrom_pos

# 3. Paste required files
echo -e "CHROM\tPOS\tknownEM_"$POP1"\tnInd_"$POP1"\tknownEM_"$POP2"\tnInd_"$POP2"" > $ANGSD_BYPOP_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF"_chrom_pos_freq_bypop.txt

paste -d "\t" $ANGSD_BYPOP_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".chrom_pos $ANGSD_BYPOP_DIR/"$POP1".knownEM_nInd $ANGSD_BYPOP_DIR/"$POP2".knownEM_nInd >> $ANGSD_BYPOP_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF"_chrom_pos_freq_bypop.txt

Rscript 01_scripts/utils/maf_nInd_fst_by_site.R $ANGSD_BYPOP_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF"_chrom_pos_freq_bypop.txt $ANGSD_FST_DIR/"$POP1"_"$POP2"/"$POP1"_"$POP2".bypos.sfs.annot.gz $ANGSD_BYPOP_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)"_site_maf_nInd_fst