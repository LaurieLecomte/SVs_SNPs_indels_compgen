#!/bin/bash

# Run RepeatMasker on SVs sequences

# manitou
# srun -p small -c 10 -J SVs_TEs -o log/SVs_TEs_%j.log /bin/sh 01_scripts/SVs_TEs.sh &

# valeria
# srun -p ibis_small -c 10 -J SVs_TEs -o log/SVs_TEs_%j.log /bin/sh 01_scripts/SVs_TEs.sh &

# VARIABLE
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

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"


POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script

TE_DIR="TE"

CPU=10

# LOAD REQUIRED MODULES
module load R/4.1
module load bcftools/1.13

module load gnu-openmpi/4.0.5
module load exonerate/2.4.0
module load RepeatMasker/4.0.8
module load ncbiblast/2.6.0
module load python/2.7
module load maker/2.31.10


# 1. Extract Required fields from VCF
## We use the original VCF that was NOT formatted for angsd
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $RAW_FST_VCF > $TE_DIR/SVs_CHR_POS_ID_REF_ALT.table

# 2. Format into a fasta for running Repeat Masker
Rscript 01_scripts/utils/vcf_to_fasta.R $TE_DIR/SVs_CHR_POS_ID_REF_ALT.table
#Rscript 01_scripts/utils/extract_SV_fasta_LL.r $TE_DIR/SVs_CHR_POS_ID_REF_ALT.table $TE_DIR/SVs_CHR_POS_ID_REF_ALT_orig.fasta

# 3. Run RepeatMasker on SV sequences
RepeatMasker $TE_DIR/SVs_CHR_POS_ID_REF_ALT.fasta -pa $CPU -lib $TE_DIR/salmo_TElib.fasta -dir $TE_DIR
#RepeatMasker $TE_DIR/SVs_CHR_POS_ID_REF_ALT_orig.fasta -pa $CPU -lib $TE_DIR/salmo_TElib.fasta -dir $TE_DIR

RepeatMasker $TE_DIR/SVs_CHR_POS_ID_REF_ALT.fasta -pa $CPU -species 'salmo salar' -dir $TE_DIR/species
#RepeatMasker $TE_DIR/SVs_CHR_POS_ID_REF_ALT_orig.fasta -pa $CPU -species 'salmo salar' -dir $TE_DIR/species

#RepeatMasker $TE_DIR/SVs_CHR_POS_ID_REF_ALT.fasta -pa $CPU -species 'salmo salar' -dir $TE_DIR/species $GENOME -gff