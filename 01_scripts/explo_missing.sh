#!/bin/bash

# 
# manitou
# srun -c 1 -p small --time=1-00:00:00 --mem=10G -J explo_missing -o log/06_filter_pop_%j.log /bin/sh ./01_scripts/explo_missing.sh &

# valeria
# srun -c 1 -p ibis_small --time=1-00:00:00 --mem=10G -J explo_missing -o log/06_filter_pop_%j.log /bin/sh ./01_scripts/explo_missing.sh &

# VARIABLES
FST_DIR="08_angsd_fst"
SV_FILT_VCF="$FST_DIR/SVs/merged_SUPP2_MAF0.05_FMISS0.5.SVsFst_RO_PU.vcf.gz"
SNP_FILT_VCF="$FST_DIR/SNPs/SNPs_MAF0.05_FMISS0.5.SNPsFst_RO_PU.vcf.gz"
INDEL_FILT_VCF="$FST_DIR/indels/indels_MAF0.05_FMISS0.5.indelsFst_RO_PU.vcf.gz"

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load htslib/1.13

# Make a list of SVs, with FST and FMISS
bcftools query -f '%CHROM\t%POS\t%ID\t%END\t%FST_RO_PU\t%F_MISSING\n' $SV_FILT_VCF > $FST_DIR/SVs/"$(basename -s .vcf.gz $SV_FILT_VCF)"_FST_FMISS.table
bcftools query -f '%CHROM\t%POS\t%ID\t%END\t%FST_RO_PU\t%F_MISSING\n' $SNP_FILT_VCF > $FST_DIR/SNPs/"$(basename -s .vcf.gz $SNP_FILT_VCF)"_FST_FMISS.table
bcftools query -f '%CHROM\t%POS\t%ID\t%END\t%FST_RO_PU\t%F_MISSING\n' $INDEL_FILT_VCF > $FST_DIR/indels/"$(basename -s .vcf.gz $INDEL_FILT_VCF)"_FST_FMISS.table