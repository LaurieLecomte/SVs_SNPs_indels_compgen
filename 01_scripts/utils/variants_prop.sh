#!/bin/bash

# Estimate fraction of genome consisting of either variants types

# srun -p small -c 1 -J variants_prop -o log/variants_prop_%j.log /bin/sh 01_scripts/utils/variants_prop.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
VCF_DIR="04_vcf"
RAW_SV_VCF="$VCF_DIR/SVs/merged_SUPP2_MAF0.05_FMISS0.5.vcf.gz"
RAW_SNPS_VCF="$VCF_DIR/SNPs/SNPs_MAF0.05_FMISS0.5.vcf.gz"
RAW_INDELS_VCF="$VCF_DIR/indels/indels_MAF0.05_FMISS0.5.vcf.gz"

SV_TABLE="$VCF_DIR/SVs/merged_SUPP2.candidates.table" # was produced in genotype_SVs_SRLR by compare_summarize_plot.sh

#ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz"


# LOAD REQUIRED MODULES
module load bcftools/1.13
module load R/4.1

# 1. Extract CHROM, POS and END for easier handling
##bcftools query -f '%CHROM\t%POS\t%END\n' $RAW_SV_VCF > $VCF_DIR/SVs/"$(basename -s .vcf.gz $RAW_SV_VCF)"_CHROM_POS_END.table

#bcftools query -f '%CHROM\t%POS\t%END\n' $RAW_SNPS_VCF > $VCF_DIR/SNPs/"$(basename -s .vcf.gz $RAW_SNPS_VCF)"_CHROM_POS_END.table

#bcftools query -f '%CHROM\t%POS\t%END\n' $RAW_INDELS_VCF > $VCF_DIR/indels/"$(basename -s .vcf.gz $RAW_INDELS_VCF)"_CHROM_POS_END.table


#less $VCF_DIR/SVs/"$(basename -s .vcf.gz $RAW_SV_VCF)"_CHROM_POS_END.table | awk '{print $3-$2}' > 
 
# 2. Count bases in fasta
GENOME_BP=$(less $GENOME | grep -E '^[ACTGN]+$' | perl -pe 's/[[:space:]]//g' | wc -c)

# 3. Compute number of bp covered by each variant type and proportion of genome covered
#Rscript 01_scripts/utils/count_bases.R $SV_TABLE $VCF_DIR/SNPs/"$(basename -s .vcf.gz $RAW_SNPS_VCF)"_CHROM_POS_END.table $VCF_DIR/indels/"$(basename -s .vcf.gz $RAW_INDELS_VCF)"_CHROM_POS_END.table $GENOME_BP
Rscript 01_scripts/utils/count_bases.R $VCF_DIR/SVs/"$(basename -s .vcf.gz $RAW_SV_VCF)"_CHROM_POS_END.table $VCF_DIR/SNPs/"$(basename -s .vcf.gz $RAW_SNPS_VCF)"_CHROM_POS_END.table $VCF_DIR/indels/"$(basename -s .vcf.gz $RAW_INDELS_VCF)"_CHROM_POS_END.table $GENOME_BP