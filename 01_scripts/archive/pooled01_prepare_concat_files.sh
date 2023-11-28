#!/bin/bash

# POOLED VARIANTS ANALYSIS : detect outliers and candidates on pooled variants set
# First prepare files and concatenante SVs, SNPs and indels

# srun -p medium -c4 --mem=50G --time=3-00:00:00 -J pooled01_prepare_concat_files -o log/pooled01_prepare_concat_files_%j.log /bin/sh 01_scripts/utils/pooled01_prepare_concat_files.sh &


# VARIABLES
SV_FST_VCF="08_angsd_fst/SVs/merged_SUPP2_MAF0.05_FMISS0.5.SVsFst_RO_PU.vcf.gz"
SNP_FST_VCF="08_angsd_fst/SNPs/SNPs_MAF0.05_FMISS0.5.SNPsFst_RO_PU.vcf.gz"
INDEL_FST_VCF="08_angsd_fst/indels/indels_MAF0.05_FMISS0.5.indelsFst_RO_PU.vcf.gz"

POP1='RO'
POP2='PU'

POOLED_DIR="pooled_analysis"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

CPU=1

# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load bedtools/2.30.0
module load R/4.1

# 1. Add TYPE tag to VCFs
# SVs
## Extract required fields
bcftools query -f '%CHROM\t%POS\tSV\n' $SV_FST_VCF > $POOLED_DIR/"$(basename -s .vcf.gz $SV_FST_VCF)".annot

## Compress 
bgzip $POOLED_DIR/"$(basename -s .vcf.gz $SV_FST_VCF)".annot -f --threads $CPU
   
## Index
tabix -s1 -b2 -e2 $POOLED_DIR/"$(basename -s .vcf.gz $SV_FST_VCF)".annot.gz -f 
   
## Prepare header : Add Fst tag
echo -e "##INFO=<ID=TYPE,Number=.,Type=String,Description=Variant type>" > $POOLED_DIR/"$(basename -s .vcf.gz $SV_FST_VCF)".annot.hdr

## Annotate and simplify  
#-a is the annotation file (tabix and bgzip, it needs at least CHROM and POS, -h are the header lines to add, -c are the meaning of the column in the annotation file
bcftools annotate -a $POOLED_DIR/"$(basename -s .vcf.gz $SV_FST_VCF)".annot.gz -h $POOLED_DIR/"$(basename -s .vcf.gz $SV_FST_VCF)".annot.hdr -c CHROM,POS,INFO/TYPE $SV_FST_VCF --threads $CPU | bcftools annotate -x ^INFO/NS,INFO/AF,INFO/MAF,INFO/FST_RO_PU,INFO/TYPE,^FORMAT/GT -Oz --threads $CPU > $POOLED_DIR/"$(basename -s .vcf.gz $SV_FST_VCF)"_annotated.vcf.gz
tabix -p vcf $POOLED_DIR/"$(basename -s .vcf.gz $SV_FST_VCF)"_annotated.vcf.gz -f
echo "done preparing SVs VCF"


# SNPs
## Extract required fields
bcftools query -f '%CHROM\t%POS\n' $SNP_FST_VCF > $POOLED_DIR/"$(basename -s .vcf.gz $SNP_FST_VCF)".CHROM_POS
Rscript 01_scripts/utils/add_unique_IDs.R $POOLED_DIR/"$(basename -s .vcf.gz $SNP_FST_VCF)".CHROM_POS "SNP"

## Compress 
bgzip $POOLED_DIR/"$(basename -s .vcf.gz $SNP_FST_VCF)".CHROM_POS_ID_TYPE.annot -f --threads $CPU
   
## Index
tabix -s1 -b2 -e2 $POOLED_DIR/"$(basename -s .vcf.gz $SNP_FST_VCF)".CHROM_POS_ID_TYPE.annot.gz -f 
   
## Prepare header : Add Fst tag
echo -e "##INFO=<ID=TYPE,Number=.,Type=String,Description=Variant type>" > $POOLED_DIR/"$(basename -s .vcf.gz $SNP_FST_VCF)".annot.hdr

## Annotate and simplify  
#-a is the annotation file (tabix and bgzip, it needs at least CHROM and POS, -h are the header lines to add, -c are the meaning of the column in the annotation file
bcftools annotate -a $POOLED_DIR/"$(basename -s .vcf.gz $SNP_FST_VCF)".CHROM_POS_ID_TYPE.annot.gz -h $POOLED_DIR/"$(basename -s .vcf.gz $SNP_FST_VCF)".annot.hdr -c CHROM,POS,ID,INFO/TYPE $SNP_FST_VCF --threads $CPU | bcftools annotate -x ^INFO/NS,INFO/AF,INFO/MAF,INFO/FST_RO_PU,INFO/TYPE,^FORMAT/GT -Oz --threads $CPU > $POOLED_DIR/"$(basename -s .vcf.gz $SNP_FST_VCF)"_annotated.vcf.gz
tabix -p vcf $POOLED_DIR/"$(basename -s .vcf.gz $SNP_FST_VCF)"_annotated.vcf.gz -f
echo "done preparing SNPs VCF"


# Indels
## Extract required fields
bcftools query -f '%CHROM\t%POS\n' $INDEL_FST_VCF > $POOLED_DIR/"$(basename -s .vcf.gz $INDEL_FST_VCF)".CHROM_POS
Rscript 01_scripts/utils/add_unique_IDs.R $POOLED_DIR/"$(basename -s .vcf.gz $INDEL_FST_VCF)".CHROM_POS "indel"

## Compress 
bgzip $POOLED_DIR/"$(basename -s .vcf.gz $INDEL_FST_VCF)".CHROM_POS_ID_TYPE.annot -f --threads $CPU
   
## Index
tabix -s1 -b2 -e2 $POOLED_DIR/"$(basename -s .vcf.gz $INDEL_FST_VCF)".CHROM_POS_ID_TYPE.annot.gz -f
   
## Prepare header : Add Fst tag
echo -e "##INFO=<ID=TYPE,Number=.,Type=String,Description=Variant type>" > $POOLED_DIR/"$(basename -s .vcf.gz $INDEL_FST_VCF)".annot.hdr

## Annotate and simplify  
#-a is the annotation file (tabix and bgzip, it needs at least CHROM and POS, -h are the header lines to add, -c are the meaning of the column in the annotation file
bcftools annotate -a $POOLED_DIR/"$(basename -s .vcf.gz $INDEL_FST_VCF)".CHROM_POS_ID_TYPE.annot.gz -h $POOLED_DIR/"$(basename -s .vcf.gz $INDEL_FST_VCF)".annot.hdr -c CHROM,POS,ID,INFO/TYPE $INDEL_FST_VCF --threads $CPU | bcftools annotate -x ^INFO/NS,INFO/AF,INFO/MAF,INFO/FST_RO_PU,INFO/TYPE,^FORMAT/GT -Oz  --threads $CPU > $POOLED_DIR/"$(basename -s .vcf.gz $INDEL_FST_VCF)"_annotated.vcf.gz
tabix -p vcf $POOLED_DIR/"$(basename -s .vcf.gz $INDEL_FST_VCF)"_annotated.vcf.gz -f
echo "done preparing indels VCF"


# 2. Concatenate VCFs
echo -e "$POOLED_DIR/"$(basename -s .vcf.gz $SV_FST_VCF)"_annotated.vcf.gz\n$POOLED_DIR/"$(basename -s .vcf.gz $SNP_FST_VCF)"_annotated.vcf.gz\n$POOLED_DIR/"$(basename -s .vcf.gz $INDEL_FST_VCF)"_annotated.vcf.gz" > 02_infos/VCFs_annot_list.txt

bcftools concat -a -f 02_infos/VCFs_annot_list.txt -Oz --threads $CPU > $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz
tabix -p vcf $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz

# 5. Convert pooled VCF to table
bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz > $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_CHROM_POS_END_ID

zless $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz | grep -v ^# | perl -pe 's/.+FST_RO_PU=([\-\.0-9a-zA-Z]+);TYPE=([A-Za-z]+).+/\1\t\2/' > $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_FST_TYPE

paste -d "\t" $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_CHROM_POS_END_ID $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_FST_TYPE > $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.table


# 3. Convert pooled VCF to genotype matrix
vcftools --gzvcf $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz --012 --out $POOLED_DIR/"$(basename -s .vcf.gz $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz)".geno_mat

# 4. Impute missing data
Rscript 01_scripts/utils/pooled_impute_missing.R $POOLED_DIR/"$(basename -s .vcf.gz $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz)".geno_mat.012 $ID_SEX_POP $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_CHROM_POS_END_ID

echo "imputation done"


# 5. Convert pooled VCF to table
#bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz > $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_CHROM_POS_END_ID

#zless $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz | grep -v ^# | perl -pe 's/.+FST_RO_PU=([\-\.0-9a-zA-Z]+);TYPE=([A-Za-z]+).+/\1\t\2/' > $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_FST_TYPE

#paste -d "\t" $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_CHROM_POS_END_ID $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_FST_TYPE > $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.table

# Clean up #
#rm $POOLED_DIR/*.annot.gz
#rm $POOLED_DIR/*.annot.hdr
##rm $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_CHROM_POS_END_ID
#rm $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_FST_TYPE