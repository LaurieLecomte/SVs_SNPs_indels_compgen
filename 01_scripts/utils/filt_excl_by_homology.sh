#!/bin/bash

# Get overlap between good (filtered) and bad (filtered out) SVs and SNPs to explore regions and genomic features leading to poor-quality SVs
# This script produces files required by the 01_scripts/utils/filt_excl_by_homology.R script

# srun -p small --time=1-00:00:00 -c 1 -J filt_excl_by_homology.sh -o log/filt_excl_by_homology_%j.log /bin/sh 01_scripts/utils/filt_excl_by_homology.sh &

# VARIABLES
SV_FILT="/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5.vcf.gz"
SV_EXCL="/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_excluded.vcf.gz"

SNP_FILT="/project/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/SNPs/SNPs_MAF0.05_FMISS0.5.vcf.gz"
SNP_EXCL="/project/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/SNPs.vcf.gz"
#INDEL_FILT
#INDEL_EXCL="/project/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/indels.vcf.gz"

HOM_BLOCKS="hom_regions/homolog_blocks_identity_75.bed" # output of 01_scripts/utils/homology_by_win.R
REPEATS="hom_regions/genome.fasta.out.table"            # output of 01_scripts/utils/RepeatMasker_manitou.sh

WIN_SIZE=1000                                           # window size used for computing variant density
CHR_BED="02_infos/chrs.bed"                             # list of chromosomes


# LOAD REQUIRED MODULES
module load bedtools
module load bcftools/1.13
module load bedops

#bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $SV_FILT > hom_regions/"$(basename -s .vcf.gz $SV_FILT)".bed #115907 SVs
#bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $SV_EXCL > hom_regions/"$(basename -s .vcf.gz $SV_EXCL)".bed #228052 SVs, a few are missing from the unfiltered 344468 SVs ?


# Check overlap between SVs and known syntenic regions, to know if each SV falls in such regions or not and the % homology
#bedtools intersect -a hom_regions/"$(basename -s .vcf.gz $SV_FILT)".bed -b $HOM_BLOCKS -wao | cut -f1-4,8,9 > hom_regions/intersect_filt_SVs_hom_regions.bed

#bedtools intersect -a hom_regions/"$(basename -s .vcf.gz $SV_FILT)".bed -b $HOM_BLOCKS -wao > hom_regions/intersect_filt_SVs_hom_regions.table

#bedtools intersect -a hom_regions/"$(basename -s .vcf.gz $SV_EXCL)".bed -b $HOM_BLOCKS -wao > hom_regions/intersect_excl_SVs_hom_regions.table

# Check overlap between SVs and known REPEATS (TEs and reps), to know if each SV falls in such regions or not
#bedtools intersect -a hom_regions/"$(basename -s .vcf.gz $SV_FILT)".bed -b $REPEATS -wao > hom_regions/intersect_filt_SVs_RM.table
#bedtools intersect -a hom_regions/"$(basename -s .vcf.gz $SV_EXCL)".bed -b $REPEATS -wao > hom_regions/intersect_excl_SVs_RM.table




# 1. Prepare bed file of filtered and excluded SVs 
bcftools query -f '%CHROM\t%POS0\t%END\t%ID\tfiltered\n' $SV_FILT > hom_regions/"$(basename -s .vcf.gz $SV_FILT)".bed
bcftools query -f "%CHROM\t%POS0\t%END\t%ID\texcluded\n" $SV_EXCL > hom_regions/"$(basename -s .vcf.gz $SV_EXCL)".bed

cat hom_regions/"$(basename -s .vcf.gz $SV_FILT)".bed hom_regions/"$(basename -s .vcf.gz $SV_EXCL)".bed > hom_regions/filtered_excluded.bed

# 2. Get overlap of good and bad SVs with known syntenic regions, to know if each SV falls in such regions or not and the % homology, allowing a 100 pb window around each SV
bedtools window -a hom_regions/filtered_excluded.bed -b $HOM_BLOCKS -w 100 > hom_regions/win100_filt_excl_SVs_hom_regions.table

# 3. Get overlap between SVs and known REPEATS (TEs and reps), to know if each SV falls in such regions or not
bedtools window -a hom_regions/filtered_excluded.bed -b $REPEATS -w 100 > hom_regions/win100_filt_excl_SVs_RM.table

# 4. Produce a bed file of WINDOW, which will be used for computing variant density in the script 
bedops --chop $WIN_SIZE -x $CHR_BED > 02_infos/chrs_win"$WIN_SIZE".bed




#bedtools intersect -a hom_regions/filtered_excluded.bed -b $HOM_BLOCKS -wa -wb -f 0.01 > hom_regions/intersect_filt_excl_SVs_hom_regions.table

#bedtools intersect -a hom_regions/filtered_excluded.bed -b "hom_regions/genome.fasta.out.table" -wa -wb -f 0.01 > hom_regions/intersect_filt_excl_SVs_RM.table


#bedtools window -a hom_regions/filtered_excluded.bed -b $HOM_BLOCKS -w 100 > hom_regions/win100_filt_excl_SVs_hom_regions.table

#bedtools window -a hom_regions/filtered_excluded.bed -b "hom_regions/genome.fasta.out.table" -w 100 > hom_regions/win100_filt_excl_SVs_RM.table


#bedops --chop $WIN_SIZE -x $CHR_BED > 02_infos/chrs_win"$WIN_SIZE".bed



# For SNPs
# 1. Prepare bed file of filtered and excluded SVs 
bcftools query -f '%CHROM\t%POS0\t%END\t%ID\tfiltered\n' $SNP_FILT > hom_regions/"$(basename -s .vcf.gz $SNP_FILT)".bed
bcftools query -f "%CHROM\t%POS0\t%END\t%ID\texcluded\n" $SNP_EXCL > hom_regions/"$(basename -s .vcf.gz $SNP_EXCL)".bed 

cat hom_regions/"$(basename -s .vcf.gz $SNP_FILT)".bed hom_regions/"$(basename -s .vcf.gz $SNP_EXCL)".bed > hom_regions/filtered_excluded_SNPs.bed

# 2. Get overlap of good and bad SNPs with known syntenic regions, to know if each SNP falls in such regions or not and the % homology, allowing a 100 pb window around each SNP
bedtools window -a hom_regions/filtered_excluded_SNPs.bed -b $HOM_BLOCKS -w 100 > hom_regions/win100_filt_excl_SNPs_hom_regions.table

# 3. Get overlap between SNPs and known REPEATS (TEs and reps), to know if each SNP falls in such regions or not
bedtools window -a hom_regions/filtered_excluded_SNPs.bed -b $REPEATS -w 100 > hom_regions/win100_filt_excl_SNPs_RM.table








# Now get density of SVs falling in a syntenic bloc, for low, elevated, high and very high % identity, and density of SVs not in a syntenic block 
# NOT USED, because we did it in R
# 2. Split genome by window
bedops --chop $WIN_SIZE -x $CHR_BED > 02_infos/chrs_win"$WIN_SIZE".bed

# 3. Get number of overlap between fixed size windows and sites 
bedtools intersect -a 02_infos/chrs_win"$WIN_SIZE".bed -b hom_regions/intersect_filt_SVs_hom_regions.bed -c > hom_regions/intersect_filt_SVs_hom_regions_density_win"$WIN_SIZE"bp.txt

bedtools intersect -a 02_infos/chrs_win"$WIN_SIZE".bed -b hom_regions/intersect_excl_SVs_hom_regions.bed -c > hom_regions/intersect_excl_SVs_hom_regions_density_win"$WIN_SIZE"bp.txt
