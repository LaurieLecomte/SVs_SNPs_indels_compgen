#!/bin/bash


# parallel -a 02_infos/chr_list.txt -j 10 srun -c 1 -p ibis_small --time=1-00:00:00 --mem=10G -J phase_SNPs_LDdecay_{} -o log/phase_SNPs_LDdecay_{}_%j.log /bin/sh 01_scripts/utils/phase_SNPs.sh {} &

# After : 
#for file in $(ls -1 LD/*_phased_RO_PU_stats*); do less $file | tail -n+2 $file >> LD/SNPs_MAF0.05_FMISS0.5_phased_RO_PU_LDstats.txt; done
#for file in $(ls -1 LD/*_phased_RO_stats*); do $file | less tail -n+2 $file >> LD/SNPs_MAF0.05_FMISS0.5_phased_RO_LDstats.txt; done
#for file in $(ls -1 LD/*_phased_PU_stats*); do $file | less tail -n+2 $file >> LD/SNPs_MAF0.05_FMISS0.5_phased_PU_LDstats.txt; done


# VARIABLES
SNPS_VCF="04_vcf/SNPs/SNPs_MAF0.05_FMISS0.5.vcf.gz"

CPU=1
CHR=$1

# LOAD REQUIRED MODULES
module load java
module load beagle
module load bcftools
module load poplddecay

# Split VCF by chrom
bcftools view -r $CHR $SNPS_VCF --threads $CPU -Oz -o 04_vcf/SNPs/"$(basename -s .vcf.gz $SNPS_VCF)"_"$CHR".vcf.gz

# Phase each chromosome
java -jar ${EBROOTBEAGLE}/beagle.22Jul22.46e.jar gt=04_vcf/SNPs/"$(basename -s .vcf.gz $SNPS_VCF)"_"$CHR".vcf.gz out=LD/"$(basename -s .vcf.gz $SNPS_VCF)"_"$CHR"_phased nthreads=$CPU

# Index each chromosome
tabix -p vcf LD/"$(basename -s .vcf.gz $SNPS_VCF)"_"$CHR"_phased.vcf.gz -f

# Run PopLDdecay
## One pop
PopLDdecay -InVCF LD/"$(basename -s .vcf.gz $SNPS_VCF)"_"$CHR"_phased.vcf.gz -OutStat LD/"$(basename -s .vcf.gz $SNPS_VCF)"_"$CHR"_phased_RO_PU_stats

## For each pop
PopLDdecay -InVCF LD/"$(basename -s .vcf.gz $SNPS_VCF)"_"$CHR"_phased.vcf.gz -SubPop 02_infos/RO_IDs.txt -OutStat LD/"$(basename -s .vcf.gz $SNPS_VCF)"_"$CHR"_phased_ROstats

PopLDdecay -InVCF LD/"$(basename -s .vcf.gz $SNPS_VCF)"_"$CHR"_phased.vcf.gz -SubPop 02_infos/PU_IDs.txt -OutStat LD/"$(basename -s .vcf.gz $SNPS_VCF)"_"$CHR"_phased_PUstats



# Split VCF by chrom
#parallel -a 02_infos/chr_list.txt -j $CPU bcftools view -r {} $SNPS_VCF --threads 1 -Oz -o 04_vcf/SNPs/SNPs_MAF0.05_FMISS0.5_"$CHR".vcf.gz

# Phase each chromosome
#parallel -a 02_infos/chr_list.txt -j $CPU java -jar ${EBROOTBEAGLE}/beagle.22Jul22.46e.jar gt=04_vcf/SNPs/SNPs_MAF0.05_FMISS0.5_"$CHR".vcf.gz out=04_vcf/SNPs/SNPs_MAF0.05_FMISS0.5_"$CHR"_phased.vcf.gz nthreads=1

# Index each chromosome
#parallel -a 02_infos/chr_list.txt -j $CPU tabix -p vcf 04_vcf/SNPs/SNPs_MAF0.05_FMISS0.5_"$CHR"_phased.vcf.gz.vcf.gz -f