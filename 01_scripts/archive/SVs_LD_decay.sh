#!/bin/bash

# Test LD decay on SVs

# Manitou
# srun -p small -c 1 --time=1:30:00 -J SVs_LD_decay -o log/SVs_LD_decay_%j.log /bin/sh 01_scripts/SVs_LD_decay.sh  &


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
#ANNOT_DIR="11_annotation/SVs"
GO_DIR="12_go/SVs"
FISHER_DIR="11_fisher_tests/SVs"

#GENOME_ANNOT="11_annotation/genome_annotation/genome_annotation_table_simplified_1.5.tsv"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

SV_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".recoded.vcf.gz"

N_IND="$(less $ID_SEX_POP | wc -l)"

MIN_MAF=0.05
MAX_MAF=0.95
MAX_MISS=0.5

#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"
#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".vcf.gz"

POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".SVsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script

LD_DIR="LD/SVs"
CHR_LIST="02_infos/chr_list.txt"

# LOAD REQUIRED MODULES
module load PopLDecay/3.41
module load vcftools/0.1.16
module load bcftools/1.13

# 0. Create output dir
if [[ ! -d $LD_DIR ]]
then
  mkdir $LD_DIR
fi


# Try on VCF that was reformated for ANGSD (e.g. SVs were recoded as SNPs)
#PopLDdecay -InVCF $SV_VCF_ANGSD -OutFilterSNP -OutStat $LD_DIR/"$(basename -s .vcf.gz $SV_VCF_ANGSD)".allsamples.stat -MAF $MIN_MAF -Miss $MAX_MISS

#vcftools --singletons --gzvcf "/project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SVs/merged_SUPP2_MAF0.05_FMISS0.5.recoded.RO.vcf.gz" --out $LD_DIR/"$(basename -s .vcf.gz $SV_VCF_ANGSD)".singletons.vcf.gz # no singletons

# First split by chromosomes
bcftools index $SV_VCF_ANGSD

less $CHR_LIST | head -n3 | while read chr; 
do
  # Subset VCF by chromosomes
  bcftools filter -r $chr $SV_VCF_ANGSD -Oz > $LD_DIR/"$(basename -s .vcf.gz $SV_VCF_ANGSD)"."$chr".vcf.gz
  # Run PopLDdecay
  PopLDdecay -OutFilterSNP -InVCF $LD_DIR/"$(basename -s .vcf.gz $SV_VCF_ANGSD)"."$chr".vcf.gz -OutStat $LD_DIR/"$(basename -s .vcf.gz $SV_VCF_ANGSD)"."$chr".RO.stat -SubPop 02_infos/RO_IDs.txt -MAF $MIN_MAF -Miss $MAX_MISS -Het 1 -OutType 6
done


# On RO only
#PopLDdecay -OutFilterSNP -InVCF $SV_VCF_ANGSD -OutStat $LD_DIR/"$(basename -s .vcf.gz $SV_VCF_ANGSD)".RO.stat -SubPop 02_infos/RO_IDs.txt -MAF $MIN_MAF -Miss $MAX_MISS -Het 1 -OutType 2

#  plot r2 as a function of distance in order to observe the decay of LD.


#Usage: PopLDdecay -InVCF  <in.vcf.gz>  -OutStat <out.stat>

#		-InVCF       <str>    Input SNP VCF Format
#		-InGenotype  <str>    Input SNP Genotype Format
#		-OutStat     <str>    OutPut Stat Dist ~ r^2 File

#		-SubPop      <str>    SubGroup SampleList of VCFFile [ALLsample]
#		-MaxDist     <int>    Max Distance (kb) between two SNP [300]
#		-MAF         <float>  Min minor allele frequency filter [0.005]
#		-Het         <float>  Max ratio of het allele filter [0.88]
#		-Miss        <float>  Max ratio of miss allele filter [0.25]
#		-EHH         <str>    To Run EHH Region decay set StartSite [NA]
#		-OutFilterSNP         OutPut the final SNP to calculate
#		-OutType     <int>    1: R^2 result 2: R^2 & D' result 3:PairWise LD Out[1]
#		                      See the Help for more OutType [1-8] details
		
#		-help                 Show more help [hewm2008 v3.42]
   
   
   
# Try on final, filtered VCF
#PopLDdecay -InVCF $RAW_FST_VCF -OutStat $LD_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".allsamples.stat