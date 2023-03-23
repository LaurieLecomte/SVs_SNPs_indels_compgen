#!/bin/bash

# Estimate genome-wide and per-site Fst by populations pairs

# manitou
# srun -p small -c 4 --mem=20G -J 02.5_SNPs_fst -o log/02.5_SNPs_fst_%j.log /bin/sh 01_scripts/02.5_SNPs_fst.sh &

# valeria
# srun -p ibis_small -c 4 --mem=20G -J 02.5_SNPs_fst -o log/02.5_SNPs_fst_%j.log /bin/sh 01_scripts/02.5_SNPs_fst.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"

RAW_VCF_DIR="04_vcf/SNPs"
RAW_SNPS_VCF="$RAW_VCF_DIR/SNPs_MAF0.05_FMISS0.5.vcf.gz"

ANGSD_INPUT_DIR="05_angsd_inputs/SNPs"
ANGSD_STATS_DIR="06_angsd_stats/SNPs"
ANGSD_BYPOP_DIR="07_angsd_bypop/SNPs"
ANGSD_FST_DIR="08_angsd_fst/SNPs"

PCA_DIR="09_pca/SNPs"
RDA_DIR="10_rda/SNPs"
ANNOT_DIR="11_annotation/SNPs"
GO_DIR="12_go/SNPs"

GENOME_ANNOT="11_annotation/genome_annotation/genome_annotation_table_simplified_1.5.tsv"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

SNPS_VCF_ANGSD=$RAW_SNPS_VCF

N_IND="$(less $ID_SEX_POP | wc -l)"

CPU=4

MIN_MAF=0.05
MAX_MAF=0.95

POPS_FILE="02_infos/pops.txt"

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

REALSFS_PATH="/prg/angsd/0.937/misc/realSFS"

# LOAD REQUIRED MODULES
module load python/3.7
module load angsd/0.937
module load pcangsd/1.10
module load bcftools/1.13

# VARIABLES FOR ANGSD
WINDOW=1000000
WIN_STEP=10000 


# 2. Calculate pairwise Fst for each populations pair, genome-wide and by site
NUM_POPS=$(less $POPS_FILE | wc -l)

# Estimate pairwise FST for all populations listed : one group = one pop pair
for i in $(seq $NUM_POPS);
do
	pop1=$(cat "$POPS_FILE" | head -"$i" | tail -1)
	for j in $(seq $[ $i + 1 ] $NUM_POPS)
	do
   
		pop2=$(cat "$POPS_FILE" | head -"$j" | tail -1)
		echo "FST between $pop1 and $pop2"
		
   # Create directory for given group in FST directory
   GROUP=$(echo "$pop1"_"$pop2")
   if [[ ! -d $ANGSD_FST_DIR/$GROUP ]]
    then
    mkdir $ANGSD_FST_DIR/$GROUP
   fi
   
   # 2.1 Compute 2dsfs priors
   echo "calcualte 2dsfs between $pop1 and $pop2"
   
   $REALSFS_PATH $ANGSD_BYPOP_DIR/"$pop1".saf.idx $ANGSD_BYPOP_DIR/"$pop2".saf.idx \
   -P $CPU -maxIter 30 > $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".prior
   
   #file="$ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".prior"
   
   Rscript 01_scripts/utils/sum_sites_2dsfs.R "$ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".prior"
   
   # 2.2 Pre-Fst steps
   echo "Pre-Fst steps"
   ## warning : realSFS outputs will be reordered by numerical order, not the order in original VCF/reference fasta
   $REALSFS_PATH fst index $ANGSD_BYPOP_DIR/"$pop1".saf.idx $ANGSD_BYPOP_DIR/"$pop2".saf.idx \
   -sfs $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".prior.2dsfs \
   -P $CPU -fstout $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2"
   
   echo "calculate SFS priori for each position"
   $REALSFS_PATH fst print $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".fst.idx -P $CPU > $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs
   
   # 2.3 Compute genome-wide Fst
   echo "get genome-wide estimate of Fst"
   $REALSFS_PATH fst stats $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".fst.idx -P $CPU > $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".SNPs.fst
   
   # 2.4 Compute Fst by sliding window
   $REALSFS_PATH fst stats2 $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".fst.idx -win $WINDOW -step $WIN_STEP -P $CPU > $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2"_win"$WINDOW"_step"$WIN_STEP".txt
   
   # 2.5 Calculate per site (= per indel) Fst
   Rscript 01_scripts/utils/per_site_Fst.R $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.bed
   
   # 2.6. Add per site Fst to input VCF
   ## Extract relevant fields from bed
   cut -f1,2,5 $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.bed > $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot
   
   ## Add field for END
   #bcftools query -f '%END\n' $FILT_VCF_FILE > $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.end
   #paste -d "\t" $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.tmp $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.end > $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot
   
   ## Compress 
   bgzip $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot -f
   
   ## Index
   tabix -s1 -b2 -e2 $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.gz -f
   
   ## Prepare header : Add Fst tag
   echo -e "##INFO=<ID=FST_"$pop1"_"$pop2",Number=.,Type=Float,Description=\"Per site Fst between "$pop1" and "$pop2"\">" > $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.hdr
   
   
   ## Annotate 
   #-a is the annotation file (tabix and bgzip, it needs at least CHROM and POS, -h are the header lines to add, -c are the meaning of the column in the annotation file
   bcftools annotate -a $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.gz -h $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.hdr -c CHROM,POS,INFO/FST_"$pop1"_"$pop2" $FILT_ANGSD_VCF -Oz --threads $CPU > "$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".SNPsFst_"$pop1"_"$pop2".vcf.gz"
   
   # 6. Add per site Fst to the input VCF that has NOT been formatted for angsd
   bcftools annotate -a $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.gz -h $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.hdr -c CHROM,POS,INFO/FST_"$pop1"_"$pop2" $RAW_SNPS_VCF -Oz --threads $CPU > "$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_SNPS_VCF)".SNPsFst_"$pop1"_"$pop2".vcf.gz"
   
   # Clean up 
   rm $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.hdr
   #rm $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.gz
   rm $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.gz.tbi
   #rm $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.end
   #rm $ANGSD_FST_DIR/$GROUP/"$pop1"_"$pop2".bypos.sfs.annot.tmp
   
   done;
   
done

