#!/bin/bash

# Format input SVs VCF for running ANGDS

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

# LOAD REQUIRED MODULES
module load bcftools/1.13


## Extract required info
bcftools query -f '%CHROM %POS %REF %ALT\n' $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.vcf > $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.variants
bcftools query -f '%CHROM\t%POS\n' $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.vcf > $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.chrpos

Rscript 01_scripts/utils/extract_SV_info_graph_LL.R $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.variants

# 2. Add flanking sequence to SV
#use Eric script (adapted by myself) to extract ref allele at position for dummy vcf - which script ? # let Claire know that it is not on GitHub

python3 01_scripts/utils/fasta_extract_flanking_regions_claire.py $GENOME $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.chrpos 1 $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.variants.ref

grep ^"#" $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.vcf > $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.header

# 4. Extract VCF content without header
grep -v ^"#" $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.vcf > $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.noheader

#edit in REF
Rscript 01_scripts/utils/make_dummy_vcf_snp.r $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.noheader $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.variants.ref $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.header

#format vcf
cat $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.header $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.noheader.withdummyREFALT > $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".forangsd.vcf.tmp
