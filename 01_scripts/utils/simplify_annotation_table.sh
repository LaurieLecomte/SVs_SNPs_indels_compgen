#!/bin/bash

# Remove header from annotation table for easier handling with bedtools

# /bin/sh 01_scripts/utils/simplify_annotation_table.sh

# VARIABLES
GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"

# Remove header and sort table by contig
less $GENOME_ANNOT | tail -n +2 | sort -k1 > $ANNOT_TABLE

# less $GENOME_ANNOT | tail -n +2 | sort -k1 > cut -f1-3,5,7,8,10-14 > $ANNOT_TABLE

