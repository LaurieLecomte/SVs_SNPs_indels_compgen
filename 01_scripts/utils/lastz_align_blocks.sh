#!/bin/bash

# srun -p medium --time=7-00:00:00 -c 1 -J lastz_align_blocks -o log/lastz_align_blocks_%j.log /bin/sh 01_scripts/lastz_align_blocks.sh &

# VARIABLES
BLOCKS="blocks"  # blocks file outputted by symap
REF="04_RepeatMasker_manitou/genome_chrs_corrected.masked_addedN.fasta"

BLOCK_FILE=blocks_different_chromosome
BLOCK_FOLDER="$BLOCK_FILE"_chain_folder


# LOAD REQUIRED MODULES
module load bedtools/2.30.0
module load samtools/1.13
module load lastz

## Adapted by Xavier Dallaire from Claire Mérot

#For each of #the 55 links between homeologous WFS blocks identified in symap,
#lastz version 1.02 (Harris 2007) was run in both directions (using the
#parameters: --gfextend --nochain --nogapped --matchcount = 100;
#similarly to Lien et al., 2016) to align these regions to one another
#and subsequently determine sequence similarity between the two.
#Following lastz alignment, matches were filtered to remove those
#with sequence similarity < 75% (in keeping with Lien et al., 2016)
#and/or smaller than 1,000 bp, and sequence similarity was then averaged
#across alignments within each block.

#All the code below is adapted from DeKayne et al 2020
#https://github.com/RishiDeKayne/WhitefishReferenceGenome/blob/master/lastz_submitscript_extend2.lsf

#we use the output from symap describing the position of synteny lockx
#01	16	1	68920136	72801554	40057108	44059771	181	0	0	0.0	0.0

#command
#lastz target[[start..end]] query[[start..end]] --gfextend --nochain --nogapped --matchcount = 100

# 1. Create fasta by chromosome
## First remove lines for intrachr blocks (first skip header)
#awk '{ if ($1 != $2) print $0 }' blocks > $BLOCK_FILE
tail -n+2 $BLOCKS | awk '{ if ($1 != $2) print $0 }' > $BLOCK_FILE

## Create directory for storing each fasta
if [[ ! -d reference ]]
then
  mkdir reference
fi

#fasta=reference/genome_hardMasked_noambiguous.fa

## Index reference
#samtools faidx $fasta
samtools faidx $REF
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $REF.fai > reference/"$(basename -s .fasta $REF)".bed


#for i in $(cat "$BLOCK_FILE" | cut -f 1 | uniq)
for i in $(cat "$BLOCK_FILE" | cut -f 1 | uniq) 
do
  ## Create a bed file for each block
  grep $(echo $i) reference/"$(basename -s .fasta $REF)".bed > reference/ref_$i.bed
  ## Extract corresponding region in reference genome
  bedtools getfasta -fi $REF -bed reference/ref_$i.bed > reference/ref_$i.fasta
done


# 2. Align syntenic blocks 
## Create output folder for results
if [[ ! -d $BLOCK_FOLDER ]]
then
  mkdir $BLOCK_FOLDER
fi

head $BLOCK_FILE
NUM_BLOCKS=$(wc -l "$BLOCK_FILE" | cut -d " " -f 1)

echo "running lastZ on "$NUM_BLOCKS "syntenic blocks"


## Remove averages.out file from previous trials if necessary
if [[ -f "$BLOCK_FOLDER"/averages.out ]]
then
  rm "$BLOCK_FOLDER"/averages.out
fi

## loop over the synteny blocks, give it the chrX.fasta for target and for query and use the start and end
for i in $(seq $NUM_BLOCKS)
do
  ChrT=$(cat "$BLOCK_FILE" | head -"$i" | tail -1 | cut -f 1)
  Tstart=$(cat "$BLOCK_FILE" | head -"$i" | tail -1 | cut -f 4)
  Tend=$(cat "$BLOCK_FILE" | head -"$i" | tail -1 | cut -f 5)
  ChrQ=$(cat "$BLOCK_FILE" | head -"$i" | tail -1 | cut -f 2)
  Qstart=$(cat "$BLOCK_FILE" | head -"$i" | tail -1 | cut -f 6)
  Qend=$(cat "$BLOCK_FILE" | head -"$i" | tail -1 | cut -f 7)
  echo Chr"$ChrT"_"$Tstart"_"$Tend"_Chr"$ChrQ"_"$Qstart"_"$Qend"

  TARGET=reference/ref_"$ChrT".fasta
  QUERY=reference/ref_"$ChrQ".fasta
  segment="$BLOCK_FOLDER"/Chr"$ChrT"_"$Tstart"_"$Tend"_Chr"$ChrQ"_"$Qstart"_"$Qend"

  lastz $TARGET["$Tstart".."$Tend"] $QUERY["$Qstart".."$Qend"] --gfextend --chain --nogapped --filter=nmatch:100 --hspthresh=6000  --format=general > "$segment".out

  #get length of mapping/overlap by splitting column 14 
  cat $segment".out" | awk '{split($14,a,"/"); print a[1]}' > $segment".awktest.out"
  paste $segment".out" $segment".awktest.out" > $segment".withcov.out"
  #head $segment".withcov.out"
  rm $segment".out"
  rm $segment".awktest.out"
  #remove % sign to filter by percentage
  sed 's/%//g' $segment".withcov.out" > $segment".withcov.nopercent.out"
  rm $segment".withcov.out"
  #head $segment".withcov.nopercent.out"
  #filter out alignments <1000 bp long
  #awk '{ if ($16 >= 1000) print $0 }' $segment".withcov.nopercent.out" > $segment".withcov.nopercent_filtered.out"

  #calculate average across whole block 
  #CAUTION this is an average unweighted by the length of the subblocks...?? Shouldn't it be?
  echo Chr"$ChrT"_"$Tstart"_"$Tend"_Chr"$ChrQ"_"$Qstart"_"$Qend" >> "$BLOCK_FOLDER"/averages.out
  cat $segment".withcov.nopercent.out" | awk '{sum+=$13} END { print "Average = ",sum/NR}' >> "$BLOCK_FOLDER"/averages.out
  #head "$BLOCK_FOLDER"/averages.out

done
