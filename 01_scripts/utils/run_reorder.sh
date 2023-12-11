#!/bin/bash

for i in $(ls /project/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SVs_vs_SNPs/SVs_SNPs_concat_100_*_SVTYPE_LEN.ld)
do
    python3 00_Script/reorder_LD.py $i 01_ld_files/$(basename $i)_ordered
done
