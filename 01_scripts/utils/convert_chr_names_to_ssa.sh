#!/bin/bash

# Get correspondance between full chromosome names (OV..) and short names (ssaX)

# VARIABLES
GENOME="03_genome/genome.fasta"


less $GENOME | grep ^'>OV' | sed -E "s/^>(OV[0-9\.]+).+chromosome:\ ([A-Za-z0-9\-]+).+/\1\t\2/" > 02_infos/OV_to_ssa.txt