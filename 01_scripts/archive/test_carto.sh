#!/bin/bash

# srun -p small --mem=200G -J map -o log/test_carto_%j.log /bin/sh 01_scripts/archive/test_carto.sh &

# LOAD REQUIRED MODULES
module load R/4.1

Rscript 01_scripts/archive/tests_carto.R