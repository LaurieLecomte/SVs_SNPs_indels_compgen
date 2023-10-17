
# 0. Access to files provided in command line arguments -------------------
argv <- commandArgs(T)
VCF <- argv[1]            # original, unformatted vcf, can be gzipped
OUT_VCF <- argv[2]  # output vcf, will be overwritten
MAX_WIDTH <- argv[3]         # minimum SVLEN for which we need to use middle position between POS and END instead of POS

# 1. Source required functions --------------------------------------------
source('01_scripts/utils/recode_SV_POS.R')

# 2. Process VCF ----------------------------------------------------

recode_SVs(VCF, 
           OUT_VCF, 
           MAX_WIDTH)