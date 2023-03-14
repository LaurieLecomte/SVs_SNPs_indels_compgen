library(data.table)
library(tidyr)
library(dplyr)

argv <- commandArgs(T)
SFS_FILE <- argv[1] # path to .sfs file
BED_OUTPUT <- argv[2] # output bed path
#TBL_OUTPUT <- argv[3] # output tbl path
#MIN_FST <- argv[4]

# 1. Import sfs
sfs <- fread(SFS_FILE, col.names = c('CHR', 'POS', 'A', 'B'), sep = "\t")

# 2. Calculate per site fst
sfs$fst <- (sfs$A / sfs$B)

# 3. Export
write.table(sfs, file = BED_OUTPUT, row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
