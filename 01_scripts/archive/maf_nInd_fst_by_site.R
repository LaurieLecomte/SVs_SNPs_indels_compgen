# Create a list of chrom - pos - freqPOP1 - nIndPOP1 - freqPOP2 - nIndPOP2 - Fst from angsd outputs


# 1. Access files from command line and import ----------------------------
argv <- commandArgs(T)
SITES <- argv[1]
#MAF_POP1 <- argv[2]
#MAF_POP2 <- argv[3]
SITES_FST <- argv[2]
OUTPUT <- argv[3]

#SITES <- '07_angsd_bypop/SVs/merged_SUPP2_MAF0.05_FMISS0.5.maf0.05_chrom_pos_freq_bypop.txt'
#SITES_FST <- "08_angsd_fst/SVs/RO_PU/RO_PU.bypos.sfs.annot.gz"

# Import
chr_pos <- read.table(SITES, header = TRUE, sep = "\t",
                      col.names = c('CHROM', 'POS', 'MAF_POP1', 'N_POP1', 'MAF_POP2', 'N_POP2'))

fst <- read.table(SITES_FST,
                  col.names = c('CHROM', 'POS', 'FST'))



# 2. Merge by CHROM and POS -----------------------------------------------
chr_pos_fst <- merge(x = chr_pos, y = fst, by = c('CHROM', 'POS'))


# 3. Export ---------------------------------------------------------------
write.table(chr_pos_fst, file = paste0(OUTPUT, '.txt'), row.names = FALSE, col.names = TRUE, quote = FALSE)

# Sort and export sorted version
chr_pos_fst <- chr_pos_fst[order(chr_pos_fst$FST, decreasing = TRUE), ]
write.table(chr_pos_fst, file = paste0(OUTPUT, '.sorted.txt'), row.names = FALSE, col.names = TRUE, quote = FALSE)