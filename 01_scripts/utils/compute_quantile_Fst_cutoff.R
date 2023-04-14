# Compute Fst quantile

#FST_TABLE <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SVs/merged_SUPP2_MAF0.05_FMISS0.5.maf0.05.SVsFst_RO_PU.SVsFst_RO_PU.table"

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
FST_TABLE <- argv[1]
CUTOFF <- as.numeric(argv[2])


Fst_table <- read.delim(FST_TABLE, header = FALSE, 
                        col.names = c('CHROM', 'POS', 'END', 'ID', 'FST'))
Fst_table$FST <- as.numeric(Fst_table$FST)

# Set negative values as 0
Fst_table$FST <- ifelse(Fst_table$FST < 0, yes = 0, no = Fst_table$FST)


# 2. Compute Fst cutoff based on quantile ---------------------------------
fst_cutoff <- (quantile(Fst_table$FST, CUTOFF, na.rm = TRUE))
print(fst_cutoff)

writeLines(text = as.character(unname(fst_cutoff)), con = paste0(unlist(strsplit(FST_TABLE, split = '.table')), '_quantile', CUTOFF, '_Fst_cutoff.txt'))

