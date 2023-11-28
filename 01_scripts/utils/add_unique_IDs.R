# Add unique variant IDs to SNPs and indels for easier processing in pooled analysis

library(data.table)

# 1. Access to files from command line and import -------------------------
argv <- commandArgs(T)
CHROM_POS <- argv[1] 
TYPE <- argv[2]

chrom_pos <- as.data.frame(fread(CHROM_POS, col.names = c('CHROM', 'POS')))


# 2. Add required unique ID -----------------------------------------------
chrom_pos$type <- TYPE
chrom_pos$ID <- paste0(TYPE, '_', seq(1, nrow(chrom_pos)))


# 3. Export ---------------------------------------------------------------
write.table(chrom_pos[, c('CHROM', 'POS', 'ID', 'type')], 
            file = paste0(CHROM_POS, '_ID_TYPE.annot'), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


