# Filter Fisher tests table

library(data.table)
library(ggplot2)

# 1. Access files in command line, import and format ----------------------
genofile <- commandArgs(trailingOnly = TRUE)[1]
popfile <- commandArgs(trailingOnly = TRUE)[2]
CHR_POS_END_ID <- commandArgs(trailingOnly = TRUE)[3]

pop1 <- commandArgs(trailingOnly = TRUE)[4]
pop2 <- commandArgs(trailingOnly = TRUE)[5]
MAX_QVAL <- commandArgs(trailingOnly = TRUE)[6]

out_file <- commandArgs(trailingOnly = TRUE)[7]

fisher_table <- paste0(out_file, '_all.txt')


# 2. Filter for Q_VAL > threshold -----------------------------------------

# Import Fisher test results table
genopos <- read.delim(fisher_table, col.names = c('CHROM', 'POS', 'END', 'ID', 'P_VAL', 'Q_VAL'))

genopos_outliers <- subset(genopos, Q_VAL < MAX_QVAL)

write.table(genopos_outliers[, c('CHROM', 'POS', 'END', 'ID', 'P_VAL', 'Q_VAL')], 
            file = paste0(out_file, '_outliers_qval', MAX_QVAL, '.txt'), sep = "\t",
            row.names = FALSE, quote = FALSE)

