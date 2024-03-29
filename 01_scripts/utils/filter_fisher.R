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

FST_TABLE <- commandArgs(trailingOnly = TRUE)[8]

fisher_table <- paste0(out_file, '_all.txt')


# 2. Filter for Q_VAL > threshold -----------------------------------------

# Import Fisher test results table
genopos <- read.delim(fisher_table, col.names = c('CHROM', 'POS', 'END', 'ID', 'P_VAL', 'Q_VAL'))

genopos_outliers <- subset(genopos, Q_VAL < MAX_QVAL)

write.table(genopos_outliers[, c('CHROM', 'POS', 'END', 'ID', 'P_VAL', 'Q_VAL')], 
            file = paste0(out_file, '_outliers_qval', MAX_QVAL, '.txt'), sep = "\t",
            row.names = FALSE, quote = FALSE)



# 3. Get min and max Fst values associated with Fisher outliers -----------
FST_table <- read.delim(FST_TABLE, col.names = c('CHROM', 'POS', 'END', 'ID', 'FST'))
FST_table$FST <- ifelse(FST_table$FST < 0, 
                        yes = 0,
                        no = FST_table$FST)
                        
genopos_outliers_Fst <- merge(x = genopos_outliers, y = FST_table,
      by = c('CHROM', 'POS', 'END', 'ID'),
      all = FALSE)

print(paste('Min Fst value associated with outliers with with Q_VAL <', MAX_QVAL, ':', min(genopos_outliers_Fst$FST)))
print(paste('Min Fst value associated with outliers with with Q_VAL <', MAX_QVAL, ':', max(genopos_outliers_Fst$FST)))
