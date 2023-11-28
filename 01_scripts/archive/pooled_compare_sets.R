library(data.table)

ALL_TABLE <- "pooled_analysis/SVs_SNPs_indels_MAF0.05_FMISS0.5.table"
#FST_TABLE <- "pooled_analysis/SVs_SNPs_indels_fisher_outliers_qval0.01.table"
FISHER_TABLE <- "pooled_analysis/SVs_SNPs_indels_fisher_genopos.txt"
RDA_TABLE <- "pooled_analysis/RDA_3sd_outliers.txt"
DIR <- 'pooled_analysis'

QUANTILE <- 0.97
MAX_QVAL <- 0.01

all_table <- as.data.frame(fread(ALL_TABLE, col.names = c('CHROM', 'POS', 'END', 'ID', 'FST', 'TYPE')))

#fst_table <- read.table(FST_TABLE, col.names = c('CHROM', 'POS', 'END', 'ID', 'FST', 'TYPE'))
RDA_cand <- read.table(RDA_TABLE, col.names = c('CHROM', 'POS', 'END', 'ID', 'loadings'), header = TRUE)
fisher_table <- as.data.frame(fread(FISHER_TABLE, col.names = c('CHROM', 'POS', 'P_VAL', 'Q_VAL'), header = TRUE))


# 0. Get proportions of each type in pooled set ---------------------------
pooled_all_prop <- table(all_table$TYPE)/nrow(all_table)

# 1. Get Fst outliers -----------------------------------------------------
# Set negative values as 0
all_table$FST <- ifelse(all_table$FST < 0, yes = 0, no = all_table$FST)

# Compute Fst cutoff based on quantile
fst_cutoff <- round(quantile(all_table$FST, QUANTILE, na.rm = TRUE), digits = 3)
fst_outliers <- subset(all_table, FST >= fst_cutoff)

# Proportions of each type ?
table(fst_outliers$TYPE)/ nrow(fst_outliers)

# 2. Get Fisher outliers --------------------------------------------------
# First match type
#fisher_table <- merge(all_table, fisher_table, by = c('CHROM', 'POS', 'END', 'ID'),
 #                     all = FALSE)

fisher_table <- cbind(fisher_table, END = all_table$END, ID = all_table$ID, TYPE = all_table$TYPE)
fisher_outliers <- subset(fisher_table, Q_VAL < MAX_QVAL)

# Proportions of each type ?
table(fisher_outliers$TYPE)/ nrow(fisher_outliers)


# 3. Get intersection outliers --------------------------------------------
intersection_outliers <- merge(fst_outliers, fisher_outliers, by = c('CHROM', 'POS', 'END', 'ID', 'TYPE'), 
                           all = FALSE, sort = FALSE)

# Proportions of each type ?
table(intersection_outliers$TYPE)/ nrow(intersection_outliers)
pooled_outliers_prop <- table(intersection_outliers$TYPE)/ nrow(intersection_outliers)
  
write.table(intersection_outliers[, c('CHROM', 'POS', 'END', 'ID', 'Q_VAL', 'FST', 'TYPE')], 
            file = paste0(DIR, '/outliers_Fst_fisher.txt'), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# 4. Get RDA candidates ---------------------------------------------------
RDA_cand <- merge(all_table, RDA_cand, by = c('CHROM', 'POS', 'ID', 'END'),
                  all = FALSE)
# Proportions of each type ?
table(RDA_cand$TYPE)/ nrow(RDA_cand)
pooled_RDA_prop <- table(RDA_cand$TYPE)/ nrow(RDA_cand)

write.table(RDA_cand[, c('CHROM', 'POS', 'END', 'ID', 'FST', 'TYPE')], 
            file = paste0(DIR, '/rda_cand_matched.txt'), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# 5. Shared between outliers and RDA candidates (non uniques variants) --------
shared_3_methods <- merge(intersection_outliers, RDA_cand, by = c('CHROM', 'POS', 'END', 'ID'),
      all = FALSE)

# Remove these shared from sets
intersection_outliers_uniques <- dplyr::anti_join(intersection_outliers, shared_3_methods, by = c('CHROM', 'POS', 'END', 'ID'))

# Proportions of each type ?
table(intersection_outliers_uniques$TYPE)/ nrow(intersection_outliers_uniques)

RDA_uniques <- dplyr::anti_join(RDA_cand, shared_3_methods, by = c('CHROM', 'POS', 'END', 'ID'))

# Proportions of each type ?
table(RDA_uniques$TYPE)/ nrow(RDA_uniques)


# 6. Compute odds ratios --------------------------------------------------

# On full intersection set
OR_outliers <- (table(intersection_outliers$TYPE)/ nrow(intersection_outliers)) / (table(all_table$TYPE)/nrow(all_table))

# On intersection uniques set
OR_outliers_uniques <- (table(intersection_outliers_uniques$TYPE)/ nrow(intersection_outliers_uniques)) / (table(all_table$TYPE)/nrow(all_table))
  (table(RDA_uniques$TYPE)/ nrow(RDA_uniques))

# On full RDA cand set
OR_RDA_cand <- (table(RDA_cand$TYPE)/ nrow(RDA_cand)) / (table(all_table$TYPE)/nrow(all_table))

# On intersection uniques set
OR_RDA_cand_uniques <- (table(RDA_uniques$TYPE)/ nrow(RDA_uniques)) / (table(all_table$TYPE)/nrow(all_table))


# 7. Chi square test on counts -----------------------------------------------
# On all types together
## Outliers
not_outliers <- table(all_table$TYPE) - table(intersection_outliers$TYPE)

mat_outliers <- matrix(data = c(table(intersection_outliers$TYPE), 
                       not_outliers), byrow = FALSE, nrow = 3)

chisq.test(mat_outliers)

## RDA candidates
not_RDA <- table(all_table$TYPE) - table(RDA_cand$TYPE)

mat_RDA <- matrix(data = c(table(RDA_cand$TYPE), 
                       not_RDA), byrow = FALSE, nrow = 3)

chisq.test(mat_RDA)


# On all types together, but not from difference between totSVs and outliersSVs
mat_all_outliers <- matrix(data = c(table(intersection_outliers$TYPE), 
                                table(all_table$TYPE)), 
                       byrow = FALSE, nrow = 3)

chisq.test(mat_all_outliers)





mat_all_type_outliers <- matrix(data = c(table(intersection_outliers$TYPE), 
                                    (nrow(all_table) - table(intersection_outliers$TYPE))
                                    ), 
                           byrow = FALSE, nrow = 3)

chisq.test(mat_all_type_outliers)

# On a per type basis

#obs_YES <- nrow(subset(intersection_outliers, TYPE == "SV")) 
#obs_NO <- nrow(all_table) - obs_YES
#exp_YES <- pooled_all_count[3]
#exp_NO <- nrow(all_table) - pooled_all_count[3]

#mat <- matrix(data = c(exp_YES, exp_NO, obs_YES, obs_NO), byrow = FALSE, nrow = 2)

#fisher.test(mat)

pooled_count <- table(all_table$TYPE)


fisher_counts <- function(type, indice, table){
  obs_YES <- nrow(subset(table, TYPE == type)) # number of outliers variants AND of given type
  obs_NO <- nrow(all_table) - obs_YES # number of variants that are NOT outlier AND of given type
  exp_YES <- pooled_count[indice] # total number of variants of given type, 
                                  # because if nul hypothesis is true,
                                  # prop of outlier of given type should be = to prop of given variant type in pooled set
  exp_NO <- nrow(all_table) - exp_YES
  
  mat <- matrix(data = c(obs_YES, obs_NO, exp_YES, exp_NO), byrow = FALSE, nrow = 2)
  print(mat)
  #fisher.test(mat)
  chisq.test(mat)
}

fisher_counts('SV', 3, intersection_outliers)
fisher_counts('SNP', 2, intersection_outliers)
fisher_counts('indel', 1, intersection_outliers)
fisher_counts('SV', 3, RDA_cand)
fisher_counts('SNP', 2, RDA_cand)
fisher_counts('indel', 1, RDA_cand)



# Chi test of adequcy
chisq.test(x = table(intersection_outliers$TYPE),
                     p = pooled_all_prop)
chisq.test(x = table(RDA_cand$TYPE),
           p = pooled_all_prop)
