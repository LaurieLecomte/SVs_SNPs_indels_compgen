# Check that outliers and RDA candidates are not missing data outliers as well

library(ggplot2)


# SVs ---------------------------------------------------------------------
ALL_SV <- "~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SVs/merged_SUPP2_MAF0.05_FMISS0.5.SVsFst_RO_PU_FST_FMISS.table"
OUTLIER_SV <- "~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SVs/SVs_RO_PU_outliers_minFst0.243_qval0.01_shared.table"
RDA_SV <- "~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/10_rda/SVs/RDA_3sd_outliers.table"

# Import and format
SV_table <- read.delim(ALL_SV, header=FALSE,
                       col.names = c('CHROM', 'POS', 'ID', 'END', 'FST', 'F_MISS'))

SV_outliers <- read.delim(OUTLIER_SV, header = FALSE, col.names = c('CHROM', 'POS', 'END', 'ID'))

## Add tag 'outlier' to outlier SVs
SV_outliers$type <- 'outliers'
SV_outliers <- merge(SV_outliers, SV_table, by = c('CHROM', 'POS', 'END', 'ID'), all.x = TRUE)

SV_RDA <- read.delim(RDA_SV, header = FALSE, col.names = c('CHROM', 'POS', 'END', 'ID', 'loadings'))

## Add tag 'RDA cand' to RDA candidate SVs
SV_RDA$type <- 'RDA candidates'
SV_RDA <- merge(SV_RDA, SV_table, by = c('CHROM', 'POS', 'END', 'ID'), all.x = TRUE)

# Merge the 3 datasets
SV_all <- merge(SV_RDA, SV_outliers, by = c('CHROM', 'POS', 'END', 'ID', 'type', 'FST', 'F_MISS'), all = TRUE)

SV_all <- merge(SV_table, SV_all, by = c('CHROM', 'POS', 'END', 'ID', 'FST', 'F_MISS'), all = TRUE)

SV_all$type <- 
ifelse(test = is.na(SV_all$type),
                    yes = 'other SVs',
       no = SV_all$type)

# Plot
ggplot(data = SV_all) +
  geom_boxplot(aes(x = type, y = F_MISS))


# SNPs --------------------------------------------------------------------
ALL_SNP <- "~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SNPs/SNPs_MAF0.05_FMISS0.5.SNPsFst_RO_PU_FST_FMISS.table"
OUTLIER_SNP <- "~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SNPs/SNPs_RO_PU_outliers_minFst0.286_qval0.01_shared.table"
RDA_SNP <- "~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/10_rda/SNPs/RDA_3sd_outliers.table"

# Import and format
SNP_table <- read.delim(ALL_SNP, header=FALSE,
                       col.names = c('CHROM', 'POS', 'ID', 'END', 'FST', 'F_MISS'))

SNP_outliers <- read.delim(OUTLIER_SNP, header = FALSE, col.names = c('CHROM', 'POS', 'END', 'ID'))

## Add tag 'outlier' to outlier SNPs
SNP_outliers$type <- 'outliers'
SNP_outliers <- merge(SNP_outliers, SNP_table, by = c('CHROM', 'POS', 'END', 'ID'), all.x = TRUE)

SNP_RDA <- read.delim(RDA_SNP, header = FALSE, col.names = c('CHROM', 'POS', 'END', 'ID', 'loadings'))

## Add tag 'RDA cand' to RDA candidate SNPs
SNP_RDA$type <- 'RDA candidates'
SNP_RDA <- merge(SNP_RDA, SNP_table, by = c('CHROM', 'POS', 'END', 'ID'), all.x = TRUE)

# Merge the 3 datasets
SNP_all <- merge(SNP_RDA, SNP_outliers, by = c('CHROM', 'POS', 'END', 'ID', 'type', 'FST', 'F_MISS'), all = TRUE)

SNP_all <- merge(SNP_table, SNP_all, by = c('CHROM', 'POS', 'END', 'ID', 'FST', 'F_MISS'), all = TRUE)

SNP_all$type <- 
  ifelse(test = is.na(SNP_all$type),
         yes = 'other SNPs',
         no = SNP_all$type)

# Plot
ggplot(data = SNP_all) +
  geom_boxplot(aes(x = type, y = F_MISS))


# indels --------------------------------------------------------------------
ALL_indel <- "~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/indels/indels_MAF0.05_FMISS0.5.indelsFst_RO_PU_FST_FMISS.table"
OUTLIER_indel <- "~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/indels/indels_RO_PU_outliers_minFst0.323_qval0.01_shared.table"
RDA_indel <- "~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/10_rda/indels/RDA_3sd_outliers.table"

# Import and format
indel_table <- read.delim(ALL_indel, header=FALSE,
                        col.names = c('CHROM', 'POS', 'ID', 'END', 'FST', 'F_MISS'))

indel_outliers <- read.delim(OUTLIER_indel, header = FALSE, col.names = c('CHROM', 'POS', 'END', 'ID'))

## Add tag 'outlier' to outlier indels
indel_outliers$type <- 'outliers'
indel_outliers <- merge(indel_outliers, indel_table, by = c('CHROM', 'POS', 'END', 'ID'), all.x = TRUE)

indel_RDA <- read.delim(RDA_indel, header = FALSE, col.names = c('CHROM', 'POS', 'END', 'ID', 'loadings'))

## Add tag 'RDA cand' to RDA candidate indels
indel_RDA$type <- 'RDA candidates'
indel_RDA <- merge(indel_RDA, indel_table, by = c('CHROM', 'POS', 'END', 'ID'), all.x = TRUE)

# Merge the 3 datasets
indel_all <- merge(indel_RDA, indel_outliers, by = c('CHROM', 'POS', 'END', 'ID', 'type', 'FST', 'F_MISS'), all = TRUE)

indel_all <- merge(indel_table, indel_all, by = c('CHROM', 'POS', 'END', 'ID', 'FST', 'F_MISS'), all = TRUE)

indel_all$type <- 
  ifelse(test = is.na(indel_all$type),
         yes = 'other indels',
         no = indel_all$type)

# Plot
ggplot(data = indel_all) +
  geom_boxplot(aes(x = type, y = F_MISS))
