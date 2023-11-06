# Check that outliers and RDA candidates are not missing data outliers as well

library(ggplot2)
library(data.table)



# All variants together, RAW variants -------------------------------------
# Import and format
RAW_SNP <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/SNPs_DP_FMISS_tag.table"

raw_SNPs <- fread(RAW_SNP, col.names = c('CHROM', 'POS', 'DP', 'F_MISS'))
raw_SNPs$var_type <- 'SNP'

RAW_INDEL <-  "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/06_merged/indels_DP_FMISS_tag.table"
raw_indels <- fread(RAW_INDEL, col.names = c('CHROM', 'POS', 'DP', 'F_MISS'))
raw_indels$var_type <- 'indel'

RAW_SV <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/08_merged/merged_SUPP2_genotyped.tagged_DP_FMISS.table"
raw_SVs <- fread(RAW_SV, col.names = c('CHROM', 'POS', 'DP', 'F_MISS'))
raw_SVs$var_type <- 'SV'

ALL_raw <- rbind(raw_SVs, raw_SNPs, raw_indels)

ALL_raw_medians <- aggregate(data = ALL_raw, F_MISS ~ var_type, FUN = median)
ALL_raw_medians$F_MISS <- round(ALL_raw_medians$F_MISS, 3)

ALL_raw_medians$y <- c(1000000, 8000000, 40000)

ggplot(data = ALL_raw) +
  facet_wrap(vars(factor(var_type, levels = c('SV', 'SNP', 'indel'))), 
             nrow = 3, scales = 'free_y', strip.position = 'right') +
  geom_histogram(aes(x = F_MISS), 
                 color = 'black', fill = 'grey70', linewidth = 0.1, binwidth = 0.02) +
  geom_vline(aes(xintercept = F_MISS), data = ALL_raw_medians, linetype = 2, color = 'firebrick') + 
  theme_bw() +
  theme(
    #panel.spacing.x = unit(0.6, 'points'),
    #panel.spacing.y = unit(3, 'points'),
    #panel.background = element_blank(),
    #panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.1),
    #panel.grid = element_blank(),
    
    strip.text.y.right = element_text(size = 8),
    strip.background.y = element_rect(color = 'black', linewidth = 0.2),
    
    #axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 4),
    #axis.title.x = element_text(size = 7),
    #axis.title.y = element_text(size = 7),
    #axis.ticks.y = element_line(linewidth = 0.3),
    #axis.line.y = element_line(linewidth = 0.1),
    #axis.line.x = element_line(linewidth = 0.1)
    
  ) +
  guides(fill = 'none') +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  #scale_fill_viridis_d(option = 'E') +
  labs(x = 'Missing genotype proportion',
       y = 'Variant count') +
  geom_text(aes(x = F_MISS + 0.08, y = y, label = sprintf("%0.3f", round(F_MISS, digits = 3))), 
            size = 2.5, color = 'firebrick', data = ALL_raw_medians)

ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/F_MISS/ALL_variants_raw_F_MISS_distrib.png',
       width = 3200,
       height = 3000,
       units = 'px',
       dpi = 700
       # device = 'pdf'
)


# SVs, SNPs and indels seperatly ------------------------------------------
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
         yes = 'other',
         no = SV_all$type)

# Plot
ggplot(data = SV_all) +
  geom_boxplot(aes(x = type, y = F_MISS))

ggplot(data = SV_all) +
  geom_histogram(aes(x = F_MISS, fill = type), color = 'black', linewidth = 0.1, binwidth = 0.02) +
  geom_vline(xintercept = median(SV_all$F_MISS))

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
         yes = 'other',
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
         yes = 'other',
         no = indel_all$type)

# Plot
ggplot(data = indel_all) +
  geom_boxplot(aes(x = type, y = F_MISS))




# All variants together, filtered and by categories -----------------------
SV_all$var_type <- 'SV'
SNP_all$var_type <- 'SNP'
indel_all$var_type <- 'indel'

ALL_filt <- rbind(SV_all, SNP_all, indel_all)

ALL_filt_medians <- aggregate(data = ALL_filt, F_MISS ~ var_type, FUN = median)
#ALL_filt_medians$F_MISS <- ALL_filt_medians$F_MISS
ALL_filt_medians$y <- c(100000, 1035000, 15500)

ggplot(data = ALL_filt) +
  facet_wrap(vars(factor(var_type, levels = c('SV', 'SNP', 'indel'))), 
             nrow = 3, scales = 'free_y', strip.position = 'right') +
  geom_histogram(aes(x = F_MISS, fill = type), 
                 color = 'black', linewidth = 0.1, binwidth = 0.02) +
  geom_vline(aes(xintercept = F_MISS), data = ALL_filt_medians, linetype = 2, color = 'firebrick') + 
  theme_bw() +
  theme(
    #panel.spacing.x = unit(0.6, 'points'),
    #panel.spacing.y = unit(3, 'points'),
    #panel.background = element_blank(),
    #panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.1),
    #panel.grid = element_blank(),
    
    strip.text.y.right = element_text(size = 8),
    strip.background.y = element_rect(color = 'black', linewidth = 0.2),
    
    #axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 4),
    #axis.title.x = element_text(size = 7),
    #axis.title.y = element_text(size = 7),
    #axis.ticks.y = element_line(linewidth = 0.3),
    #axis.line.y = element_line(linewidth = 0.1),
    #axis.line.x = element_line(linewidth = 0.1)
    
  ) +
  #guides(fill = 'none') +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  scale_fill_manual(values = c('gray70', 'yellow', 'purple4')) +
  labs(x = 'Missing genotype proportion',
       y = 'Variant count') +
  geom_text(aes(x = F_MISS + 0.08, y = y, label = sprintf("%0.3f", round(F_MISS, digits = 3))), size = 2.5, color = 'firebrick', data = ALL_filt_medians)


ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/F_MISS/ALL_variants_filt_F_MISS_distrib.png',
       width = 3200,
       height = 3000,
       units = 'px',
       dpi = 700
       # device = 'pdf'
)
