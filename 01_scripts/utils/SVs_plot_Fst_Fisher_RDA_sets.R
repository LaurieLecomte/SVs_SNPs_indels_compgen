library(ggplot2)
library(VennDiagram)
library(dplyr)
#library(ggVennDiagram)

INPUT_TABLE <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SVs/merged_SUPP2_MAF0.05_FMISS0.5.SVsFst_RO_PU.table"

FST_OUTLIERS <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SVs/SVs_RO_PU_outliers_minFst0.278.table"
FISHER_OUTLIERS <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/11_fisher_tests/SVs/SVs_fisher_RO_PU_outliers_qval0.01.table"
  
SHARED_OUTLIERS <- "08_angsd_fst/SVs/SVs_RO_PU_outliers_minFst0.278_qval0.01_shared.table"

RDA_CAND <- "10_rda/SVs/RDA_3sd_outliers.table"

OV_2_SSA <- "02_infos/OV_to_ssa.txt"

MIN_FST <- 0.278

MEAN_FST <- "08_angsd_fst/SVs/RO_PU/RO_PU.SVs.fst"


SD <- 3
MAX_QVAL <- 0.01


# Per site Fst
FST_table <- read.table(INPUT_TABLE, 
                        col.names = c('CHROM', 'POS', 'END', 'ID', 'FST'),
                        na.strings = '.'
)


FST_outliers <- read.table(FST_OUTLIERS, 
                           col.names = c('CHROM', 'POS', 'END', 'ID', 'FST'),
                           na.strings = '.'
)

FST_outliers$FST_outlier <- 1
FST_outliers$code <- paste0(FST_outliers$CHROM, '_', FST_outliers$POS, '_', FST_outliers$END)
  
  
Fisher_outliers <- read.table(FISHER_OUTLIERS, 
                           col.names = c('CHROM', 'POS', 'END', 'ID', 'PVAL', 'QVAL'),
                           na.strings = '.'
)
Fisher_outliers$Fisher_outlier <- 1
Fisher_outliers$code <- paste0(Fisher_outliers$CHROM, '_', Fisher_outliers$POS, '_', Fisher_outliers$END)

shared_outliers <- read.table(SHARED_OUTLIERS, 
                        col.names = c('CHROM', 'POS', 'END', 'ID'),
                        na.strings = '.'
)


#outliers$type <- 'Fst/Fisher outlier'

RDA_cand <- read.table(RDA_CAND, 
                       col.names = c('CHROM', 'POS', 'END', 'ID', 'loadings'),
                       na.strings = '.'
)
RDA_cand$code <- paste0(RDA_cand$CHROM, '_', RDA_cand$POS, '_', RDA_cand$END)
RDA_cand$RDA_cand <- 1



mean_fst <- read.table(MEAN_FST, col.names = c('unweighted_FST', 'weighted_FST'))

## Convert OV chromosomes to ssa notation
OV_2_sasa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))

FST_table <- merge(x = FST_table, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')

#outliers <- merge(x = outliers, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')

#RDA_cand <- merge(x = RDA_cand, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')

## Convert negative Fst values to 0
FST_table$FST <- ifelse(FST_table$FST < 0, 
                        yes = 0,
                        no = FST_table$FST)




FST_Fisher_merge <- merge(x = select(FST_outliers, -FST), 
                          y = select(Fisher_outliers, -c(PVAL, QVAL)), 
                          by = c('CHROM', 'POS', 'END', 'ID', 'code'),
                          all = TRUE)
FST_Fisher_RDA_merge <- merge(x = FST_Fisher_merge, 
                          y = select(RDA_cand, -loadings), 
                          by = c('CHROM', 'POS', 'END', 'ID', 'code'),
                          all = TRUE)

FST_Fisher_RDA_merge[is.na(FST_Fisher_RDA_merge)] <- 0

#FST_Fisher_RDA_merge <- merge(x = FST_Fisher_RDA_merge, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')

 
for (i in 1:nrow(FST_Fisher_RDA_merge)){
  combi <- paste0(FST_Fisher_RDA_merge$FST_outlier[i], FST_Fisher_RDA_merge$Fisher_outlier[i], FST_Fisher_RDA_merge$RDA_cand[i])
  FST_Fisher_RDA_merge$group[i] <-
  switch(combi,
         '111' = 'Fst + Fisher + RDA',
         '100' = 'Fst',
         '010' = 'Fisher',
         '001' = 'RDA',
         '110' = 'Fst + Fisher',
         '011' = 'RDA',
         '101' = 'RDA')
}

for (i in 1:nrow(FST_Fisher_RDA_merge)){
  combi <- paste0(FST_Fisher_RDA_merge$FST_outlier[i], FST_Fisher_RDA_merge$Fisher_outlier[i], FST_Fisher_RDA_merge$RDA_cand[i])
  FST_Fisher_RDA_merge$cand_group[i] <-
    switch(combi,
           '111' = 'Shared (Fst + Fisher + RDA)',
           '100' = 'Fst',
           '010' = 'Fisher',
           '001' = 'RDA candidate',
           '110' = 'Outlier (Fst AND Fisher)',
           '011' = 'RDA candidate',
           '101' = 'RDA candidate')
}

FST_Fisher_RDA_all_sites <- merge(x = FST_Fisher_RDA_merge, y = FST_table, by = c('CHROM', 'POS', 'END', 'ID'), 
                                  all.x = TRUE)


ggplot(data = FST_table) +
  facet_wrap(~ CHROM_SSA, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = POS, y = FST), alpha = 2 * FST_table$FST, size = 2*FST_table$FST, col = 'grey40') +
  theme(panel.spacing = unit(0.1, 'points'),
        strip.text.x = element_text(size = 6),
        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        panel.background = element_rect(color = "gray70"),
        strip.placement = "inside",
        strip.background = element_rect(colour = 'gray70')
  ) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) + 
  ylim(0, 1) + 
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = 'Fst') +
  #geom_hline(yintercept = mean_fst$weighted_FST, col = 'red') +
 # geom_hline(yintercept = MIN_FST, col = 'red', linetype  = 2) +
  geom_point(data = FST_Fisher_RDA_all_sites, 
             aes(x = POS, y = FST, col = group, shape = cand_group), size = 0.8) + 
  #geom_point(data = RDA_cand, aes(x = POS, y = FST), col = 'green', size = 0.8) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) 




# best one :
ggplot(data = FST_table) +
  facet_wrap(~ CHROM_SSA, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = POS, y = FST), alpha = 2 * FST_table$FST, size = 4*FST_table$FST, col = 'grey40') +
  theme(panel.spacing = unit(0.5, 'points'),
        strip.text.x = element_text(size = 6),
        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        panel.background = element_rect(color = "gray90"),
        strip.placement = "inside",
        strip.background = element_rect(colour = 'gray90')
  ) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) + 
  ylim(0, 0.8) + 
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = bquote(F[ST])) +
  geom_hline(yintercept = MIN_FST, col = 'red', linetype  = 2) +
  geom_point(data = subset(FST_Fisher_RDA_all_sites, group != 'Fst' & group != 'Fisher'), 
             aes(x = POS, y = FST, fill = cand_group), shape = 21, col = 'black') + 
  #geom_point(data = RDA_cand, aes(x = POS, y = FST), col = 'green', size = 0.8) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) +
  scale_fill_viridis_d(option = 'D') 



# Venn Diagram

## Get hex codes
hex_cols <- viridisLite::viridis(n = 3, option = 'D')
scales:::show_col(hex_cols)


display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(
  x = list(FST_outliers$code, Fisher_outliers$code, RDA_cand$code),
  category.names = c("Fst" , "Fisher" , "RDA"),
  fill = hex_cols
  #filename = paste0(unlist(strsplit(SHARED_OUTLIERS, split = '_shared.table')), '.venn.png')
  #output=TRUE
)






## Match Fst values with sets
#outliers <- merge(x = outliers, y = fst_table, by = c('CHROM', 'CHROM_SSA', 'POS', 'END', 'ID'))
#RDA_cand <- merge(x = RDA_cand, y = fst_table, by = c('CHROM', 'CHROM_SSA', 'POS', 'END', 'ID'))

ggplot(data = fst_table) +
  facet_wrap(~ CHROM_SSA, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = POS, y = FST), alpha = 1.5 * fst_table$FST, size = 1.2*fst_table$FST, col = 'darkblue') +
  theme(panel.spacing = unit(0.1, 'points'),
        strip.text.x = element_text(size = 6),
        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        panel.background = element_rect(color = "gray70"),
        strip.placement = "inside",
        strip.background = element_rect(colour = 'gray70')
  ) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) + 
  ylim(0, 1) + 
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = 'Fst') +
  geom_hline(yintercept = mean_fst$weighted_FST, col = 'red') +
  geom_hline(yintercept = MIN_FST, col = 'red', linetype  = 2) +
  geom_point(data = outliers, aes(x = POS, y = FST), col = 'orange', size = 0.8) + 
  geom_point(data = RDA_cand, aes(x = POS, y = FST), col = 'green', size = 0.8) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) 
