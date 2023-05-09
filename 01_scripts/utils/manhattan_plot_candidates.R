# Produce per-site manhattan plot for Fst values, including info on candidate sets

library(ggplot2)
library(dplyr)

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
#ALL_SITES <- "08_angsd_fst/SVs/merged_SUPP2_MAF0.05_FMISS0.5.SVsFst_RO_PU.table"
ALL_SITES <- argv[1]

#FST_OUTLIERS <- "08_angsd_fst/SVs/SVs_RO_PU_outliers_minFst0.243.table"
FST_OUTLIERS <- argv[2]

#RDA_CAND <- "10_rda/SVs/RDA_3sd_outliers.table"
RDA_CAND <- argv[3]

#FISHER_OUTLIERS <- "11_fisher_tests/SVs/SVs_fisher_RO_PU_outliers_qval0.01.table"
FISHER_OUTLIERS <- argv[4]

#INTERSECT_OUTLIERS <- "08_angsd_fst/SVs/SVs_RO_PU_outliers_minFst0.243_qval0.01_shared.table"
INTERSECT_OUTLIERS <- argv[5]

#SHARED_3_METHODS <- "08_angsd_fst/SVs/SVs_RO_PU_outliers_minFst0.243_qval0.01_RDA_3sd_shared.table"
SHARED_3_METHODS <- argv[6]

QUANTILE <- argv[7]
MIN_FST <- as.numeric(argv[8])
SD <- argv[9]
MAX_QVAL <- argv[10]

OVERLAP_WIN <- argv[11]

#OV_2_SSA <- "02_infos/OV_to_ssa.txt"
OV_2_SSA <- argv[12]


# Import all known variant sites
all_sites <- read.table(ALL_SITES, 
                        col.names = c('CHROM', 'POS', 'END', 'ID', 'FST'),
                        na.strings = '.'
)

## Convert CHROM OV to SSA
OV_2_sasa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))
OV_2_sasa$CHROM_SIMPL <- sapply(X = OV_2_sasa$CHROM_SSA, FUN = function(x){unlist(strsplit(x = x, split = 'ssa'))[2]})


all_sites <- merge(x = all_sites, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')

## Convert negative Fst values to 0
all_sites$FST <- ifelse(all_sites$FST < 0, 
                        yes = 0,
                        no = all_sites$FST)



# Import Fst outliers
FST_outliers <- read.table(FST_OUTLIERS, 
                           col.names = c('CHROM', 'POS', 'END', 'ID', 'FST'),
                           na.strings = '.'
)

# Import Fisher outliers
Fisher_outliers <- read.table(FISHER_OUTLIERS, 
                              col.names = c('CHROM', 'POS', 'END', 'ID', 'PVAL', 'QVAL'),
                              na.strings = '.'
)


# Import outlier intersection set
intersect_outliers <- read.table(INTERSECT_OUTLIERS, 
           col.names = c('CHROM', 'POS', 'END', 'ID'),
           na.strings = '.'
)

intersect_outliers <- merge(x = intersect_outliers, y = all_sites, 
                            by = c('CHROM', 'POS', 'END', 'ID'),
                            all.x = TRUE)

# Import RDA candidates
RDA_cand <- read.table(RDA_CAND, 
           col.names = c('CHROM', 'POS', 'END', 'ID', 'loadings'),
           na.strings = '.'
)


RDA_cand <- merge(x = RDA_cand, y = all_sites,
                  by = c('CHROM', 'POS', 'END', 'ID'),
                  all.x = TRUE)


# Import shared across 3 methods set
shared_3 <- read.table(SHARED_3_METHODS, 
                       col.names = c('CHROM', 'POS', 'END', 'ID'),
                       na.strings = '.'
                       )
shared_3 <- merge(x = shared_3, y = all_sites,
                  by = c('CHROM', 'POS', 'END', 'ID'),
                  all.x = TRUE)

print('done import')

# 2. Plot -----------------------------------------------------------------

# Get colors 
hex_cols <- viridisLite::viridis(n = 3, option = 'D')
scales:::show_col(hex_cols)
print('done colors')

manh_plot <- 
ggplot(data = all_sites) +
  facet_wrap(~ CHROM_SIMPL, nrow = 1, scales = 'free_x') +
  geom_point(aes(x = POS, y = FST), alpha = 2 * all_sites$FST, size = 3*all_sites$FST, col = 'grey30') +
  theme(panel.spacing = unit(1, 'points'),
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
       #y = bquote(F[ST]), # somehow this throws a Error: Discrete value supplied to continuous scale
       y = expression(F[ST]),
       fill = 'Candidate variants group') +
  geom_hline(yintercept = MIN_FST, col = 'red', linetype  = 2) +
  ## Outliers set layer
  geom_point(data = intersect_outliers, 
             aes(x = POS, y = FST, fill = 'Outlier intersection'), shape = 21, size = 1.5) +
  ## RDA candidates set layer
  geom_point(data = RDA_cand, 
             aes(x = POS, y = FST, fill = 'RDA candidate'), shape = 21, size = 1.5) +
  ## Shared across 3 methods set
  geom_point(data = shared_3, 
             aes(x = POS, y = FST, fill = 'Shared (Fst + Fisher + RDA)'), shape = 23, size = 2) +
  scale_fill_manual(values = hex_cols) 

manh_plot

ggsave(manh_plot, filename = paste0(unlist(strsplit(ALL_SITES, split = '.table'))[1], '.manh_plot.jpg'))

saveRDS(manh_plot, file = paste0(unlist(strsplit(ALL_SITES, split = '.table'))[1], '.manh_plot.rds')
        )
#save(manh_plot, all_sites, intersect_outliers, RDA_cand, shared_3, MIN_FST, file = paste0(unlist(strsplit(ALL_SITES, split = '.table'))[1], '.manh_plot.Rdata'))

#outliers_candidates <- merge(x = intersect_outliers, y = RDA_cand,
#      by = c('CHROM', 'POS', 'END', 'ID', 'FST', 'CHROM_SSA', 'CHROM_SIMPL'), all = TRUE)  
