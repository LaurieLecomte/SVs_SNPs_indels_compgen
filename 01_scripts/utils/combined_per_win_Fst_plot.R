# Plot per window Fst for SVs, SNPs and indels together

library(ggplot2)

# 1. Import and format ----------------------------------------------------
SV_FST_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SVs/RO_PU/RO_PU_win100000_step10000.txt"
SNP_FST_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SNPs/RO_PU/RO_PU_win100000_step10000.txt"
INDEL_FST_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/indels/RO_PU/RO_PU_win100000_step10000.txt"


SVs_Fst_win <- read.table(SV_FST_WIN, skip = 1, sep = "\t",
                      col.names = c('REGION', 'CHROM', 'midPOS', 'Nsites', 'FST')
                      )
## Add a column for type, for easier plotting
SVs_Fst_win$TYPE <- 'SVs'

SNPs_Fst_win <- read.table(SNP_FST_WIN, skip = 1, sep = "\t",
                      col.names = c('REGION', 'CHROM', 'midPOS', 'Nsites', 'FST')
)
SNPs_Fst_win$TYPE <- 'SNPs'

indels_Fst_win <- read.table(INDEL_FST_WIN, skip = 1, sep = "\t",
                           col.names = c('REGION', 'CHROM', 'midPOS', 'Nsites', 'FST')
)
indels_Fst_win$TYPE <- 'indels'

# Bind these 3 datasets together
ALL_Fst_win <- do.call("rbind", list(SVs_Fst_win, SNPs_Fst_win, indels_Fst_win))

# Convert OV to SSA chrom
OV_2_SSA <- '02_infos/OV_to_ssa.txt'
OV_2_ssa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))
## simplify for plotting
OV_2_ssa$CHROM_NUM <- sapply(X = OV_2_ssa$CHROM_SSA, FUN = function(x){
  unlist(strsplit(x, split = 'ssa'))[2]}
)

## merge with full dataset
ALL_Fst_win <- merge(x = ALL_Fst_win, y = OV_2_ssa, 
                     by.x = 'CHROM', by.y = 'CHROM_OV',
                     sort = FALSE)


# 2. Plot  ----------------------------------------------------------------
ggplot(data = ALL_Fst_win) +
  #facet_wrap(~ CHROM_SSA, nrow = 2, scales = 'free_x') +
  facet_grid(factor(TYPE, levels = c('SVs', 'SNPs', 'indels')) ~ CHROM_NUM, 
             scales = 'free_x', space = 'free_x') +
  geom_point(aes(x = midPOS, y = FST, col = TYPE),
            alpha = 1.5 * ALL_Fst_win$FST, 
            size = 1.2*ALL_Fst_win$FST, 
            col = 'darkblue') +
  theme(panel.spacing = unit(0.2, 'points'),
        strip.text.x = element_text(size = 5),
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
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = expression(paste('Per window ', italic('F'['ST'])))
                      ) + 
  guides(color = FALSE)
# +
  #geom_hline(yintercept = mean_fst$weighted_FST, col = 'red')

#ggsave()
