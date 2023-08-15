# Plot per window Fst for SVs, SNPs and indels together

library(ggplot2)

# 1. Import and format ----------------------------------------------------
SV_FST_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SVs/RO_PU/RO_PU_win100000_step10000.txt"
SNP_FST_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SNPs/RO_PU/RO_PU_win100000_step10000.txt"
INDEL_FST_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/indels/RO_PU/RO_PU_win100000_step10000.txt"

SVs_mean_Fst <- 0.044
SNPs_mean_Fst <- 0.065
indels_mean_Fst <- 0.079


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
## Add mean Fst for each type
mean_Fst <- data.frame(TYPE = c("SVs", "SNPs", "indels"), 
                     mean_Fst = c(SVs_mean_Fst,
                                  SNPs_mean_Fst,
                                  indels_mean_Fst))
ALL_Fst_win <- merge(ALL_Fst_win, mean_Fst, by = c('TYPE'))

## Plot
ggplot(data = ALL_Fst_win) +
  #facet_wrap(~ CHROM_SSA, nrow = 2, scales = 'free_x') +
  facet_grid(factor(TYPE, levels = c('SVs', 'SNPs', 'indels')) ~ CHROM_NUM, 
             scales = 'free_x', space = 'free_x') +
  geom_point(aes(x = midPOS, y = FST, col = TYPE),
            alpha = 1.5 * ALL_Fst_win$FST, 
            size = 1.2*ALL_Fst_win$FST, 
            col = 'midnightblue') +
  theme(panel.spacing.x = unit(0.2, 'points'),
        panel.spacing.y = unit(2, 'points'),
        strip.text.x.top = element_text(size = 4,
                                        margin = margin(3,0,3,0, 'pt')),
        strip.text.y.right = element_text(size = 5,
                                          margin = margin(0,3,0,3, 'pt')),
        strip.placement = "inside",
        strip.background = element_rect(colour = 'gray70'),
        
        axis.text.x = element_text(angle = 45, size = 3, hjust = 1),
        axis.text.y = element_text(size = 4),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.2),
        panel.grid.minor.y = element_line(linewidth = 0.2),
        panel.grid.major.y = element_line(linewidth = 0.3),
        panel.background = element_rect(color = "gray70")
        ) +
  scale_x_continuous(
    breaks = seq(0, 1.6*(10^8), by = (0.4*10^8)),
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) + 
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = expression(italic('F'['ST']))
                      ) + 
  guides(color = 'none') +
  geom_hline(aes(yintercept = mean_Fst), 
             color = 'white', linewidth = 0.5, linetype = 'dotdash') +
  geom_text(data = subset(ALL_Fst_win, CHROM_NUM == '01-23'),
            aes(y = mean_Fst + 0.03, 
                x = 45000000, 
            label = mean_Fst),
            size = 1.4, color = 'white')
  
  
# Save to external file
ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/combined_per_win_Fst.png',
       width = 2800,
       height = 3100,
       units = 'px',
       dpi = 700
)
