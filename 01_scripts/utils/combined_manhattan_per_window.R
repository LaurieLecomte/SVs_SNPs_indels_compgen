# Plot Fst per window for SVs, SNPs and indels together

SV_FST_WIN <- "08_angsd_fst/SVs/RO_PU/RO_PU_win1000000_step10000.txt"
SNP_FST_WIN <- "08_angsd_fst/SNPs/RO_PU/RO_PU_win1000000_step10000.txt"
INDEL_FST_WIN <- "08_angsd_fst/indels/RO_PU/RO_PU_win1000000_step10000.txt"
OV_2_SSA <- '02_infos/OV_to_ssa.txt'
SV_MEAN_FST <- "08_angsd_fst/SVs/RO_PU/RO_PU.SVs.fst"
SNP_MEAN_FST <- "08_angsd_fst/SNPs/RO_PU/RO_PU.SNPs.fst"
INDEL_MEAN_FST <- "08_angsd_fst/indels/RO_PU/RO_PU.indels.fst"

# Genome wide Fst
SVs_mean_fst <- read.table(SV_MEAN_FST, col.names = c('unweighted_FST', 'weighted_FST'))
SNPs_mean_fst <- read.table(SNP_MEAN_FST, col.names = c('unweighted_FST', 'weighted_FST'))
indels_mean_fst <- read.table(INDEL_MEAN_FST, col.names = c('unweighted_FST', 'weighted_FST'))

# Per window Fst 
SVs_Fst_win <- read.table(SV_FST_WIN, skip = 1, sep = "\t",
                      col.names = c('REGION', 'CHROM', 'midPOS', 'Nsites', 'FST')
)

SNPs_Fst_win <- read.table(SNP_FST_WIN, skip = 1, sep = "\t",
                           col.names = c('REGION', 'CHROM', 'midPOS', 'Nsites', 'FST')
)

indels_Fst_win <- read.table(INDEL_FST_WIN, skip = 1, sep = "\t",
                             col.names = c('REGION', 'CHROM', 'midPOS', 'Nsites', 'FST')
)

## Convert negative Fst values to 0
SVs_Fst_win$FST <- ifelse(SVs_Fst_win$FST < 0, 
                        yes = 0,
                        no = SVs_Fst_win$FST)

SNPs_Fst_win$FST <- ifelse(SNPs_Fst_win$FST < 0, 
                          yes = 0,
                          no = SNPs_Fst_win$FST)

indels_Fst_win$FST <- ifelse(indels_Fst_win$FST < 0, 
                           yes = 0,
                           no = indels_Fst_win$FST)

## Convert OV chromosomes to ssa notation
OV_2_sasa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))
OV_2_sasa$CHROM_SIMPL <- sapply(X = OV_2_sasa$CHROM_SSA, FUN = function(x){unlist(strsplit(x = x, split = 'ssa'))[2]})

SVs_Fst_win <- merge(x = SVs_Fst_win, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')
SNPs_Fst_win <- merge(x = SNPs_Fst_win, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')
indels_Fst_win <- merge(x = indels_Fst_win, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')


# Plot
SVs_Fst_win_plot <- 
ggplot(data = SVs_Fst_win) +
  facet_wrap(~ CHROM_SIMPL, nrow = 1, scales = 'free_x') +
  geom_point(aes(x = midPOS, y = FST), alpha = SVs_Fst_win$FST, size = 3*SVs_Fst_win$FST, 
             col = 'darkblue') +
  theme(panel.spacing = unit(1.5, 'points'),
        strip.text.x = element_text(size = 6),
        #axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        panel.background = element_rect(color = "gray70"),
        strip.placement = "inside",
        strip.background = element_rect(colour = 'gray70'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank()
  ) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) + 
  ylim(0, 0.4) +
  geom_hline(yintercept = SVs_mean_fst$weighted_FST, col = 'red')

SNPs_Fst_win_plot <- 
  ggplot(data = SNPs_Fst_win) +
  facet_wrap(~ CHROM_SIMPL, nrow = 1, scales = 'free_x') +
  geom_point(aes(x = midPOS, y = FST), alpha = SNPs_Fst_win$FST, size = 3*SNPs_Fst_win$FST, 
             col = 'darkblue') +
  theme(panel.spacing = unit(1.5, 'points'),
        strip.text.x = element_text(size = 6),
        #axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        panel.background = element_rect(color = "gray70"),
        strip.placement = "inside",
        strip.background = element_rect(colour = 'gray70'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()
  ) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) + 
  ylim(0, 0.4) + 
  #labs(
   #    y = expression(Per~window~F[ST])) +
  geom_hline(yintercept = SNPs_mean_fst$weighted_FST, col = 'red')


indels_Fst_win_plot <- 
  ggplot(data = indels_Fst_win) +
  facet_wrap(~ CHROM_SIMPL, nrow = 1, scales = 'free_x') +
  geom_point(aes(x = midPOS, y = FST), alpha = indels_Fst_win$FST, size = 3*indels_Fst_win$FST, 
             col = 'darkblue') +
  theme(panel.spacing = unit(1.5, 'points'),
        strip.text.x = element_text(size = 6),
        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        panel.background = element_rect(color = "gray70"),
        strip.placement = "inside",
        strip.background = element_rect(colour = 'gray70'),
        axis.title.y = element_blank()
        #axis.title.x = element_blank()
  ) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) + 
  ylim(0, 0.4) + 
  labs(x = expression(paste('Position (', 10^8, ' bp)' ))) +
  geom_hline(yintercept = indels_mean_fst$weighted_FST, col = 'red')


combined <- 
  ggarrange(SVs_Fst_win_plot, 
            SNPs_Fst_win_plot, 
            indels_Fst_win_plot, 
            ncol = 1, nrow = 3)

combined_annot <- 
annotate_figure(combined,
                left = text_grob(expression(Per~window~F[ST]), rot = 90),
)


ggsave(combined_annot, filename = 'SVs_SNPs_indels_win1000000_step10000.jpg', width = 10, height = 10)
pdf(file = 'SVs_SNPs_indels_win1000000_step10000.pdf')
jpeg(file = 'SVs_SNPs_indels_win1000000_step10000_annot.jpg')
combined_annot
dev.off()

dev.print(file="test.png", device=png, width=1500)
