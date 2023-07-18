# Plot Fst distribution and per site Fst

library(ggplot2)

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
INPUT_TABLE <- argv[1]
OV_2_SSA <- argv[2]
MEAN_FST <- argv[3]
WIN_FST <- argv[4]
WIN_SIZE <- as.numeric(argv[5])

# Per site Fst
fst_table <- read.table(INPUT_TABLE, 
                        col.names = c('CHROM', 'POS', 'END', 'ID', 'FST'),
                        na.strings = '.'
                        )

# Per window Fst 
fst_win <- read.table(WIN_FST, skip = 1, sep = "\t",
                      col.names = c('REGION', 'CHROM', 'midPOS', 'Nsites', 'FST')
)
                

## Convert negative Fst values to 0
fst_table$FST <- ifelse(fst_table$FST < 0, 
                        yes = 0,
                        no = fst_table$FST)

fst_win$FST <- ifelse(fst_win$FST < 0, 
                        yes = 0,
                        no = fst_win$FST)

## Get min and max
print(paste('Min Fst :', min(fst_table$FST)))
print(paste('Max Fst :', max(fst_table$FST)))

## Convert OV chromosomes to ssa notation
OV_2_sasa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))
## simplify for plotting
OV_2_sasa$CHROM_NUM <- sapply(X = OV_2_sasa$CHROM_SSA, FUN = function(x){
  unlist(strsplit(x, split = 'ssa'))[2]}
  )

## Add info on chrom type ### NOT USED, we do not know each chrom's type
#meta <- c()
#OV_2_ssa$CHROM_TYPE <- ifelse(OV_2_ssa$CHROM_SSA %in% meta,
#                               yes = 'meta',
#                               no = 'acro')


fst_table <- merge(x = fst_table, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')

fst_win <- merge(x = fst_win, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')

# Genome wide Fst
mean_fst <- read.table(MEAN_FST, col.names = c('unweighted_FST', 'weighted_FST'))


# 2. Plot Fst distribution ------------------------------------------------
# Split Fst by bins
fst_bins <- seq(0, 1, by = 0.02)

fst_bins_names <- vector(mode = 'character', length = length(fst_bins))
for(i in 1:length(fst_bins)){
  if (i < length(fst_bins)){
    fst_bins_names[i] <- (paste0('[', fst_bins[i], '-', fst_bins[i+1], '['))
  } else {
    fst_bins_names[i] <- ('1')
  }
  
}

fst_bins <- c( fst_bins, Inf)
fst_table$Fst_bin <-
  cut(abs(fst_table$FST), breaks = fst_bins, labels = fst_bins_names, right = FALSE)

# Get category for weighted genome wide Fst
wFst <-  as.vector(
  cut(abs(mean_fst$weighted_FST), breaks = fst_bins, 
      labels = fst_bins_names, right = FALSE)
  )


# Plot
## Note : NA are variants present in table, but that were filtered out prior to Fst calculation
fst_distrib_plot <- 
ggplot(data = fst_table) + 
  geom_bar(aes(x = Fst_bin)) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      size = 4,
      hjust = 1
    ),
    plot.title = element_text(size = 11)
    ) +
  geom_vline(aes(xintercept = wFst), col = 'red') +
  annotate('text', x = levels(fst_table$Fst_bin)[(which(levels(fst_table$Fst_bin) == wFst)) + 3], 
           y = 0.75 * max(table(fst_table$Fst_bin)), 
           label= paste("weighted Fst :\n", round(mean_fst$weighted_FST, 3)),
           col ='red', size = 3) +
  labs(x = 'Fst', y = 'Count', title = 'Fst distribution') +
  theme(
    plot.title = element_text(size = 10, face = 'bold', hjust = 0.5)
  )

## save both plot and wFst
datalist <- list(wFst, fst_distrib_plot)
saveRDS(datalist, file = paste0(strsplit(INPUT_TABLE, split = '.table')[[1]], '_fst_distrib.rds'))

jpeg(file = paste0(strsplit(INPUT_TABLE, split = '.table')[[1]], '_fst_distrib.jpg'))
fst_distrib_plot
dev.off()

# 3. Per site Fst ---------------------------------------------------------
fst_per_site_plot <- 
ggplot(data = fst_table) +
  facet_grid(.~CHROM_NUM, scales = 'free_x') +
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
  geom_hline(yintercept = mean_fst$weighted_FST, col = 'red')

jpeg(file = paste0(strsplit(INPUT_TABLE, split = '.table')[[1]], '_fst_per_site.jpg'))
fst_per_site_plot

saveRDS(fst_per_site_plot, file = paste0(strsplit(INPUT_TABLE, split = '.table')[[1]], '_fst_per_site.rds'))
dev.off()


# 4. Per window Fst -------------------------------------------------------
fst_win_plot <- 
  ggplot(data = fst_win) +
  #facet_wrap(~ CHROM_SSA, nrow = 2, scales = 'free_x') +
  facet_grid(.~CHROM_NUM, scales = 'free_x', space = 'free_x') +
  geom_point(aes(x = midPOS, y = FST), alpha = 1.5 * fst_win$FST, size = 1.2*fst_win$FST, col = 'darkblue') +
  theme(panel.spacing = unit(0.8, 'points'),
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
  #ylim(0, 1) + 
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = 'Per window Fst') +
  geom_hline(yintercept = mean_fst$weighted_FST, col = 'red')

jpeg(file = paste0(strsplit(WIN_FST, split = '.txt')[[1]], '_fst_per_win.jpg'))
fst_win_plot

saveRDS(fst_win_plot, file = paste0(strsplit(WIN_FST, split = '.txt')[[1]], '_fst_per_win.rds'))
dev.off()


# meta vs centro
#acro <- subset(fst_win, CHROM_TYPE == 'acro')
#fst_win_chr_acro <-
#  ggplot(data = acro) + 
#  facet_grid(.~CHROM_NUM, scales = 'free_x', space = 'free_x') +
#  geom_point(aes(x = midPOS, y = FST), size = 0.8, alpha = 0.6,
#             col = 'firebrick') +
#  theme(panel.spacing = unit(0.1, 'points'),
#        strip.text.x = element_text(size = 6),
#        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
#        panel.background = element_rect(color = "gray70"),
#        strip.placement = "inside",
#        strip.background = element_rect(colour = 'gray70'),
#        legend.position = 'bottom'
#  ) + 
#  scale_x_continuous(
#    labels = function(x) {
#      round(x/10^8, 1)
#    }
#  ) + 
#  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
#       y = paste('Fst by', WIN_SIZE/1000000, 'Mb window')) 

#meta <- subset(fst_win, CHROM_TYPE == 'meta')
#fst_win_chr_meta <-
#  ggplot(data = meta) + 
#  facet_grid(.~CHROM_NUM, scales = 'free_x', space = 'free_x') +
#  geom_point(aes(x = midPOS, y = FST), size = 0.8, alpha = 0.6,
#             col = 'blue') +
#  theme(panel.spacing = unit(0.1, 'points'),
#        strip.text.x = element_text(size = 6),
#        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
#        panel.background = element_rect(color = "gray70"),
#        strip.placement = "inside",
#        strip.background = element_rect(colour = 'gray70'),
#        legend.position = 'bottom'
#  ) + 
#  scale_x_continuous(
#    labels = function(x) {
#      round(x/10^8, 1)
#    }
#  ) + 
#  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
#       y = paste('Fst by', WIN_SIZE/1000000, 'Mb window')) 


#library(ggpubr)
#fst_win_plot_chr_type <- 
#  ggarrange(fst_win_chr_acro, 
#            fst_win_chr_meta, 
#            ncol = 1, nrow = 2, common.legend = TRUE, legend = 'bottom',
#            label.x = 0.05,
#            label.y = 0.8,
#            labels = list('acro', 'meta'))

#ggsave(fst_win_plot_chr_type, filename = paste0(unlist(strsplit(WIN_FST, split = '.txt')), '_chr_type_plot.jpg'), width = 10, height = 6)


#saveRDS(fst_win_plot_chr_type, file = paste0(unlist(strsplit(WIN_FST, split = '.txt')), '_chr_type_plot.rds'))


