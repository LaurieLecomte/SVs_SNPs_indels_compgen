# Compute mean variant density by chromosome and plot

options(scipen=999)

library(ggplot2)

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
DENSITY_TABLE <- argv[1] 
OV_2_SSA <- argv[2]
WIN_SIZE <- as.numeric(argv[3])

var_density <-  read.delim(DENSITY_TABLE, header = FALSE,
                           col.names = c('CHROM', 'BIN_START', 'BIN_STOP', 'DENSITY'))


## Convert OV chromosomes to ssa notation
OV_2_sasa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))
## simplify for plotting
OV_2_sasa$CHROM_NUM <- sapply(X = OV_2_sasa$CHROM_SSA, FUN = function(x){
  unlist(strsplit(x, split = 'ssa'))[2]}
)
## add info on acro vs metacentric chromosomes
OV_2_sasa$CHROM_TYPE <- ifelse(OV_2_sasa$CHROM_SSA %in% c('ssa01-23', paste0('ssa0', seq(2,8))),
                               yes = 'meta',
                               no = 'acro')

var_density_OV <- merge(x = var_density, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV', sort = FALSE)


# 2. Compute mean density by CHROM ----------------------------------------

density_by_chrom <- aggregate(DENSITY ~ CHROM_SSA, data = var_density_OV, FUN = "mean")

write.table(density_by_chrom, file = paste0(unlist(strsplit(DENSITY_TABLE, split = '.txt')), '_byCHROM.txt'),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


# 3. Plot density by bins -------------------------------------------------
## Compute bin middle 
var_density_OV$BIN_MIDDLE <- var_density_OV$BIN_START + WIN_SIZE/2

density_plot <-
ggplot(data = var_density_OV) + 
  facet_grid(.~CHROM_NUM, scales = 'free_x', space = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = DENSITY), size = 0.8, alpha = 0.5) +
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
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = paste('Density by', WIN_SIZE/1000000, 'Mb window'))

saveRDS(density_plot, file = paste0(unlist(strsplit(DENSITY_TABLE, split = '.txt')), '_plot.rds'))

jpeg(file = paste0(unlist(strsplit(DENSITY_TABLE, split = '.txt')), '_plot.jpg'))
density_plot
dev.off()

# meta vs centro
acro <- subset(var_density_OV, CHROM_TYPE == 'acro')
density_plot_chr_acro <-
  ggplot(data = acro) + 
  facet_grid(.~CHROM_NUM, scales = 'free_x', space = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = DENSITY), size = 0.8, alpha = 0.8,
             col = 'firebrick') +
  theme(panel.spacing = unit(0.1, 'points'),
        strip.text.x = element_text(size = 6),
        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        panel.background = element_rect(color = "gray70"),
        strip.placement = "inside",
        strip.background = element_rect(colour = 'gray70'),
        legend.position = 'bottom'
  ) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) + 
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = paste('Density by', WIN_SIZE/1000000, 'Mb window')) 

meta <- subset(var_density_OV, CHROM_TYPE == 'meta')
density_plot_chr_meta <-
  ggplot(data = meta) + 
  facet_grid(.~CHROM_NUM, scales = 'free_x', space = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = DENSITY), size = 0.8, alpha = 0.8,
             col = 'blue') +
  theme(panel.spacing = unit(0.1, 'points'),
        strip.text.x = element_text(size = 6),
        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        panel.background = element_rect(color = "gray70"),
        strip.placement = "inside",
        strip.background = element_rect(colour = 'gray70'),
        legend.position = 'bottom'
  ) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) + 
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = paste('Density by', WIN_SIZE/1000000, 'Mb window')) 


library(ggpubr)
density_plot_chr_type <- 
  ggarrange(density_plot_chr_acro, 
            density_plot_chr_meta, 
            ncol = 1, nrow = 2, common.legend = TRUE, legend = 'bottom',
            label.x = 0.05,
            label.y = 0.8,
            labels = list('acro', 'meta'))

ggsave(density_plot_chr_type, filename = paste0(unlist(strsplit(DENSITY_TABLE, split = '.txt')), '_chr_type_plot.jpg'), width = 10, height = 6)


saveRDS(density_plot_chr_type, file = paste0(unlist(strsplit(DENSITY_TABLE, split = '.txt')), '_chr_type_plot.rds'))
