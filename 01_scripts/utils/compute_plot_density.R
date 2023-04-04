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
OV_2_ssa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))

var_density_OV <- merge(x = var_density, y = OV_2_ssa, by.x = 'CHROM', by.y = 'CHROM_OV', sort = FALSE)


# 2. Compute mean density by CHROM ----------------------------------------

density_by_chrom <- aggregate(DENSITY ~ CHROM_SSA, data = var_density_OV, FUN = "mean")

write.table(density_by_chrom, file = paste0(unlist(strsplit(DENSITY_TABLE, split = '.txt')), '_byCHROM.txt'),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


# 3. Plot density by bins -------------------------------------------------
## Compute bin middle 
var_density_OV$BIN_MIDDLE <- var_density_OV$BIN_START + WIN_SIZE/2

density_plot <-
ggplot(data = var_density_OV) + 
  facet_wrap(~ CHROM_SSA, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = DENSITY), size = 1) +
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
