# Plot SVs, SNPs and indels density

options(scipen=999)

library(ggplot2)

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
SV_DENSITY <- argv[1] 
SNP_DENSITY <- argv[2] 
INDEL_DENSITY <- argv[3] 
OV_2_SSA <- argv[4]
WIN_SIZE <- as.numeric(argv[5])

SV_density <-  read.delim(SV_DENSITY, header = FALSE,
                           col.names = c('CHROM', 'BIN_START', 'BIN_STOP', 'DENSITY'))

SNP_density <-  read.delim(SNP_DENSITY, header = FALSE,
                          col.names = c('CHROM', 'BIN_START', 'BIN_STOP', 'DENSITY'))

indel_density <-  read.delim(INDEL_DENSITY, header = FALSE,
                           col.names = c('CHROM', 'BIN_START', 'BIN_STOP', 'DENSITY'))

OV_to_ssa <- read.table(OV_2_SSA, col.names = c('OV', 'SSA'))


# Add variable for type
SV_density$TYPE <- 'SV'
SNP_density$TYPE <- 'SNP'
indel_density$TYPE <- 'indel'


# Bind the 3 dataframes together

var_density <- rbind(SV_density, SNP_density, indel_density)

# Convert chr names to ssa notation
var_density <- merge(x = var_density, y = OV_to_ssa, 
                    by.x = 'CHROM', by.y = 'OV', sort = FALSE)


# 2. Plot -----------------------------------------------------------------

## Compute bin middle 
var_density$BIN_MIDDLE <- var_density$BIN_START + WIN_SIZE/2

density_plot <-
  ggplot(data = var_density) + 
  facet_wrap(~ SSA, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = BIN_MIDDLE, y = DENSITY, col = TYPE), size = 1) +
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

