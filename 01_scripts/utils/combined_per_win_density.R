# Plot per window density for SVs, SNPs and indels together

library(ggplot2)

# 1. Import and format ----------------------------------------------------
SV_DENS_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/SVs/SVs_density_win100000bp.txt"
SNP_DENS_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/SNPs/SNPs_density_win100000bp.txt"
INDEL_DENS_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/indels/indels_density_win100000bp.txt"

WIN_SIZE <- 100000


SVs_dens_win <-  read.delim(SV_DENS_WIN, header = FALSE,
                           col.names = c('CHROM', 'BIN_START', 'BIN_STOP', 'DENSITY'))

SVs_dens_win$TYPE <- 'SVs'

SNPs_dens_win <-  read.delim(SNP_DENS_WIN, header = FALSE,
                            col.names = c('CHROM', 'BIN_START', 'BIN_STOP', 'DENSITY'))
SNPs_dens_win$TYPE <- 'SNPs'

indels_dens_win <-  read.delim(INDEL_DENS_WIN, header = FALSE,
                             col.names = c('CHROM', 'BIN_START', 'BIN_STOP', 'DENSITY'))
indels_dens_win$TYPE <- 'indels'


# Bind these 3 datasets together
ALL_dens_win <- do.call("rbind", list(SVs_dens_win, SNPs_dens_win, indels_dens_win))

# Convert OV to SSA chrom
OV_2_SSA <- '02_infos/OV_to_ssa.txt'
OV_2_ssa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))
## simplify for plotting
OV_2_ssa$CHROM_NUM <- sapply(X = OV_2_ssa$CHROM_SSA, FUN = function(x){
  unlist(strsplit(x, split = 'ssa'))[2]}
)

## add info on acro vs metacentric chromosomes
OV_2_ssa$CHROM_TYPE <- ifelse(OV_2_ssa$CHROM_SSA %in% c('ssa01-23', paste0('ssa0', seq(2,8))),
                               yes = 'meta',
                               no = 'acro')

## merge with full dataset
ALL_dens_win <- merge(x = ALL_dens_win, y = OV_2_ssa, 
                     by.x = 'CHROM', by.y = 'CHROM_OV',
                     sort = FALSE)

# 3. Plot density by bins -------------------------------------------------
## Compute bin middle 
ALL_dens_win$midBIN <- ALL_dens_win$BIN_START + WIN_SIZE/2

#density_plot <-
ggplot(data = ALL_dens_win) +
  facet_grid(factor(TYPE, levels = c('SVs', 'SNPs', 'indels')) ~ CHROM_NUM,
             scales = 'free',
             space = 'free_x') +
  geom_point(
    aes(x = midBIN, y = DENSITY),
    size = 0.8,
    alpha = 0.5,
    col = 'darkblue'
  ) +
  theme(
    panel.spacing = unit(0.1, 'points'),
    strip.text.x = element_text(size = 6),
    axis.text.x = element_text(
      angle = 45,
      size = 4,
      hjust = 1
    ),
    panel.background = element_rect(color = "gray70"),
    strip.placement = "inside",
    strip.background = element_rect(colour = 'gray70')
  ) +
  scale_x_continuous(
    labels = function(x) {
      round(x / 10 ^ 8, 1)
    }
  ) +
  labs(x = expression(paste('Position (', 10 ^ 8, ' bp)')),
       y = paste('Density by', WIN_SIZE / 1000000, 'Mb window'))

  
  # meta vs centro
acro <- subset(ALL_dens_win, CHROM_TYPE == 'acro')
#density_plot_chr_acro <-
  ggplot(data = acro) +
  facet_grid(factor(TYPE, levels = c('SVs', 'SNPs', 'indels')) ~ CHROM_NUM, 
             scales = 'free', space = 'free_x') +
  geom_point(
    aes(x = midBIN, y = DENSITY),
    size = 0.8,
    alpha = 0.8,
    col = 'darkorange3'
  ) +
  theme(
    panel.spacing = unit(0.1, 'points'),
    strip.text.x = element_text(size = 6),
    axis.text.x = element_text(
      angle = 45,
      size = 4,
      hjust = 1
    ),
    panel.background = element_rect(color = "gray70"),
    strip.placement = "inside",
    strip.background = element_rect(colour = 'gray70'),
    legend.position = 'bottom'
  ) +
  scale_x_continuous(
    labels = function(x) {
      round(x / 10 ^ 8, 1)
    }
  ) +
  labs(x = expression(paste('Position (', 10 ^ 8, ' bp)')),
       y = paste('Density by', WIN_SIZE / 1000000, 'Mb window'))

meta <- subset(ALL_dens_win, CHROM_TYPE == 'meta')
#density_plot_chr_meta <-
  ggplot(data = meta) +
  facet_grid(. ~ CHROM_NUM, scales = 'free_x', space = 'free_x') +
  geom_point(
    aes(x = midBIN, y = DENSITY),
    size = 0.8,
    alpha = 0.8,
    col = 'magenta4'
  ) +
  theme(
    panel.spacing = unit(0.1, 'points'),
    strip.text.x = element_text(size = 6),
    axis.text.x = element_text(
      angle = 45,
      size = 4,
      hjust = 1
    ),
    panel.background = element_rect(color = "gray70"),
    strip.placement = "inside",
    strip.background = element_rect(colour = 'gray70'),
    legend.position = 'bottom'
  ) +
  scale_x_continuous(
    labels = function(x) {
      round(x / 10 ^ 8, 1)
    }
  ) +
  labs(x = expression(paste('Position (', 10 ^ 8, ' bp)')),
       y = paste('Density by', WIN_SIZE / 1000000, 'Mb window')) 
  
  
