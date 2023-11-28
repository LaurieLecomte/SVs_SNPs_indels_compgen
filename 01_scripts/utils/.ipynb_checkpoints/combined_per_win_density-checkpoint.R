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

## Add info on chrom type ### NOT USED, we do not know each chrom's type
#meta <- c()
#OV_2_ssa$CHROM_TYPE <- ifelse(OV_2_ssa$CHROM_SSA %in% meta,
#                               yes = 'meta',
#                               no = 'acro')

## merge with full dataset
ALL_dens_win <- merge(x = ALL_dens_win, y = OV_2_ssa, 
                     by.x = 'CHROM', by.y = 'CHROM_OV',
                     sort = FALSE)

# 3. Plot density by bins -------------------------------------------------
## Compute bin middle 
ALL_dens_win$midBIN <- ALL_dens_win$BIN_START + WIN_SIZE/2

#density_plot <-
#ggplot(data = ALL_dens_win) +
#  facet_grid(factor(TYPE, levels = c('SVs', 'SNPs', 'indels')) ~ CHROM_NUM,
#             scales = 'free',
#             space = 'free_x') +
#  geom_point(
#    aes(x = midBIN, y = DENSITY),
#    size = 0.3,
#    col = 'midnightblue'
#  ) +
#  theme(panel.spacing.x = unit(0, 'points'),
#        panel.spacing.y = unit(2, 'points'),
#        strip.text.x.top = element_text(size = 4,
#                                        margin = margin(3,0,3,0, 'pt')),
#        strip.text.y.right = element_text(size = 5,
#                                          margin = margin(0,3,0,3, 'pt')),
#        strip.background.y = element_rect(color = 'grey80', linewidth = 0.1),
#        
#        #strip.placement = "inside",
#        strip.background.x = element_rect(colour = 'grey80', linewidth = 0.1),
#        
#        axis.text.x = element_text(angle = 45, size = 3, hjust = 1),
#        axis.text.y = element_text(size = 4),
#        axis.title.x = element_text(size = 7),
#        axis.title.y = element_text(size = 7),
#        axis.ticks.x = element_line(linewidth = 0.2),
#        axis.ticks.y = element_line(linewidth = 0.3),
#        #panel.grid.minor.x = element_blank(),
#        #panel.grid.major.x = element_line(linewidth = 0.2),
#        #panel.grid.minor.y = element_line(linewidth = 0.2),
#        panel.grid.major.y = element_line(linewidth = 0.1, color = 'grey80'),
#        
#        ## Background
#        panel.background = element_blank(),
#        panel.border = element_rect(color = 'grey80', fill = NA, linewidth = 0.1),
#        panel.grid = element_blank(),
#        #panel.grid.major.y = element_line(linewidth = 0.1, color = "black" )
#        
#  ) +
#  scale_x_continuous(
#    breaks = seq(0, 1.6*(10^8), by = (0.4*10^8)),
#    labels = function(x) {
#      round(x / 10^8, 1)
#    }
#  ) +
#  labs(x = expression(paste('Position (', 10 ^ 8, ' bp)')),
#       #y = paste('Density by', WIN_SIZE / 1000000, 'Mb window'))
#       y = 'Variant density'
#  )

# Save to external file
#ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/combined_per_win_density.png',
#       width = 2800,
#       height = 3100,
#       units = 'px',
#       dpi = 700
#)

ggplot(data = ALL_bp_win) +
  facet_grid(factor(TYPE, levels = c('SVs', 'SNPs', 'indels')) ~ CHROM_NUM, 
             scales = 'free', space = 'free_x') +
  geom_point(aes(x = midBIN, y = prop, color = CHROM_NUM),
             #alpha = 1.5 * ALL_bp_win$prop, 
             size = 0.04) +
  theme(
    # Panels and background
    panel.spacing.x = unit(0.6, 'points'),
    panel.spacing.y = unit(3, 'points'),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    # Strips
    strip.text.x.top = element_text(size = 3, 
                                    margin = margin(3,0,3,0, 'pt')),
    strip.text.y.right = element_text(size = 4,
                                      margin = margin(0,1,0,1, 'pt')),
    strip.background.y = element_rect(color = 'black', linewidth = 0.1),
    strip.background.x = element_rect(colour = 'black', linewidth = 0.1),
    
    # Axis
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.line.y = element_line(linewidth = 0.08),
    axis.line.x = element_line(linewidth = 0.08)
    
  ) + 
  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  
  guides(color = 'none') +
  
  scale_color_manual(values = rep(c('black', 'grey60'), 
                                  length(unique(ALL_bp_win$CHROM_NUM))/2)) +
  labs(x = 'Position along each chromosome',
       y = 'Variant density')

  
# Save to external file
#ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/combined_per_win_density.png',
#       width = 2800,
#       height = 3100,
#       units = 'px',
#       dpi = 700,
#        device = 'pdf'
#)

ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/combined_per_win_density_updated.png',
       width = 2800,
       height = 3100,
       units = 'px',
       dpi = 700,
       #device = 'pdf'
)
# meta vs centro
#acro <- subset(ALL_dens_win, CHROM_TYPE == 'acro')
##density_plot_chr_acro <-
#  ggplot(data = acro) +
#  facet_grid(factor(TYPE, levels = c('SVs', 'SNPs', 'indels')) ~ CHROM_NUM, 
#             scales = 'free', space = 'free_x') +
#  geom_point(
#    aes(x = midBIN, y = DENSITY),
#    size = 0.8,
#    alpha = 0.8,
#    col = 'darkorange3'
#  ) +
#  theme(
#    panel.spacing = unit(0.1, 'points'),
#    strip.text.x = element_text(size = 6),
#    axis.text.x = element_text(
#      angle = 45,
#      size = 4,
#      hjust = 1
#    ),
#    panel.background = element_rect(color = "gray70"),
#    strip.placement = "inside",
#    strip.background = element_rect(colour = 'gray70'),
#    legend.position = 'bottom'
#  ) +
#  scale_x_continuous(
#    labels = function(x) {
#      round(x / 10 ^ 8, 1)
#    }
#  ) +
#  labs(x = expression(paste('Position (', 10 ^ 8, ' bp)')),
#       y = paste('Density by', WIN_SIZE / 1000000, 'Mb window'))

#meta <- subset(ALL_dens_win, CHROM_TYPE == 'meta')
##density_plot_chr_meta <-
#  ggplot(data = meta) +
#  facet_grid(. ~ CHROM_NUM, scales = 'free_x', space = 'free_x') +
#  geom_point(
#    aes(x = midBIN, y = DENSITY),
#    size = 0.8,
#    alpha = 0.8,
#    col = 'magenta4'
#  ) +
#  theme(
#    panel.spacing = unit(0.1, 'points'),
#    strip.text.x = element_text(size = 6),
#    axis.text.x = element_text(
#      angle = 45,
#      size = 4,
#      hjust = 1
#    ),
#    panel.background = element_rect(color = "gray70"),
#    strip.placement = "inside",
#    strip.background = element_rect(colour = 'gray70'),
#    legend.position = 'bottom'
#  ) +
#  scale_x_continuous(
#    labels = function(x) {
#      round(x / 10 ^ 8, 1)
#    }
#  ) +
#  labs(x = expression(paste('Position (', 10 ^ 8, ' bp)')),
#       y = paste('Density by', WIN_SIZE / 1000000, 'Mb window')) 
  
  
