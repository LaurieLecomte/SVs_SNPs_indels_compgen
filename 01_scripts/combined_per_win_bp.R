# Plot per window proportion of bp covered, for SVs, SNPs and indels together

library(ggplot2)
library(scales)

# 1. Import and format ----------------------------------------------------
## These files were produced by scripts XX_per_win_bp_prop.sh
SV_BP_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/SVs/SVs_intersect_win100000bp_bp_prop.txt"
SNP_BP_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/SNPs/SNPs_intersect_win100000bp_bp_prop.txt"
INDEL_BP_WIN <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/indels/indels_intersect_win100000bp_bp_prop.txt"


SVs_bp_win <- read.table(SV_BP_WIN, skip = 1, sep = "\t",
                          col.names = c('CHROM_NUM', 'midBIN', 'LEN', 'prop')
)
## Add a column for type, for easier plotting
SVs_bp_win$TYPE <- 'SVs'

SNPs_bp_win <- read.table(SNP_BP_WIN, skip = 1, sep = "\t",
                          col.names = c('CHROM_NUM', 'midBIN', 'LEN', 'prop')
)

SNPs_bp_win$TYPE <- 'SNPs'

indels_bp_win <- read.table(INDEL_BP_WIN, skip = 1, sep = "\t",
                          col.names = c('CHROM_NUM', 'midBIN', 'LEN', 'prop')
)

indels_bp_win$TYPE <- 'indels'

# Bind these 3 datasets together
ALL_bp_win <- do.call("rbind", list(SVs_bp_win, SNPs_bp_win, indels_bp_win))



# 2. Plot  ----------------------------------------------------------------
#ggplot(data = ALL_bp_win) +
#  #facet_wrap(~ CHROM_SSA, nrow = 2, scales = 'free_x') +
#  facet_grid(factor(TYPE, levels = c('SVs', 'SNPs', 'indels')) ~ CHROM_NUM, 
#             scales = 'free', space = 'free_x') +
#  geom_point(aes(x = midBIN, y = prop),
#             #alpha = 1.5 * ALL_bp_win$prop, 
#             size = 0.3, 
#             col = 'midnightblue') +
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
#  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
#  scale_x_continuous(
#    breaks = seq(0, 1.6*(10^8), by = (0.4*10^8)),
#    labels = function(x) {
#      round(x/10^8, 1)
#    }
#  ) + 
#  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
#       y = 'Proportion of base pairs covered') + 
#  guides(color = 'none')

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
    strip.text.x.top = element_text(size = 3.5, 
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
       y = 'Proportion of base pairs covered') 
  

# Save to external file
#ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/combined_per_win_bp_covered.png',
#       width = 2800,
#       height = 3100,
#       units = 'px',
#       dpi = 700
#)

ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/density/combined_per_win_bp_covered_updated.png',
       width = 2800,
       height = 3100,
       units = 'px',
       dpi = 700
      # device = 'pdf'
)

# 3. Compute per win bp ratio between SVs and SNPs ------------------------
# First merge both datasets together to get stats per window
SVs_SNPs_bp_win <- merge(SVs_bp_win, SNPs_bp_win, by = c('CHROM_NUM', 'midBIN'), all = TRUE)

# Compute ratio
SVs_SNPs_bp_win$SVs_bp_to_SNPs_bp <- (SVs_SNPs_bp_win$prop.x / SVs_SNPs_bp_win$prop.y)

# Compute mean ratio
mean(SVs_SNPs_bp_win$SVs_bp_to_SNPs_bp, na.rm = TRUE)
max(SVs_SNPs_bp_win$SVs_bp_to_SNPs_bp, na.rm = TRUE)
min(SVs_SNPs_bp_win$SVs_bp_to_SNPs_bp, na.rm = TRUE)


SVs_indels_bp_win <- merge(SVs_bp_win, indels_bp_win, by = c('CHROM_NUM', 'midBIN'), all = TRUE)
SVs_indels_bp_win$SVs_bp_to_indels_bp <- (SVs_indels_bp_win$prop.x / SVs_indels_bp_win$prop.y)

mean(SVs_indels_bp_win$SVs_bp_to_indels_bp, na.rm = TRUE)
max(SVs_indels_bp_win$SVs_bp_to_indels_bp, na.rm = TRUE)
min(SVs_indels_bp_win$SVs_bp_to_indels_bp, na.rm = TRUE)
