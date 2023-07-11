# Plot per window proportion of bp covered, for SVs, SNPs and indels together

library(ggplot2)

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
ggplot(data = ALL_bp_win) +
  #facet_wrap(~ CHROM_SSA, nrow = 2, scales = 'free_x') +
  facet_grid(factor(TYPE, levels = c('SVs', 'SNPs', 'indels')) ~ CHROM_NUM, 
             scales = 'free', space = 'free_x') +
  geom_point(aes(x = midBIN, y = prop),
             #alpha = 1.5 * ALL_bp_win$prop, 
             size = 1.2*ALL_bp_win$prop, 
             col = 'darkblue') +
  theme(panel.spacing = unit(0.2, 'points'),
        strip.text.x = element_text(size = 5),
        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        axis.text.y = element_text(size = 8),
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
       y = 'Proportion of base pairs covered') + 
  guides(color = FALSE)
# +
#geom_hline(yintercept = mean_fst$weighted_FST, col = 'red')


# 3. Compute per win bp ratio between SVs and SNPs ------------------------
# First merge both datasets together to get stats per window
SVs_SNPs_bp_win <- merge(SVs_bp_win, SNPs_bp_win, by = c('CHROM_NUM', 'midBIN'), all = TRUE)

# Compute ratio
SVs_SNPs_bp_win$SVs_bp_to_SNPs_bp <- (SVs_SNPs_bp_win$prop.x / SVs_SNPs_bp_win$prop.y)

# Compute mean ratio
mean(SVs_SNPs_bp_win$SVs_bp_to_SNPs_bp, na.rm = TRUE)
max(SVs_SNPs_bp_win$SVs_bp_to_SNPs_bp, na.rm = TRUE)
min(SVs_SNPs_bp_win$SVs_bp_to_SNPs_bp, na.rm = TRUE)
