# Plot Fst distribution

library(ggplot2)

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
INPUT_TABLE <- argv[1]
OV_2_SASA <- argv[2]

fst_table <- read.table(INPUT_TABLE, 
                        col.names = c('CHROM', 'POS', 'ID', 'FST'),
                        na.strings = '.'
                        )

fst_table$FST <- abs(fst_table$FST)

# Split Fst by bins
## Convert candidate SVs length to num and bins

fst_bins <- seq(0, 1, by = 0.01)

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

ggplot(data = fst_table) + 
  geom_bar(aes(x = Fst_bin)) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      size = 4,
      hjust = 1
    ),
    plot.title = element_text(size = 11)
    )

# By site
# OV to sasa
OV_2_sasa <- read.table(OV_2_SASA, col.names = c('OV', 'SSA'))

fst_table_ssa <- merge(x = fst_table, y = OV_2_sasa, by.x = 'CHROM', by.y = 'OV')


ggplot(data = fst_table_ssa) +
  facet_wrap(~ SSA, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = POS, y = FST), alpha = 0.7, size = 2*fst_table_ssa$FST) +
  theme(panel.spacing = unit(0.1, 'points'),
        strip.text.x = element_text(size = 6),
        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        panel.background = element_rect(color = "gray60"),
        strip.placement = "inside",
        strip.background = element_rect(colour = 'gray60')
  ) + 
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
  scale_x_continuous(
    labels = function(x) {
      round(x/10^8, 1)
    }
  ) + 
  scale_color_manual(values = c('blue', 'red')) + 
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = 'Fst') 
