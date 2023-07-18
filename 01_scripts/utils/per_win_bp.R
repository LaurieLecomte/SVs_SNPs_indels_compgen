# Compute proportion of bp covered by variant by window

library(ggplot2)

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
WIN_OVERLAP <- argv[1]
OV_2_SSA <- argv[2]
WIN_SIZE <- as.numeric(argv[3])

# Import
win_table <- read.table(WIN_OVERLAP, col.names = c('WIN_CHROM', 'BIN_START', 'BIN_END', 'CHROM', 'POS', 'END', 'ID'))

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

win_table <- merge(x = win_table, y = OV_2_sasa, by.x = 'CHROM', by.y = 'CHROM_OV')


# 2. Compute variant size and middle position -----------------------------
win_table$LEN <- ifelse(win_table$POS != win_table$END,
                        yes = win_table$END - win_table$POS,
                        no = 1)

win_table$midBIN <- win_table$BIN_START + WIN_SIZE/2



# 3. Compute number of bp covered by window -------------------------------
bp_per_win <- aggregate(data = win_table, 
          LEN ~ CHROM_NUM + midBIN, 
          FUN = sum)

bp_per_win$prop <- bp_per_win$LEN/WIN_SIZE

# Export
write.table(bp_per_win, file = paste0(unlist(strsplit(WIN_OVERLAP, split = '.txt')), '_bp_prop.txt'),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

# Compute mean
print(paste('mean number of bp covered by',  WIN_SIZE, 'bp window :', mean(bp_per_win$LEN)))
print(paste('mean proportion of bp covered by',  WIN_SIZE, 'bp window :', mean(bp_per_win$prop)))


# 4. Plot -----------------------------------------------------------------
prop_plot <- 
ggplot(data = bp_per_win) +
  facet_grid(.~CHROM_NUM, scales = 'free_x', space = 'free_x') +
  geom_point(aes(x = midBIN, y = prop), alpha = 0.6, col = 'darkblue') +
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
  ylim(0, ifelse(test = max(bp_per_win$prop) > 0.2,
                 yes = 0.75,
                 no = 0.2)) + 
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = paste('Base pairs covered by', WIN_SIZE/1000000, 'Mb window'))

ggsave(prop_plot, filename = paste0(unlist(strsplit(WIN_OVERLAP, split = '.txt')), '_bp_prop_plot.jpg'), 
       width = 10, height = 4)
saveRDS(prop_plot, file = paste0(unlist(strsplit(WIN_OVERLAP, split = '.txt')), '_bp_prop_plot.rds'))
