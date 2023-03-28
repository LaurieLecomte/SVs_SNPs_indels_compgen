# Plot per site Fst for SVs, SNPs and indels at the same time

library(data.table)
library(ggplot2)
library(dplyr)


# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)

SV_TABLE <- argv[1] 
SNP_TABLE <- argv[2]
INDEL_TABLE <- argv[3]
SV_GW_FST <- argv[4] 
SNP_GW_FST <- argv[5]
INDEL_GW_FST <- argv[6]

OV_2_SSA <- argv[7]



SV_TABLE <- "08_angsd_fst/SVs/merged_SUPP2_MAF0.05_FMISS0.5.SVsFst_RO_PU.table"
SNP_TABLE="08_angsd_fst/SNPs/SNPs_MAF0.05_FMISS0.5.SNPsFst_RO_PU.table"
INDEL_TABLE="08_angsd_fst/indels/indels_MAF0.05_FMISS0.5.indelsFst_RO_PU.table"

SV_GW_FST="08_angsd_fst/SVs/RO_PU/RO_PU.SVs.fst"
SNP_GW_FST="08_angsd_fst/SNPs/RO_PU/RO_PU.SNPs.fst"
INDEL_GW_FST="08_angsd_fst/indels/RO_PU/RO_PU.indels.fst"
OV_2_SSA <- '02_infos/OV_to_ssa.txt'

# Convertion table between OV and SSA chromosomes
OV_2_ssa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))


# SVs
SVs_table <- fread(SV_TABLE, 
                        col.names = c('CHROM', 'POS', 'ID', 'FST'),
                        na.strings = '.'
)

## Convert OV to ssa
SVs_table <- merge(x = SVs_table, y = OV_2_ssa, by.x = 'CHROM', by.y = 'CHROM_OV')

## Mean genome wide
SVs_mean_fst <- read.table(SV_GW_FST, col.names = c('unweighted_FST', 'weighted_FST'))

## Add value for type
SVs_table$TYPE <- 'SV'


# SNPs
SNPs_table <- fread(SNP_TABLE, 
                   col.names = c('CHROM', 'POS', 'ID', 'FST'),
                   na.strings = '.'
)
## Convert OV to ssa
SNPs_table <- merge(x = SNPs_table, y = OV_2_ssa, by.x = 'CHROM', by.y = 'CHROM_OV')

## Add value for type
SNPs_table$TYPE <- 'SNP/indel'

## Mean genome wide
SNPs_mean_fst <- read.table(SNP_GW_FST, col.names = c('unweighted_FST', 'weighted_FST'))


# Indels
indels_table <- fread(INDEL_TABLE, 
                    col.names = c('CHROM', 'POS', 'ID', 'FST'),
                    na.strings = '.'
)

## Convert OV to ssa
indels_table <- merge(x = indels_table, y = OV_2_ssa, by.x = 'CHROM', by.y = 'CHROM_OV')

## Mean genome wide
indels_mean_fst <- read.table(INDEL_GW_FST, col.names = c('unweighted_FST', 'weighted_FST'))

## Add value for type
indels_table$TYPE <- 'SNP/indel'

# Bind 3 datasets together
ALL_types <- rbind(SNPs_table, indels_table, SVs_table)

# Convert negative Fst to 0
ALL_types$FST <- ifelse(ALL_types$FST < 0, 
                        yes = 0, 
                        no = ALL_types$FST)

# 2. Plot -----------------------------------------------------------------
# Assign colors
types <- c('SNP/indel', 'SV')

#ALL_types$TYPE_reordered <- factor(ALL_types$TYPE, 
#                                    levels = types)


## Get hex code for as many colors as svtypes for a given viridis palette
hex_types <- viridisLite::viridis(n = length(types), option = 'D')
scales::show_col(hex_types)

## Assign a color to each svtype in a named vector
cols_types <- vector(mode = 'character', length = length(types))

for (i in 1:length(types)) {
  names(cols_types)[i] <- types[i]
  cols_types[i] <- hex_types[i]
}



# Plot
fst_per_site_plot <- 
  ggplot(data = ALL_types) +
  facet_wrap(~ CHROM_SSA, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = POS, y = FST, col = TYPE, shape = TYPE), alpha = 1.5 * ALL_types$FST, size = 1.2*ALL_types$FST) +
  theme(panel.spacing = unit(0.1, 'points'),
        strip.text.x = element_text(size = 6),
        axis.text.x = element_text(angle = 45, size = 4, hjust = 1),
        panel.background = element_rect(color = "gray60"),
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
  scale_color_manual(values = cols_types) #+ 
  #geom_hline(yintercept = SVs_mean_fst$weighted_FST, col = cols_types[which(names(cols_types[1]) == 'SV')]) +
  #geom_hline(yintercept = SNPs_mean_fst$weighted_FST, col = cols_types[which(names(cols_types[1]) == 'SNP')]) +
  #geom_hline(yintercept = indels_mean_fst$weighted_FST, col = cols_types[which(names(cols_types[1]) == 'indel')])
