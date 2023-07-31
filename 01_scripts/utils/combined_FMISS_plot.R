# Plot distrib of missing data, for SVs, SNPs and indels together

library(ggplot2)
library(ggpp)
library(data.table)

# 1. Import and format ----------------------------------------------------
## These files were produced by first generating a table of F_MISS per site
## with bcftools +fill-tags -- -t F_MISSING, 
## then running compute_DP_FMISS_stats.R scripts 
## in the genotype_SVs_SRLR and SNPs_indels_SR pipelines

SV_FMISS <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5_DP_FMISS_FMISS_stats.table'
SNP_FMISS <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/SNPs/SNPs_MAF0.05_FMISS0.5_DP_FMISS_FMISS_stats.table'
INDEL_FMISS <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/03_SNP/SNPs_indels_SR/07_filtered/indels/indels_MAF0.05_FMISS0.5_DP_FMISS_FMISS_stats.table'

SVs_FMISS <- fread(SV_FMISS, sep = "\t",
                        col.names = c('CHROM', 'POS', 'DP', 'F_MISSING',
                                      'F_MISS_bins', 'F_MISS_groups', 'TYPE'))

SNPs_FMISS <- fread(SNP_FMISS, sep = "\t",
                        col.names = c('CHROM', 'POS', 'DP', 'F_MISSING',
                                      'F_MISS_bins', 'F_MISS_groups', 'TYPE'))

indels_FMISS <- fread(INDEL_FMISS, sep = "\t",
                         col.names = c('CHROM', 'POS', 'DP', 'F_MISSING',
                                       'F_MISS_bins', 'F_MISS_groups', 'TYPE'))


# 2. Bind these 3 datasets together and format ----------------------------
ALL_FMISS <- rbind(SVs_FMISS, SNPs_FMISS, indels_FMISS)

ALL_FMISS$F_MISS_bins <- cut_interval(ALL_FMISS$F_MISSING, length = 0.1, right = FALSE)


# 3. Plot  ----------------------------------------------------------------
ggplot(data = ALL_FMISS) +
  #facet_grid(.~factor(TYPE, levels = c('SVs', 'SNPs', 'indels')), 
  facet_grid(factor(TYPE, levels = c('SVs', 'SNPs', 'Indels')) ~ ., 
             scales = 'free_y') +
  geom_bar(aes(F_MISS_bins, fill = F_MISS_groups), color = 'black', linewidth = 0.1) + 
  theme(
    ## Axis
    axis.text.x = element_text(angle = 45, size = 6, hjust = 1),
    axis.text.y = element_text(size = 6, hjust = 1),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    strip.text.y.right = element_text(size = 7),
    ## Legend
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 7),
    legend.key.size = unit(5, 'mm')
    
  ) + scale_fill_manual(values = c('red', 'grey60')) +
  labs(x = 'Proportion of missing genotypes',
       y = 'Count',
       fill = 'F_MISSING filter') +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

# Save to external file
ggsave(filename = "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/FMISS_plot_by_type.png", 
       width = 2600,
       height = 2800,
       units = 'px',
       dpi = 600
)

