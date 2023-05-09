# Combine Manhattan plots for SVs, SNPs and indels together

library(ggplot2)
library(ggpubr)

# 1. Access files in command line, import and format ----------------------
SV_PLOT <- "08_angsd_fst/SVs/merged_SUPP2_MAF0.05_FMISS0.5.SVsFst_RO_PU.manh_plot.rds"
SNP_PLOT <- "08_angsd_fst/SNPs/SNPs_MAF0.05_FMISS0.5.SNPsFst_RO_PU.manh_plot.rds"
INDEL_PLOT <- "08_angsd_fst/indels/indels_MAF0.05_FMISS0.5.indelsFst_RO_PU.manh_plot.rds"
     
SV_plot <- readRDS(SV_PLOT)
SV_plot_updated <- SV_plot + theme(axis.title.x = element_blank(), 
                                   axis.text.x = element_blank(),
                                   axis.title.y = element_blank())

SNP_plot <- readRDS(SNP_PLOT)
SNP_plot_updated <- SNP_plot + theme(axis.title.x = element_blank(), 
                                     axis.text.x = element_blank())
indel_plot <- readRDS(INDEL_PLOT)
indel_plot_updated <- indel_plot + theme(axis.title.y = element_blank(),
                                         axis.text.x = element_text(angle = 45, size = 6, hjust = 1))
 

# 2. Combine plots --------------------------------------------------------
combined <- 
ggarrange(SV_plot_updated, 
          SNP_plot_updated, 
          indel_plot_updated, 
          ncol = 1, nrow = 3, common.legend = TRUE, legend = 'bottom')

ggsave(combined, filename = 'SVs_SNPs_indels_outliers_Fst0.97_fisher0.01_RDA3sd.jpg', width = 8, height = 10)
