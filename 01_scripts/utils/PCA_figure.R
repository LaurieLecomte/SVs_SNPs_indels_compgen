# Produce 3 panels PCA figure

library(ggpubr)
library(ggplot2)

# 1. Access files in command line, import and format ----------------------
#argv <- commandArgs(T)
#SVs_PCA <- argv[1]
#SNPs_PCA <- argv[2]
#INDELs_PCA <- argv[3]

SVs_PCA <- readRDS("09_pca/SVs/merged_SUPP2_MAF0.05_FMISS0.5_PC1_PC2_formatted.rds")
indels_PCA <- readRDS("09_pca/indels/indels_MAF0.05_FMISS0.5_PC1_PC2_formatted.rds")
SNPs_PCA <- readRDS("09_pca/SNPs/SNPs_MAF0.05_FMISS0.5_PC1_PC2_formatted.rds")

ggarrange(SVs_PCA,
          SNPs_PCA,
          indels_PCA,
          ncol = 3, nrow = 1, common.legend = TRUE, legend = 'right')

ggarrange(merged_SUPP2_MAF0.05_FMISS0.5_PC1_PC2,
          SNPs_MAF0.05_FMISS0.5.maf0.05_PC1_PC2,
          indels_MAF0.05_FMISS0.5.maf0.05_PC1_PC2,
          ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

SNPs_MAF0.05_FMISS0.5.maf0.05_PC1_PC2 + 
  ylim(-0.5, 0.5) + xlim(-0.5, 0.5) +
  stat_ellipse(linewidth = 0.5, aes(group = POP, col = POP), show.legend = FALSE) +
  annotate("text", x = -0.2, y = 0.4, label = "A")



# 2. PC3 ------------------------------------------------------------------
SVs_PCA <- readRDS("09_pca/SVs/merged_SUPP2_MAF0.05_FMISS0.5_PC1_PC3_formatted.rds")
indels_PCA <- readRDS("09_pca/indels/indels_MAF0.05_FMISS0.5_PC1_PC3_formatted.rds")
SNPs_PCA <- readRDS("09_pca/SNPs/SNPs_MAF0.05_FMISS0.5_PC1_PC3_formatted.rds")

