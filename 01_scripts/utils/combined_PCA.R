# Plot PCA of SVs, SNPs and small indels together

library(ggplot2)
library(dplyr)
library(ggpubr)

# 1. Import and format ----------------------------------------------------
# PCA files
SV_PCA <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/09_pca/SVs/merged_SUPP2_MAF0.05_FMISS0.5.cov.pca"
SNP_PCA <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/09_pca/SNPs/SNPs_MAF0.05_FMISS0.5.cov.pca"
INDEL_PCA <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/09_pca/indels/indels_MAF0.05_FMISS0.5.cov.pca"

ID_SEX_POP <- '02_infos/ID_Sex_Pop_updated.txt'

PCA_SVs <- read.table(SV_PCA, 
                      row.names = 1,
                      header = TRUE)
PCA_SNPs <- read.table(SNP_PCA, 
                       row.names = 1,
                       header = TRUE)

PCA_indels <- read.table(INDEL_PCA, 
                         row.names = 1,
                         header = TRUE)

## Add ID sex pop info for coloring points
ID_sex_pop <- read.table(ID_SEX_POP, sep = "\t", header = FALSE) [, 1:3]
colnames(ID_sex_pop) <- c('ID', 'sex', 'pop')


## Assign pop and sex values according to IDs 
POP <- c()
SEX <- c()
ID <- c()

for (i in rownames(PCA_SVs)) {
  # extract this ID in ID_SEX_POP
  sample_df <- subset(ID_sex_pop, ID == i)
  ID[i] <- sample_df$ID
  POP[i] <- sample_df$pop
  SEX[i] <- sample_df$sex
}

PCA_SVs <- cbind(PCA_SVs, ID, POP, SEX)

for (i in rownames(PCA_SNPs)) {
  # extract this ID in ID_SEX_POP
  sample_df <- subset(ID_sex_pop, ID == i)
  ID[i] <- sample_df$ID
  POP[i] <- sample_df$pop
  SEX[i] <- sample_df$sex
}
PCA_SNPs <- cbind(PCA_SNPs, ID, POP, SEX)

for (i in rownames(PCA_indels)) {
  # extract this ID in ID_SEX_POP
  sample_df <- subset(ID_sex_pop, ID == i)
  ID[i] <- sample_df$ID
  POP[i] <- sample_df$pop
  SEX[i] <- sample_df$sex
}
PCA_indels <- cbind(PCA_indels, ID, POP, SEX)

# Import variance (for axis labels)
SV_var <- as.data.frame(t(read.table(file = paste0(unlist(strsplit(SV_PCA, spli = '.pca', fixed = TRUE))[1], '.eig'),
                       row.names = 1, header = TRUE)))
colnames(SV_var) <- paste0('PC', colnames(SV_var))

SNP_var <- as.data.frame(t(read.table(file = paste0(unlist(strsplit(SNP_PCA, spli = '.pca', fixed = TRUE))[1], '.eig'),
                      row.names = 1, header = TRUE)))
colnames(SNP_var) <- paste0('PC', colnames(SNP_var))
SNP_var <- format(round(SNP_var, 2), nsmall=2)

indel_var <- as.data.frame(t(read.table(file = paste0(unlist(strsplit(INDEL_PCA, spli = '.pca', fixed = TRUE))[1], '.eig'),
                        row.names = 1, header = TRUE)))
colnames(indel_var) <- paste0('PC', colnames(indel_var))



# 2. Individual plots -----------------------------------------------------
SVs_PC1_PC2 <- 
  ggplot(data = PCA_SVs, aes(x = round(PC1, 2), y = round(PC2, 2), label = ID)) +
  geom_point(aes(col = POP, shape = POP), size = 0.8) +
  stat_ellipse(linewidth = 0.5, aes(group = POP, col = POP), show.legend = FALSE) + 
  labs(y = paste0("PC2 (", SV_var$PC2, ' %)'), x = paste0("PC1 (", SV_var$PC1, ' %)'),
       color = 'Population', shape = 'Population') +
  scale_color_manual(values = c("red", "blue")) +
  geom_text(
    mapping = aes(x = -Inf, y = Inf, label = 'A'),
    vjust   = 1.2, hjust = -0.5,
    size = 4) + 
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.01)
  ) + 
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.01)
  ) +
  theme(
    ## Axis
    axis.text.x = element_text(angle = 45, size = 5, hjust = 1),
    axis.text.y = element_text(size = 5, hjust = 1),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    ## Legend
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(5, 'mm'),
    ## Background
    panel.background = element_rect(fill = 'gray95'),
    panel.grid.minor.x = element_line(linewidth = 0.3),
    panel.grid.major.x = element_line(linewidth = 0.4),
    panel.grid.minor.y = element_line(linewidth = 0.3),
    panel.grid.major.y = element_line(linewidth = 0.4)
  )

SNPs_PC1_PC2 <- 
  ggplot(data = PCA_SNPs, aes(x = round(PC1, 2), y = round(PC2, 2), label = ID)) +
  geom_point(aes(col = POP, shape = POP), size = 0.8) +
  stat_ellipse(linewidth = 0.5, aes(group = POP, col = POP), show.legend = FALSE) + 
  labs(y = paste0("PC2 (", SNP_var$PC2, ' %)'), x = paste0("PC1 (", SNP_var$PC1, ' %)'),
       color = 'Population', shape = 'Population') +
  scale_color_manual(values = c("red", "blue")) +
  geom_text(
    mapping = aes(x = -Inf, y = Inf, label = 'B'),
    vjust   = 1.2, hjust = -0.5,
    size = 4) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.01)
  ) + 
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.01)
  ) +
  theme(
    ## Axis
    axis.text.x = element_text(angle = 45, size = 5, hjust = 1),
    axis.text.y = element_text(size = 5, hjust = 1),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    ## Legend
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(5, 'mm'),
    ## Background
    panel.background = element_rect(fill = 'gray95'),
    panel.grid.minor.x = element_line(linewidth = 0.3),
    panel.grid.major.x = element_line(linewidth = 0.4),
    panel.grid.minor.y = element_line(linewidth = 0.3),
    panel.grid.major.y = element_line(linewidth = 0.4)
  )


indels_PC1_PC2 <- 
  ggplot(data = PCA_indels, aes(x = round(PC1, 2), y = round(PC2, 2), label = ID)) +
  geom_point(aes(col = POP, shape = POP), size = 0.8) +
  stat_ellipse(linewidth = 0.5, aes(group = POP, col = POP), show.legend = FALSE) + 
  labs(y = paste0("PC2 (", indel_var$PC2, ' %)'), x = paste0("PC1 (", indel_var$PC1, ' %)'),
       color = 'Population', shape = 'Population') +
  scale_color_manual(values = c("red", "blue")) +
  geom_text(
    mapping = aes(x = -Inf, y = Inf, label = 'C'),
    vjust   = 1.2, hjust = -0.5,
    size = 4) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.01)
  ) + 
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.01)
  ) +
  theme(
    ## Axis
    axis.text.x = element_text(angle = 45, size = 5, hjust = 1),
    axis.text.y = element_text(size = 5, hjust = 1),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    ## Legend
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(5, 'mm'),
    ## Background
    panel.background = element_rect(fill = 'gray95'),
    panel.grid.minor.x = element_line(linewidth = 0.3),
    panel.grid.major.x = element_line(linewidth = 0.4),
    panel.grid.minor.y = element_line(linewidth = 0.3),
    panel.grid.major.y = element_line(linewidth = 0.4)
  )


# 3. Combined plot --------------------------------------------------------

ggarrange(SVs_PC1_PC2,
          SNPs_PC1_PC2,
          indels_PC1_PC2,
          ncol = 3, nrow = 1, common.legend = TRUE, legend = 'bottom')

# Save to external file
ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/09_pca/combined_PCA.png',
       width = 3000,
       height = 1400,
       units = 'px',
       dpi = 600,
       bg = "white"
)
