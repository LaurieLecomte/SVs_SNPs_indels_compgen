library(data.table)
library(ggplot2)

setwd("~/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen")

LD_TABLE <- fread('~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SVs_vs_SNPs/LD_SVs_SNPs_OV354430.1.txt',
                  col.names = c('CHROM', 'POS', 'END', 'ID', 'SNP_POS', 'D', 'Rsq', 'dist'))

# Subset for testing
LD_table_small <- head(LD_TABLE, 50000)

# Summarize by SV



ggplot(data = LD_TABLE) +
  geom_point(aes(x = dist, y = Rsq))


LD_TABLE_SNPs <- '~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SVs_vs_SNPs/LD_SNPs_vs_SNPs_test_OV354430.1.txt'

LD_TABLE <- fread(LD_TABLE_SNPs,
                  col.names = c('CHROM', 'POS', 'END', 'ID', 'SNP_POS', 'D', 'Rsq', 'dist'))

PLINK_TABLE <- '~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SVs_vs_SNPs/LD_SVs_SNPs_OV354430.1.ld'
LD_TABLE <- fread(PLINK_TABLE)

# Remove SNPs-SNPs and SVs-SVs comparisons
LD_TABLE_filt <- subset(LD_TABLE, (SNP_A == '.' & SNP_B != '.')|(SNP_A != '.' & SNP_B == '.'))

LD_TABLE_filt$dist <- LD_TABLE_filt$BP_B-LD_TABLE_filt$BP_A
ggplot(data = LD_TABLE_filt) +
  geom_point(aes(x = dist, y = R2))

# Extract SV from each comparison
LD_TABLE_filt$SV_ID <- ifelse(LD_TABLE_filt$SNP_A != '.',
                              yes = LD_TABLE_filt$SNP_A,
                              no = LD_TABLE_filt$SNP_B)

aggregate(R2 ~ SV_ID, data = LD_TABLE_filt, FUN = var)

LD_TABLE_SNPs <- subset(LD_TABLE, (SNP_A == '.' & SNP_B == '.'))
LD_TABLE_SNPs$dist <- LD_TABLE_SNPs$BP_B-LD_TABLE_SNPs$BP_A
ggplot(data = LD_TABLE_SNPs) +
  geom_point(aes(x = dist, y = R2))







PLINK_TABLE <- '~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SVs_vs_SNPs/SVs_SNPs_concat_100_OV354430.1.ld'
LD_TABLE <- data.table::fread(PLINK_TABLE)

# Remove SNPs-SNPs and SVs-SVs comparisons
LD_TABLE_filt <- subset(LD_TABLE, (SNP_A == '.' & SNP_B != '.')|(SNP_A != '.' & SNP_B == '.'))


# Compute dist in each pair
LD_TABLE_filt$dist <- LD_TABLE_filt$BP_B-LD_TABLE_filt$BP_A


LD_decay_plot <- 
ggplot(data = LD_TABLE_filt) + xlim(c(0, 50000))+
  geom_point(aes(x = dist, y = R2, color = CHR_A)) 



mean_r2 <- aggregate(R2 ~ dist + CHR_A, data = LD_TABLE_filt, FUN = mean)
ggplot(data = mean_r2) +
  geom_point(aes(x = dist, y = R2, color = CHR_A))


# Bootstrap 
n_boot_repeats <- 20 


CHR <- 'OV354429.1'




# For each chromosome :

# Extract SVs-SNPs pairs for given CHR
chr_sub <- subset(LD_TABLE_filt, CHR_A == CHR)

# Initialize list
halfR2 <- vector(mode = 'numeric', length = n_boot_repeats) # we store one mean LD50 dist per bootstrap iteration

# Bootstrap 
for (i in 1:n_boot_repeats){
  # Subsamble CHR dataset for about 75% of dataset size
  lines_to_sample <- sort(sample(1:nrow(chr_sub), 
                                 size = ceiling(0.75 * nrow(chr_sub)), 
                                 replace = TRUE))
  boot_table <- chr_sub[lines_to_sample, ]
  
  # Compute r^2 mean per dist value 
  boot_mean <- aggregate(R2 ~ dist + CHR_A, data = boot_table, FUN = mean)
  
  # Compute max R2 at dist = 0
  #max_R2 <- boot_mean$R2[boot_mean$dist == 0]
  max_R2 <- boot_mean$R2[which.min(boot_mean$dist)]
  
  # Compute half decay dist : get dist for R2 closest to half of max R2
  halfR2_dist <- boot_mean$dist[which.min(abs(boot_mean$R2 - max_R2/2))]

  halfR2[i] <- halfR2_dist
}

# Compute confidence intervals

conf_int_min <- mean(halfR2) - 1.96*(sd(halfR2)/sqrt(n_boot_repeats))
conf_int_max <- mean(halfR2) + 1.96*(sd(halfR2)/sqrt(n_boot_repeats))

# travailler sur vecteur de valeurs de dist obtenues après bootstrap ()






# Get LD decay distance at 50% of max r^2 ---------------------------------

# Import and filter
PLINK_TABLE <- '~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SVs_vs_SNPs/SVs_SNPs_concat_100_OV354430.1.ld'
LD_TABLE <- fread(PLINK_TABLE)

## Remove SNPs-SNPs and SVs-SVs comparisons
LD_TABLE_filt <- subset(LD_TABLE, (SNP_A == '.' & SNP_B != '.')|(SNP_A != '.' & SNP_B == '.'))

## Compute dist between markers in each pair
LD_TABLE_filt$dist <- LD_TABLE_filt$BP_B-LD_TABLE_filt$BP_A


# Calculate LD50 based on 100 bootstrap iterations for each chromosome
n_boot_repeats <- 50 # number of bootstrap iterations

## Extract chromosomes
chroms <- unique(LD_TABLE_filt$CHR_A)

## Initialize empty dataframe for output
dist_decay <- data.frame('CHROM' = chroms,               # chromosome
                 mean_LD_dist = rep(NA, length(chroms)), # mean LD50
                 CI_min = rep(NA, length(chroms)),       # lower CI for mean LD50
                 CI_max = rep(NA, length(chroms)),       # upper CI for mean LD50
                 stringsAsFactors=FALSE) 

#### From per dist means, Compute mean R2 at smallest distance
min_dist_pairs <- subset(LD_TABLE_filt, dist == min(LD_TABLE_filt$dist))

max_R2 <- mean(min_dist_pairs$R2)

## Loop over chromosomes
for (CHR in chroms){
  
  ### Extract SVs-SNPs pairs for given CHR
  chr_sub <- subset(LD_TABLE_filt, CHR_A == CHR)
  
  ### Initialize list for storing each iteration's results
  halfR2 <- vector(mode = 'numeric', length = n_boot_repeats)
  
  ### Bootstrap iterations 
  for (i in 1:n_boot_repeats){
    #### Subsamble CHR dataset for about 90% of dataset size
    lines_to_sample <- sort(sample(1:nrow(chr_sub), 
                                   size = ceiling(0.9 * nrow(chr_sub)), 
                                   replace = TRUE))
    boot_table <- chr_sub[lines_to_sample, ]
    
    #### Compute r^2 mean per dist value 
    boot_mean <- aggregate(R2 ~ dist + CHR_A, data = boot_table, FUN = mean)
    
    #### From per dist means, Compute mean R2 at smallest distance
    #max_R2 <- boot_mean$R2[which.min(boot_mean$dist)]
    #max_R2 <- boot_mean$R2[boot_mean$dist == 2]
    
    
    #### Compute half decay dist : get dist for which R2 closest to half of max R2
    halfR2_dist <- boot_mean$dist[which.min(abs(boot_mean$R2 - max_R2/2))]
    
    halfR2[i] <- halfR2_dist
  }
  
  ### Compute confidence intervals
  dist_decay$mean_LD_dist[which(dist_decay$CHROM == CHR)] <- mean(halfR2) 
  dist_decay$CI_min[which(dist_decay$CHROM == CHR)] <- mean(halfR2) - 1.96*(sd(halfR2)/sqrt(n_boot_repeats))
  dist_decay$CI_max[which(dist_decay$CHROM == CHR)] <- mean(halfR2) + 1.96*(sd(halfR2)/sqrt(n_boot_repeats))
}

LD50_plot <- 
ggplot(data = LD_TABLE_filt) +
  geom_point(aes(x = dist, y = R2), color = 'grey50', size = 0.5, alpha = 0.2) +
  geom_point(data = boot_mean, aes(x = dist, y = R2), color = 'black', size = 0.8) +
  geom_smooth(aes(x = dist, y = R2)) +
  geom_vline(xintercept = dist_decay$mean_LD_dist, linewidth = 1, color = 'firebrick')+
  geom_vline(xintercept = dist_decay$CI_min, linewidth = 0.5, linetype = 2, color = 'firebrick') +
  geom_vline(xintercept = dist_decay$CI_max, linewidth = 0.5, linetype = 2, color = 'firebrick')


# Import and filter
PLINK_SNP <- '~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SNPs_vs_SNPs/SNPs_SNPs_no_overlap_100_OV354430.1.ld'
LD_SNPs <- fread(PLINK_SNP)

## Remove SNPs-SNPs and SVs-SVs comparisons
#LD_SNPs_filt <- subset(LD_SNPs, (SNP_A == '.' & SNP_B != '.')|(SNP_A != '.' & SNP_B == '.'))

## Compute dist between markers in each pair
LD_SNPs$dist <- LD_SNPs$BP_B-LD_SNPs$BP_A


# Calculate LD50 based on 100 bootstrap iterations for each chromosome
n_boot_repeats <- 20 # number of bootstrap iterations

## Extract chromosomes
chroms <- unique(LD_SNPs$CHR_A)

## Initialize empty dataframe for output
dist_decay <- data.frame('CHROM' = chroms,               # chromosome
                         mean_LD_dist = rep(NA, length(chroms)), # mean LD50
                         CI_min = rep(NA, length(chroms)),       # lower CI for mean LD50
                         CI_max = rep(NA, length(chroms)),       # upper CI for mean LD50
                         stringsAsFactors=FALSE) 

#### From per dist means, Compute mean R2 at smallest distance
min_dist_pairs <- subset(LD_SNPs, dist == min(LD_SNPs$dist))

max_R2 <- mean(min_dist_pairs$R2)

## Loop over chromosomes
for (CHR in chroms){
  
  ### Extract SVs-SNPs pairs for given CHR
  chr_sub <- subset(LD_SNPs, CHR_A == CHR)
  
  ### Initialize list for storing each iteration's results
  halfR2 <- vector(mode = 'numeric', length = n_boot_repeats)
  
  ### Bootstrap iterations 
  for (i in 1:n_boot_repeats){
    #### Subsamble CHR dataset for about 90% of dataset size
    lines_to_sample <- sort(sample(1:nrow(chr_sub), 
                                   size = ceiling(0.9 * nrow(chr_sub)), 
                                   replace = TRUE))
    boot_table <- chr_sub[lines_to_sample, ]
    
    #### Compute r^2 mean per dist value 
    boot_mean <- aggregate(R2 ~ dist + CHR_A, data = boot_table, FUN = mean)
    
    #### From per dist means, Compute max R2 at smallest distance
    max_R2 <- boot_mean$R2[which.min(boot_mean$dist)]
    #max_R2 <- boot_mean$R2[boot_mean$dist == 2]
    
    #### Compute half decay dist : get dist for which R2 closest to half of max R2
    halfR2_dist <- boot_mean$dist[which.min(abs(boot_mean$R2 - max_R2/2))]
    
    halfR2[i] <- halfR2_dist
  }
  
  ### Compute confidence intervals
  dist_decay$mean_LD_dist[which(dist_decay$CHROM == CHR)] <- mean(halfR2) 
  dist_decay$CI_min[which(dist_decay$CHROM == CHR)] <- mean(halfR2) - 1.96*(sd(halfR2)/sqrt(n_boot_repeats))
  dist_decay$CI_max[which(dist_decay$CHROM == CHR)] <- mean(halfR2) + 1.96*(sd(halfR2)/sqrt(n_boot_repeats))
}

#LD50_plot <- 
  #ggplot(data = LD_SNPs) +
  #geom_point(aes(x = dist, y = R2), color = 'grey', fill = 'white', size = 0.5, alpha = 0.2, shape = 21) +
  ggplot(data = boot_mean)+
  geom_point(aes(x = dist, y = R2), color = 'black', size = 0.8) +
  geom_smooth(aes(x = dist, y = R2)) +
  geom_vline(xintercept = dist_decay$mean_LD_dist, linewidth = 1, color = 'firebrick')+
  geom_vline(xintercept = dist_decay$CI_min, linewidth = 0.5, linetype = 2, color = 'firebrick') +
  geom_vline(xintercept = dist_decay$CI_max, linewidth = 0.5, linetype = 2, color = 'firebrick')

ggsave(filename = '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SNPs_vs_SNPs/SNPs_vs_SNPs.png',
       width = 2800,
       height = 2800,
       units = 'px',
       dpi = 700
       # device = 'pdf'
)















mean_r2 <- aggregate(R2 ~ dist + CHR_A, data = LD_TABLE_filt, FUN = mean)

LD_TABLE_filt_chr1 <- subset(LD_TABLE_filt, CHR_A == 'OV354430.1')
mean_r2_chr1 <- subset(mean_r2, CHR_A == 'OV354430.1')

ggplot(data = LD_TABLE_filt_chr1) +
 # facet_wrap(~CHR_A, nrow = 4) +
  geom_point(aes(x = dist, y = R2), color = 'grey', size = 0.5) +
  #geom_point(data = mean_r2_chr1, aes(x = dist, y = R2)) +
  #stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE) 
  geom_smooth(method = "nls",
              formula = y ~ a*x^(-b),
              method.args = list(start=c(a=20, b=0.01)),  # 
              se = F,
              data = mean_r2_chr1, aes(x = dist, y = R2))


  geom_smooth(data = mean_r2_chr1, aes(x = dist, y = R2), method = "glm",
              method.args = list(family = gaussian(link = 'log'), start = c(-0.5, 1)))
              
              
  geom_smooth(method ="lm", formula = y ~ poly(x,2), data = mean_r2_chr1, aes(x = dist, y = R2))


  geom_smooth(data = mean_r2_chr1, aes(x = dist, y = R2), method = "glm", family = gaussian(lin="log"))
              #formula = y ~ a * x^b, 
              #se = FALSE, method.args = list(start = c(a = 1, b = 1)))

 stat_smooth(data = mean_r2_chr1, aes(x = dist, y = R2), method = "nls")
 
 
 PLINK_TABLE <- '~/projects/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SVs_vs_SNPs/SVs_SNPs_concat_1000_OV354430.1.ld'
 
