
library(data.table)
library(ggplot2)


# SV - SNP pairs ----------------------------------------------------------
# Import and filter
PLINK_TABLE <- 'SVs_SNPs_concat_100_OV354430.1.ld'
LD_TABLE <- fread(PLINK_TABLE)

## Keep SV-SNP pairs only 
LD_TABLE_filt <- subset(LD_TABLE, (SNP_A == '.' & SNP_B != '.')|(SNP_A != '.' & SNP_B == '.'))


# Compute dist between markers in each pair
LD_TABLE_filt$dist <- LD_TABLE_filt$BP_B - LD_TABLE_filt$BP_A


# Perform bootstrap to get LD50 and CI on each chromosome
n_boot_repeats <- 50 # number of bootstrap iterations

## Initialize empty dataframe for output
dist_decay <- data.frame('CHROM' = chroms,               # chromosome
                         mean_LD_dist = rep(NA, length(chroms)), # mean LD50
                         CI_min = rep(NA, length(chroms)),       # lower CI for mean LD50
                         CI_max = rep(NA, length(chroms)),       # upper CI for mean LD50
                         stringsAsFactors=FALSE) 

## Get max R2 value
### Extract pairs where dist is the smallest
min_dist_pairs <- subset(LD_TABLE_filt, dist == min(LD_TABLE_filt$dist))
### Get mean R2 at min dist
max_R2 <- mean(min_dist_pairs$R2)


## Extract chromosomes
chroms <- unique(LD_TABLE_filt$CHR_A)

## Bootstrap on each chromosome
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
    
    #### Compute LD50 dist : get dist for which R2 is closest to 50% of max R2
    halfR2_dist <- boot_mean$dist[which.min(abs(boot_mean$R2 - max_R2/2))]
    
    halfR2[i] <- halfR2_dist
  }
  
  ### Compute confidence intervals
  dist_decay$mean_LD_dist[which(dist_decay$CHROM == CHR)] <- mean(halfR2) 
  dist_decay$CI_min[which(dist_decay$CHROM == CHR)] <- mean(halfR2) - 1.96*(sd(halfR2)/sqrt(n_boot_repeats))
  dist_decay$CI_max[which(dist_decay$CHROM == CHR)] <- mean(halfR2) + 1.96*(sd(halfR2)/sqrt(n_boot_repeats))
}


# Plot
ggplot(data = LD_TABLE_filt) +
  geom_point(aes(x = dist, y = R2), color = 'grey50', size = 0.5, alpha = 0.2) +
  geom_point(data = boot_mean, aes(x = dist, y = R2), color = 'black', size = 0.8) +
  geom_smooth(aes(x = dist, y = R2)) +
  geom_vline(xintercept = dist_decay$mean_LD_dist, linewidth = 1, color = 'firebrick')+
  geom_vline(xintercept = dist_decay$CI_min, linewidth = 0.5, linetype = 2, color = 'firebrick') +
  geom_vline(xintercept = dist_decay$CI_max, linewidth = 0.5, linetype = 2, color = 'firebrick')



# SNP - SNP pairs ---------------------------------------------------------
# Import and filter
PLINK_SNP <- 'SNPs_vs_SNPs/SNPs_SNPs_no_overlap_100_OV354430.1.ld'
LD_SNPs <- fread(PLINK_SNP)

# Compute dist between markers in each pair
LD_SNPs$dist <- LD_SNPs$BP_B - LD_SNPs$BP_A


# Perform bootstrap to get LD50 and CI on each chromosome
n_boot_repeats <- 50 # number of bootstrap iterations

## Initialize empty dataframe for output
dist_decay <- data.frame('CHROM' = chroms,               # chromosome
                         mean_LD_dist = rep(NA, length(chroms)), # mean LD50
                         CI_min = rep(NA, length(chroms)),       # lower CI for mean LD50
                         CI_max = rep(NA, length(chroms)),       # upper CI for mean LD50
                         stringsAsFactors=FALSE) 

## Get max R2 value
### Extract pairs where dist is the smallest
min_dist_pairs <- subset(LD_SNPs, dist == min(LD_SNPs$dist))
### Get mean R2 at min dist
max_R2 <- mean(min_dist_pairs$R2)


## Extract chromosomes
chroms <- unique(LD_SNPs$CHR_A)

## Bootstrap on each chromosome
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
    
    #### Compute LD50 dist : get dist for which R2 is closest to 50% of max R2
    halfR2_dist <- boot_mean$dist[which.min(abs(boot_mean$R2 - max_R2/2))]
    
    halfR2[i] <- halfR2_dist
  }
  
  ### Compute confidence intervals
  dist_decay$mean_LD_dist[which(dist_decay$CHROM == CHR)] <- mean(halfR2) 
  dist_decay$CI_min[which(dist_decay$CHROM == CHR)] <- mean(halfR2) - 1.96*(sd(halfR2)/sqrt(n_boot_repeats))
  dist_decay$CI_max[which(dist_decay$CHROM == CHR)] <- mean(halfR2) + 1.96*(sd(halfR2)/sqrt(n_boot_repeats))
}


# Plot
ggplot(data = LD_SNPs) +
  geom_point(aes(x = dist, y = R2), color = 'grey50', size = 0.5, alpha = 0.2) +
  geom_point(data = boot_mean, aes(x = dist, y = R2), color = 'black', size = 0.8) +
  geom_smooth(aes(x = dist, y = R2)) +
  geom_vline(xintercept = dist_decay$mean_LD_dist, linewidth = 1, color = 'firebrick')+
  geom_vline(xintercept = dist_decay$CI_min, linewidth = 0.5, linetype = 2, color = 'firebrick') +
  geom_vline(xintercept = dist_decay$CI_max, linewidth = 0.5, linetype = 2, color = 'firebrick')
