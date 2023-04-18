# Perform RDA on a imputed 012 genotype matrix

library(dplyr)
library(data.table)
library(vegan)
library(ggplot2)

# 1. Access to files from command line and import -------------------------
argv <- commandArgs(T)
FILE <- argv[1] # "$MAT_file"
ID_SEX_POP <- argv[2]
CHR_POS_END_ID <- argv[3]
SD <- as.numeric(argv[4])
OUT_DIR <- argv[5]

GEN_MAT_PATH <- paste0(FILE, ".imp")
INDV_FILE <- paste0(FILE, '.indv')

# Imputed 012 matrix 
geno.012.imp <- fread(GEN_MAT_PATH, header = TRUE, showProgress = TRUE) 
geno.012.imp <- as.data.frame(geno.012.imp)
print('done importing genotype matrix')

## Check missing data (should be 0)
sum(geno.012.imp == '-1')

# Import 012 matrix samples IDs
gen_mat_ID <- fread(INDV_FILE, header = FALSE, col.names = 'ID')


# Predictors dataframe : ID, sex and pop
ID_sex_pop <- read.table(ID_SEX_POP, header = FALSE,
                         col.names = c("ID", "sex", "pop"))

#merge(gen_mat_ID, ID_sex_pop, by = 'ID')

## Reorder predictors IDs according to gen_mat_ID
ID_sex_pop <- ID_sex_pop[match(gen_mat_ID$ID, ID_sex_pop$ID), ]

## Add IDs as row names for genotype matrix
rownames(geno.012.imp) <- gen_mat_ID$ID



# 2. Run RDA --------------------------------------------------------------

print('starting RDA')
pop.rda <- rda(geno.012.imp ~ ID_sex_pop$pop, scale = TRUE)

print(paste('completed RDA. saving to ', OUT_DIR))
saveRDS(pop.rda, file = paste0(OUT_DIR, '/pop.rda.rds'))

###pop.rda <- readRDS('rda/pop_rda.rds')

# Calculate R^2
print('Calculating R^2')
RsquareAdj(pop.rda) # r.squared=0.0181 # adj.r.squared=0.00118
summary(pop.rda)$concont

## export table
write.table(RsquareAdj(pop.rda), paste0(OUT_DIR, "/adjR2_pop.txt"), row.names = FALSE, quote = FALSE)
#pdf(file=paste0(OUT_DIR, "/screeplot.pdf"))

#dev.off()


# 3. Identify candidate sites ---------------------------------------------

# Extract loadings
load.rda <- summary(pop.rda)$species[, 1:3]
#load.rda <- scores(pop.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
print('done extracting species')

pdf(file = paste0(OUT_DIR, "/hist_loadings_pop.pdf"))
par(mfrow = c(1, 3))
hist(load.rda[, 1], main = "Loadings on RDA1")
hist(load.rda[, 2], main = "Loadings on RDA2")
hist(load.rda[, 3], main = "Loadings on RDA3")
dev.off()
## export table
write.table(load.rda[, 1:3], paste0(OUT_DIR, "/rda_loading_pop.txt"), quote=FALSE)

# Outlier detection function
print('starting outlier detection')
outliers <- function(x, z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     ## find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               ## locus names in these tails
}

# Detect outliers for RDA1 and number of SDs specified in args
## 3 sd is standard, 2 sd is less conservative
cand1 <- outliers(load.rda[, 1], SD) # candidates SNPs for RDA1
cand2 <- outliers(load.rda[, 2], SD) # candidates SNPs for RDA2
cand3 <- outliers(load.rda[, 3], SD) 

## Get number of candidate outlier SNPs
ncand <- length(cand1) + length(cand2) + length(cand3)
cat(paste('number of candidate outlier sites : ', ncand, "\n"))

#lim <- 3
#if (ncand == 0){
#  cand1 <- outliers(load.rda[, 1], 2.5) 
#  cand2 <- outliers(load.rda[, 2], 2.5) 
#  cand3 <- outliers(load.rda[, 3], 2.5) 
#  ncand<-length(cand1) + length (cand2) + length ( cand3) # total nb of snps
#  lim <- 2.5
#  print("no outlier site with 3 sd, outlier detection re-done with 2.5 sd")
#}

## Store candidate variants with axis, variant  name, loading (=RDA1), & correlation with each predictor
### For each axis/RDA
cand1 <- cbind.data.frame(axis = rep(1, times = length(cand1)), snp = names(cand1), loading = unname(cand1))
cand2 <- cbind.data.frame(axis = rep(2, times = length(cand2)), snp = names(cand2), loading = unname(cand2))
cand3 <- cbind.data.frame(axis = rep(3, times = length(cand3)), snp = names(cand3), loading = unname(cand3))
#cand3 <- cbind.data.frame(rep(3,times = length(cand3)), names(cand3), unname(cand3))
#colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis", "snp", "loading")

## Bind all together
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)
cand$CHROM <- sapply(X = cand$snp, FUN = function(x) unlist(strsplit(x, split = '_'))[1])
cand$POS <- sapply(X = cand$snp, FUN = function(x) as.numeric(unlist(strsplit(x, split = '_'))[2]))
#cand$END <- cand$POS + 1

## Add variant ID (relevant for SVs only)
chr_pos_end_ID <- read.table(CHR_POS_END_ID, header = FALSE, col.names = c('CHROM', 'POS', 'END', 'ID'))

cand <- merge(cand, chr_pos_end_ID, by = c('CHROM', 'POS'))

write.table(cand[, c('CHROM', 'POS', 'END', 'ID', 'loading')], file = paste0(OUT_DIR, '/RDA_', SD, 'sd_outliers.txt'), sep = "\t", quote = FALSE, row.names = FALSE)

## Check for duplicates
length(cand$snp[duplicated(cand$snp)])
cand <- cand[!duplicated(cand$snp), ]


# 4. Plot RDA outlier sites -----------------------------------------------

# sites at middle (red), samples are black, vectors are predictors
print('plotting')
#jpeg(file = paste0(OUT_DIR, "/rda_1-2_1-3_pop.jpg"))
par(mfrow = c(1, 2)) 
plot(pop.rda, scaling = 3)
text(pop.rda, display = "sites", col = 1, scaling = 3)
plot(pop.rda, choices = c(1,3), scaling = 3) # axis 1 and 3
text(pop.rda, display = "sites", col = 1, scaling = 3)
#dev.off()

print('done')

# All sites : load.rda
all_sites <- as.data.frame(load.rda[, 1:2])
all_sites$site <- row.names(load.rda)

all_sites$CHROM <- sapply(X = all_sites$site, FUN = function(x) unlist(strsplit(x, split = '_'))[1])
all_sites$POS <- sapply(X = all_sites$site, FUN = function(x) as.numeric(unlist(strsplit(x, split = '_'))[2]))

all_sites <- merge(all_sites, chr_pos_end_ID, by = c('CHROM', 'POS'))
write.table(all_sites[, c('CHROM', 'POS', 'END', 'ID', 'RDA1', 'PC1')], 
            file = paste0(OUT_DIR, '/RDA_', SD, 'sd_all_sites.txt'), sep = "\t", quote = FALSE, row.names = FALSE)


all_sites_RDA1 <- 
ggplot(data = all_sites) +
  facet_wrap(~ CHROM, nrow = 2, scales = 'free_x') +
  geom_point(aes(x = POS, y = abs(RDA1)), alpha = 0.5, size = 0.5) +
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
  labs(x = expression(paste('Position (', 10^8, ' bp)' )),
       y = 'abs(RDA1)') 


#jpeg(file = paste0(OUT_DIR, "/plot_RDA1.jpg"))
all_sites_RDA1

saveRDS(all_sites_RDA1, file = paste0(OUT_DIR, "/plot_RDA1_", SD, "sd.rds"))
#dev.off()

#jpeg(file = paste0(OUT_DIR, "/plot_RDA1_outliers.jpg"))

all_sites_RDA1_outliers <- all_sites_RDA1 + 
  geom_point(data = cand, aes(x = POS, y = abs(loading)), col = 'red', size = 0.5) 
all_sites_RDA1_outliers

saveRDS(all_sites_RDA1_outliers, file = paste0(OUT_DIR, "/plot_RDA1_", SD, "sd_outliers.rds"))
dev.off()

# 5. Check RDA signifiance using ANOVA ------------------------------------
print('checking sgnifiance using ANOVA')

signif.full <- anova.cca(pop.rda, parallel = getOption("mc.cores")) # default is permutation=999
#signif.full <- anova.cca(pop.rda, parallel = CORES)
#signif.full <- anova.cca(pop.rda)
signif.full
saveRDS(signif.full, file = paste0(OUT_DIR, '/signif_full_', SD, 'sd_pop.rds'))

signif.axis <- anova.cca(pop.rda, by = "axis", parallel = getOption("mc.cores"))
#signif.axis <- anova.cca(pop.rda, by = "axis", parallel = CORES)
#signif.axis <- anova.cca(pop.rda, by = "axis")
signif.axis
saveRDS(signif.axis, file = paste0(OUT_DIR, '/signif_axis_', SD, 'sd_pop.rds'))



# 6. Add correlation with each predictor ----------------------------------
# Add correlation with each predictor : NOT RELEVANT FOR A SINGLE CATEGORICAL PREDICTOR
# n <- dim(ID_sex_pop)[2]
#foo <- matrix(nrow = (ncand), ncol = 1)  # ncol = number of predictors, nrow = nrow(cand)
#colnames(foo) <- c("sex")

#for (i in 1:length(cand$snp)) {
#  nam <- cand[i, 2] # SNP's name
#  snp.gen <- geno.012.imp[, ..nam] # # extract all GT for a given SNP
#  #foo[i, ] <- apply(ID_sex_pop, 2, function(x) cor(x, snp.gen)) # calculate correlation for each ID_SEX_POP column
#  foo[i, ] <- cor(ID_sex_pop$sex, snp.gen)
#}
#table of candidate snp with loading on each axis and correlation with env predictors
#cand <- cbind.data.frame(cand,foo)  
#head(cand)

#cand <- cbind.data.frame(cand,foo)  
#head(cand)

# 8. Investigate candidate SNPs
# Any SNPs associated with several axis? if yes remove them
#n_dupli <- length(cand$snp[duplicated(cand$snp)])
#n_dupli
#if (n_dupli >= 1){cand <- cand[!duplicated(cand$snp)]}

# Find the most stronly correlated predictor
#n < -dim(cand)[2]
#for (i in 1:length(cand$snp)) {
#  bar <- cand[i,]
#  cand[i, n+1] <- names(which.max(abs(bar[4:n]))) # gives the variable
#  cand[i, n+2] <- max(abs(bar[4:n]))              # gives the correlation
#}

#colnames(cand)[n+1] <- "predictor"
#colnames(cand)[n+2] <- "correlation"
#head(cand)
#table(cand$predictor) 
#write.table(cand, paste0(OUTPUT_FOLDER,ENV,"_candidate_SNP_",lim,"_sd.txt"), quote=FALSE, sep=" ", row.names=FALSE)
