#Impute

library(dplyr)
library(data.table)
library(vegan)

# Access to files provided in command line arguments
argv <- commandArgs(T)
FILE <- argv[1] # "$MAT_file"
ID_SEX_POP <- argv[2]
CHROM_POS_END_ID <- argv[3]
#pca_path <- argv[2] # $pca_path


# 1. Load SNPs info
## Create a new column for SNP info name (CHR + position AND unique ID)
geno.012.pos <- as.data.frame(
  fread(paste0(FILE, ".pos"),
        col.names = c('CHROM', 'POS')))

print(nrow(geno.012.pos))

## match with IDs
chr_pos_end_ID <- as.data.frame(fread(CHROM_POS_END_ID,
                             col.names = c('CHROM', 'POS', 'END', 'ID')))

geno.012.pos <- cbind(geno.012.pos, ID = chr_pos_end_ID$ID)

geno.012.pos$locus <- paste0(geno.012.pos$CHROM, '_', geno.012.pos$POS, '-', geno.012.pos$ID)

#geno.012.pos <- read.table(paste0(FILE, ".pos")) %>% #read genotype_matrix.pos; 1st column = CHR, 2nd column = position
#  mutate(., locus = paste(V1, V2, sep = '_')) #paste CHR (V1) and position (V2) together and add this new "locus" column to genotype_matrix.pos


# 2. Load individuals info
geno.012.indv <- read.table(paste0(FILE, ".indv"), col.names = c("ID")) #read genotype_matrix.indv
#geno.012.indv <- read.table(paste0(FILE, ".pos")) 
N_SAMPLES <- nrow(geno.012.indv)

# 3. Load genotype matrix 
geno.012 <- fread(FILE, showProgress = TRUE)[, -1] #read genotype_matrix and exclude 1st column from .012 matrix file 

## Set rownames and colnames for the genotype matrix
colnames(geno.012) <- geno.012.pos$locus #each column is a SNP
rownames(geno.012) <- geno.012.indv$ID #each row is a different indv

## Inspect matrix
geno.012[1:12, 1:20] #show 12 first indv and 20 first SNPs 


# 4. Imputation of missing genotypes : replace missing geno (-1) with the most frequent genotype (0, 1, 2) for given SNP
# For ALL samples, not by pop
geno.012.imp <- apply(
  geno.012, 2, function(x){
    replace(
      x, 
      which(x == "-1"), 
      max(
        0, 
        as.numeric(
          names(
            which.max(
              table(x)     # genotypes frequencies
            )
          )
        )
      )
    )
  }
)


## Save as table to .imp file
write.table(geno.012.imp, paste0(FILE, ".imp"), sep = "\t", row.names = F, quote = F)
