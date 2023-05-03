# Adapted from original work of Florent Sylvestre

#perform per snps intersex fisher test on allele count to detect signifcantly differentiated snps between male and female
# need a genotype file format 0 1 2 (from vcftools for exemple) and a popmap
#popmap is a tab delimieted file indicating samples of first collumn and sex (M or F) on second
#Usage : Rscript Fishert_test.R genofile popmap, no header


library(data.table)
library(ggplot2)

# 1. Access files in command line, import and format ----------------------
genofile <- commandArgs(trailingOnly = TRUE)[1]
popfile <- commandArgs(trailingOnly = TRUE)[2]
CHR_POS_END_ID <- commandArgs(trailingOnly = TRUE)[3]

pop1 <- commandArgs(trailingOnly = TRUE)[4]
pop2 <- commandArgs(trailingOnly = TRUE)[5]
MAX_QVAL <- commandArgs(trailingOnly = TRUE)[6]

out_file <- commandArgs(trailingOnly = TRUE)[7]

genotype <- data.frame(fread(genofile, stringsAsFactors = FALSE, header = FALSE))[,-1]
genopos <- data.frame(fread(paste(genofile, ".pos", sep=""), stringsAsFactors = FALSE, header = FALSE, 
                            col.names = c('CHROM', 'POS')))
genoind <- data.frame(fread(paste(genofile, ".indv", sep=""), stringsAsFactors = FALSE, header = FALSE,
                            col.names = c('ID')))
popmap <- data.frame(fread(popfile, stringsAsFactors = FALSE,header = FALSE,
                           col.names = c('ID', 'POP')))

chr_pos_end_ID <- data.frame(fread(CHR_POS_END_ID, stringsAsFactors = FALSE,header = FALSE, 
                                   col.names = c('CHROM', 'POS', 'END', 'ID') ))



# 2. Perform Fisher test --------------------------------------------------
# Define function
fisher_function <- function(genotype, pop){
  genocount <- matrix(data = rep(0,4),ncol=2)
  colnames(genocount) <- c(pop1, pop2)
  rownames(genocount) <- c("A","R")  
  
  for(i in 1:length(genotype)){
    if(genotype[i] == -1){
      next()
    } else if(genotype[i] == 0){
      genocount["R", pop[i] ] <-  genocount["R", pop[i] ] + 2
      
    }else if(genotype[i] == 2){
      genocount["A", pop[i] ] <-  genocount["A", pop[i] ] + 2
      
    }else if(genotype[i] == 1){
      genocount["A", pop[i] ] <-  genocount["A", pop[i] ] + 1
      genocount["R", pop[i] ] <-  genocount["R", pop[i] ] + 1
      
    }
  }
  if(sum(genocount[1,]) == 0) { return(NA)}else{
    return(fisher.test(genocount)$p.value)}
  
}

# Assign pop based on popmap
genoind$pop <- merge(genoind, popmap, sort = FALSE)[, 2]

# Get pval for each site
genopos$P_VAL <- unlist(apply(genotype, 2 , fisher_function, pop = genoind$pop))
genopos$Q_VAL <- p.adjust(genopos$P_VAL, method = "BH")

# Merge with original sites list, including END field
genopos <- merge(genopos, chr_pos_end_ID, by = c('CHROM', 'POS'), sort = FALSE)

#write.table(genopos[, c('CHROM', 'POS', 'END', 'ID', 'P_VAL', 'Q_VAL')], 
#            file = paste0(out_dir, '/', unlist(strsplit(basename(genofile), split = '.geno_mat.012'))[1], ".fisher.txt"), sep = "\t",
#            row.names = FALSE, quote = FALSE)


write.table(genopos[, c('CHROM', 'POS', 'END', 'ID', 'P_VAL', 'Q_VAL')], 
                        file = paste0(out_file, '_all.txt'), sep = "\t",
                        row.names = FALSE, quote = FALSE)

# Filter for Q_VAL > threshold
#genopos_outliers <- subset(genopos, Q_VAL < MAX_QVAL)

#write.table(genopos_outliers[, c('CHROM', 'POS', 'END', 'ID', 'P_VAL', 'Q_VAL')], 
#            file = paste0(out_file, '_outliers_qval', MAX_QVAL, '.txt'), sep = "\t",
#            row.names = FALSE, quote = FALSE)

# Plot Q_VAL distribution
q_val_distrib <- ggplot(data = genopos) + 
  geom_histogram(aes(x = Q_VAL), binwidth = 0.025) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      size = 8,
      hjust = 1
    ),
    plot.title = element_text(size = 11)
  ) +
  labs(x = 'Q_VAL', y = 'Count', title = 'Q_VAL distribution') +
  theme(
    plot.title = element_text(size = 10, face = 'bold', hjust = 0.5)
  )

jpeg(file = paste0(out_file, '_qval_distrib.jpg'))
q_val_distrib
dev.off()
