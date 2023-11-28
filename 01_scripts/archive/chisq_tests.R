library(dplyr)
outliers <- c(2079, 71953, 6223)

RDA_cand <- c(1280, 114637, 13871)
shared_3 <- c(43, 1231, 100)
RDA_cand_no_shared <- RDA_cand - shared_3

total_cands <- outliers + RDA_cand_no_shared

sets_2_groups <- data.frame('outliers_not_shared' = outliers - shared_3, 'RDA_not_shared' = RDA_cand - shared_3, 'TOT' = total_cands - shared_3)

#sets_2_all <- data.frame('outliers' = outliers, 'RDA cand' = RDA_cand)


# chi square test WITHOUT shared across 3 methods
#sets_2_uniques <- sets_2_all - shared_3
chisq.test(as.matrix(select(sets_2_groups, -TOT)))





# chi square test WITH shared across 3 methods
sets_3 <- data.frame('outliers_not_shared' = (outliers - shared_3), 'RDA_not_shared' = RDA_cand_no_shared, 'shared_3' = shared_3)
sets_3$TOT <- rowSums(sets_3)
chisq.test(as.matrix(select(sets_3, -TOT)))


# Variants vs genes distribution
variants_all <- c(115907, 8777832, 1089321)
variants_near_genes <- c(76241, 5527078, 740372)
 
variants_nogenes <- variants_all - variants_near_genes
var_genes <- data.frame('near_gene' = variants_near_genes, 'no_genes' = variants_nogenes,
                        row.names = c('SVs', 'SNPs', 'indels'))
                        
chisq.test(as.matrix(var_genes))


# chi square test of adequacy outliers vs all genes
## Get proportion of variants located within 10 kb of genes
expected_near_genes <- var_genes$near_gene/variants_all

outliers_all <- (outliers - shared_3)
#outliers_near_genes <- c(1364, 47164, 4331)
outliers_near_genes <- c(1334, 46337, 4263)
outliers_nogenes <- outliers_all - outliers_near_genes

chi_outliers <- vector(mode = 'numeric', length = 3)
pval_outliers <- vector(mode = 'numeric', length = 3)

for (i in seq(1,3)){
  chi <- 
    chisq.test(
      c(outliers_near_genes[i], outliers_nogenes[i]),
      p = c(expected_near_genes[i], 1 - expected_near_genes[i])
      )
  chi_outliers[i] <- chi$statistic
  pval_outliers[i] <- chi$p.value
}


#chisq.test(
#  c(outliers_near_genes[1], outliers_nogenes[1]),
#  p = c(expected_near_genes[1], 1 - expected_near_genes[1])
#)
#chisq.test(
#  c(outliers_near_genes[2], outliers_nogenes[2]),
#  p = c(expected_near_genes[2], 1 - expected_near_genes[2])
#)
#chisq.test(
#  c(outliers_near_genes[3], outliers_nogenes[3]),
#  p = c(expected_near_genes[3], 1 - expected_near_genes[3])
#)



RDA_cand_all <- RDA_cand - shared_3
#RDA_cand_near_genes <- c(848, 76080, 9803)
RDA_cand_near_genes <- c(818, 75253, 9735)
RDA_cand_nogenes <- RDA_cand_all - RDA_cand_near_genes

chi_cand <- vector(mode = 'numeric', length = 3)
pval_cand <- vector(mode = 'numeric', length = 3)

for (i in seq(1,3)){
  chi <- 
    chisq.test(
      c(RDA_cand_near_genes[i], RDA_cand_nogenes[i]),
      p = c(expected_near_genes[i], 1 - expected_near_genes[i])
    )
  chi_cand[i] <- chi$statistic
  pval_cand[i] <- chi$p.value
}


#chisq.test(
#  c(RDA_cand_near_genes[1], RDA_cand_nogenes[1]),
#  p = c(expected_near_genes[1], 1 - expected_near_genes[1])
#)
#chisq.test(
#  c(RDA_cand_near_genes[2], RDA_cand_nogenes[2]),
#  p = c(expected_near_genes[2], 1 - expected_near_genes[2])
#)
#chisq.test(
#  c(RDA_cand_near_genes[3], RDA_cand_nogenes[3]),
#  p = c(expected_near_genes[3], 1 - expected_near_genes[3])
#)

p.adjust(c(pval_outliers, pval_cand))
