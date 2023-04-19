

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
FST_OUTLIERS <- argv[1] 
RDA_OUTLIERS <- argv[2] 
FISHER_OUTLIERS <- argv[3] 

#RDA_OUTLIERS <- "12_go/SVs/SVs_RO_PU_outliers_RDA_overlap10000bp_outlierIDs.txt"
RDA_OUTLIERS <- "12_go/SVs/SVs_RO_PU_outliers_RDA_3sd_overlap10000bp_outlierIDs.txt"
FST_OUTLIERS <- "12_go/SVs/SVs_RO_PU_outliers_minFst0.28_overlap10000bp_outlierIDs.txt"
FISHER_OUTLIERS <- "12_go/SVs/SVs_fisher_RO_PU_outliers_qval0.01_overlap10000bp_outlierIDs.txt"

outliers_FST <- read.delim(FST_OUTLIERS, header = FALSE, col.names = c('FST_genes'))
outliers_RDA <- read.delim(RDA_OUTLIERS, header = FALSE, col.names = c('RDA_genes'))
outliers_Fisher <- read.delim(FISHER_OUTLIERS, header = FALSE, col.names = c('FISHER_genes'))

# 2. Compare outlier genes sets between Fst, RDA and Fisher ---------------
# Union of all 3 methods
union_3sets <- Reduce(union, c(outliers_FST$FST_genes, outliers_RDA$RDA_genes, outliers_Fisher))

# Intersection of all 3 methods
intersect_3sets <- Reduce(intersect, c(outliers_FST$FST_genes, outliers_RDA$RDA_genes, outliers_Fisher))

# Shared between Fst and RDA
intersect_Fst_RDA <- intersect(outliers_RDA$RDA_genes, outliers_FST$FST_genes)
print(paste(round(length(intersect_Fst_RDA)/nrow(outliers_RDA), 3), 'of RDA outliers are shared with Fst'))
print(paste(round(length(intersect_Fst_RDA)/nrow(outliers_FST), 3), 'of FST outliers are shared with RDA'))

# Unique to RDA
RDA_uniques <- setdiff(outliers_RDA$RDA_genes, outliers_FST$FST_genes)
print(paste(round(length(RDA_uniques)/nrow(outliers_RDA), 3), 'of RDA outliers are unique'))
writeLines(text = RDA_uniques, con = paste0(unlist(strsplit(RDA_OUTLIERS, split = '.txt')), '_uniques.txt'))


# Shared between Fst and Fisher
intersect_Fst_Fisher <- intersect(outliers_FST$FST_genes, outliers_Fisher$FISHER_genes)
print(paste(round(length(intersect_Fst_Fisher)/nrow(outliers_FST), 3), 'of Fst outliers are shared with Fisher'))
print(paste(round(length(intersect_Fst_Fisher)/nrow(outliers_Fisher), 3), 'of Fisher outliers are shared with Fst'))

writeLines(text = intersect_Fst_Fisher, con = paste0(unlist(strsplit(FST_OUTLIERS, split = '.txt')), 
                                                     '_shared_Fisher_qval', MAX_QVAL, '.txt'))

# Shared between RDA and Fisher
intersect_RDA_Fisher <- intersect(outliers_RDA$RDA_genes, outliers_Fisher$FISHER_genes)
print(paste(round(length(intersect_RDA_Fisher)/nrow(outliers_RDA), 3), 'of RDA outliers are shared with Fisher'))
print(paste(round(length(intersect_RDA_Fisher)/nrow(outliers_Fisher), 3), 'of Fisher outliers are shared with RDA'))

