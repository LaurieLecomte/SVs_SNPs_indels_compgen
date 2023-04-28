# Get intersection of Fst and Fisher outliers to produce a single "outliers" set, 
# and extract candidates unique to RDA (for further validation)

# 1. Access files in command line, import and format ----------------------
## Files
argv <- commandArgs(T)
FST_OUTLIERS <- argv[1] 
RDA_CAND <- argv[2] 
FISHER_OUTLIERS <- argv[3] 

## Cutoffs and values
OVERLAP_WIN <- as.numeric(argv[4])
QUANTILE <- as.numeric(argv[5])
MIN_FST <- as.numeric(argv[6])
SD <- as.numeric(argv[7])
MAX_QVAL <- as.numeric(argv[8])

FST_FISHER_OUTPUT <- argv[9]
RDA_UNIQUE_OUTPUT <- argv[10]

outliers_Fst <- read.table(FST_OUTLIERS, col.names = c('CHROM', 'POS', 'END', 'ID', 'FST'))
cand_RDA <- read.table(RDA_CAND, col.names = c('CHROM', 'POS', 'END', 'ID', 'loadings'))
outliers_Fisher <- read.table(FISHER_OUTLIERS, col.names = c('CHROM', 'POS', 'END', 'ID', 'PVAL', 'QVAL'))


# 2. Fst vs Fisher outliers -----------------------------------------------
# Merge by CHROM and POS 
## We cannot only do intersect on IDs only, because it would not work with SNPs and indels that do not have a unique ID
shared_Fst_Fisher <- merge(outliers_Fst, outliers_Fisher, by = c('CHROM', 'POS', 'END', 'ID'), all = FALSE, sort = FALSE)

print(paste(nrow(outliers_Fst), 'Fst outliers with min Fst =', MIN_FST, 'for quantile', QUANTILE))
print(paste(nrow(outliers_Fisher), 'Fisher outliers with max qval =', MAX_QVAL))

print(paste(
  round( (nrow(shared_Fst_Fisher) / nrow(outliers_Fst)), 2),
  'of Fst outliers are also Fisher outliers'
  )
)

print(paste(
  round( (nrow(shared_Fst_Fisher) / nrow(outliers_Fisher)), 2),
  'of Fisher outliers are also Fst outliers'
  )
)

print(paste(
  nrow(shared_Fst_Fisher), ' shared outliers between Fst and Fisher = high confidence outliers set'
)
)


# Export
write.table(shared_Fst_Fisher[, c('CHROM', 'POS', 'END', 'ID')], 
            file = FST_FISHER_OUTPUT,
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")



# 3. Get candidates unique to RDA -----------------------------------------
# Get RDA candidates shared with Fst outliers
shared_Fst_RDA <- merge(outliers_Fst, cand_RDA, by = c('CHROM', 'POS', 'END', 'ID'), all = FALSE, sort = FALSE)

print(paste(nrow(cand_RDA), 'RDA candidates with', SD, 'SDs'))

print(paste(
  round( (nrow(shared_Fst_RDA) / nrow(outliers_Fst)), 2),
  'of Fst outliers are also RDA candidates'
)
)
print(paste(
  round( (nrow(shared_Fst_RDA) / nrow(cand_RDA)), 2),
  'of RDA candidates are also Fst outliers'
)
)

# Get RDA candidates shared with Fisher outliers
shared_Fisher_RDA <- merge(outliers_Fisher, cand_RDA, by = c('CHROM', 'POS', 'END', 'ID'), all = FALSE, sort = FALSE)
print(paste(
  round( (nrow(shared_Fisher_RDA) / nrow(outliers_Fisher)), 2),
  'of Fisher outliers are also RDA candidates'
)
)

print(paste(
  round( (nrow(shared_Fisher_RDA) / nrow(cand_RDA)), 2),
  'of RDA candidates are also Fisher outliers'
)
)

RDA_shared <- merge(shared_Fisher_RDA, shared_Fst_RDA, by = c('CHROM', 'POS', 'END', 'ID'), all = TRUE )

# How many RDA candidates are also Fisher AND Fst outliers ?
print(paste(nrow(RDA_shared[! is.na(RDA_shared$QVAL) & ! is.na(RDA_shared$FST), ]), 
            'RDA candidates are also Fst AND Fisher outliers'))

# Remove shared RDA candidates from RDA candidates list by performing anti-join
RDA_uniques <- dplyr::anti_join(cand_RDA, RDA_shared, by = c('CHROM', 'POS', 'END', 'ID'))

print(paste(
  round( (nrow(RDA_uniques) / nrow(cand_RDA)), 3),
  'of RDA candidates are uniques'
)
)

print(paste0(
  round( (nrow(RDA_shared) / nrow(cand_RDA)), 3), ' (', nrow(RDA_shared), ') ',
  'of RDA candidates are also either a Fst outlier or a Fisher outlier'
)
)

# Export
write.table(RDA_uniques[, c('CHROM', 'POS', 'END', 'ID')], 
            file = RDA_UNIQUE_OUTPUT,
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
