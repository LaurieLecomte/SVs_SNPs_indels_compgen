# Calculate proportion of genome covered by SVs, indels and SNPs

library(data.table)

options(scipen=999)

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
SV_TABLE <- argv[1]
SNP_TABLE <- argv[2]
INDEL_TABLE <- argv[3]
GENOME_BP <- as.numeric(argv[4])
SV_MATCHED <- argv[5]

#SV_TABLE <- "04_vcf/SVs/merged_SUPP2.candidates.table"
#SV_TABLE <- "04_vcf/SVs/merged_SUPP2_MAF0.05_FMISS0.5_CHROM_POS_END_ID.table"
#SNP_TABLE <- "04_vcf/SNPs/SNPs_MAF0.05_FMISS0.5_CHROM_POS_END.table"
#INDEL_TABLE <- "04_vcf/indels/indels_MAF0.05_FMISS0.5_CHROM_POS_END.table"
#SV_MATCHED <- "04_vcf/SVs/merged_SUPP2_MAF0.05_FMISS0.5_matched_offset5bp.txt"

#SVs <- fread(SV_TABLE, header = TRUE)[, c(1,2,6:8)]
#colnames(SVs) <- c('CHROM', 'POS', 'SVTYPE', 'SVLEN', 'END')

SVs <- fread(SV_TABLE, col.names = c('CHROM', 'POS', 'END', 'ID'))

SNPs <- fread(SNP_TABLE, col.names = c('CHROM', 'POS', 'END'))

indels <- fread(INDEL_TABLE, col.names = c('CHROM', 'POS', 'END'))

# For SVs only, we'll assign a type by merging with the list 
# of genotyped SVs matched with a known candidate
# This is done for taking INS length into account when possible
matched_SVs <- fread(SV_MATCHED, header = TRUE)

SVs <- merge(SVs, matched_SVs, by = c('CHROM', 'POS', 'ID'), all.x = TRUE)


# 2. Compute variants length for indels and SVs ---------------------------
indels$LEN <- indels$END - indels$POS
SVs$LEN <- SVs$END - SVs$POS

# 3. Compute number of bases covered by each type -------------------------
## By SVs, all types
#bp_SVs <- sum(abs(SVs$SVLEN))

## For each SV type 
#SV_types <- unique(SVs$SVTYPE)
#print('bp py SV type')
#sapply(X = SV_types, FUN = function(x){
#  sum(abs(SVs$SVLEN[SVs$SVTYPE == x]))
#})

bp_SVs_noINS <- sum(abs(SVs$LEN)) # this excludes INSs as their END and POS is identical
bp_SVs <- sum(abs(SVs$CAND_SVLEN), na.rm = TRUE)

bp_SNPs <- nrow(SNPs) # 1 bp per SNP

bp_indels <- sum(abs(indels$LEN))


# 4. Compute proportion of genome consisting of each type -----------------
prop_SVs <- bp_SVs/GENOME_BP
prop_SNPs <- bp_SNPs/GENOME_BP
prop_indels <- bp_indels/GENOME_BP

print(paste('Genome size :', GENOME_BP))

print(paste('SVs : ', bp_SVs, 'bp,', round(prop_SVs, 5), 'of genome'))
#print('Prop by SV type :')
#sapply(X = SV_types, FUN = function(x){
#  sum(abs(SVs$SVLEN[SVs$SVTYPE == x]))/GENOME_BP
#})

print(paste('SNPs : ', bp_SNPs, 'bp,', round(prop_SNPs, 5), 'of genome'))
print(paste('Indels : ', bp_indels, 'bp,', round(prop_indels, 5), 'of genome'))
