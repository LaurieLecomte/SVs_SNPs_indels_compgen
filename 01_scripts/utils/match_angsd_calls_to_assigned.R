ANGSD_SV <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SVs/merged_SUPP2_MAF0.05_FMISS0.5.maf0.05.SVsFst_RO_PU.table"

MATCHED_CAND <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/genotype_SVs_SRLR/09_filtered/merged_SUPP2_MAF0.05_FMISS0.5_matched_offset5bp.txt"

SVs_set <- read.delim(ANGSD_SV, header=FALSE,
                      col.names = c('CHROM', 'POS', 'END', 'VG_ID', 'FST'))
matched_SVs <- read.delim(MATCHED_CAND, header = TRUE)

matched_angsd_SVs <- merge(x = SVs_set, y = matched_SVs, by = c('CHROM', 'POS'),
                           all.x = TRUE)