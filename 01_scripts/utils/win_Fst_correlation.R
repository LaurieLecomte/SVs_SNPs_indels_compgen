
# 1. Access files, import and format ----------------------
SV_WIN_FST <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SVs/RO_PU/RO_PU_win1000000_step10000.txt"
SNP_WIN_FST <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/SNPs/RO_PU/RO_PU_win1000000_step10000.txt"
INDEL_WIN_FST <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/08_angsd_fst/indels/RO_PU/RO_PU_win1000000_step10000.txt"

SVs_fst_win <- read.table(SV_WIN_FST, skip = 1, sep = "\t",
                      col.names = c('REGION', 'CHROM', 'midPOS', 'Nsites', 'FST')
)
SNPs_fst_win <- read.table(SNP_WIN_FST, skip = 1, sep = "\t",
                          col.names = c('REGION', 'CHROM', 'midPOS', 'Nsites', 'FST')
)
indels_fst_win <- read.table(INDEL_WIN_FST, skip = 1, sep = "\t",
                          col.names = c('REGION', 'CHROM', 'midPOS', 'Nsites', 'FST')
)


# 2. Merge all datasets pairs to have the same sites ----------------------
SVs_SNPs_win <- merge(SVs_fst_win, SNPs_fst_win, by = c('CHROM', 'midPOS'), all.x = TRUE)
SVs_indels_win <- merge(SVs_fst_win, indels_fst_win, by = c('CHROM', 'midPOS'), all.x = TRUE)
indels_SNPs_win <- merge(indels_fst_win, SNPs_fst_win, by = c('CHROM', 'midPOS'), all.x = TRUE)


# 3. Perform correlation tests --------------------------------------------
cor.test(SVs_SNPs_win$FST.y, SVs_SNPs_win$FST.x,  method = "pearson", use = "complete.obs")
cor.test(SVs_indels_win$FST.y, SVs_indels_win$FST.x,  method = "pearson", use = "complete.obs")
cor.test(indels_SNPs_win$FST.y, indels_SNPs_win$FST.x,  method = "pearson", use = "complete.obs")

