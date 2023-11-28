# Edit PLINK output to add info on SV type

library(data.table)

# Get list of chromosomes
CHR_LIST <- '/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/02_infos/chr_list.txt'

chrs <- readLines(CHR_LIST)

# Get list of matched SVs, with info on candidate SVs' type and size
MATCHED_SV <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/04_vcf/SVs/merged_SUPP2_MAF0.05_FMISS0.5_matched_offset5bp.txt"

matched_SVs <- fread(MATCHED_SV, header = TRUE)

# Get list of genotyped + filtered SVs to have END position
GENO_FILT_SV <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/04_vcf/SVs/merged_SUPP2_MAF0.05_FMISS0.5.bed"
geno_SVs <- fread(GENO_FILT_SV, header = FALSE, col.names = c('CHROM', 'POS', 'END', 'ID'))

matched_SVs <- merge(matched_SVs, geno_SVs[, c('END', 'ID')], by = c('ID'), all.x = TRUE)


# Loop over chromosomes ---------------------------------------------------
for (CHR in chrs) {
  # Import
  PLINK_OUTPUT <- paste0('/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SVs_vs_SNPs/SVs_SNPs_concat_100_',
                         CHR, '.ld')
  LD_table <- fread(PLINK_OUTPUT, header = TRUE)
  
  # For SV-SV pairs
  ## Keep SV-SNP pairs only 
  LD_table <- subset(LD_table, (SNP_A == '.' & SNP_B != '.')|(SNP_A != '.' & SNP_B == '.'))
  
  
  ## Match each SV in PLINK output with a known candidate based on CHROM and POS

  ### Extract the SV's genotyped ID for each pair
  LD_table$SV_ID <- 
    ifelse(test = LD_table$SNP_A == '.',
           yes = LD_table$SNP_B,
           no = LD_table$SNP_A)
  
  
  ## Match by SV ID
  LD_table_matched <- merge(x = LD_table, y = matched_SVs,
                            by.x = 'SV_ID', by.y = 'ID', all.x = TRUE)
  
  
  
  LD_table_matched <- as.data.frame(LD_table_matched)
  print(paste0(unlist(strsplit(PLINK_OUTPUT, split = '.ld')), '_SVTYPE_LEN.ld'))
  
  ## Export
  write.table(file = paste0(unlist(strsplit(PLINK_OUTPUT, split = '.ld')), '_SVTYPE_LEN.ld'),
              LD_table_matched[, c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2", "POS", "END", "CAND_SVTYPE", 'CAND_SVLEN')],
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}



