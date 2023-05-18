ANNOT <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/03_genome/annotation/genome_annotation_table_simplified_1.5.table"

GO <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/12_go/SVs/SVs_RO_PU_RDA_3sd_outliers_overlap10000bp_GO.fdr0.1_depth1.csv"

CHR_SIZES <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/03_genome/genome.fasta.bed"

OV_2_SSA <- "02_infos/OV_to_ssa.txt"
ov_2_ssa <- read.table(OV_2_SSA, col.names = c('CHROM_OV', 'CHROM_SSA'))

chr_sizes <- read.table(CHR_SIZES, col.names = c('CHROM', 'START', 'STOP'))
chr_sizes <- chr_sizes[!grepl(chr_sizes$CHROM, pattern = 'CAK'), c('CHROM', 'STOP')]
chr_sizes <- merge(x = chr_sizes, y = ov_2_ssa, by.x = 'CHROM', by.y = 'CHROM_OV')

chr_end <- data.frame(TranscriptName = chr_sizes$CHROM_SSA, CHROM_SSA = chr_sizes$CHROM_SSA, 
                      FromPosition = chr_sizes$STOP, ToPosition = chr_sizes$STOP)
chr_start <- data.frame(TranscriptName = chr_sizes$CHROM_SSA, CHROM_SSA = chr_sizes$CHROM_SSA, 
                        FromPosition = 0, ToPosition = chr_sizes$STOP)

chr <- rbind(chr_start, chr_end)




annot <- read.delim(ANNOT, header = FALSE,
                    col.names = c('ScaffoldName', 'FromPosition', 'ToPosition', 'Sense', 'TranscriptName', 'TranscriptPath',  'GeneAccession', 'GeneName',
                                  'GeneAltNames', 'GenePfam', 'GeneGo', 'CellularComponent', 'MolecularFunction', 'BiologicalProcess'))
annot <- merge(x = annot, y = ov_2_ssa, 
               by.x = 'ScaffoldName', by.y = 'CHROM_OV', all = TRUE)

go_terms <- read.delim(GO, header = TRUE)



#transcripts <- data.frame(transcript_ID = unlist(strsplit(x = go_terms$study_items[1], split = ', ')))

#transcripts_infos <- merge(x = transcripts, y = annot[,  c('ScaffoldName', 'FromPosition', 'ToPosition', 'TranscriptName')],
#      by.x = 'transcript_ID', by.y = 'TranscriptName', all.y = FALSE)
#table(transcripts_infos$ScaffoldName)


for (i in 1:nrow(go_terms)){
  # Split transcripts list into a vector to a data frame (for easier merging)
  transcripts <- data.frame(TranscriptName = unlist(strsplit(x = go_terms$study_items[i], split = ', ')))
  # Get CHROM and POS info for each transcript
  transcripts_infos <- merge(x = transcripts, y = annot[,  c('CHROM_SSA', 'FromPosition', 'ToPosition', 'TranscriptName')],
                             by = 'TranscriptName', all.y = FALSE)
  
  # 3. Add a dummy SNP to set chromosome length
  dummy <- transcripts_infos[nrow(transcripts_infos), ]
  transcripts_infos <- rbind(transcripts_infos, chr)
  # 4. Sort by ScaffoldName
  transcripts_infos <- transcripts_infos[order(transcripts_infos$CHROM_SSA), ]
  #transcripts_infos <- merge(x = transcripts_infos, y = ov_2_ssa,
  #                           by.x = 'ScaffoldName', by.y = 'CHROM_OV',
  #                           all.x = TRUE)
  
  CMplot(transcripts_infos, type = "p", plot.type = "d", 
         chr.pos.max = TRUE, 
         bin.size = 100, chr.den.col = c("#440154FF", "#21908CFF", "#FDE725FF"), 
         #memo = go_terms$name[i], 
         #memo = paste0('SVs_', key, '_', gsub(x = df$name[p], ' ', '_')), dpi = 300, file = 'jpg',
         file.output = FALSE,
         #main = paste('SVs -', key, ':', eval(df$name[p])), main.cex = 1.2, 
         main = go_terms$name[i],
         width = 9, height = 6
  )
}
