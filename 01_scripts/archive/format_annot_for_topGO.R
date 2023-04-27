genome_annotation_table_simplified_1.5 <- read.delim("03_genome/annotation/genome_annotation_table_simplified_1.5.tsv")

#genome_annotation_table_simplified_1.5$taxid <- 'salmo_salar'
#genome_annotation_table_simplified_1.5$gene_id <- genome_annotation_table_simplified_1.5$TranscriptName
#genome_annotation_table_simplified_1.5$gene_symbol <- genome_annotation_table_simplified_1.5$GeneName
#genome_annotation_table_simplified_1.5$evidence <- 'IEA'


GO_list <- vector(mode = 'list', length = nrow(genome_annotation_table_simplified_1.5))


for (i in 1:nrow(genome_annotation_table_simplified_1.5)){
  vec <- unlist(strsplit(x = sub(genome_annotation_table_simplified_1.5$GeneGo[i], pattern = ' ', replacement = ''), 
                         split = ';'))
  names(vec) <- rep(genome_annotation_table_simplified_1.5$TranscriptName[i], time = length(vec))
  GO_list[[i]] <- as.list(vec)
}

full_list <- do.call(c, GO_list)
list_df <- stack(full_list)
colnames(list_df) <- c('GOID', 'gene_id')


merged_GOs <- 
merge(genome_annotation_table_simplified_1.5, list_df, 
      by.x = 'TranscriptName', by.y = 'gene_id')

merged_GOs$taxid <- 'salmo_salar'
merged_GOs$evidence <- 'IEA'

merged_GOs <- merged_GOs[, c('taxid', 'TranscriptName', 'TranscriptName', 'GOID', 'evidence')]
colnames(merged_GOs) <- c('taxid', 'gene_id', 'gene_symbol', 'GOID', 'evidence')

write.table(merged_GOs, file = '03_genome/annotation/custom_annot_viseago.txt',
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

not_in_db <- read.table("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/03_genome/annotation/not_in_db.txt", quote="\"", comment.char="")

not_in_db <- as.vector(not_in_db$V1)

merged_GOs_corrected <- merged_GOs[ !merged_GOs$GOID %in% not_in_db, ]
write.table(merged_GOs_corrected, file = '03_genome/annotation/custom_annot_viseago_corrected.txt',
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
