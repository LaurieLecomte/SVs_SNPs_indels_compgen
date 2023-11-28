# Create a list of exons, CDS, mRNA, etc, only for filtered genes 


library(data.table)


# 1. Import ---------------------------------------------------------------

# Complete, unfiltered table (produced with gff2bed from bedops and genome.gff3 from GAWN pipeline)
full_table <- as.data.frame(fread("03_genome/annotation/genome.bed",
                            col.names = c('CHROM', 'POS', 'END', 'V4', 'V5', 'STRAND',
                                          'SOURCE', 'TYPE', 'V9', 'ATT')))

## Remove unplaced contigs to reduce data weight
full_table <- full_table[! grepl(full_table$CHROM, pattern = 'CAK'), ]

## Extract parent gene ID for each element
full_table$GENE_ID <- sub(".*Name=([-._A-Za-z0-9]+);.*", "\\1", full_table$ATT)

# Filtered table, outputted by GAWN
simpl_table <- as.data.frame(fread("03_genome/annotation/genome_annotation_table_simplified_1.5.table",
                                   col.names = c('CHROM', 'POS', 'END', 'STRAND', 'ID', 'PATH',
                                                 'GeneAccession', 'GeneName', 'GeneAltNames',
                                                 'GenePfam', 'GeneGo', 'CellularComponent', 
                                                 'MolecularFunction', 'BiologicalProcess')))

## Remove unplaced contigs to reduce data weight
simpl_table <- simpl_table[! grepl(simpl_table$CHROM, pattern = 'CAK'), ]



# 2. Merge both tables ----------------------------------------------------
merged_table <- 
merge(x = simpl_table, y = full_table,
      by.x = c('CHROM', 'ID'), by.y = c('CHROM', 'GENE_ID'), all = FALSE)

merged_table_simpl <- 
merged_table[, c('CHROM', 'POS.y', 'END.y', 'ID', 'TYPE')]

# Export
write.table(merged_table_simpl, file = '03_genome/annotation/filtered_elements.table',
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
