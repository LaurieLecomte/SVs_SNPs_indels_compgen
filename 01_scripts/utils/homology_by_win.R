# Get homology by windows 
setwd("~/projects/lalec31/RDC_Romaine/RepeatMasker_202303")

library(ggplot2)
library(magrittr)
library(dplyr)

options(scipen=999)

PCT_IDENTITY <- 85

# Mapped blocks
files <- list.files(path = "blocks_different_chromosome_chain_folder/", 
                   pattern = ".withcov.nopercent.out")
#files = sub("NC_0", "NC0", files)
#files =	sub("NC_0", "NC0", files)



files.dd = do.call(rbind.data.frame, strsplit(files, "\\_|\\.withcov.nopercent.out"))
colnames(files.dd) = c("Chr1", "start1", "end1",
                       "Chr2", "start2", "end2")
mapped_blocks <- files.dd

# Format chromosomes names
mapped_blocks$Chr1 <- substr(mapped_blocks$Chr1, 4, nchar(mapped_blocks$Chr1))
mapped_blocks$Chr2 <- substr(mapped_blocks$Chr2, 4, nchar(mapped_blocks$Chr2))

# Convert positions to numeric
mapped_blocks$start1 <- as.numeric(mapped_blocks$start1)
mapped_blocks$end1 <- as.numeric(mapped_blocks$end1)
mapped_blocks$start2 <- as.numeric(mapped_blocks$start2)
#mapped_blocks$end2 <- sapply(strsplit(mapped_blocks$end2, ".withcov.nopercent.out"), "[[", 1)
mapped_blocks$end2 <- as.numeric(mapped_blocks$end2)


#mapped_blocks <- mapped_blocks[-which(mapped_blocks$Chr1 == "grp1"),]

# List of blocks outputted by symap
symap_blocks <- read.delim("blocks_different_chromosome", header = FALSE,
                    col.names = c('grp1', 'grp2', 'block',   
                                  'start1', 'end1', 'start2', 'end2',
                                  'hits', 'gene1', 'gene2', '%gene1', '%gene2'))

list.chr <- unique(symap_blocks$grp1)
win.size <- 1000000

hom_blocks <- data.frame(chr = NA, start = NA, end = NA,
                identity = NA, identity.weighted = NA,
                pct_cov = NA, n_anchor = NA,
                mean_hsp = NA)


# Loop by window
#for(i in 1:length(list.chr)){
##for(i in 7){
#  chr = list.chr[i]
#  
#  row_in_dd = nrow(dd)
#  
#  for(j in 1:ceiling(max(symap_blocks$V5[symap_blocks$V1 == chr])/win.size)){
#    j2 = j + row_in_dd
#    start = (j-1)*win.size +1
#    end = j*win.size
#    
#    dd[j2, "chr"] = chr
#    dd$start[j2] = start
#    dd$end[j2] = end
#    
#    mapped_blocks.sub = which(mapped_blocks$Chr1 == chr &
#                              ((mapped_blocks$start1 > start & mapped_blocks$start1 < end) |
#                              (mapped_blocks$end1 > start & mapped_blocks$end1 < end) |
#                              (mapped_blocks$start1 < start & mapped_blocks$end1 > end)))
#    
#    for(k in mapped_blocks.sub){
#      
#      anchor = read.table(files[k])
#      anchor.win = anchor[anchor$V5 > start & anchor$V6 < end & anchor$V13 > 75,]
#      dd$identity[j2] = mean(anchor.win$V13)
#      dd$identity.weighted[j2] = sum(anchor.win$V13 * anchor.win$V16) / sum(anchor.win$V16)
#      dd$pct_cov[j2] = 
#    }
#  }
#}

#Loop by mapped blocks
for(i in 1:nrow(mapped_blocks)){
  # Extract chr and positions
  chr <- mapped_blocks$Chr1[i]
  # Get anchors from corresponding lastz output file
  anchor <- read.delim(paste0("blocks_different_chromosome_chain_folder/", files[i])
                      )
  
  ## Get positions in Mb
  first <- floor(mapped_blocks$start1[i]/win.size) + 1
  last <- ceiling(mapped_blocks$end1[i]/win.size) + 1
  
  #row_in_hom_blocks <- nrow(hom_blocks)

  
  for(j in first:last){
    start <- (j - 1) * win.size
    end <- j * win.size
    
    k <- nrow(hom_blocks) + 1
    
    hom_blocks[k, "chr"] <- chr
    hom_blocks$start[k] <- start
    hom_blocks$end[k] <- end
    
    # Extract anchors that fall inside block and that have % identity over 75 %
    anchor.win <- anchor[anchor$zstart1 > start & anchor$end1 < end & anchor$idPct > PCT_IDENTITY, ]
    
    if (nrow(anchor.win) == 0) print('no anchors')
    
    hom_blocks$identity[k] <- mean(anchor.win$idPct)
    hom_blocks$identity.weighted[k] <- sum(anchor.win$idPct * anchor.win$coverage.1) / sum(anchor.win$coverage.1)
    hom_blocks$pct_cov[k] <- sum(anchor.win$coverage.1) / win.size * 100
    hom_blocks$n_anchor[k] <- nrow(anchor.win)
    hom_blocks$mean_length[k] <- mean(anchor.win$coverage.1)
  }
}


# Remove blank entries
hom_blocks <- hom_blocks[-1,]
hom_blocks <- hom_blocks[hom_blocks$n_anchor > 0,]
hom_blocks <- hom_blocks[order(hom_blocks$start),]
hom_blocks <- hom_blocks[order(hom_blocks$chr),]

hom_blocks %<>% group_by(chr, start, end) %>%
  summarise(identity = weighted.mean(identity, n_anchor),
            identity.weighted = weighted.mean(identity.weighted, n_anchor),
            pct_cov = sum(pct_cov),
            n_anchor = sum(n_anchor))

# Correct chromosome name
hom_blocks$chr <- paste0(hom_blocks$chr, '.1')
#ggplot(dd[dd$chr == "HG993270",], aes(x = start, y = identity.weighted)) + geom_col() +
#  coord_cartesian(ylim = c(70,100)) + theme(legend.position = "none") +
#  geom_segment(data = block[block$V1 == "HG993270",], aes(y = 100, yend = 100, x = V4, xend = V5, col = V2))

# Export as bed file
write.table(x = hom_blocks[,1:5], paste0('homolog_blocks_identity_', PCT_IDENTITY, '.bed'), row.names = FALSE, col.names = FALSE, 
            sep = "\t", quote = FALSE)

