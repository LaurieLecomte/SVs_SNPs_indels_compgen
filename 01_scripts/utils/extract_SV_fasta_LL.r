library(data.table)
library(tidyr)
library(dplyr)
argv <- commandArgs(T)
BED <- argv[1] # path to VCF contents (without header)
OUTPUT <- argv[2] # output

#input file is 4 columns extracted from a vcf output by vg call with CHR start, reference sequence in the sv bubble, alt sequence#
#ref seq of 1 and alt seq >1 is an insertion relatively to ref genome
#ref seq >1 and alt seq=1 is a deletion
#it is yet still unclear how we can get inversions
# I believe we also have SV with simply alternative path which are neither inversion nor deletions, nor insertions

BED <- 'TE/SVs_CHR_POS_ID_REF_ALT.table'

# read the vcf information file
bed <- fread(BED, sep = "\t")
colnames(bed) <- c("CHR", "POS", "ID", "REF", "ALT")


# count the number of base in ref and alt
bed$ref_length <- nchar(bed$REF)
bed$alt_length <- nchar(bed$ALT)

# calculate length
bed$LEN <- bed$alt_length - bed$ref_length

# categorize as INS or DEL relatively to reference
#vcf_info$TYPE <- "INS"
#vcf_info$TYPE[vcf_info$LEN < 0] <- "DEL"

bed$start_bed <- bed$POS - 1 #because bedfiles are 0-index
bed$stop_bed <- bed$start_bed + bed$ref_length
bed$chr_pos <- paste(bed$CHR, bed$POS, sep = "_")
bed$max_len <- pmax(bed$alt_length, bed$ref_length) -1

# Create empty vector for storing sequences
nSV <- dim(bed)[1]
fasta_vec <- vector(length = 4*nSV)

# Loop over SVs
for (i in 1:nSV) {
  if (bed$ref_length[i] > 1) { # means that SV is likely a DEL, so we annotate the REF 
    fasta_vec[4*i-3] <- paste0(">", bed$chr_pos[i],  '_', bed$ID[i], '_', bed$ref_length[i], "-REF")
    fasta_vec[4*i-2] <- bed$REF[i]
  }
  if (bed$alt_length[i] > 1) { # means that SV is an INS, so we annotate the ALT sequence 
    fasta_vec[4*i-1] <- paste0(">", bed$chr_pos[i],  '_', bed$ID[i], '_', bed$alt_length[i], "-ALT")
    fasta_vec[4*i] <- bed$ALT[i]
  }
}

fasta_vec2 <- fasta_vec[-which(fasta_vec == FALSE)]


#write.table(vcf_info[, c(1,9,10,11)], paste0(VCF_INFO, ".bed"), sep = "\t", row.names = F, quote = F, col.names = F)
write.table(bed, paste0(BED, ".table"), sep = "\t", row.names = F, quote = F, col.names = F)

#writeLines(fasta_vec2, paste0(strsplit(BED, split = ".bed" )[[1]][1], '.fasta'), sep = "\n")
writeLines(fasta_vec2, paste0(OUTPUT), sep = "\n")

#write.table(cbind(fasta_vec2), paste0(strsplit(BED, split = "/|.bed" )[[1]] %>% .[length(.)], ".fasta"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)