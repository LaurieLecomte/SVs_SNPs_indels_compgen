# Format SV VCF into a fasta for running RepeatMasker
# Based on previous work by Claire Merot (https://github.com/clairemerot/genotyping_SV/blob/main/01_scripts/Rscripts/extract_SV_fasta.r)

# 1. Access files from command line and import ----------------------------
argv <- commandArgs(T)
SV_SEQS <- argv[1]

SV_infos <- read.table(SV_SEQS, header = FALSE, 
                       col.names = c('CHROM', 'POS', 'ID', 'REF', 'ALT'))


# 2. Infer SV type from ALT and REF alleles length ------------------------
# Count the number of base in ref and alt
SV_infos$ref_length <- nchar(SV_infos$REF)
SV_infos$alt_length <- nchar(SV_infos$ALT)

# Calculate length
SV_infos$LEN <- SV_infos$alt_length-SV_infos$ref_length

# Categorize as INS or DEL relatively to reference
SV_infos$TYPE <- "INS"
SV_infos$TYPE[SV_infos$LEN < 0] <- "DEL"

# Convert POS and END to 0 based coordinates
SV_infos$start_bed <- SV_infos$POS - 1 ## because bedfiles are 0-index
SV_infos$stop_bed <- SV_infos$start_bed + SV_infos$ref_length
#SV_infos$chr_pos_id <- paste(SV_infos$CHR, SV_infos$POS, SV_infos$ID, sep = "_") # we can't use ID because it contains > and RepeatMasker interprets this as multiple sequences
SV_infos$chr_pos <- paste(SV_infos$CHR, SV_infos$POS, sep = '_')
SV_infos$max_len <- pmax(SV_infos$alt_length, SV_infos$ref_length) - 1

head(SV_infos)
tail(SV_infos)



# 3. Extract ALT or REF seq based on inferred SV type ---------------------
# Initialize container for storing results
nSV <- dim(SV_infos)[1]
fasta_vec <- vector(length = 4 * nSV)

# Loop over SVs
for (i in 1 : nSV ) {
  if (SV_infos$ref_length[i] > 1) { # likely a DEL, so we use REF seq
    fasta_vec[4*i - 3] <- paste0(">", SV_infos$chr_pos[i], '_', SV_infos$ref_length[i], "-REF")
    fasta_vec[4*i - 2] <- SV_infos$REF[i]
    } 
  if (SV_infos$alt_length[i] > 1) { # likely an INS, so we use ALT seq
    fasta_vec[4*i - 1] <- paste0(">", SV_infos$chr_pos[i], '_', SV_infos$alt_length[i], "-ALT")
    fasta_vec[4*i] <- SV_infos$ALT[i]
    }
  }
# we use if instead of if else because it allows to keep both ALT and REF in cases where REF nor ALT has length 1
fasta_vec2 <- fasta_vec[-which(fasta_vec == FALSE)]


# 4. Export ---------------------------------------------------------------
# Export to bed
write.table(SV_infos[, c('CHROM', 'start_bed', 'stop_bed', 'chr_pos')], 
            file = paste0(strsplit(SV_SEQS, split = '.table')[[1]], ".bed"), 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

#write.table(SV_infos[, c(1,9,10,11)], paste0(strsplit(SV_SEQS, split = '.table')[[1]], ".bed"), 
#            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# Export to bed (why ?)
write.table(SV_infos, 
            file = paste0(strsplit(SV_SEQS, split = '.table')[[1]], ".info"), 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# Save fasta
#write.table(cbind(fasta_vec2), 
#            file = paste0(strsplit(SV_SEQS, split = '.table')[[1]], ".fasta"), 
#            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

writeLines(fasta_vec2, paste0(strsplit(SV_SEQS, split = '.table')[[1]], ".fasta"), sep = "\n")