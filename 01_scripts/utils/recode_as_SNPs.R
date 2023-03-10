recode_as_SNPs <- function(input_vcf, genome, output_vcf) {
  
  # Check reference genome files
  if (is.null(genome) || !file.exists(genome)) {
    stop("Reference genome fasta file missing")
  }
  
  # Check reference genome index
  if (!file.exists(paste0(genome, ".fai"))) {
    stop("Reference genome .fai index missing")
  }
  
  # Open connection for reading input VCF
  input_con <- file(input_vcf, open = "rt")
  on.exit(close(input_con), add = TRUE)
  
  # Open connection for writing to output VCF
  output_con <- file(output_vcf, open = "wt")
  on.exit(close(output_con), add = TRUE)
  
  # Lines are then used one by one to output the header
  #while(grepl("^#",
  #            cur_line <- scan(input_con, what = character(), sep = "\n", n = 1, quiet = TRUE))
  #) {
  #  if (grepl('^#CHROM', cur_line)) {
  #    cat(paste0(cur_line, "\n"), file = output_con)
  #  }
  #}
  
  # Write header lines to output
  while (grepl("^#",
               cur_line <-
               scan(
                 input_con,
                 what = character(),
                 sep = "\n",
                 n = 1,
                 quiet = TRUE
               ))) {
    cat(paste0(cur_line, "\n"), file = output_con)
  }
  
  # Read input VCF, skipping header lines
  vcf <- read.table(input_vcf, comment.char = "#", stringsAsFactors = FALSE)
  
  # Recode variants as SNPs
  
  ## Get position range on reference
  ref_range <- GenomicRanges::GRanges(seqnames = vcf$V1,
                                      ranges = IRanges::IRanges(start = vcf$V2, end = vcf$V2))
  
  ## Extract nt at given pos from reference sequence
  ref <- Rsamtools::scanFa(genome, ref_range)
  vcf$V4 <- unname(as.character(ref))
  
  ## Code dummy ALT allele for each REF
  #vcf$V5 <- unname(sapply(
  #  X = vcf$V4,
  #  FUN = function(x) {
  #    switch(
  #      x,
  #      'A' = sample(x = c('T', 'G', 'C'), size = 1),
  #      'T' = sample(x = c('A', 'G', 'C'), size = 1),
  #      'G' = sample(x = c('T', 'A', 'C'), size = 1),
  #      'C' = sample(x = c('T', 'G', 'A'), size = 1)
  #    )
  #  }
  #)
  #)
  vcf$V5 <- unname(sapply(
    X = vcf$V4,
    FUN = function(x) {
      switch(
        x,
        'A' = 'C',
        'T' = 'A',
        'G' = 'T',
        'C' = 'G'
      )
    }
  )
  )
  
  
  # Write to output
  write.table(
    vcf,
    file = output_con,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )

  return(invisible(NULL))
  
}