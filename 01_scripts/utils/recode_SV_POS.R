recode_SVs <- function(input_vcf, output_vcf, max_width) {
  
  # Disable scientific notation
  options(scipen=999)
  
  # Opening a connection to read the vcf file
  input_con <- file(input_vcf, open = "rt")
  on.exit(close(input_con), add = TRUE)
  
  # Opening another connection to the output file
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
  while(grepl("^#", cur_line <- scan(input_con, what = character(), sep = "\n", n = 1, quiet = TRUE))) {
    cat(paste0(cur_line, "\n"), file = output_con)
  }
  
  # Reading the vcf from file wnad ignoring header lines
  vcf <- read.table(input_vcf, comment.char = "#", stringsAsFactors = FALSE)
  
  # Extracting some useful information for each variant
  chrs   <- vcf[[1]]
  starts <- vcf[[2]]
  # END field can be at the beginning or in the middle of INFO fields, so we need to extract accordingly
  ends <- ifelse(test = grepl("^END=", x = vcf[[8]]),
                 yes = as.integer(sub("^END=(-?[0-9]+).*", "\\1", vcf[[8]])),
                 no = as.integer(sub(".*;END=(-?[0-9]+).*", "\\1", vcf[[8]]))
  )

  widths <- ends - starts ######### addition by me
  


  # Recode position if POS and END are not equal and if difference is greater than min_width
  new_pos <- ifelse(widths >= max_width,
                    yes = ceiling(starts + (widths/2)),
                    no = starts
                    )
  
  # Assigning the results to the right columns of the vcf file
  vcf[[2]] <- new_pos

  
  
  # Writing the data.frame to the output file
  write.table(vcf, file = output_con, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # Writing
  return(invisible(NULL))
}

