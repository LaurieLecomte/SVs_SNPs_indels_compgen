library(data.table)
library(tidyr)
library(dplyr)
argv <- commandArgs(T)
INPUT <- argv[1] # "$FILT_POP_NAME"_NS30_2all.noheader  # VCF contents
SITES_FILE <- argv[2] # "$FILT_POP_NAME"_NS30_2all.variants.ref contains sites with flanking sequences
HEADER <- argv[3] # "$FILT_POP_NAME"_NS30_2all.variants.header contains VCF header


# 0. Extract fields names from header file
## Import header lines from header file
header_lines <- readLines(HEADER)
## Extract the line that contains field names in the VCF : the one that starts with #CHROM
VCF_fields_names <- header_lines[grep(header_lines, pattern = "^#CHROM", perl = TRUE)]
## Format this line as a character vector
input_columns <- unlist(strsplit(VCF_fields_names, split = "\t"))
input_columns[1] <- 'CHROM'
input_columns


# 1. Import VCF contents (without header)
#vcf <- fread(INPUT)
vcf <- fread(INPUT, sep = "\t") # the fields we need are separated from each other with tabs
#colnames(vcf)[1:2] <- c("CHR", "POS")
#print(paste('VCF has', ncol(vcf), 'columns'))

colnames(vcf) <- input_columns
#ncol <- dim(vcf)[2]

## remove REF and ALT column : these will be replaces later on by dummy allele
vcf <- select(vcf, -c('REF', 'ALT'))


# 2. Import sites file with flanking regions
sites <- fread(SITES_FILE)
colnames(sites) <- c("CHR_POS", "REF")
head(sites)

# 3. Split sites' CHR_POS column into CHR and POS
## Here, the separator is _. Beware if CHR names contains this character.
sites2 <- separate(sites,  CHR_POS, into = c("CHROM", "POS"), sep = '_') # rename CHR to CHROM to make merging step easier
head(sites2) # sites2 will have CHROM, POS, REF

# 3. Create column for storing dummy ALT as factors 
# make a dummy alt allele diff from REF
## To assign a dummy allele, we'll store the REF allele as a factor, which will be encoded as an integer 1,2,...
sites2$ALT <- as.numeric(as.factor(sites$REF))
## Assign ATCG to factor levels in sites2
sites2$ALT[sites2$ALT=="1"] <- "C"
sites2$ALT[sites2$ALT=="2"] <- "G"
sites2$ALT[sites2$ALT=="3"] <- "T"
sites2$ALT[sites2$ALT=="4"] <- "A"

## Convert POS to integers
sites2$POS <- as.integer(sites2$POS)

head(sites2)

# 4. Combine VCF contents with dummy data
vcf_sites <- left_join(vcf, sites2) # will merge based on CHROM and POS

# 5. Assign N and A to missing values
vcf_sites$REF[which(is.na(vcf_sites$REF))] <- "N"
vcf_sites$ALT[which(is.na(vcf_sites$ALT))] <- "A"

head(vcf_sites)
colnames(vcf_sites)


# 6. Format output
## Main fields column names
fields <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
## Genotype columns : columns that are NOT main fields
samples_GT <- setdiff(colnames(vcf_sites), fields)

output_columns <- c(fields, samples_GT)
output_columns

output <- vcf_sites[, ..output_columns]

head(output)

write.table(output,  paste0(INPUT,  ".withdummyREFALT"),  sep = "\t",  col.names = FALSE,  row.names = FALSE,  quote = FALSE)