argv <- commandArgs(T)
VCF_INFO <- argv[1]
#input file is 4 columns extracted from a vcf output by vg call with CHR start, reference sequence in the sv bubble, alt sequence#
#ref seq of 1 and alt seq >1 is an insertion relatively to ref genome
#ref seq >1 and alt seq=1 is a deletion
#it is yet still unclear how we can get inversions
# I believe we also have SV with simply alternative path which are neither inversion nor deletions, nor insertions


#read the vcf information file
vcf_info <- read.table(VCF_INFO, header = FALSE, stringsAsFactor = FALSE)
colnames(vcf_info) <- c("CHR","start","REF","ALT")
#vcf_info[1:2,]

#count the number of base in ref and alt
vcf_info$ref_length <- nchar(vcf_info$REF)
vcf_info$alt_length <- nchar(vcf_info$ALT)
#calculate length
vcf_info$LEN <- vcf_info$alt_length - vcf_info$ref_length
#categorize as INS or DEL relatively to reference
vcf_info$TYPE <- "INS"
vcf_info$TYPE[vcf_info$LEN<0] <- "DEL"

vcf_info$start_bed <- vcf_info$start-1 #because bedfiles are 0-index
vcf_info$stop_bed <- vcf_info$start_bed + vcf_info$ref_length
vcf_info$chr_pos <- paste(vcf_info$CHR, vcf_info$start, sep="_")

head(vcf_info)

vcf_info$max_LEN <- pmax(vcf_info$alt_length, vcf_info$ref_length) - 1

nSV = dim(vcf_info)[1]
fasta_vec <- vector(length = 4*nSV)

#for (i in 1 : nSV )
#{
#if (vcf_info$ref_length[i]>1)
#{
#fasta_vec[4*i-3]<-paste0(">",vcf_info$chr_pos[i],"-REF")
#fasta_vec[4*i-2]<-vcf_info$REF[i]
#}
#if (vcf_info$alt_length[i]>1)
#{
#fasta_vec[4*i-1]<-paste0(">",vcf_info$chr_pos[i],"-ALT")
#fasta_vec[4*i]<-vcf_info$ALT[i]
#}
#}

# Loop over SVs
for (i in 1:nSV) {
  if (vcf_info$ref_length[i] > 1) { # means that SV is likely a DEL, so we annotate the REF 
    fasta_vec[4*i-3] <- paste0(">", vcf_info$chr_pos[i], '_', vcf_info$SVTYPE[i], "-REF")
    fasta_vec[4*i-2] <- vcf_info$REF[i]
  }
  if (vcf_info$alt_length[i] > 1) { # means that SV is an INS, so we annotate the ALT sequence 
    fasta_vec[4*i-1] <- paste0(">", vcf_info$chr_pos[i], '_', vcf_info$SVTYPE[i], "-ALT")
    fasta_vec[4*i] <- vcf_info$ALT[i]
  }
}

fasta_vec2 <- fasta_vec[-which(fasta_vec == FALSE)]


write.table(vcf_info[,c("CHR", "start",  "start_bed", "stop_bed")], paste0(VCF_INFO, ".bed"), sep = "\t", row.names = F, quote = F, col.names = F)

#write.table(vcf_info[,c(1,10,11,12)],paste0(VCF_INFO,".bed"),sep = "\t", row.names = F, quote = F, col.names = F)

write.table(vcf_info, paste0(VCF_INFO, ".info"),sep = "\t", row.names = F, quote = F, col.names = F)

write.table(cbind(fasta_vec2), paste0(VCF_INFO, ".fasta"), sep = "\t", row.names = F, quote = F, col.names = F)