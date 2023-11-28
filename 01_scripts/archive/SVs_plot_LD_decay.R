# Plot LDdecay

library(ggplot2)

# 1. Import and format ----------------------------------------------------
SV_LD_CHR1 <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/LD/SVs/merged_SUPP2_MAF0.05_FMISS0.5.recoded.OV354430.1.RO.stat.gz"

SVs_LD_chr1 <- read.table(SV_LD_CHR1, header = FALSE, 
                     col.names = c('Dist', 'Mean_r^2', 'Mean_D', 'Sum_r^2', 'Sum_D', 'NumberPairs')
                    # colClasses = c('character')
)

SVs_LD_chr1$CHR <- 'OV354430.1'

#SVs_LD_chr1 <- SVs_LD_chr1[order(SVs_LD_chr1$subjects, decreasing = FALSE), ]   

ggplot(data = SVs_LD_chr1) +
  geom_point(aes(x = Dist, y = Mean_r.2))
