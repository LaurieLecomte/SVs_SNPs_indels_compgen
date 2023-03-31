# Filter GO enrichment results on FDR and level

# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
GO_TABLE <- argv[1] 
MAX_FDR <- as.numeric(argv[2])
MIN_LEVEL <- as.numeric(argv[3])
#MAX_LEVEL <- as.numeric(argv[4])


# Import 
GO_table <- read.delim(GO_TABLE, header = TRUE)

#GO_table$ratio_obs <- sapply(X = GO_table$ratio_in_study,
#                             FUN = function(x){
#                               eval(parse(text = x))
#                             })

# 2. Filter on FDR and depth (level) --------------------------------------
GO_table_filt <- subset(GO_table, 
                        p_fdr_bh <= MAX_FDR & depth >= MIN_LEVEL)


# 3. Export to table ------------------------------------------------------
write.table(GO_table_filt,
            file = paste0(unlist(strsplit(GO_TABLE, split = '.csv'))[1], 
                          '.fdr', MAX_FDR, '_depth', MIN_LEVEL, '.csv'),
            col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = "\t"
            )

