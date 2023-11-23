# Plot kept/excluded variant density in repeats and syntenic regions
# Required files are produced by the 01_scripts/utils/filt_excl_by_homology.sh script

library(data.table)
library(ggplot2)
library(dplyr)
options(scipen = 999)

MIN_IDY_PCT <- 85 # min % of identity of syntenic regions to consider
DENS_WIN <- 10000  # size of window for computing variant density

# SVs ---------------------------------------------------------------------

# SV intersection with syntenic and repeats -------------------------------

# Intersection with syntenic regions
INTERSECT_HOM <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/win100_filt_excl_SVs_hom_regions.table"

intersect_hom <- read.table(INTERSECT_HOM, 
           col.names = c('CHROM', 'POS', 'END', 'ID', 'type', 'chrom', 'start', 'stop', 'hom_pct', 'wg_hom_pct'))

## add unique suffixe to IDs, because some SVs were genotyped more than once and 
## a given ID can be bin both filt and excl SVs
intersect_hom$ID <- ifelse(intersect_hom$type == 'excluded',
                           yes = paste0(intersect_hom$ID, '_excl'),
                           no = intersect_hom$ID)


# Intersection with repeats (RepeatMasker output)
INTERSECT_RM <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/win100_filt_excl_SVs_RM.table"

intersect_rm <- read.table(INTERSECT_RM, 
           col.names = c('CHROM', 'POS', 'END', 'ID', 'type', 'chrom', 'start', 'stop', 'element'))

## add unique suffixe to IDs, because some SVs were genotyped more than once and 
## a given ID can be bin both filt and excl SVs
intersect_rm$ID <- ifelse(intersect_rm$type == 'excluded',
                           yes = paste0(intersect_rm$ID, '_excl'),
                           no = intersect_rm$ID)


# Combine both datasets by merging on SV CHROM, POS, END, ID and type (filtered or excluded)
intersect_both <- 
merge(x = intersect_hom[ c('CHROM', 'POS', 'END', 'ID', 'type', 'hom_pct', 'wg_hom_pct')], 
      y = intersect_rm[, c('CHROM', 'POS', 'END', 'ID', 'type', 'element')],
      by = c('CHROM', 'POS', 'END', 'ID', 'type'), all = TRUE)

# Remove non-unique rows
intersect_both <- distinct(intersect_both)


# Classify as in repetitive or non repetitive 
intersect_both$repetitive <- ifelse(is.na(intersect_both$element),
                                    yes = 'non-repetitive',
                                    no = 'repetitive')


# Add variable for homology level for SVs overlapping syntenic regions
intersect_both$wg_hom_pct <- round(intersect_both$wg_hom_pct, digits = 0)

# Remove syntenic regions with % identity lower than MIN_IDY_PCT
intersect_both <- subset(intersect_both, wg_hom_pct > MIN_IDY_PCT)

intersect_both$homology <- as.character(
  cut(intersect_both$wg_hom_pct, 
      breaks = c(85, 90, 95, 97.5, 100), 
      labels = c('Low (85 - 90%)', 'Elevated (90 - 95%)', 
                 'High (95 - 97.5%)', 'Very high (97.5 - 100%)'), 
      right = FALSE)  )


intersect_both$homology <- factor(intersect_both$homology, 
                                       levels=c('Low (85 - 90%)', 'Elevated (90 - 95%)', 
                                                'High (95 - 97.5%)', 'Very high (97.5 - 100%)')
)

# Export result table to file : will be used later for plotting an histogram 
write.table(file = "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/SVs_intersect_both.table",
            intersect_both,
            row.names = FALSE,
            sep = "\t")


ggplot(data = intersect_both) + 
  geom_histogram(aes(wg_hom_pct, fill = homology), binwidth = 1) +
  theme_bw() + 
  labs(x = 'Homology (% identity)') +
  guides(fill = 'none') +
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))


# SV density by window ----------------------------------------------------

# Get windows
WIN_CHUNKS <- paste0("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/02_infos/chrs_win", DENS_WIN, '.bed')

win_chunks <- as.data.table(read.table(WIN_CHUNKS,
                                       col.names = c('CHROM', 'START', 'STOP')))

# For each SV, get corresponding window
intersect_both_dt <- as.data.table(intersect_both)

intersect_both_win <- 
  intersect_both_dt[win_chunks, 
                    on = .(CHROM, 
                           POS >= START, END <= STOP), # POS (from intersect_both_dt) must be between START and END from win_chunks
                    .(x.CHROM, x.POS, x.END, x.ID, 
                      x.type, x.repetitive, x.homology, x.element,
                      i.START, i.STOP 
                    )] 

## remove windows where there is nothing
intersect_both_win <- subset(intersect_both_win, !is.na(x.CHROM))



## Get density : count by group  
by_homgroup <-
  intersect_both_win %>% count(x.CHROM, i.START, i.STOP, 
                               x.type, x.homology, x.repetitive, sort = TRUE)

by_homgroup$x.homology_group <- as.character(by_homgroup$x.homology)

by_homgroup$x.homology_group <- ifelse(is.na(by_homgroup$x.homology_group),
                                                          yes = 'Non-syntenic',
                                                          no = by_homgroup$x.homology_group)

by_homgroup$x.homology_group <- factor(by_homgroup$x.homology_group , 
                                                          levels=c('Non-syntenic',
                                                                   'Low (85 - 90%)', 'Elevated (90 - 95%)', 
                                                                   'High (95 - 97.5%)', 'Very high (97.5 - 100%)')
)

by_homgroup$x.type <- ifelse(by_homgroup$x.type == 'excluded',
                                                          yes = 'filtered out',
                                                          no = 'kept')

# Export result table to file, for plotting on a local machine where ggpattern could be installed
write.table(file = paste0("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/SVs_by_homgroup_win", DENS_WIN/1000, "kb.table"),
            by_homgroup,
            row.names = FALSE,
            sep = "\t")

# Plot
ggplot(data = by_homgroup) + 
  facet_grid(~ x.type, scales = 'free_y') +
  geom_boxplot(aes(y = n, fill = x.homology_group, x = x.repetitive), 
               outlier.size = 0.5, linewidth = 0.2, size = 0.5) +
  labs(y = paste0('SV count per ', DENS_WIN/1000, ' kb window'), x = 'SV class', fill = 'Homology (% identity)') +
  theme_bw() #+
  #scale_fill_viridis_d(option = 'H') #+
  #theme(axis.text.x = element_blank(),
   #     axis.ticks.x = element_blank())



ggplot(data = by_homgroup) + 
  facet_grid(~ x.type, scales = 'free_y') +
  geom_boxplot(aes(y = n, fill = x.homology_group, x = x.repetitive), 
               outlier.size = 0.5, linewidth = 0.2, size = 0.5) +
  labs(y = 'SV count per 10 kb', x = 'SV class', fill = 'Homology (% identity)') +
  theme_bw() 



ggplot(data = by_homgroup) + 
  #facet_grid(~ x.type, scales = 'free_y') +
  geom_boxplot(aes(y = n, fill = x.homology, color = x.repetitive, x = x.type), linewidth = 0.2) +
  labs(y = 'SV count per 10 kb', x = 'SV class', fill = 'Homology (% identity)') +
  theme_bw() +
  scale_color_manual(values = c('black', 'grey50')) 
  #theme(axis.text.x = element_blank(),
  #           axis.ticks.x = element_blank())


ggplot(data = intersect_both) +
  geom_histogram(aes(x = wg_hom_pct, fill = homology), binwidth = 1)
table(intersect_both$wg_hom_pct, intersect_both$homology)








# SNPs ---------------------------------------------------------------------

# SNP intersection with syntenic and repeats -------------------------------
INTERSECT_HOM_SNP <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/win100_filt_excl_SNPs_hom_regions.table"

intersect_hom_SNPs <- read.table(INTERSECT_HOM_SNP, 
                            col.names = c('CHROM', 'POS', 'END', 'ID', 'type', 'chrom', 'start', 'stop', 'hom_pct', 'wg_hom_pct'))

intersect_hom_SNPs$ID <- ifelse(intersect_hom_SNPs$type == 'excluded',
                           yes = paste0(intersect_hom_SNPs$ID, '_excl'),
                           no = intersect_hom_SNPs$ID)

INTERSECT_RM_SNP <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/win100_filt_excl_SNPs_RM.table"

intersect_rm_SNPs <- read.table(INTERSECT_RM_SNP, 
                           col.names = c('CHROM', 'POS', 'END', 'ID', 'type', 'chrom', 'start', 'stop', 'element'))
intersect_rm_SNPs$ID <- ifelse(intersect_rm_SNPs$type == 'excluded',
                          yes = paste0(intersect_rm_SNPs$ID, '_excl'),
                          no = intersect_rm_SNPs$ID)

intersect_both_SNPs <- 
  merge(x = intersect_hom_SNPs[ c('CHROM', 'POS', 'END', 'ID', 'type', 'hom_pct', 'wg_hom_pct')], 
        y = intersect_rm_SNPs[, c('CHROM', 'POS', 'END', 'ID', 'type', 'element')],
        by = c('CHROM', 'POS', 'END', 'ID', 'type'), all = TRUE)

intersect_both_SNPs <- distinct(intersect_both_SNPs)


# Classify as in repetitive or non repetitive 
intersect_both_SNPs$repetitive <- ifelse(is.na(intersect_both_SNPs$element),
                                    yes = 'non-repetitive',
                                    no = 'repetitive')


# Add variable for homology level for SVs overlapping syntenic regions
intersect_both_SNPs$wg_hom_pct <- round(intersect_both_SNPs$wg_hom_pct, digits = 0)

# Remove syntenic regions with % identity lower than MIN_IDY_PCT
intersect_both_SNPs <- subset(intersect_both_SNPs, wg_hom_pct > MIN_IDY_PCT)

intersect_both_SNPs$homology <- as.character(
  cut(intersect_both_SNPs$wg_hom_pct, 
      breaks = c(85, 90, 95, 97.5, 100), 
      labels = c('Low (85 - 90%)', 'Elevated (90 - 95%)', 
                 'High (95 - 97.5%)', 'Very high (97.5 - 100%)'), 
      right = FALSE)  )

intersect_both_SNPs$homology <- factor(intersect_both_SNPs$homology, 
                                  levels=c('Low (85 - 90%)', 'Elevated (90 - 95%)', 
                                           'High (95 - 97.5%)', 'Very high (97.5 - 100%)')
)

# Export result table to file : will be used later for plotting an histogram 
write.table(file = "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/SNPs_intersect_both.table",
            intersect_both_SNPs,
            row.names = FALSE,
            sep = "\t")


ggplot(data = intersect_both_SNPs) + 
  geom_histogram(aes(wg_hom_pct, fill = homology), binwidth = 1) +
  theme_bw() + 
  labs(x = 'Homology (% identity)') +
  guides(fill = 'none') +
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))



# SNP density by window ---------------------------------------------------
WIN_CHUNKS <- paste0("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/02_infos/chrs_win", DENS_WIN, '.bed')


win_chunks <- as.data.table(read.table(WIN_CHUNKS,
                                       col.names = c('CHROM', 'START', 'STOP')))


# Get corresponding chunk for each SV
intersect_both_dt_SNPs <- as.data.table(intersect_both_SNPs)

intersect_both_win_SNPs <- 
  intersect_both_dt_SNPs[win_chunks, 
                    on = .(CHROM, 
                           POS >= START, END <= STOP),
                    .(x.CHROM, x.POS, x.END, x.ID, 
                      x.type, x.repetitive, x.homology, x.element,
                      i.START, i.STOP # variables I need in merged output (i = overlap_all_bed) (x = win_chunks)
                    )] 

# remove windows where there is nothing
intersect_both_win_SNPs <- subset(intersect_both_win_SNPs, !is.na(x.CHROM))



## Get density : count by group  
by_homgroup_SNPs <-
  intersect_both_win_SNPs %>% count(x.CHROM, i.START, i.STOP, 
                               x.type, x.homology, x.repetitive, sort = TRUE)

by_homgroup_SNPs$x.homology_group <- as.character(by_homgroup_SNPs$x.homology)

by_homgroup_SNPs$x.homology_group <- ifelse(is.na(by_homgroup_SNPs$x.homology_group),
                                                          yes = 'Non-syntenic',
                                                          no = by_homgroup_SNPs$x.homology_group)

by_homgroup_SNPs$x.homology_group <- factor(by_homgroup_SNPs$x.homology_group , 
                                                          levels=c('Non-syntenic',
                                                                   'Low (85 - 90%)', 'Elevated (90 - 95%)', 
                                                                   'High (95 - 97.5%)', 'Very high (97.5 - 100%)')
)

by_homgroup_SNPs$x.type <- ifelse(by_homgroup_SNPs$x.type == 'excluded',
                                                yes = 'filtered out',
                                                no = 'kept')


# Export result table to file, for plotting on a local machine where ggpattern could be installed
write.table(file = paste0("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/SNPs_by_homgroup_win", DENS_WIN/1000, "kb.table"),
            by_homgroup,
            row.names = FALSE,
            sep = "\t")

# Plot
ggplot(data = by_homgroup_SNPs) + 
  facet_grid(~ x.type, scales = 'free_y') +
  geom_boxplot(aes(y = n, fill = x.homology_group, x = x.repetitive), 
               outlier.size = 0.5, linewidth = 0.2, size = 0.5) +
  labs(y = 'SNP count per 10 kb', x = 'SNP class', fill = 'Homology (% identity)') +
  theme_bw() #+
#theme(axis.text.x = element_blank(),
#     axis.ticks.x = element_blank())

ggplot(data = by_homgroup_SNPs) + 
  #facet_grid(~ x.type, scales = 'free_y') +
  geom_boxplot(aes(y = n, fill = x.homology, color = x.repetitive, x = x.type), linewidth = 0.2) +
  labs(y = 'SNP count per 10 kb', x = 'SNP class', fill = 'Homology (% identity)') +
  theme_bw() +
  scale_color_manual(values = c('black', 'grey50')) 
#theme(axis.text.x = element_blank(),
#           axis.ticks.x = element_blank())























# Poubelle ----------------------------------------------------------------

 # Split in repetitive and non-repetitive region SVs


overlap_all$repetitive <- ifelse(overlap_all$hom_pct == '.',
                                  'non-repetitive',
                                  'repetitive')

overlap_all$hom_pct <- ifelse(overlap_all$hom_pct == '.',
                               0, 
                               overlap_all$hom_pct)

overlap_all$hom_pct <- as.numeric(overlap_all$hom_pct)






overlap_all$homology <- 
  cut(overlap_all$hom_pct, 
      breaks = c(0, 85, 90, 95, 97.5, 100), 
      labels = c('Non syntenic', 'Low (85 - 90 %)', 'Elevated (90 - 95 %)', 
                 'High (95 - 97.5 %)', 'Very high (97.5 - 100 %)'), 
      right = FALSE)


WIN_CHUNKS <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/02_infos/chrs_win1000.bed"


win_chunks <- as.data.table(read.table(WIN_CHUNKS,
                         col.names = c('CHROM', 'START', 'STOP')))


# Get corresponding chunk for each SV
overlap_all_bed <- as.data.table(overlap_all[, c('CHROM', 'POS', 'END', 'ID', 'BP_OVERLAP', 'group', 'repetitive', 'homology')])


library(data.table)



overlap_all_by_win <- 
overlap_all_bed[win_chunks, 
     on = .(CHROM, 
            POS >= START, END <= STOP),
     .(x.CHROM, x.POS, x.END, x.ID, 
       x.group, x.repetitive, x.homology,
       i.START, i.STOP # variables I need in merged output (i = overlap_all_bed) (x = win_chunks)
       )] 



overlap_all_by_win_by_homgroup <-
aggregate(data = overlap_all_by_win, 
          . ~ x.CHROM + i.START + i.STOP + x.group + x.homology + x.repetitive, FUN = length)



ggplot(data = overlap_all_by_win_by_homgroup) + 
  #facet_grid(~ x.group) +
  geom_boxplot(aes(y = x.ID, x = x.group, fill = x.homology)) +
  labs(y = 'SV count per 10 kb', x = 'SV class', fill = 'Homology (% identity)') +
  theme_bw()

OVERLAP_FILT <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/intersect_filt_SVs_hom_regions.table" 

overlap_filt <- read.table(OVERLAP_FILT, 
                           col.names = c('CHROM', 'POS', 'END', 'ID', 'chrom', 'start', 'stop', 'hom_pct', 'wg_hom_pct', 'BP_OVERLAP'))
overlap_filt$group <- 'filtered'

OVERLAP_EXCL <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/intersect_excl_SVs_hom_regions.table" 

overlap_excl <- read.table(OVERLAP_EXCL, 
                           col.names = c('CHROM', 'POS', 'END', 'ID', 'chrom', 'start', 'stop', 'hom_pct', 'wg_hom_pct', 'BP_OVERLAP'))
overlap_excl$group <- 'excluded'

## add unique suffixe to IDs, because some SVs were genotyped more than once and 
## a given ID can be bin both filt and excl SVs
overlap_excl$ID <- paste0(overlap_excl$ID, '_excl')

#excl <- read.table("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/merged_SUPP2_excluded.bed", 
#col.names = c('CHROM', 'POS', 'END', 'ID'))
#filt <- read.table("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/merged_SUPP2_MAF0.05_FMISS0.5.bed", 
#col.names = c('CHROM', 'POS', 'END', 'ID'))


# Bind together filtered and excluded SVs
overlap_all <- rbind(overlap_filt, overlap_excl)




OVERLAP_FILT_RM <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/intersect_filt_SVs_RM.table"

overlap_filt_RM <- read.table(OVERLAP_FILT_RM, 
                              col.names = c('CHROM', 'POS', 'END', 'ID', 'chrom', 'start', 'stop', 'element', 'BP_OVERLAP'))
overlap_filt_RM$group <- 'filtered'

OVERLAP_EXCL_RM <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/hom_regions/intersect_excl_SVs_RM.table"

overlap_excl_RM <- read.table(OVERLAP_FILT_RM, 
                              col.names = c('CHROM', 'POS', 'END', 'ID', 'chrom', 'start', 'stop', 'element', 'BP_OVERLAP'))
overlap_excl_RM$group <- 'excluded'
overlap_excl_RM$ID <- paste0(overlap_excl_RM$ID, '_excl')



overlap_all_RM <- rbind(overlap_filt_RM, overlap_excl_RM)




# Merge hom blocks and repeats intersect SV sets
all_SVs_groups <- 
  merge(x = overlap_all[ c('CHROM', 'POS', 'END', 'ID', 'group', 'hom_pct', 'wg_hom_pct')], 
        y = overlap_all_RM[, c('CHROM', 'POS', 'END', 'ID', 'group', 'element')],
        by = c('CHROM', 'POS', 'END', 'ID', 'group'), all = TRUE)


#all_SVs_groups <- 
#merge(x = overlap_all[ c(1:4,8,9,11)], y = overlap_all_RM[, c(1:4,8,10:12)],
#      by = c('CHROM', 'POS', 'END', 'ID', 'group'), all = TRUE)

# Classify as syntenic or non syntenic
all_SVs_groups$syntenic <- ifelse(is.na(all_SVs_groups$hom_pct),
                                  yes = 'non-syntenic',
                                  no = 'syntenic'
) 

# Classify as in repetitive or non repetitive 
all_SVs_groups$repetitive <- ifelse(all_SVs_groups$element == '.' | is.na(all_SVs_groups$element),
                                    yes = 'non-repetitive',
                                    no = 'repetitive')


# Classify as repeat or TE
all_SVs_groups$rep_type <-
  ifelse(all_SVs_groups$repetitive == 'repetitive' & grepl(all_SVs_groups$element, pattern = ')n'),
         yes = 'rep',
         no = ifelse(test = all_SVs_groups$repetitive == 'repetitive',
                     yes = 'TE',
                     no = 'non-repetitive'))

# Classify as homology levels
all_SVs_groups$homology <- as.character(
  cut(all_SVs_groups$wg_hom_pct, 
      breaks = c(85, 90, 95, 97.5, 100), 
      labels = c('Low (85 - 90 %)', 'Elevated (90 - 95 %)', 
                 'High (95 - 97.5 %)', 'Very high (97.5 - 100 %)'), 
      right = FALSE)  )

all_SVs_groups$homology <- ifelse(is.na(all_SVs_groups$homology),
                                  yes = 'Non-syntenic',
                                  no = all_SVs_groups$homology)

#sapply(all_SVs_groups$homology, FUN = function(x) ifelse(is.na(x), yes = 'NS', no = x))




WIN_CHUNKS <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/03_SR_LR/SVs_SNPs_indels_compgen/02_infos/chrs_win100000.bed"


win_chunks <- as.data.table(read.table(WIN_CHUNKS,
                                       col.names = c('CHROM', 'START', 'STOP')))


# Get corresponding chunk for each SV
all_SVs_groups_dt <- as.data.table(all_SVs_groups)

all_SVs_groups_win <- 
  all_SVs_groups_dt[win_chunks, 
                    on = .(CHROM, 
                           POS >= START, END <= STOP),
                    .(x.CHROM, x.POS, x.END, x.ID, 
                      x.group, x.repetitive, x.homology, x.syntenic, x.element,
                      i.START, i.STOP # variables I need in merged output (i = overlap_all_bed) (x = win_chunks)
                    )] 

# remove windows where there is nothing
all_SVs_groups_win <- subset(all_SVs_groups_win, !is.na(x.CHROM))

# Count by group  
all_SVs_groups_win_by_homgroup <-
  all_SVs_groups_win %>% count(x.CHROM, i.START, i.STOP, 
                               x.group, x.syntenic, x.homology, x.repetitive, sort = TRUE)




ggplot(data = all_SVs_groups_win_by_homgroup) + 
  #facet_grid(~ x.group) +
  geom_boxplot(aes(y = n, x = x.group, fill = x.homology, color = x.repetitive)) +
  labs(y = 'SV count per 10 kb', x = 'SV class', fill = 'Homology (% identity)') +
  theme_bw()


ggplot(data = all_SVs_groups_win) + 
  #facet_grid(~ x.group) +
  geom_bar(aes(x = x.group, fill = x.homology, color = x.repetitive)) +
  labs(y = 'SV count per kb', x = 'SV class', fill = 'Homology (% identity)') +
  theme_bw()


