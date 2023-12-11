# SVs_SNPs_indels_compgen

Scripts for performing  various population genomics analyses for the paper Investigating structural variant, indel and single nucleotide polymorphism differentiation between locally adapted Atlantic salmon populations using whole genome sequencing and a hybrid genomic polymorphism detection approach by Lecomte et al. (2023) ([See on bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.12.557169v1).

All analyses are performed seperatly for SVs, SNPs and indels. 

This pipeline uses sets of genotyped variants outputted by 3 pipelines :

1. SVs : [genotype_SVs_SRLR](https://github.com/LaurieLecomte/genotype_SVs_SRLR)
2. SNPs : [SNPs_indels_SR](https://github.com/LaurieLecomte/SNPs_indels_SR)
3. Indels : [SNPs_indels_SR](https://github.com/LaurieLecomte/SNPs_indels_SR)


## Pipeline overview

1. Format VCF files for subsequent steps : scripts `0X.1_format_VCF.sh` (not required for SNP VCF)
2. Convert VCFs to beagle : scripts `0X.2_X_vcf_to_begale.sh`
3. Perform **PCA** by population : scripts `0X.3_X_pca.sh`

4. Compute population and per-site _F<sub>ST</sub>_ : scripts `0X.4_X_saf_maf_bypop.sh` and `0X.5_X_fst.sh`

5. Get _F<sub>ST</sub>_ outliers : scripts `0X.6_X_fst_outliers.sh`

6. Perform RDA and get RDA candidate variants : scripts `0X.7_X_rda.sh`

7. Perform GO enrichment analysis on RDA candidates : scripts `0X.8_X_rda_go.sh`

8. Perform Fisher's tests : scripts `0X.9_X_fisher_tests.sh` 

10. Get intersection of _F<sub>ST</sub>_ and Fisher outliers, and perform GO enrichment analysis on these : script `0X.10_X_compare_outliers.sh`

 
* Variant density : script `X_density.sh` and script `combined_per_win_density.R` for plotting 
* Base pairs covered by variants by window : script `X_per_win_bp_prop.sh` and script `combined_per_win_bp.R` for plotting
* _F<sub>ST</sub>_ plot : script `combined_per_win_Fst_plot.R`
 

## Prerequisites

### Files

#### For almost everything
* A reference genome (`.fasta`) and its index (`.fai`) in `03_genome`

* A VCF of genotyped and filtered SVs from the [genotype_SVs_SRLR pipeline](https://github.com/LaurieLecomte/genotype_SVs_SRLR)
* A VCF of genotyped and filtered SNPs from the [SNPs_indels_SR pipeline](https://github.com/LaurieLecomte/SNPs_indels_SR)
* A VCF of genotyped and filtered indels from the [SNPs_indels_SR pipeline](https://github.com/LaurieLecomte/SNPs_indels_SR)

* A chromosomes list (or contigs, or sites) in `02_infos`. This list is not currently used in the pipeline, but may be required for future improvements. It can be produced from the indexed genome file (`"$GENOME".fai`) : `less "$GENOME".fai | cut -f1 > 02_infos/chr.txt`. 

* A tab-seperated list of sample IDs with their source population (for performing coloring PCA plots and _F<sub>ST</sub>_ calculation)

#### For performing GO enrichment analysis
* A GO database (`go-basic.obo`) in `12_go/go_db`
* A table of GO annotations for known genes (`all_go_annotations.csv`) in `12_go/go_db`, outputted by the [GAWN pipeline](https://github.com/enormandeau/gawn)


## Aditionnal analyses

### Exploration of genomic features or regions prone to low-quality SVs

#### Procedure overview
1. Identify **repeated regions** : `01_scripts/utils/RepeatMasker_manitou.sh` (this script is provided for reference and is not intended to be used in its current state, as it comes from another local private repo)
2. Identify **syntenic regions** due to whole genome duplication by performing whole-genome self alignment : `01_scripts/utils/nucmer_symap_addedN.sh` (this script also comes from another local and private repo). This script first "patches" a problematic region preventing successful alignment.
3. Get the % identity (homology level) of syntenic regions by mapping them to each other : `01_scripts/utils/lastz_align_blocks.sh` and `01_scripts/utils/homology_by_win.R` (these scripts also come from another local and private repo and are based on previous work by [Xavier Dallaire et al. (2023)](https://www.biorxiv.org/content/10.1101/2023.07.27.550877v1))
4. Get the overlap between good and bad SVs with syntenic regions and/or repeats, compute variant density by filtered/excluded status and by window, then plot : `01_scripts/utils/filt_excl_by_homology.sh` and `01_scripts/utils/filt_excl_by_homology.R`

### Compare missing data distribution variant groups

See script `01_scripts/utils/F_MISS_by_type_vartype.R` 


### Get LD decay between SVs and SNPs
1. Calculate _r^2_ between SVs and SNPs and between SNPs : `01_scripts/utils/LD_SVs_SNPs_plink.sh`
2. Add SV info to the table outputted by `plink` :  `01_scripts/utils/LD_add_SV_info.R`
3. Format and reorder output table and perform bootstrap on per-window median _r^2_ to get confidence intervals, then plot : `01_scripts/utils/reorder_LD.py`, `01_scripts/utils/LD_analysis.R` (these scripts were produced by Florent Sylvestre)

  
