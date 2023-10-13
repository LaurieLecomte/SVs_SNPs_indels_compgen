#!/usr/bin/env python3
"""Compute LD between SVs and SNPs, using 2 tab-delimited files containing CHROM, POS, END, ID and GT fields for a given chromosome.
 
First produce required files using bcftools in script 01_scripts/utils/LD_SVs_SNPs.sh

Usage:
    <program> SV_FILE SNP_FILE MAX_DIST OUTFILE

Parallel:
    
"""

# Modules
import re
import sys


# Parsing user input
try:
    SV_FILE = sys.argv[1]
    SNP_FILE = sys.argv[2]
    MAX_DIST = int(sys.argv[3])
    OUTFILE = sys.argv[4]
    #flanking_size = int(sys.argv[3])
    #output_genes = sys.argv[4]
except:
    print(__doc__)
    sys.exit(1)

# Functions
## Count genotypes at given site
def count_gt(gt):
    ## Collapse genotypes
    gt_coll = []
    for x in gt:
        gt_coll.append(re.sub("\/", '', x))
    ## Count occurences of each genotype   
    gen_patterns = ["00", "01", "10", "11", ".."]
    gt_counts = [gt_coll.count(pat) for pat in gen_patterns]
    ## Combine het counts together
    AA_count = gt_counts[0]
    Aa_count = gt_counts[1]+gt_counts[2]
    aa_count = gt_counts[3]
    NN_count = gt_counts[4]
    
    final_counts = [AA_count, Aa_count, aa_count, NN_count]
    
    return(final_counts)

## Get allelic frequencies from genotype count at each site
def allele_freqs(gt_counts):
    n_alleles = sum(gt_counts) - gt_counts[3]
    p = ((2*gt_counts[0]) + gt_counts[1])/(2*n_alleles)
    q = (1 - p) 
    freqs = [p, q]
    
    return(freqs)
    
## Get possible haplotypes for the given pair of sites
def get_haplo(gt1, gt2):
    ## Collapse genotypes into a single string
    gt1_coll = ''.join([re.sub("\|", '', x) for x in gt1])
    gt2_coll = ''.join([re.sub("\|", '', x) for x in gt2])
    
    ## Generate haplotypes
    #haps = [a+b for a in gt1_coll for b in gt2_coll]
    haps = [(gt1_coll[x]+gt2_coll[x]) for x in range(0, len(gt1_coll))]
    
    return(haps)
    
# Main function
with open(OUTFILE, "wt") as outfile:
    with open(SV_FILE) as SVfile:
        for line in SVfile:
            SV = line.strip().split("\t")
        
            # Extract SV CHROM, POS and ID
            SV_CHROM = SV[0]
            SV_POS = int(SV[1])
        
            SV_END = int(SV[2])
            SV_ID = SV[3]
            
            MID_POS = ((SV_END - SV_POS) /2) + SV_POS

            # Extract GT fields (4th field to last field of row) 
            SV_GTs = SV[5:-1] 
        
            # Loop over each SNP
            with open(SNP_FILE) as SNPfile:
                for line in SNPfile:
                    SNP = line.strip().split("\t")
                
                    # Extract CHROM, POS and ID
                    SNP_CHROM = SNP[0]
                    SNP_POS = int(SNP[1])
                    SNP_END = int(SNP[2])
                    SNP_ID = SNP[3]
                    
                    # Check if this SNPs is located within MAX_DIST bp of SV
                    ## Calculate distance between both sites
                    ### Decide which value to use for SV position
                    if SV_END != SV_POS:
                        SV_SITE = MID_POS
                    elif SV_END == SV_POS:
                        SV_SITE = SV_POS
                    
                    ### Calculate physical distance     
                    if SNP_POS < SV_POS:
                        dist = SV_SITE - SNP_POS
                    elif SNP_POS > SV_POS:
                        dist = SNP_POS - SV_SITE
                        
                    # Perform LD calculation if dist criteria is met                   
                    if dist < MAX_DIST:
                
                        # Extract GT fields (4th field to last field of row)
                        SNP_GTs = SNP[5:-1]
        
                        # Get allele frequencies from genotype counts at both sites
                        pA = allele_freqs(count_gt(SV_GTs))[0]
                        pB = allele_freqs(count_gt(SNP_GTs))[0]
                
                        # Calculate pAB:
                        haplotypes = get_haplo(SV_GTs, SNP_GTs)
                        ### Now, just count the fraction of times "00" is the hap:
                        pAB = haplotypes.count('00')/len(haplotypes)
                
                        # Calculate D and r-squared
                        D = pAB - (pA * pB)
                        rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))

                        # Write output
                        out_line = [SV_CHROM, SV_POS, SV_END, SV_ID, SNP_POS, D, rsq, dist]
                        outfile.write('\t'.join([str(x) for x in out_line]) + "\n")
                        #outfile.write('\t'.join([str(x) for x in out_line]))