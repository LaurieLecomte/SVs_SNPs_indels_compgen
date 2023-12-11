#!/bon/env python3

"""
Script produced by Florent Sylvestre

This script uses a modified pairwise SNP-SV LD file (PLINK output : CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2, with manually added columns for END, CAND_SVTYPE, CAND_SVLEN) 
as input and reorders it so that the SV is always marker B and the SNP is always marker A
"""


## Libraries
import sys

## Parsing input
input_path = sys.argv[1]
output_path = sys.argv[2]

## Function
def reorder(ld_line):
    l = ld_line.split("\t")
    infos = {'chra': l[0],
             'posa': l[1],
             'bpa': l[2],
             'chrb': l[3],
             'posb':l[4],
             'bpb': l[5],
             'other': "\t".join(l[6:])}

    if infos["bpa"] == ".":
        return ld_line

    return f"{infos['chrb']}\t{infos['posb']}\t{infos['bpb']}\t{infos['chra']}\t{infos['posa']}\t{infos['bpa']}\t{infos['other']}"


## Main
with open(input_path, "r") as inp:
    with open(output_path, "w") as oup:
        oup.write(inp.readline())
        for line in inp:
            oup.write(reorder(line))
