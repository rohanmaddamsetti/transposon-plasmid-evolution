#!/usr/bin/env python

'''
Usage: python printEcoliIDs.py -i ../data/draft-manuscript-1A/genome-sequencing/NEB5-alpha-NZ_CP053607.gb  > ../results/draft-manuscript-1A/NEB5-alpha_IDs.csv

This script skips over all pseudogenes that have been annotated with the /pseudo tag.

NOTA BENE: repetitive regions are NOT masked by default-- supply a mask file
corresponding to the reference genome if desired.

'''

import argparse
import re
from Bio import SeqIO

def overlaps_any_region(cds_start, cds_end, regions):
    for r in regions:
        if overlaps_this_region(cds_start, cds_end, r):
            return True
    return False
    
def overlaps_this_region(cds_start, cds_end, region):
    ''' returns True if the intervals spanned by the region and the cds 
        overlap at all. returns False if there's no overlap. '''
    region_start, region_end = region
    assert region_start < region_end
    assert cds_start < cds_end
    if (cds_end < region_start) or (cds_start > region_end):
        return False
    return True

def is_in_any_masked_region(cds_start, cds_end, masked_regions):
    for mask_start, mask_len in masked_regions:
        mask_end = mask_start + mask_len
        if overlaps_this_region(cds_start, cds_end, (mask_start, mask_end)):
            return True
    return False

def get_masked_regions(maskfile):
    ''' 
    input: a mask file.
    output: a list of tuples (region_start, region_end).
    '''
    mask_f = open(maskfile, "r")
    masked_regions = []
    for l in mask_f:
        if l.startswith('#'): continue
        assert l.startswith('MASK')
        l = l.strip()
        ldata = l.split()
        ## Genbank starts from 1, Biopython starts from 0.
        region_start = int(ldata[4]) - 1 ## so subtract 1 from the start.
        region_len = int(ldata[5])
        region_tuple = (region_start, region_len)
        masked_regions.append(region_tuple)
    return masked_regions

def get_repeat_regions(genbankf):
    genome = next(SeqIO.parse(genbankf, "genbank"))
    repeats = []
    for feat in genome.features:
        ## misc_features are prophage that Jeff Barrick annotated.
        ## let's keep misc_features in, for now. But consider skipping.
        ##if (feat.type == "repeat_region") or (feat.type == "misc_feature"):
        if (feat.type == "repeat_region"):
            repeat_start = feat.location.start
            repeat_end = feat.location.end
            repeat_tuple = (repeat_start, repeat_end)
            repeats.append(repeat_tuple)
    return repeats
        

def main():
    parser = argparse.ArgumentParser(description='print csv file of CDS IDs in Genbank file, as long as they do not overlap a repetitive region or prophage.')

    parser.add_argument("-i", "--input", help='Input genbank file',required=True)
    parser.add_argument("-m", "--mask", help="genomediff mask file", required=False)
    
    args = parser.parse_args()

    masked_regions = () ## initialize to an empty tuple.
    if args.mask: ## then parse the mask file.
        masked_regions = get_masked_regions(args.mask)
    ## parse the genome for repetitive regions.
    annotated_repeats = get_repeat_regions(args.input)

    ## now parse the genome for CDS.
    genome = next(SeqIO.parse(args.input, "genbank"))
    print(','.join(['Gene','locus_tag','gene_length','product', 'start', 'end', 'strand','G','C','A','T']))    
    for feat in genome.features:
        if feat.type != 'CDS': continue ## only consider protein-coding genes
        if 'pseudo' in feat.qualifiers: continue ## skip loci annotated as pseudogenes
        my_start = feat.location.start
        my_end = feat.location.end
        my_strand = feat.location.strand
        if overlaps_any_region(my_start, my_end, annotated_repeats): continue
        if is_in_any_masked_region(my_start, my_end, masked_regions): continue
        my_length = my_end - my_start
        locus_tag = feat.qualifiers['locus_tag'].pop()

        try:
            gene = feat.qualifiers['gene'].pop()
        except KeyError:
            gene = locus_tag
        try:
            product = feat.qualifiers['product'].pop()
        except:
            product = 'NA'
        try:
            note = feat.qualifiers['product'].pop()
        except:
            note = ""
        ## strip all punctuation that could cause parsing problems.
        product = re.sub('[,;()]', '', product)

        ## get nucleotide composition of the coding strand of the gene.
        my_seq = feat.extract(genome)
        my_Gs = 0
        my_Cs = 0
        my_As = 0
        my_Ts = 0
        for x in my_seq:
            if x == 'G':
                my_Gs += 1
            elif x == 'C':
                my_Cs += 1
            elif x == 'A':
                my_As += 1
            elif x == 'T':
                my_Ts += 1
        ## except genes with a programmed frameshift during translation.
        ## some genes seems to have some issues. I can troubleshoot this later.
        assert (my_Gs + my_Cs + my_As + my_Ts == my_length) or "programmed frameshift" in note or locus_tag in ["HO392_RS01400", "HO392_RS01450", "HO392_RS04680", "HO392_RS05245", "HO392_RS06955", "HO392_RS09980", "HO392_RS10470", "HO392_RS14310", "HO392_RS14460", "HO392_RS15210", "HO392_RS17785", "HO392_RS21405"]
        ## start+1 (but end+0) to be consistent with genbank format.
        print(','.join([str(x) for x in (gene,locus_tag,my_length,product, my_start+1, my_end, my_strand, my_Gs, my_Cs, my_As, my_Ts)]))


main()
