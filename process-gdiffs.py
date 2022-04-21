#!/usr/bin/env python

"""
process-gdiffs.py by Rohan Maddamsetti.

I take the evolved genomes, and write out evolved-mutations.csv for
downstream analysis in R.

"""

from os.path import join, dirname, realpath
from os import listdir, walk
import pandas as pd

## add genomediff parser to path.
import sys
sys.path.append("genomediff-python")
import genomediff

def write_evolved_pop_mutations(evol_pop_labels, mixed_pop_paths, outf):    
    outfh = open(outf,"w")
    outfh.write("Sample,Transposon,Plasmid,Population,Mutation,Mutation_Category,Gene,Position,Allele,Frequency\n")
    
    for input_dir in mixed_pop_paths:
        for path, subdirs, files in walk(input_dir):
            for f in [x for x in files if x.endswith('annotated.gd')]:
                full_f = join(path, f)
                infh = open(full_f, 'r', encoding='utf-8')
                sample = full_f.split("mixed-pops/")[1].split("/output/evidence/")[0]
                my_row = evol_pop_labels[evol_pop_labels.Sample == sample].iloc[0]
                transposon = my_row["Transposon"]
                plasmid = my_row["Plasmid"]
                population = str(int(my_row["Population"]))

                gd = genomediff.GenomeDiff.read(infh)
                muts = gd.mutations
                muttype = ""
                allele = ""
                for rec in muts:
                    pos = str(rec.attributes['position'])
                    mut_category = rec.attributes['mutation_category']
                    gene = rec.attributes['gene_name']
                    ## handle transpositions of the synthetic Tn5 transposon.
                    if "repeat_name" in rec.attributes and rec.attributes["repeat_name"] == " Synthetic Tn5": ## the leading space is important.
                        gene = "tetA-Tn5-" + gene
                    freq = str(rec.attributes['frequency'])
                    if 'new_seq' in rec.attributes:
                        allele = rec.attributes['ref_seq'] + "->" + rec.attributes['new_seq']
                    else:
                        allele = rec.type
                    if rec.type == 'SNP':
                        muttype = rec.attributes['snp_type']
                    else:
                        muttype = rec.type
                    mutation_row_data = ','.join([sample, transposon, plasmid, population, muttype, mut_category, gene, pos, allele, freq])
                    outfh.write(mutation_row_data + "\n")


def main():
    srcdir = dirname(realpath(__file__))
    assert srcdir.endswith("src")
    projdir = dirname(srcdir)
    assert projdir.endswith("transposon-plasmid-evolution")

    genome_results_dir = join(projdir,"results","draft-manuscript-1A","genome-analysis")
    breseq_pops_results_dir = join(genome_results_dir,"mixed-pops")

    pop_clone_label_f = join(projdir,"data","draft-manuscript-1A","populations-and-clones.csv")
    
    all_pop_clone_labels = pd.read_csv(pop_clone_label_f)

    mixed_pops = [x for x in listdir(breseq_pops_results_dir)]
    mixed_pop_paths = [join(breseq_pops_results_dir, x, "output", "evidence") for x in mixed_pops if x.startswith('RM')]
    
    ''' tabulate the mutations in the evolved populations.'''
    evolved_pop_labels = all_pop_clone_labels[all_pop_clone_labels['Sample'].isin(mixed_pops)]
    mut_table_outf = join(genome_results_dir,"evolved_mutations.csv")
    write_evolved_pop_mutations(evolved_pop_labels, mixed_pop_paths, mut_table_outf)

    
main()
    
