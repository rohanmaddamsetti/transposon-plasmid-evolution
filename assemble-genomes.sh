#!/usr/bin/env bash

## assemble-genomes.sh by Rohan Maddamsetti.
## This shell script does a bunch of genome assembly tasks using breseq and gdtools.
## COMMENT OUT all tasks that don't need to be done.
## Potential TODO: rewrite this task and re-script these tasks using SnakeMake.

## first, assemble the ancestral genomes using the NEB5-alpha-NZ_CP053607.gb reference genome.

## assemble RM6-176-18 using the NEB5-alpha-NZ_CP053607.gb reference genome.
##breseq -o ../results/genome-analysis/RM6-176-18 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP053607.gb -r ../data/genome-sequencing/B30-miniTn5-TetA.gb -r ../data/genome-sequencing/a18-puc-synoriginal.gb ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_18/*.fastq.gz

## assemble RM6-200-6 using the NEB5-alpha-NZ_CP053607.gb reference genome.
##breseq -o ../results/genome-analysis/RM6-200-6 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP053607.gb -r ../data/genome-sequencing/B30-miniTn5-TetA.gb ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_6/*.fastq.gz

## apply mutations in the ancestral B30+A18 genome, RM6-176-18, to the NEB5-alpha-NZ_CP053607.gb reference genome.

## NOTE: the output GFF3 files are combined references for the chromosome, the transposon, and the plasmid.

##gdtools APPLY -o ../results/genome-analysis/RM6-176-18.gff3 -f GFF3 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP053607.gb -r ../data/genome-sequencing/B30-miniTn5-TetA.gb -r ../data/genome-sequencing/a18-puc-synoriginal.gb  ../results/genome-analysis/RM6-176-18/output/output.gd

## apply mutations in the ancestral B30+noplasmid genome, RM6-200-6, to the NEB5-alpha-NZ_CP053607.gb reference genome.

##gdtools APPLY -o ../results/genome-analysis/RM6-200-6.gff3 -f GFF3 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP053607.gb -r ../data/genome-sequencing/B30-miniTn5-TetA.gb ../results/genome-analysis/RM6-200-6/output/output.gd

## now, test the new references, by re-mapping reads.

##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/remapped-RM6-176-18 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_18/*.fastq.gz"

##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/remapped-RM6-200-6 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_6/*.fastq.gz"

## now assemble evolved mixed population samples with respect to their ancestor.
## set the minimum Phred score threshold for base quality to be 10, and polymorphic sites must have at least 4 reads on each strand supporting the mutation.

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/RM6-176-13 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_13/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/RM6-176-14 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_14/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/RM6-176-15 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_15/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/RM6-176-16 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_16/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/RM6-176-17 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_17/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/RM6-200-1 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_1/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/RM6-200-2 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_2/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/RM6-200-3 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_3/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/RM6-200-4 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_4/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/RM6-200-5 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_5/*.fastq.gz"

## now assemble the two clonal genomes from the RM6-200-3 mixed pop. sample.

##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/RM6-147-1 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210630/063021_150/*.fastq.gz"

##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/RM6-147-2 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210630/063021_151/*.fastq.gz"
