#!/usr/bin/env bash

## assemble-genomes.sh by Rohan Maddamsetti.
## This shell script does a bunch of genome assembly tasks using breseq and gdtools.
## COMMENT OUT all tasks that don't need to be done.
## Potential TODO: rewrite this task and re-script these tasks using SnakeMake.


################################################################################
## first, assemble the ancestral genomes using the NEB5-alpha-NZ_CP017100.gb reference genome.

## assemble the B30 + pUC ancestral strain RM6-176-18 using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/RM6-176-18 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B30-miniTn5-TetA.gb -r ../data/genome-sequencing/A18-pUC.gb ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_18/*.fastq.gz"

## assemble the B30 ancestral strain RM6-200-6 using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/RM6-200-6 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B30-miniTn5-TetA.gb ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_6/*.fastq.gz"

## assemble the B30 + p15A ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/RM7-72-1 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B30-miniTn5-TetA.gb  -r ../data/genome-sequencing/A31-p15A.gb ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_1/*.fastq.gz"

## assemble the B20 + pUC ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/RM7-72-2 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B20-miniTn5-TetA.gb  -r ../data/genome-sequencing/A18-pUC.gb ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_2/*.fastq.gz"

## assemble the B20 ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/RM7-72-3 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B20-miniTn5-TetA.gb ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_3/*.fastq.gz"

## assemble the B20 + p15A ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/RM7-72-4 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B20-miniTn5-TetA.gb  -r ../data/genome-sequencing/A31-p15A.gb ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_4/*.fastq.gz"

################################################################################
## apply mutations in the ancestral B30+A18 genome, RM6-176-18, to the NEB5-alpha-NZ_CP017100.gb reference genome.

## NOTE: the output GFF3 files are combined references for the chromosome, the transposon, and the plasmid.

##gdtools APPLY -o ../results/genome-analysis/RM6-176-18.gff3 -f GFF3 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B30-miniTn5-TetA.gb -r ../data/genome-sequencing/A18-pUC.gb  ../results/genome-analysis/RM6-176-18/output/output.gd

## apply mutations in the ancestral B30+noplasmid genome, RM6-200-6, to the NEB5-alpha-NZ_CP017100.gb reference genome.

##gdtools APPLY -o ../results/genome-analysis/RM6-200-6.gff3 -f GFF3 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B30-miniTn5-TetA.gb ../results/genome-analysis/RM6-200-6/output/output.gd

## apply mutations in the ancestral B30+A31 genome, RM7-72-1, to the NEB5-alpha-NZ_CP017100.gb reference genome.
##gdtools APPLY -o ../results/genome-analysis/RM7-72-1.gff3 -f GFF3 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B30-miniTn5-TetA.gb -r ../data/genome-sequencing/A31-p15A.gb  ../results/genome-analysis/RM7-72-1/output/output.gd

## apply mutations in the ancestral B20+pUC genome, RM7-72-2, to the NEB5-alpha-NZ_CP017100.gb reference genome.
##gdtools APPLY -o ../results/genome-analysis/RM7-72-2.gff3 -f GFF3 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B20-miniTn5-TetA.gb -r ../data/genome-sequencing/A18-pUC.gb  ../results/genome-analysis/RM7-72-2/output/output.gd

## apply mutations in the ancestral B20+noplasmid genome, RM7-72-3, to the NEB5-alpha-NZ_CP017100.gb reference genome.
##gdtools APPLY -o ../results/genome-analysis/RM7-72-3.gff3 -f GFF3 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B20-miniTn5-TetA.gb ../results/genome-analysis/RM7-72-3/output/output.gd

## apply mutations in the ancestral B20+p15A genome, RM7-72-4, to the NEB5-alpha-NZ_CP017100.gb reference genome.
##gdtools APPLY -o ../results/genome-analysis/RM7-72-4.gff3 -f GFF3 -r ../data/genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/genome-sequencing/B20-miniTn5-TetA.gb -r ../data/genome-sequencing/A31-p15A.gb  ../results/genome-analysis/RM7-72-4/output/output.gd

################################################################################
## now, test the new references, by re-mapping reads.

sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/clones/remapped-RM6-176-18 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_18/*.fastq.gz"

sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/clones/remapped-RM6-200-6 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_6/*.fastq.gz"

sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/clones/remapped-RM7-72-1 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_1/*.fastq.gz"

sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/clones/remapped-RM7-72-2 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_2/*.fastq.gz"

sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/clones/remapped-RM7-72-3 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_3/*.fastq.gz"

sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/clones/remapped-RM7-72-4 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_4/*.fastq.gz"


################################################################################
## now assemble evolved mixed population samples with respect to their ancestor.
## set the minimum Phred score threshold for base quality to be 10, and polymorphic sites must have at least 4 reads on each strand supporting the mutation.

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM6-176-13 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_13/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM6-176-14 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_14/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM6-176-15 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_15/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM6-176-16 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_16/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM6-176-17 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_17/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM6-200-1 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_1/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM6-200-2 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_2/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM6-200-3 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_3/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM6-200-4 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_4/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM6-200-5 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_5/*.fastq.gz"

################################################################################
## now assemble the two clonal genomes from the RM6-200-3 mixed pop. sample.

sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/clones/RM6-147-1 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210630/063021_150/*.fastq.gz"

sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/clones/RM6-147-2 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti210630/063021_151/*.fastq.gz"

################################################################################
## Let's assemble the DH5a + B30 + A18 pUC ancestor to see whether it has any secondary mutations (say the one mutation in all the evolved strains).
## D'oh! I already did this control! Let's compare this sample to that one, as an extra check.

sbatch --mem=2G -c 1 --wrap="breseq -o ../results/genome-analysis/clones/RM7-60-5 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_60_5/*.fastq.gz"

################################################################################
## Now let's assemble the DH5a + B20 evolved populations.

## RM7.4.31-35 are DH5a + B20 evolved pops 1-5.
sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-31 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_31/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-32 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_32/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-33 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_33/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-34 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_34/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-35 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_35/*.fastq.gz"

## RM7.4.36-40 are DH5a + B20 + A18 evolved pops 1-5.
sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-36 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_36/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-37 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_37/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-38 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_38/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-39 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_39/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-40 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_40/*.fastq.gz"


## RM7.4.41-45 are DH5a + B30 + A31 evolved pops 1-5.
sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-41 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_41/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-42 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_42/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-43 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_43/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-44 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_44/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-4-45 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_45/*.fastq.gz"


################################################################################
## RM7.72.5-9 is B20 + A31 (Tet 0 control populations).
sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-5 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_5/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-6 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_6/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-7 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_7/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-8 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_8/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-9 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_9/*.fastq.gz"

## RM7.72.10-14 is B20 + A31 (Tet 50 populations).
sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-10 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_10/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-11 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_11/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-12 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_12/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-13 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_13/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-14 -r ../results/genome-analysis/RM7-72-4.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_14/*.fastq.gz"

## RM7.72.15-19 is B30+noplasmid (Tet 0 control populations).
sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-15 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_15/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-16 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_16/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-17 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_17/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-18 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_18/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-19 -r ../results/genome-analysis/RM6-200-6.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_19/*.fastq.gz"

## RM7.72.20-24 is B30+A18 (Tet 0 control populations).
sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-20 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_20/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-21 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_21/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-22 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_22/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-23 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_23/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-24 -r ../results/genome-analysis/RM6-176-18.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_24/*.fastq.gz"

## RM7.72.25-29 is B20 + A18 (Tet 0 control populations).
sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-25 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_25/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-26 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_26/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-27 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_27/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-28 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_28/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-29 -r ../results/genome-analysis/RM7-72-2.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_29/*.fastq.gz"

## RM7.72.30-34 is B20+noplasmid (Tet 0 control populations).
sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-30 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_30/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-31 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_31/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-32 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_32/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-33 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_33/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-34 -r ../results/genome-analysis/RM7-72-3.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_34/*.fastq.gz"

## RM7.72.35-39 is B30 + A31 (Tet 0 control populations).
sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-35 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_35/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-36 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_36/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-37 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_37/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-38 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_38/*.fastq.gz"

sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 --base-quality-cutoff 10 -o ../results/genome-analysis/mixed-pops/RM7-72-39 -r ../results/genome-analysis/RM7-72-1.gff3 ../data/genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_39/*.fastq.gz"
