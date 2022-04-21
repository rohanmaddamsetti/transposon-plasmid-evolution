""" 

specificity-analysis.jl by Rohan Maddamsetti.

I use Fisher's exact test to test genes for the clustering of occurrences of mutations
in that gene within treatments.

I also use a randomization test to test for similarity of metagenomic evolution within
treatments (similar to analyses that use Dice coefficient as a similarity metric),

and use a randomization test to compare the similarity of metagenomic evolution
to external genomic datasets from evolution experiments with Tetracycline selection.

"""

using DataFrames, DataFramesMeta, CSV, HypothesisTests, StatsBase, Random, FLoops


function GeneSpecificityFisherTest(gene_specificity_summary, gene)
    """
    The idea here is the following contingency table:

                                                   pUC  | noPlasmid
    Number of metagenomes w/  mutations in gene x | A   |     B
    Number of metagenomes w/o mutations in gene x | C   |     D

    """
    
    POPS_PER_TREATMENT = 5
        
    cur_summary = @subset(gene_specificity_summary, :Gene .== gene)
    A = 0
    try ## when A == 0, there is no corresponding row in the df.
        A = @subset(cur_summary, :Treatment .== "B30_pUC").NumPops[1]
    catch BoundsError end

    C = 0
    try ## when C == 0, there is no corresponding row in the df.
        C = @subset(cur_summary, :Treatment .== "B30_None").NumPops[1]
    catch BoundsError end
    
    B = POPS_PER_TREATMENT - A
    D = POPS_PER_TREATMENT - C
    gene_test = FisherExactTest(A,B,C,D)
    pval = pvalue(gene_test)
    if pval < 0.05
        println(gene)
        println(gene_test)
        println()
    end
end


function RunGeneSpecificityAnalysis(evolved_mutations)
    ## for documentation of this syntax, see:
    ## https://juliadata.github.io/DataFramesMeta.jl/stable/#Comparison-with-dplyr-and-LINQ
    gene_specificity_summary = @chain evolved_mutations begin
        @rtransform(:Treatment = :Transposon * "_" * :Plasmid)
        groupby([:Gene, :Treatment])
        @combine(:MutationCount = length(:Mutation),
                 :SummedAlleleFrequency = sum(:Frequency),
                 :NumPops = length(unique(:Sample)))
        @orderby(-:MutationCount, :Gene, :Treatment)
    end

    for gene in unique(gene_specificity_summary.Gene)
        GeneSpecificityFisherTest(gene_specificity_summary, gene)
    end
end


function GeneSampleTotalAlleleFreqMatrix(evolved_mutations, treatment)
    genes_for_rows = unique(evolved_mutations.Gene)
    gene_to_matrix_idx = Dict(gene => i for (i,gene) in enumerate(genes_for_rows))

    treatment_data = @rsubset(evolved_mutations, :Treatment == treatment)
    samples_for_cols = unique(treatment_data.Sample)

    gene_sample_matrix = zeros(length(genes_for_rows),length(samples_for_cols))
    ## fill in the values of the matrix with the summed allele frequency at the given locus.
    for (j, cur_sample) in enumerate(samples_for_cols)
        pop_data = @rsubset(treatment_data, :Sample == cur_sample)
        for x = 1:nrow(pop_data)
            cur_gene = pop_data.Gene[x]
            cur_frequency = pop_data.Frequency[x]
            i = gene_to_matrix_idx[cur_gene]
            gene_sample_matrix[i,j] += cur_frequency
        end
    end
    return gene_sample_matrix
end


function WeightedJaccardSimilarity(vec1, vec2)
    @assert(length(vec1) == length(vec2))
    sum_min = 0
    sum_max = 0
    for i = 1:length(vec1)
        cur_min = minimum([vec1[i], vec2[i]])
        cur_max = maximum([vec1[i], vec2[i]])
        sum_min += cur_min
        sum_max += cur_max
    end
    return (sum_min/sum_max)
end


function MeanWeightedJaccardSimilarityForTreatment(M)
    Sum_Similarity = 0.0
    ncols_of_M = size(M,2)
    for i = 1:ncols_of_M
        for j = 1:ncols_of_M
            Sum_Similarity += WeightedJaccardSimilarity(M[:,i], M[:,j])
        end
    end
    Mean_Similarity = Sum_Similarity/(ncols_of_M^2)
    return Mean_Similarity
end


function MeanWeightedJaccardSimilarityBetweenTreatments(M1, M2)
    Sum_Similarity = 0.0
    ncols_of_M1 = size(M1,2)
    ncols_of_M2 = size(M2,2)
    for i = 1:ncols_of_M1
        for j = 1:ncols_of_M2
            Sum_Similarity += WeightedJaccardSimilarity(M1[:,i], M2[:,j])
        end
    end
    Mean_Similarity = Sum_Similarity/(ncols_of_M1 * ncols_of_M2)
    return Mean_Similarity
end


function WeightedJaccardRandomizationTest(M1, M2, N = 100_000, ncores = 4)

    glued_matrix = hcat(M1,M2)
    Grand_Mean = MeanWeightedJaccardSimilarityForTreatment(glued_matrix)

    M1_mean = MeanWeightedJaccardSimilarityForTreatment(M1)
    M2_mean = MeanWeightedJaccardSimilarityForTreatment(M2)

    ## following notation in Deatherage et al. (2014).
    Sw = mean([M1_mean,M2_mean]) ## within treatment similarity
    Sb = MeanWeightedJaccardSimilarityBetweenTreatments(M1,M2)

    JaccardDifference_in_Data = Sw - Sb
    
    ncols_of_M1 = size(M1,2)
    ncols_of_M2 = size(M2,2)

    num_random_draws_more_similar_than_data = 0
    @floop ThreadedEx(basesize = N ÷ ncores) for _ in 1:N
        ## shuffle columns of the glued_matrix. OK to modify the original object.
        rand_glued_matrix = glued_matrix[:, shuffle(1:end)]
        ## split into two randomized treatments
        rand_M1 = rand_glued_matrix[:, 1:ncols_of_M1]
        rand_M2 = rand_glued_matrix[:, ncols_of_M1+1:(ncols_of_M1+ncols_of_M2)]
        
        rand_M1_mean = MeanWeightedJaccardSimilarityForTreatment(rand_M1)
        rand_M2_mean = MeanWeightedJaccardSimilarityForTreatment(rand_M2)
        
        rand_Sw = mean([rand_M1_mean,rand_M2_mean]) ## within treatment similarity
        rand_Sb = MeanWeightedJaccardSimilarityBetweenTreatments(rand_M1, rand_M2)
        
        Randomized_JaccardDifference = rand_Sw - rand_Sb
        greater_than_data = (Randomized_JaccardDifference > JaccardDifference_in_Data)

        @reduce() do (num_random_draws_more_similar_than_data = 0; greater_than_data)
            num_random_draws_more_similar_than_data += greater_than_data
        end
    end

    p_val = num_random_draws_more_similar_than_data / N
    return p_val
end


function filter_Card2021_mutations(card2021_mutations, evolved_mutations)
    ## have to get a little fancy in order to parse the genes in these data ...
    
    ## get all genes in evolved_mutations, for pattern matching.
    loci_of_interest = unique(evolved_mutations.Gene)
    loci_to_match = []
    for s in loci_of_interest
        gene_list = [String(g) for g in split(s, "/")]
        for g in gene_list
            push!(loci_to_match, g)
        end
    end
    loci_regex = Regex(join(loci_to_match,"|"))
    
    gene_in_gene_column = [occursin(loci_regex, x) ? true : false for x in card2021_mutations.gene]
    gene_in_description_column = [occursin(loci_regex, x) ? true : false for x in card2021_mutations.description]

    card2021_mutations.matches_gene = gene_in_gene_column
    card2021_mutations.matches_description = gene_in_description_column
    
    filtered_card2021_mutations = @rsubset(card2021_mutations, :matches_gene || :matches_description)
    return filtered_card2021_mutations
end


################################################################################
## import the data.
## NOTE: This currently has all mutations >5%, should probably filter by
## Frequency > 10%.
## But filtering by >10% is way too conservative for the no Plasmid data,
## which has great coverage-- 200X-300X. 5% in those data is like
## 20% in the pUC plasmid treatment.

## for instance, I find parallel evolution at nucleotide level at ilvC/ppiC,
## which is 7-9% in 3 of the no plasmid populations.

evolved_mutations = CSV.read(
    "../results/draft-manuscript-1A/genome-analysis/evolved_mutations.csv", DataFrame)
## Add a Treatment Column.
@rtransform!(evolved_mutations, :Treatment = :Transposon * "_" * :Plasmid)
## for now, I keep mutations > 5% in the no plasmid data, and keep
## mutations >10% in the pUC plasmid data.
@rsubset!(evolved_mutations, :Plasmid == "None" || :Frequency > 0.1)

pUC_treatment_mutations = @rsubset(evolved_mutations, :Treatment == "B30_pUC")
noPlasmid_treatment_mutations = @rsubset(evolved_mutations, :Treatment == "B30_None")
###############################################################################
## Gene-specific specificity analyses. Use Fisher's exact test, and extensions
## to metagenomic data.

## Fisher's exact test is really too conservative but it gets the job done.
## TODO: find a better way that doesn't needlessly throw out data.

## TODO: generalize this function to other pairs of treatments, when the data arrives.
## (treatments are currently hardcoded).

RunGeneSpecificityAnalysis(evolved_mutations)

## Based on this population-level presence/absence analysis,
## following genes show evidence of specificity:
## lysO/aqpZ
## yohP/dusC
## tetA
## KanR (tet-transposon insertions on the plasmid)
## frmR/yaiO
## yeaD/yeaE

################################################################################
## I use the weighted Jaccard index to calculate the specificity of evolutionary
## paths within experimental treatments.
## See: https://en.wikipedia.org/wiki/Jaccard_index#Weighted_Jaccard_similarity_and_distance
## https://en.wikipedia.org/wiki/Fuzzy_set
## https://en.wikipedia.org/wiki/T-norm

## First, make matrices for the metagenomes in each treatment. Each row is a different
## gene or intergenic region that is mutated in at least one mixed pop. sample.

noPlasmid_matrix = GeneSampleTotalAlleleFreqMatrix(evolved_mutations, "B30_None")
pUC_matrix = GeneSampleTotalAlleleFreqMatrix(evolved_mutations, "B30_pUC")

## more parallelism in the no Plasmid treatment than in the pUC treatment
## (given differences in chromosomal coverage)
noPlasmid_Mean = MeanWeightedJaccardSimilarityForTreatment(noPlasmid_matrix)
pUC_Mean = MeanWeightedJaccardSimilarityForTreatment(pUC_matrix)

## I follow the logic of the calculation in Deatherage et al. (2014) in PNAS.
glued_matrix = hcat(noPlasmid_matrix,pUC_matrix)

Grand_Mean = MeanWeightedJaccardSimilarityForTreatment(glued_matrix)
withinTreatment_Mean = mean([noPlasmid_Mean,pUC_Mean])
betweenTreatment_Mean = MeanWeightedJaccardSimilarityBetweenTreatments(pUC_matrix, noPlasmid_matrix)
                                                                       
## now, do the randomization test for statistical significance.
@time WeightedJaccardRandomizationTest(pUC_matrix, noPlasmid_matrix, 100_000)
## p = 0.00026 with 100,000 bootstraps.

################################################################################
## comparison of evolutionary paths to other published evolution experiments.

## Bottery 2017 in Nature Ecology and Evolution.
## from the Results:
## Thirty independent isogenic populations of E. coli MG1655 carrying the MDR
## plasmid RK214, which encodes resistances to TET and AMP, were experimentally
## evolved for ~530 generations (80 days), under five antibiotic treatments
## (six independently evolving lines per treatment):
## no antibiotic (N); AMP (A); TET (T); AMP plus TET (AT); and 24 h cycling
## between AMP and TET (A/T) (see Methods). Six control populations of plasmid-free
## MG1655 were similarly experimentally evolved with no antibiotic (C).
## Plasmids remained at high frequency in all populations for the duration of the
## selection experiment. Plasmid-free segregants were observed only at very low
## frequency in two of the six populations from treatment N (Supplementary Fig. 1),
## whereas transposition of resistance genes from RK2 onto the host’s chromosome was
## never observed.

## import the Bottery 2017 data.
bottery_mutations = CSV.read("../data/draft-manuscript-1A/Bottery2017-Figure_2_S2_S3__data.csv", DataFrame)
## replace "-" with "/" in the gene_name column for compatibility with my data.
## This goes into a new Gene column.
@rtransform!(bottery_mutations, :Gene = replace(:gene_name, r"-" => "/"))

relevant_bottery_muts = @rsubset(bottery_mutations, :Gene in evolved_mutations.Gene)


pUC_and_bottery_muts = @rsubset(bottery_mutations, :Gene in pUC_treatment_mutations.Gene)
noPlasmid_and_bottery_muts = @rsubset(bottery_mutations, :Gene in noPlasmid_treatment_mutations.Gene)

## there is much more gene-level parallelism between the no Plasmid treatment
## and the Bottery experiment than between the pUC treatment and the Bottery experiment.


bottery_frmR_yaiO_muts = @rsubset(bottery_mutations, :Gene == "frmR/yaiO")
my_frmR_yaiO_muts = @rsubset(evolved_mutations, :Gene == "frmR/yaiO")
## by examining the GenBank: U00096.3 K-12 MG1655 reference used by Bottery et al. (2017),
## the 380012 +G frmR/yaiO K-12 MG1655 mutation 
## 286010, -131/+56 frmR/yaiO in DH5-alpha are homologous!
## They affect the same G(9) or G(10) repeat.
## NOTA BENE: this mutation was seen in 1/6 Tet lines, and 1/6 no plasmid no antibiotic
## lines in the Bottery dataset. Another mutation at 379915, or at the -34 position
## relative to frmR was observed in 1/6 of the +plasmid -antibiotic treatement in
## Bottery et al. (2017).

## phase variation in the same G(8) repeat is in the LTEE Genomes (see Shiny app).
## only in populations Ara-2, Ara-4, Ara+3, Ara-3. However, these mutations
## occur in the hypermutator LTEE populations (but not Ara-1 nor Ara+6).

## This paper reports that frmR represses conjugation of IncP1-alpha plasmids:
## Isolation and Analysis of Donor Chromosomal Genes Whose Deficiency Is Responsible
## for Accelerating Bacterial and Trans-Kingdom Conjugations by IncP1 T4SS Machinery
## the frm operon is induced in response to formaldehyde

## idea for how to do this analysis, from Deatherage et al. (2014):
## We measured the similarity of genomic evolution in the TEE to the LTEE by calculating
## a recapitulation index (RI) equal to the percentage of qualifying mutations in a
## particular thermal group that affected genes that subsequently acquired a qualifying
## mutation in an LTEE clone from the source population. For example, 10 of the 22
## qualifying mutations in the 37 °C TEE group are found in hslU (five lines), ybaL (two
## lines), iclR, insJ-5, and mrdA, and all five of these genes also harbor qualifying
## mutations in the 20K clone in the LTEE population from which the progenitor was
## sampled. The 20K RI is therefore 45% (10/22) for the set of TEE lines that
## evolved at 37 °C. 

####################################

## Let's compare to mutations in Kyle Card's 2021 PNAS paper.
## nothing really to write home about-- one deletion removing frmR and yaiO
## in Tet treatment, and some mutations affecting envZ, which only occurs
## at 5% frequency in RM6.200.1.
card2021_mutations = CSV.read("../data/draft-manuscript-1A/Card2021_Mutations_curated.csv", DataFrame)
relevant_card2021_mutations = filter_Card2021_mutations(card2021_mutations, evolved_mutations)

####################################
## Lukacisinova 2020 in Nature Communications.
## the main finding here is the lack of parallelism between these data
## and the data in my work. Probably the big difference is that
## there's no tet gene in these strains, so one of the main pathways
## in this paper is gene duplication/amplification of native efflux pumps.

lukacisinova2020_tet_mutations = CSV.read("../data/draft-manuscript-1A/Lukacisinova2020-TET-S1Data.csv", DataFrame)

lukacisinova_genes = unique(lukacisinova2020_tet_mutations.gene)

lukacisinova_regex = Regex(join(lukacisinova_genes,"|"))
    
gene_in_lukacisinova_gene_column = [occursin(lukacisinova_regex, x) ? true : false for x in evolved_mutations.Gene]

evolved_mutations.matches_lukacisinova = gene_in_lukacisinova_gene_column

## minimal overlap: three mutations in two matching loci: envZ and acrR.
evolved_mutations_in_lukacisinova_genes = @rsubset(evolved_mutations, :matches_lukacisinova)

