""" 
specificity-analysis.jl by Rohan Maddamsetti.

I use Fisher's exact test to test genes for the clustering of occurrences of mutations
in that gene within treatments.

I also use a randomization test to test for similarity of metagenomic evolution within
treatments (similar to analyses that use Dice coefficient as a similarity metric),

and use a randomization test to compare the similarity of metagenomic evolution
to external genomic datasets from evolution experiments with tetracycline selection.

This analysis builds upon the analysis framework in:
Deatherage et al. (2017) "Specificity of genome evolution in experimental populations of
Escherichia coli evolved at different temperatures"

TODO: correct p-values for multiple testing as needed in the future.
"""

using DataFrames, DataFramesMeta, CSV, StatsBase, HypothesisTests, MultipleTesting, Random, FLoops


function GeneSpecificityFisherTest(specificity_summary, gene, treatment1, treatment2)
    """
    The idea here is the following contingency table:

                                                                                             pUC  | noPlasmid
    Number of metagenomes w/  mutations in gene x   |  A  |  B
    Number of metagenomes w/o mutations in gene x  |  C |  D

Possible values for the treatment_strs:
B20_None_50
B20_p15A_50
B20_pUC_50
B30_None_50
B30_p15A_50
B30_pUC_50
    """
    
    POPS_PER_TREATMENT = 5
        
    cur_summary = @subset(specificity_summary, :Gene .== gene)
    A = 0
    try ## when A == 0, there is no corresponding row in the df.
        A = @subset(cur_summary, :Treatment .== treatment1).NumPops[1]
    catch BoundsError end

    C = 0
    try ## when C == 0, there is no corresponding row in the df.
        C = @subset(cur_summary, :Treatment .== treatment2).NumPops[1]
    catch BoundsError end
    
    B = POPS_PER_TREATMENT - A
    D = POPS_PER_TREATMENT - C
    try
        gene_test = FisherExactTest(A,B,C,D)
        pval = pvalue(gene_test)
        if pval < 0.05
            println(gene)
            println(gene_test)
            println()
        end
    catch ArgumentError end ## If the FisherExactTest fails, return silently.
end


function RunGeneSpecificityAnalysis(evolved_mutations)
    ## for documentation of this syntax, see:
    ## https://juliadata.github.io/DataFramesMeta.jl/stable/#Comparison-with-dplyr-and-LINQ
    gene_specificity_summary = @chain evolved_mutations begin
        @rtransform(:Treatment = :Transposon * "_" * :Plasmid * "_" * string(:Tet))
        groupby([:Gene, :Treatment])
        @combine(:MutationCount = length(:Mutation),
                 :SummedAlleleFrequency = sum(:Frequency),
                 :NumPops = length(unique(:Sample)))
        @orderby(-:MutationCount, :Gene, :Treatment)
    end

    ## let's examine the following pairs of treatments.
    treatments_to_compare = [("B30_None_50", "B30_pUC_50"), ("B20_None_50", "B20_pUC_50"),
                             ("B20_None_50", "B20_p15A_50"), ("B20_pUC_50", "B20_p15A_50"),
                             ("B20_None_50", "B30_None_50"), ("B20_p15A_50", "B30_p15A_50"), ("B20_pUC_50", "B30_pUC_50")]
        
    ## iterate over the pairs of treatments to compare.
    for treatment_tuple in treatments_to_compare
        treatment1, treatment2 = treatment_tuple
        ## critical step: subset the dataframe by the treatments of interest.
        ## we only want to compare genes found in at least one of these treatments.
        sliced_specificity_summary = @rsubset(gene_specificity_summary, :Treatment in treatment_tuple)

        println("Gene Specificity Comparison between " * treatment1 * " and " * treatment2)
        println("***********************************")
        for gene in unique(sliced_specificity_summary.Gene)
            GeneSpecificityFisherTest(gene_specificity_summary, gene, treatment1, treatment2)
        end
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

    ## following notation in Deatherage et al. (2017).
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


## IMPORTANT TODO: I may need to apply a similar reweighting of allele frequency for B20 pUC, and the p15A
## transposition events, since there are multiple transpositions that sometimes combine to have an allele frequency > 100%.
## This is technically possible if some plasmids have double insertions, but if each plasmid has a single insertion,
## then the allele frequencies for each transposition need to be renormalized to sum to 100%.
function reweight_tetA_Tn5_KanR_plasmid_frequencies(df)
    ## For the purpose of this analysis, the allele frequencies of the tetA-Tn5-KanR in the B30 treatment
    ## need to be set to 1.0, to treat  the tetA-plasmid balanced polymorphism as a fixed genotype.
    ## This is a hack... 
    ## See: https://stackoverflow.com/questions/66586623/julia-dataframe-preferred-way-to-update-values-in-one-column-based-on-the-valu/66587724#66587724
    updated_df = df
    updated_df.Frequency = @. ifelse(df.Treatment == "B30_pUC_50" && df.Gene == "tetA-Tn5-KanR", 1.0, df.Frequency)
    return updated_df
end


function read_evolved_mutations()
    ## I wrapped up this code into a function so that refactoring doesn't accidentally break this logic.
    evolved_mutations = CSV.read(
        "../results/genome-analysis/evolved_mutations.csv", DataFrame)
    ## Add a Treatment Column.
    @rtransform!(evolved_mutations, :Treatment = :Transposon * "_" * :Plasmid * "_" * string(:Tet))
    evolved_mutations = reweight_tetA_Tn5_KanR_plasmid_frequencies(evolved_mutations)    
    return evolved_mutations
end


################################################################################
## import the data. All mutation filtering criteria is done when breseq is used to call variants,
## so no additional filters are applied at this stage of the data analysis.
evolved_mutations = read_evolved_mutations() ## wrapping up some ugly code that needs to be run together.

B20_noPlasmid_50_treatment_mutations = @rsubset(evolved_mutations, :Treatment == "B20_None_50")
B20_p15A_50_treatment_mutations = @rsubset(evolved_mutations, :Treatment == "B20_p15A_50")
B20_pUC_50_treatment_mutations = @rsubset(evolved_mutations, :Treatment == "B20_pUC_50")

B30_noPlasmid_50_treatment_mutations = @rsubset(evolved_mutations, :Treatment == "B30_None_50")
B30_p15A_50_treatment_mutations = @rsubset(evolved_mutations, :Treatment == "B30_p15A_50")
B30_pUC_50_treatment_mutations = @rsubset(evolved_mutations, :Treatment == "B30_pUC_50")

## Consider the analysis when we restrict to genes showing parallel evolution.
parallel_genes = @chain evolved_mutations begin
    @rtransform(:Treatment = :Transposon * "_" * :Plasmid)
    groupby(:Gene)
    @combine(:MutationCount = length(:Mutation))
    @rsubset(:MutationCount > 1)
    innerjoin(evolved_mutations, on = :Gene)
end

parallel_evolved_mutations = @rsubset(evolved_mutations, :Gene in parallel_genes.Gene)

###############################################################################
## Gene-specific specificity analyses. Use Fisher's exact test, and extensions
## to metagenomic data.

## Fisher's exact test is really too conservative but it gets the job done.
## TODO: find a better way that doesn't needlessly throw out data;
## this is a population-level presence/absence analysis,.

# the treatments are currently hardcoded in this function.
RunGeneSpecificityAnalysis(evolved_mutations)
## the following genes show evidence of specificity between B30_None_50 and B30_pUC_50.
## tetA
## tetA-Tn5-KanR (tet-transposon insertions on the plasmid)
## –/KanR (repeat expansion on the plasmid, probably not adaptive because happens in Tet 0 treatment).

## the following genes show evidence of specificity between B20_None_50 and B20_pUC_50
## tetA-Tn5-KanR (tet-transposon insertions on the plasmid)
## tetA-Tn5-NEB5A_RS07165 (tet-transposon insertions on the chromosome with B20)

## the following genes show evidence of specificity between B20_None_50 and B20_p15A_50:
## –/KanR (repeat expansion on the plasmid, probably not adaptive because happens in Tet 0 treatment).
## tetA-Tn5-KanR (tet-transposon insertions on the plasmid)
## tetA-Tn5-NEB5A_RS07165 (tet-transposon insertions on the chromosome)

## NO  genes show evidence of specificity between B20_p15A_50 and B20_pUC_50

## the following genes show evidence of specificity between B20_None_50 and B30_None_50.
## tetA-Tn5-NEB5A_RS07165 (tet-transposon insertions on the chromosome with B20, not B30).

## No genes show evidence of specificity between B20_pUC and B30_pUC.

################################################################################
## IMPORTANT TODO: generalize this to multiple treatments, contrasting one treatment against all others.
function WeightedJaccardIndexSpecificityAnalysis(mutations, treatment_str1, treatment_str2, nBootstraps=100_000)
    ## I use the weighted Jaccard index to calculate the specificity of evolutionary
    ## paths within experimental treatments.
    ## See: https://en.wikipedia.org/wiki/Jaccard_index#Weighted_Jaccard_similarity_and_distance
    ## https://en.wikipedia.org/wiki/Fuzzy_set
    ## https://en.wikipedia.org/wiki/T-norm

    println("Treatment 1: " * treatment_str1)
    println("Treatment 2: " * treatment_str2)
    ## First, make matrices for the metagenomes in each treatment. Each row is a different
    ## gene or intergenic region that is mutated in at least one mixed pop. sample.
    treatment1_matrix = GeneSampleTotalAlleleFreqMatrix(mutations, treatment_str1)
    treatment2_matrix = GeneSampleTotalAlleleFreqMatrix(mutations, treatment_str2)
    
    ## Calculate the mean similarity within each treatment and print.
    treatment1_Mean = MeanWeightedJaccardSimilarityForTreatment(treatment1_matrix)
    println("Treatment 1 Mean: ", treatment1_Mean)
    treatment2_Mean = MeanWeightedJaccardSimilarityForTreatment(treatment2_matrix)
    println("Treatment 2 Mean: ", treatment2_Mean)
    
    ## I follow the logic of the calculation in Deatherage et al. (2017) in PNAS.
    glued_matrix = hcat(treatment1_matrix, treatment2_matrix)

    Grand_Mean = MeanWeightedJaccardSimilarityForTreatment(glued_matrix)
    println("Grand Mean: ", Grand_Mean)
    withinTreatment_Mean = mean([treatment1_Mean, treatment2_Mean])
    println("Within Treatment Mean: ", withinTreatment_Mean)
    betweenTreatment_Mean = MeanWeightedJaccardSimilarityBetweenTreatments(treatment1_matrix, treatment2_matrix)
    println("Between Treatment Mean: ", betweenTreatment_Mean)
    ## now, do the randomization test for statistical significance.
    println("Time to run and p-value:")
    @time WeightedJaccardRandomizationTest(treatment1_matrix, treatment2_matrix, nBootstraps)
end

## This is just for checking out the frequency reweighting to 1.0 I applied to pUC tetA-Tn5-KanR mutations.
@rsubset(evolved_mutations, :Treatment == "B30_pUC_50").Gene
@rsubset(evolved_mutations, :Treatment == "B30_pUC_50").Frequency

## Within treatment populations are much more similar than between treatment populations
## (perhaps trivial), but no compelling evidence that plasmid copy number INCREASES parallel evolution.
## There's strongly parallel evolution (evolutionary convergence) in the no plasmid treatments!

## Comparing B30 None to B30 pUC.
WeightedJaccardIndexSpecificityAnalysis(evolved_mutations, "B30_None_50", "B30_pUC_50")
## Comparing B30 None to B30 p15A.
WeightedJaccardIndexSpecificityAnalysis(evolved_mutations, "B30_None_50", "B30_p15A_50")
## Comparing B30 p15A to B30 pUC.
WeightedJaccardIndexSpecificityAnalysis(evolved_mutations, "B30_p15A_50", "B30_pUC_50")

## re-run comparisons, with B20.
## more specificity with greater plasmid copy number with B20.
WeightedJaccardIndexSpecificityAnalysis(evolved_mutations, "B20_None_50", "B20_p15A_50")
WeightedJaccardIndexSpecificityAnalysis(evolved_mutations, "B20_p15A_50", "B20_pUC_50")


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

pUC_treatment_mutations = @rsubset(evolved_mutations, :Plasmid == "pUC")
p15A_treatment_mutations = @rsubset(evolved_mutations, :Plasmid == "p15A")
noPlasmid_treatment_mutations = @rsubset(evolved_mutations, :Plasmid == "None")

pUC_and_bottery_muts = @rsubset(bottery_mutations, :Gene in pUC_treatment_mutations.Gene)
p15A_and_bottery_muts = @rsubset(bottery_mutations, :Gene in p15A_treatment_mutations.Gene)
noPlasmid_and_bottery_muts = @rsubset(bottery_mutations, :Gene in noPlasmid_treatment_mutations.Gene)

## there is much more gene-level parallelism between the no Plasmid treatment
## and the Bottery experiment than between the pUC treatment and the Bottery experiment.

## idea for how to do this analysis, from Deatherage et al. (2017):
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
## in Tet treatment, and some mutations affecting envZ and marR.
card2021_mutations = CSV.read("../data/draft-manuscript-1A/Card2021_Mutations_curated.csv", DataFrame)
relevant_card2021_mutations = filter_Card2021_mutations(card2021_mutations, evolved_mutations)
## show all columns.
show(relevant_card2021_mutations, allcols=true)

####################################
## Lukacisinova 2020 in Nature Communications.
## the main finding here is the lack of parallelism between these data
## and the data in my work. Probably the big difference is that
## there's no tet gene in these strains, so one of the main pathways
## in this paper is gene duplication/amplification of native efflux pumps.

lukacisinova2020_tet_mutations = CSV.read("../data/draft-manuscript-1A/Lukacisinova2020-TET-S1Data.tsv", DataFrame)

lukacisinova_genes = unique(lukacisinova2020_tet_mutations.gene)

lukacisinova_regex = Regex(join(lukacisinova_genes,"|"))
    
gene_in_lukacisinova_gene_column = [occursin(lukacisinova_regex, x) ? true : false for x in evolved_mutations.Gene]

evolved_mutations.matches_lukacisinova = gene_in_lukacisinova_gene_column

## minimal overlap: 11 mutations in 5 matching loci: acrA/acrR, marR, soxR, lon, acrR (omit Tet0 match to argA).
evolved_mutations_in_lukacisinova_genes = @rsubset(evolved_mutations, :matches_lukacisinova)

