## transposon-to-plasmid-metagenomics.R by Rohan Maddamsetti.

## IMPORTANT TODO:
## re-run breseq using the most applicable reference genome listed here:
## https://international.neb.com/tools-and-resources/usage-guidelines/competent-e-coli-genome-sequences-tool

## I believe I am inadvertently using an NEB50alpha F’ (lacIq) reference genome,
## rather than the vanilla NEB 5-alpha reference genome.
## I doubt this will make a substantial difference, but good to get it right.

library(tidyverse)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)

## colorblind-friendly palette.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## assert that we are in the src directory, such that
## projdir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("transposon-plasmid-evolution","src")))
projdir <- file.path("..")


pop.clone.labels <- read.csv(
  file.path(projdir,
            "data/draft-manuscript-1A/populations-and-clones.csv"),
  stringsAsFactors=FALSE)


##############
## TEMPORARY HACK: remove "fixed" mutations with Frequency == 1.0
## that are found in the ancestral DH5a+B20 strain.
ancestral.mut.positions <- c(897327,1105685,1105791,1105812,1105910,1106014,1106021,
                             1106558,1441958,1589337,2140692,2257581,3189165,3804752,4510097,
                             4549970,4550084,
                             ## this next position could be real??? CmR/TcR repeat expansion mutation in all B20 pops.
                             ## let's just remove the CmR/TcR fixed mutations for now.
                             1438,3020,3024,3043,3047)

## best guess is that there are a number of point mutations that need to be added to the B20 transposon reference
## sequence.

##############

## there's the issue of parallelism of high frequency mutations,
## and parallelism of low frequency mutations. Both are important!
## TODO: Dissect the two cases carefully.

## I require a minimum of 4 reads per strand (8 total) to support any variant call.
evolved.mutations <- read.csv(
    file.path(projdir,
              "results/draft-manuscript-1A/genome-analysis/evolved_mutations.csv"),
    stringsAsFactors=FALSE) %>%
    ## CRITICAL TODO: rotate genome coordinates based on oriC in DH5-alpha.
##    mutate(Mbp.coordinate=oriC.coordinate/1000000) %>%
    mutate(Mbp.coordinate=Position/1000000) %>%
    ## TEMPORARY HACK: remove "fixed" mutations with Frequency == 1.0
    ## that are found in the ancestral DH5a+B20 strain.
    filter(!(Position %in% ancestral.mut.positions))

high.freq.evolved.mutations <- evolved.mutations %>%
    filter(Frequency > 0.5)

B30.evolved.mutations <- evolved.mutations %>%
    filter(Transposon == "B30")

B20.evolved.mutations <- evolved.mutations %>%
    filter(Transposon == "B20")

## THIS IS CRITICAL TO LOOK AT-- for checking and filtering ancestral mutations.
fixed.mutations <- evolved.mutations %>%
    filter(Frequency == 1.0) %>%
    arrange(Position,Sample)

###############################################

## Figure 2: Plot the distribution of measured allele frequencies in each population.

Fig2A <- ggplot(B30.evolved.mutations, aes(x=Frequency)) +
    geom_histogram(bins = 100) +
    geom_text_repel(data=filter(B30.evolved.mutations, Frequency>0.2),
                    aes(x=Frequency,y=1,label=Gene),
                    fontface = "italic",size=3,show.legend=FALSE,inherit.aes=FALSE) +
    theme_classic() +
    ylab("Count") +
    xlab("Allele Frequency") +
    facet_grid(Population~Plasmid,scales="free") 

Fig2B <- ggplot(B20.evolved.mutations, aes(x=Frequency)) +
    geom_histogram(bins = 100) +
    geom_text_repel(data=filter(B20.evolved.mutations, Frequency>0.2),
                    aes(x=Frequency,y=1,label=Gene),
                    fontface = "italic",size=3,show.legend=FALSE,inherit.aes=FALSE) +
    theme_classic() +
    ylab("Count") +
    xlab("Allele Frequency") +
    facet_grid(Population~Plasmid,scales="free") 

## This is a very important plot: what does this distribution say about
## the possibility of false positives? how can I interpret this with
## reference to population genetic theory?

## a couple ideas for empirical controls:

## 1) downsample reads from the treatment without plasmid, and re-run breseq
## to see if false positive mutation calls arise when coverage is ~40X rather than
## 300X.

## 2) do pop gen. simulations, and compare allele frequency spectrum.

fig2A.output <- "../results/draft-manuscript-1A/B30-allele-frequency-spectrum.pdf"
fig2B.output <- "../results/draft-manuscript-1A/B20-allele-frequency-spectrum.pdf"
ggsave(Fig2A, file=fig2A.output,width=8,height=6)
ggsave(Fig2B, file=fig2B.output,width=8,height=6)

###############################################
## I also make versions of Figures 3, 4, 5, in which
## the counts are weighted by allele frequency.

## Figure 3: make a stacked bar plot of the kinds of mutations in each clone.


plot.mutations.stackbar <- function(fig.df,panel,leg=FALSE, weight.by.freq=FALSE) {
    if (panel == "No plasmid") {
        muts <- filter(fig.df, Plasmid == 'None')
        my.title <- 'Tet50, No plasmid'
    } else if (panel == "pUC") {
        muts <- filter(fig.df, Plasmid == 'pUC')
        my.title <- 'Tet50, pUC plasmid'
    }  else if (panel == "p15A") {
        muts <- filter(fig.df, Plasmid == 'p15A')
        my.title <- 'Tet50, p15A plasmid'
    }  else {
        stopifnot(TRUE == FALSE) ## panic if we get here.
    }
    
    if (weight.by.freq) {
        fig <- ggplot(muts,aes(x=Population, y=WeightedCount, fill=Mutation)) +
            ylab("Summed Allele Frequency")
    } else {
        fig <- ggplot(muts,aes(x=Population, y=Count, fill=Mutation)) +
            ylab("Count")
    }
    fig <- fig +
        geom_bar(stat='identity') +
        ##ylim(c(0,300)) +
        scale_fill_brewer(palette = "RdYlBu", direction=-1,drop=FALSE) +        
        theme_classic(base_family='Helvetica') +
        theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
              axis.text.y=element_text(size=12),
              panel.border=element_blank(),
              plot.title=element_text(hjust=0, size = 12, face = "bold"),
              strip.text=element_text(hjust=0,size=12),
              panel.spacing.x=unit(1, "cm"),
              panel.spacing.y=unit(0.5, "cm")) +
    ggtitle(my.title)

    if (leg == TRUE) {
        fig <- fig +
            theme(legend.title=element_text(size=12, face="bold"),
                  legend.title.align=0.5,
                  legend.text=element_text(size=12))
    } else {
        fig <- fig + guides(fill = "none")
    }
    
    return(fig)
}

B30.mutation.class.df <- B30.evolved.mutations %>%
    ## give nicer names for mutation classes.
    mutate(Mutation=recode(Mutation,
                  MOB = "Mobile element transposition",
                  DEL = "Indel",
                  INS = "Indel",
                  SUB = "Multiple-base substitution",
                  nonsynonymous = "Nonsynonymous",
                  synonymous = "Synonymous",
                  nonsense = "Nonsense",
                  pseudogene = "Pseudogene",
                  intergenic = "Intergenic",
                  )) %>%
    group_by(Sample, Transposon, Plasmid, Population, Mutation) %>%
    summarize(Count=n(),WeightedCount = sum(Frequency)) %>%
    ungroup() %>%
    data.frame() %>%
    mutate(Mutation=as.factor(as.character(Mutation)))

B20.mutation.class.df <- B20.evolved.mutations %>%
    ## give nicer names for mutation classes.
    mutate(Mutation=recode(Mutation,
                  MOB = "Mobile element transposition",
                  DEL = "Indel",
                  INS = "Indel",
                  SUB = "Multiple-base substitution",
                  nonsynonymous = "Nonsynonymous",
                  synonymous = "Synonymous",
                  nonsense = "Nonsense",
                  pseudogene = "Pseudogene",
                  intergenic = "Intergenic",
                  )) %>%
    group_by(Sample, Transposon, Plasmid, Population, Mutation) %>%
    summarize(Count=n(),WeightedCount = sum(Frequency)) %>%
    ungroup() %>%
    data.frame() %>%
    mutate(Mutation=as.factor(as.character(Mutation)))

## make sure colors are the same across plots by setting levels.
fig3A <- plot.mutations.stackbar(B30.mutation.class.df, panel="No plasmid") + ylim(0,25)
fig3B <- plot.mutations.stackbar(B30.mutation.class.df, panel="pUC",leg=TRUE) + ylim(0,25)
## pull the legend from fig3B.
legend <- cowplot::get_legend(fig3B)
## remove the legend from fig3B.
fig3B <- fig3B + guides(fill = "none")

## IMPORTANT NOTE: This may be an overestimate, since I need to subtract mutations in the ancestor.
fig3C <- plot.mutations.stackbar(B20.mutation.class.df, panel="No plasmid") + ylim(0,25)
fig3D <- plot.mutations.stackbar(B20.mutation.class.df, panel="p15A") + ylim(0,25)
fig3E <- plot.mutations.stackbar(B20.mutation.class.df, panel="pUC") + ylim(0,25)

## plot B30 on the top, B20 on the bottom.
fig3 <- plot_grid(
    fig3A,fig3B,legend,
    fig3C, fig3D, fig3E,              
    labels=c('A','B',"",
             'C', 'D', 'E'),
                  ncol = 3,
                  nrow=2)

fig3.output <- "../results/draft-manuscript-1A/mutation-classes.pdf"
ggsave(fig3, file=fig3.output,width=8,height=8)

## Repeat, but weight by allele frequency.

weighted.fig3A <- plot.mutations.stackbar(B30.mutation.class.df, panel="No plasmid",weight.by.freq=TRUE) + ylim(0,9)
weighted.fig3B <- plot.mutations.stackbar(B30.mutation.class.df, panel="pUC",leg=TRUE,weight.by.freq=TRUE) +
    ylim(0,9)
## remove the legend from weighted.fig3B.
weighted.fig3B <- weighted.fig3B + guides(fill = "none") 

## IMPORTANT NOTE: This may be an overestimate, since I need to subtract mutations in the ancestor.
weighted.fig3C <- plot.mutations.stackbar(B20.mutation.class.df, panel="No plasmid", weight.by.freq=TRUE) + ylim(0,9)
weighted.fig3D <- plot.mutations.stackbar(B20.mutation.class.df, panel="p15A", weight.by.freq=TRUE) + ylim(0,9)
weighted.fig3E <- plot.mutations.stackbar(B20.mutation.class.df, panel="pUC", weight.by.freq=TRUE) + ylim(0,9)


weighted.fig3 <- plot_grid(weighted.fig3A,weighted.fig3B,legend,
                           weighted.fig3C, weighted.fig3D, weighted.fig3E,
                           labels=c('A','B',"",
                                    "C","D","E"),
                           ncol=3,
                           nrow=2)

weighted.fig3.output <- "../results/draft-manuscript-1A/weighted-mutation-classes.pdf"
ggsave(weighted.fig3, file=weighted.fig3.output,width=8,height=7)

    

## let's examine the mutation spectrum, and take a look at the distribution
## of mutations over the genome.
## CRITICAL TODO: find the oriC for the DH5-alpha genome, and rotate coordinates.


## helper to map Allele to class of point mutation.
SNPToSpectrumMap <- function(SNP) {
    if (SNP == 'A->G') {
        return("A:T→G:C")
    } else if (SNP == "T->C") {
        return("A:T→G:C")
    } else if (SNP == "G->A") {
        return("G:C→A:T")
    } else if (SNP == "C->T") {
        return("G:C→A:T")
    } else if (SNP == "A->T") {
        return("A:T→T:A")
    } else if (SNP == "T->A") {
        return("A:T→T:A")
    } else if (SNP == "G->T") {
        return("G:C→T:A")
    } else if (SNP == "C->A") {
        return("G:C→T:A")
    } else if (SNP == "A->C") {
        return("A:T→C:G")
    } else if (SNP == "T->G") {
        return("A:T→C:G")
    } else if (SNP == "G->C") {
        return("G:C→C:G")
    } else if (SNP == "C->G") {
        return("G:C→C:G")
    }
}


SNPSpectrumToClassMap <- function(Spec) {
    if (Spec == 'A:T→G:C') {
        return("Transition")
    } else if (Spec == 'G:C→A:T') {
        return("Transition")
    } else {
        return("Transversion")
    }
}


make.summed.plot <- function(df, number.of.bins = 46) {
    ggplot(df, aes(x=Mbp.coordinate, fill=Spectrum)) +
        geom_histogram(bins = number.of.bins) + 
        theme_classic() +
        ylab("Count") +
        xlab("Genomic position (Mb)") +
        theme(legend.position="bottom")
}


make.facet.mut.plot <- function(df) {
    make.summed.plot(df) + facet_wrap(Plasmid~Transposon,scales="free",nrow=4) 
}


make.spectrum.plot <- function(spectrum.df,weight.by.freq=FALSE) {
    if (weight.by.freq) {
        p <- ggplot(spectrum.df, aes(x=Population, y=WeightedCount, fill=Spectrum)) +
            ylab("Summed Allele Frequency")
    } else {
        p <- ggplot(spectrum.df, aes(x=Population, y=Count, fill=Spectrum)) +
            ylab("Count")
    }
    
    p <- p +
        geom_bar(stat='identity') +
        theme_classic() +
        facet_wrap(Plasmid~Transposon,scales="free") + 
        xlab("Population") +
        theme(legend.position="bottom")
    return(p)
}


make.spectrum.class.plot <- function(spectrum.df,weight.by.freq=FALSE) {

    if (weight.by.freq) {
        p <- ggplot(spectrum.df, aes(x=Population, y=WeightedCount, fill=Spectrum.Class)) +
            ylab("Summed Allele Frequency")
    } else {
        p <- ggplot(spectrum.df, aes(x=Population, y=Count, fill=Spectrum.Class)) +
            ylab("Count") 
    }
    
    p <- p + geom_bar(stat='identity') +
        theme_classic() +
        facet_wrap(Plasmid~Transposon,scales="free") + 
        xlab("Population") +
        theme(legend.position="bottom")
    return(p)
}


## we need a consistent color scale for all 6 classes of mutations in all plots.
## let's use the viridis color scheme.
## https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin

evolved.snps <- evolved.mutations %>%
    filter(str_detect(Mutation_Category, "^snp_")) %>%
    mutate(Spectrum = sapply(Allele, SNPToSpectrumMap)) %>%
    mutate(Spectrum = as.factor(Spectrum)) %>%
    mutate(Spectrum.Class = sapply(Spectrum, SNPSpectrumToClassMap)) %>%
    mutate(Spectrum.Class = as.factor(Spectrum.Class))

spectrum.level.vec <- levels(evolved.snps$Spectrum)
pal <- viridisLite::viridis(length(spectrum.level.vec))
names(pal) <- spectrum.level.vec
COL_SCALE <- scale_fill_manual(name = "Spectrum", values = pal)

Fig4 <- make.facet.mut.plot(evolved.snps) + COL_SCALE
ggsave("../results/draft-manuscript-1A/genomic-mutations.pdf",
       Fig4, width=7, height=5)


point.mut.spectrum.df <- evolved.snps %>%
    group_by(Sample, Transposon, Plasmid, Population, Spectrum) %>%
    summarize(Count=n(), WeightedCount = sum(Frequency)) %>%
    ungroup() %>%
    data.frame()

Fig5A <- make.spectrum.plot(point.mut.spectrum.df) + COL_SCALE + ylim(0, 15)

point.mut.spectrum.class.df <- evolved.snps %>%
    group_by(Sample, Transposon, Plasmid, Population, Spectrum.Class) %>%
    summarize(Count=n(), WeightedCount = sum(Frequency)) %>%
    ungroup() %>%
    data.frame()

Fig5B <- make.spectrum.class.plot(point.mut.spectrum.class.df) + ylim(0, 15)

Fig5 <- plot_grid(Fig5A, Fig5B, labels = c('A','B'),nrow=1)

ggsave("../results/draft-manuscript-1A/mutation-dynamics.pdf",
       Fig5, width=10, height=10)

## repeat, but weigh by allele frequency.

weighted.Fig5A <- make.spectrum.plot(point.mut.spectrum.df,weight.by.freq=TRUE) +
    COL_SCALE + ylim(0,4)

weighted.Fig5B <- make.spectrum.class.plot(point.mut.spectrum.class.df,
                                           weight.by.freq=TRUE) + ylim(0,4)

weighted.Fig5 <- plot_grid(weighted.Fig5A, weighted.Fig5B,
                           labels = c('A','B'),nrow=1) 

ggsave("../results/draft-manuscript-1A/weighted-mutation-dynamics.pdf",
       weighted.Fig5, width=10, height=10)

#####################################################################################

## let's take a close look at the different kinds of evolved mutations.

## TODO: Make a figure to show the granular changes, bp-level parallelism here!

## the frmR/yaiO mutation is a +G mutation, causing a (G)9 repeat to become a (G)10 repeat.
## This is a likely candidate for a hypermutable contingency locus.
evolved.INDEL <- evolved.mutations %>% filter(Mutation == "INS" | Mutation == "DEL")

## ALL of these MOB insertiosn are miniTn5-Tet insertions, either into the KanR gene
## on the plasmid, or into chromosomal genes in the no plasmid treatment.
evolved.MOB <- evolved.mutations %>% filter(Mutation == "MOB")

## very strong evidence of parallel evolution in three intergenic regions:
## lysO/aqpZ, yeaD/yeaE, and yohP/dusC.
## These show extremely strong parallel evolution at the base-pair level as well.

## TODO: search for unannotated sRNAs and very short proteins in the following
## regions: lysO/aqpZ, yohP/dusC, and for 3'-UTR sRNA just downstream of both
## yeaD/yeaE.

## As a quick check, I BLASTed these regions in the CP053607.1 reference genome
## again the well-annotated NC_000913.3 K-12 MG1655 reference genome,
## since there are often noncoding RNAs annotated in K-12 that are not annotated
## in other E. coli reference genomes. I did not find anything this way.

evolved.intergenic <- evolved.mutations %>% filter(Mutation == "intergenic") %>%
    arrange(Gene, Position, Plasmid, Sample) %>%
    select(-Transposon, -Mutation, -Mutation_Category, -Population)

## extremely strong parallel evolution in tetA-- seen only in no plasmid treatment.
## the parallel evolution is at the nucleotide level in the no plasmid treatment:
## (3100, 3102) in the tetA promoter, (3134 has two independent dS mutations, 3135
## has a nonsynonymous mutation,), and (4285, 4285, 4286 are independent nonsense mutations
## that truncate the very end of the protein.)
evolved.tetA <- evolved.mutations %>% filter(str_detect(Gene, "tetA")) %>%
    arrange(Sample, Position)

## most synonymous mutations in pops 3,4,5  of the pUC plasmid treatment.
## This will certainly be a significant association. perhaps some cryptic hypermutator
## clades here?
evolved.synonymous <- evolved.mutations %>% filter(Mutation == "synonymous") %>%
    arrange(Gene, Position, Sample)

## parallel evolution of robA in no plasmid treatment.
evolved.nonsynonymous <- evolved.mutations %>% filter(Mutation == "nonsynonymous") %>%
    arrange(Gene, Position, Sample)


#####################################################################################
## examine DNA repair and DNA polymerase/replication genes for mutator and anti-mutator
## candidates.

DNA.repair.loci <- read.csv("../data/draft-manuscript-1A/DNA-repair-and-replication.csv",header=TRUE,as.is=TRUE)

DNA.repair.muts <- filter(evolved.mutations, Gene %in% unique(DNA.repair.loci$Gene))

## no parallelism at nucleotide level.
DNA.repair.parallel.nuc <- DNA.repair.muts %>%
    group_by(Gene, Position) %>% summarize(count=n()) %>%
    arrange(desc(count)) %>% filter(count>1)

## no evidence of parallel evolution.
DNA.repair.parallel.gene <- DNA.repair.muts %>%
    group_by(Gene) %>% summarize(count=n()) %>%
    arrange(desc(count)) %>% filter(count>1)

DNA.repair.parallel.loci.muts <- DNA.repair.muts %>%
    filter(Gene %in% DNA.repair.parallel.gene$Gene)

#################################################################################
## analysis of parallel evolution at the same nucleotide.
## discuss numbers and finding in the text (no figure.).
## This could be a Supplementary Table.

bp.parallel.mutations <- evolved.mutations %>% group_by(Position) %>%
summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.MOB <- filter(bp.parallel.mutations,Mutation=='MOB')
## no parallel DEL muts at bp level.
parallel.INS <- filter(bp.parallel.mutations,Mutation=='INS')
parallel.dN <- filter(bp.parallel.mutations,Mutation=='nonsynonymous')
parallel.dS <- filter(bp.parallel.mutations,Mutation=='synonymous')

## examine parallel evolution at amino acid level (only one case, in robA).
parallel.AA.dN <- evolved.mutations %>% filter(Mutation=='nonsynonymous') %>% group_by(Position) %>% summarize(count=n()) %>% filter(count > 1)
parallel.dN.Table <- filter(evolved.mutations, Position %in% parallel.AA.dN$Position) %>% arrange(Position)


##################################################################################
## analysis of parallel evolution at the gene level (including intergenic regions).

gene.level.parallel.mutations <- evolved.mutations %>% group_by(Gene) %>%
summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.genes <- gene.level.parallel.mutations %>%
    select(Gene, count, Plasmid, Transposon) %>%
    distinct() %>%
    arrange(desc(count))

################################################################################
### Figure 6: make a matrix plot of genes with mutations in two or more clones.
################################################################################

MakeMutCountMatrixFigure <- function(evolved.mutations, pop.clone.labels, show.all=FALSE) {

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.mutations %>%
        select(Gene, Sample) %>%
        group_by(Sample, Gene) %>%
        summarize(mutation.count = n()) %>%
        left_join(pop.clone.labels) %>%
        unite("Treatment", Transposon:Plasmid, remove = FALSE) %>%
        select(Gene, Sample, mutation.count, Transposon, Plasmid, Treatment)
    
    total.muts <- matrix.data %>%
        group_by(Gene) %>%
        summarize(total.mutation.count = sum(mutation.count))
    
    matrix.data <- left_join(matrix.data, total.muts)

    if (!show.all) ## then filter out genes that are only hit in one sample.
        matrix.data <- matrix.data %>% filter(total.mutation.count > 1)
    
    ## sort genes by number of mutations in each row.
    ## this is overwritten by the alternate sorting method that follows,
    ## but keeping this snippet because it's useful to have.
    gene.hit.sort <- matrix.data %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(hits=sum(mutation.count)) %>%
        arrange(desc(hits))
    matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.hit.sort$Gene))
    
    ## alternate sorting method: difference in hits between environments.
    pUC.mut.count <- filter(matrix.data,Plasmid=="pUC") %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(pUC.mut.count=sum(mutation.count))
    
    noPlasmid.mut.count <- filter(matrix.data,Plasmid=="None") %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(noPlasmid.mut.count=sum(mutation.count))
    
    treatment.mut.sort <- inner_join(pUC.mut.count, noPlasmid.mut.count) %>%
        mutate(mut.diff = abs(pUC.mut.count - noPlasmid.mut.count)) %>%
        arrange(desc(mut.diff))
    
    ## now use these calculations to sort genes by the absolute value of the
    ## difference in number of mutations between the pUC and No Plasmid  treatments.
    matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(treatment.mut.sort$Gene))
    
    ## cast mutation.count into a factor for plotting.
    matrix.data$mutation.count <- factor(matrix.data$mutation.count)
    
    make.matrix.panel <- function(mdata, treatment, leg=FALSE) {
        fig <- ggplot(filter(mdata,Treatment==treatment),
                      aes(x=Sample,
                          y=Gene,
                          fill=mutation.count,
                          frame=Treatment)
                      ) +
            geom_tile(color="black",size=0.1) +
            ylab("Gene") +
            xlab("Metagenome") +
            ggtitle(paste(treatment,'metagenomes')) +
            theme_tufte(base_family='Helvetica') +
            theme(axis.ticks = element_blank(),
                  axis.text.x = element_text(size=10,angle=45,hjust=1),
                  axis.text.y = element_text(size=10,hjust=1,face="italic"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  ) +
            scale_y_discrete(drop=FALSE) + ## don't drop missing genes.
            scale_fill_manual(name="Mutations",
                              values = c("#ffdf00", "#bebada", ##c("white", "#ffdf00", "#bebada",
                                         "#fb8072", "#80b1d3", "#fdb462"))
        
        if (leg == FALSE) {
            fig <- fig + guides(fill= "none")
        }
        return(fig)
    }

    B20.noPlasmid.matrix.panel <- make.matrix.panel(matrix.data,"B20_None")
    B20.A31.matrix.panel <- make.matrix.panel(matrix.data,"B20_p15A")
    B20.A18.matrix.panel <- make.matrix.panel(matrix.data,"B20_pUC")

    B30.noPlasmid.matrix.panel <- make.matrix.panel(matrix.data,"B30_None")
    B30.A18.matrix.panel <- make.matrix.panel(matrix.data,"B30_pUC")
    
    matrix.figure <- plot_grid(
        B20.noPlasmid.matrix.panel,
        B20.A31.matrix.panel,    
        B20.A18.matrix.panel,
        B30.noPlasmid.matrix.panel,
        B30.A18.matrix.panel,
        nrow=1,
        align = 'vh')
    return(matrix.figure)
}

Fig6 <- MakeMutCountMatrixFigure(evolved.mutations, pop.clone.labels)
matrix.outf <- "../results/draft-manuscript-1A/Fig6.pdf"
ggsave(matrix.outf, Fig6,width=20)

Fig6singles <- MakeMutCountMatrixFigure(evolved.mutations, pop.clone.labels, show.all=T)
ggsave("../results/draft-manuscript-1A/Fig6-singles.pdf", Fig6singles,width=20, height=20)


## Redo, but used summed allele frequency for the heatmap.
MakeSummedAlleleFrequencyMatrixFigure <- function(evolved.mutations, pop.clone.labels,
                                                  allele.freq.threshold = 0.2,
                                                  show.all=FALSE) {

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.mutations %>%
        select(Gene, Sample, Frequency) %>%
        group_by(Sample, Gene) %>%
        summarize(mutation.count = n(), summed.Allele.Frequency = sum(Frequency)) %>%
        left_join(pop.clone.labels) %>%
        unite("Treatment", Transposon:Plasmid, remove = FALSE) %>%
        select(Gene, Sample, mutation.count, summed.Allele.Frequency,
               Transposon, Plasmid, Treatment) %>%
        ungroup()
    
    total.allele.freqs <- matrix.data %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(total.Allele.Frequency = sum(summed.Allele.Frequency))

    ## filter matrix.data for genes that pass the allele frequency threshold,
    ## based on total allele frequency summed across all pops.
    if (!show.all) 
        matrix.data <- left_join(matrix.data, total.allele.freqs) %>%
            filter(total.Allele.Frequency > allele.freq.threshold)
        
    ## sorting method: difference in allele frequency between environments.
    pUC.allele.freq <- filter(matrix.data,Plasmid=="pUC") %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(pUC.allele.frequency=sum(summed.Allele.Frequency))
    
    noPlasmid.allele.freq <- filter(matrix.data,Plasmid=="None") %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(noPlasmid.allele.frequency=sum(summed.Allele.Frequency))
    
    treatment.freq.sort <- full_join(pUC.allele.freq, noPlasmid.allele.freq) %>%
        replace_na(list(pUC.allele.frequency = 0, noPlasmid.allele.frequency = 0)) %>%
        mutate(allele.diff = abs(pUC.allele.frequency - noPlasmid.allele.frequency)) %>%
        arrange(desc(allele.diff)) %>%
        mutate(Gene = factor(Gene))
    
    ## now use these calculations to sort genes by the absolute value of the
    ## difference in allele frequency between the pUC and No Plasmid treatments.
    matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(treatment.freq.sort$Gene))

    make.allele.freq.matrix.panel <- function(mdata, treatment, leg=FALSE) {
        fig <- ggplot(filter(mdata,Treatment==treatment),
                      aes(x=Sample,
                          y=Gene,
                          fill=summed.Allele.Frequency,
                          frame=Treatment)
                      ) +
            geom_tile(color="black",size=0.1) +
            ylab("Gene") +
            xlab("Metagenome") +
            ggtitle(paste(treatment,'metagenomes')) +
            theme_tufte(base_family='Helvetica') +
            theme(axis.ticks = element_blank(),
                  axis.text.x = element_text(size=10,angle=45,hjust=1),
                  axis.text.y = element_text(size=10,hjust=1,face="italic"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  ) +
            scale_y_discrete(drop=FALSE) + ## don't drop missing genes.
            scale_fill_viridis_c(option = "inferno")


        
        if (leg == FALSE) {
            fig <- fig + guides(fill= "none")
        }
        return(fig)
    }
    

    B20.noPlasmid.matrix.panel <- make.allele.freq.matrix.panel(matrix.data,"B20_None")
    B20.A31.matrix.panel <- make.allele.freq.matrix.panel(matrix.data,"B20_p15A")
    B20.A18.matrix.panel <- make.allele.freq.matrix.panel(matrix.data,"B20_pUC")

    
    B30.noPlasmid.matrix.panel <- make.allele.freq.matrix.panel(matrix.data,"B30_None")
    B30.A18.matrix.panel <- make.allele.freq.matrix.panel(matrix.data,"B30_pUC")

    matrix.figure <- plot_grid(
        B20.noPlasmid.matrix.panel,
        B20.A31.matrix.panel,    
        B20.A18.matrix.panel,
        B30.noPlasmid.matrix.panel,
        B30.A18.matrix.panel,
        nrow=1,
        align = 'vh')
    return(matrix.figure)

}

Fig7 <- MakeSummedAlleleFrequencyMatrixFigure(evolved.mutations, pop.clone.labels)
allele.freq.matrix.outf <- "../results/draft-manuscript-1A/Fig7.pdf"
ggsave(allele.freq.matrix.outf, Fig7, width=20)


Fig7singles <- MakeSummedAlleleFrequencyMatrixFigure(evolved.mutations, pop.clone.labels, show.all=T)
ggsave("../results/draft-manuscript-1A/Fig7-singles.pdf", Fig7singles,height=20, width=20)

###########################################################

## Evolutionary rate comparisons across treatments,
## following Deatherage et al. (2017)
## and Blount, Maddamsetti, Grant et al. (2020).

## This section will be expanded, as more genome data comes in.

## Sum the allele frequencies and number of mutations across genes in each sample.
total.allele.frequency.summary <- evolved.mutations %>%
    group_by(Sample, Transposon, Plasmid, Population) %>%
    summarize(mutation.count = n(), summed.Allele.Frequency = sum(Frequency)) %>%
    left_join(pop.clone.labels) %>%
    unite("Treatment", Transposon:Plasmid, remove = FALSE) %>%
    select(Sample, mutation.count, summed.Allele.Frequency,
           Transposon, Plasmid, Treatment) %>%
    ungroup()

pUC.treatment.total.allele.frequencies <- total.allele.frequency.summary %>%
    filter(Plasmid == "pUC")

noPlasmid.treatment.total.allele.frequencies <- total.allele.frequency.summary %>%
    filter(Plasmid == "None")

## two-tailed Mann-Whitney U-test: no difference in numbers of mutations
## or summed allele frequencies between treatments. BUT note that
## RM6-176-16 and RM6-176-17 samples have more mutations, probably
## due to some cryptic hypermutator phenotype.
wilcox.test(pUC.treatment.total.allele.frequencies$mutation.count,
            noPlasmid.treatment.total.allele.frequencies$mutation.count)

wilcox.test(pUC.treatment.total.allele.frequencies$summed.Allele.Frequency,
            noPlasmid.treatment.total.allele.frequencies$summed.Allele.Frequency)

## we can also use the Kruskal-Wallis test to show no difference between treatments
## in number of mutations or total allele frequency.
## NOTE: the df is different between these two calls-- might this be a bug?
kruskal.test(pUC.treatment.total.allele.frequencies$mutation.count,
             noPlasmid.treatment.total.allele.frequencies$mutation.count)

kruskal.test(pUC.treatment.total.allele.frequencies$summed.Allele.Frequency,
             noPlasmid.treatment.total.allele.frequencies$summed.Allele.Frequency)
