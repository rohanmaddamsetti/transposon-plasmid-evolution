## transposon-to-plasmid-metagenomics.R by Rohan Maddamsetti.

## IMPORTANT TODO:
## re-run breseq using the most applicable reference genome listed here:
## https://international.neb.com/tools-and-resources/usage-guidelines/competent-e-coli-genome-sequences-tool

## I believe I am inadvertently using an NEB50alpha F’ (lacIq) reference genome,
## rather than the vanilla NEB 5-alpha reference genome.
## I doubt this will make a substantial difference, but good to get it right.

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggthemes)
library(viridis)
library(ggrepel)

## colorblind-friendly palette.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## assert that we are in the src directory, such that
## projdir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("transposon-plasmid-evolution","src")))
projdir <- file.path("..")

## This file *just* has the evolved populations (Day 9).
pop.clone.labels <- read.csv(
  file.path(projdir,
            "data/draft-manuscript-1A/evolved-populations-and-clones.csv"),
  stringsAsFactors=FALSE)


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
    mutate(Mbp.coordinate=Position/1000000)


high.freq.evolved.mutations <- evolved.mutations %>%
    filter(Frequency > 0.5)

B30.evolved.mutations <- evolved.mutations %>%
    filter(Transposon == "B30")

B20.evolved.mutations <- evolved.mutations %>%
    filter(Transposon == "B20")

B30.Tet0.evolved.mutations <- B30.evolved.mutations %>%
    filter(Tet==0)

B30.Tet50.evolved.mutations <- B30.evolved.mutations %>%
    filter(Tet==50)

B20.Tet0.evolved.mutations <- B20.evolved.mutations %>%
    filter(Tet==0)

B20.Tet50.evolved.mutations <- B20.evolved.mutations %>%
    filter(Tet==50)


## THIS IS CRITICAL TO LOOK AT-- for checking and filtering ancestral mutations.
fixed.mutations <- evolved.mutations %>%
    filter(Frequency == 1.0) %>%
    arrange(Position,Sample)

###############################################

## Figure 2: Plot the distribution of measured allele frequencies in each population.

make.allele.freq.histogram <- function(evolved.mutations.df, my.title,annotate=FALSE) {
    p <- ggplot(evolved.mutations.df, aes(x=Frequency)) +
        geom_histogram(bins = 100) +
        theme_classic() +
        ylab("Count") +
        xlab("Allele Frequency") +
        scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1), limits = c(0,1.1)) +
        ggtitle(my.title) +
        facet_grid(Plasmid~.) +
    geom_vline(xintercept=0.25,color="red",linetype="dashed",size=0.2)

    muts.to.label <- filter(evolved.mutations.df, Frequency>0.25)
    if (annotate && nrow(muts.to.label) > 0) {
        p <- p +
            geom_text_repel(
                ## label mutations with > 25% allele frequency
                data= muts.to.label,
                aes(x=Frequency,y=1,label=Gene),
                fontface = "italic",size=1.5,show.legend=FALSE,inherit.aes=FALSE)
        }
    return(p)
}

Fig2A <- make.allele.freq.histogram(B30.Tet0.evolved.mutations, "B30, Tet 0 populations",TRUE)
Fig2B <- make.allele.freq.histogram(B30.Tet50.evolved.mutations, "B30, Tet 50 populations",TRUE)
Fig2C <- make.allele.freq.histogram(B20.Tet0.evolved.mutations, "B20, Tet 0 populations",TRUE)
Fig2D <- make.allele.freq.histogram(B20.Tet50.evolved.mutations, "B20, Tet 50 populations",TRUE)
## This is a very important plot: what does this distribution say about
## the possibility of false positives? how can I interpret this with
## reference to population genetic theory?

## a couple ideas for empirical controls:

## 1) downsample reads from the treatment without plasmid, and re-run breseq
## to see if false positive mutation calls arise when coverage is ~40X rather than
## 300X.

## 2) do pop gen. simulations, and compare allele frequency spectrum.

Fig2 <- Fig2A + Fig2B + Fig2C + Fig2D
fig2.output <- "../results/draft-manuscript-1A/Fig2.pdf"
ggsave(Fig2, file=fig2.output,width=10,height=8)

###############################################
## I also make versions of Figures 3, 4, 5, in which
## the counts are weighted by allele frequency.
## Figure 3: make a stacked bar plot of the kinds of mutations in each clone.

plot.mutations.stackbar <- function(fig.df, my.title, leg=FALSE, weight.by.freq=FALSE) {

    if (str_detect(my.title, "no plasmid")) {
        muts <- filter(fig.df, Plasmid == 'None')
    } else if (str_detect(my.title, "pUC")) {
        muts <- filter(fig.df, Plasmid == 'pUC')
    }  else if (str_detect(my.title, "p15A")) {
        muts <- filter(fig.df, Plasmid == 'p15A')
    } else {
        stopifnot(TRUE == FALSE) ## panic if we get here.
    }

    if (weight.by.freq) {
        fig <- ggplot(muts,aes(x=Population, y=WeightedCount, fill=Mutation)) +
            ylim(0,7.5) +
            ylab("Summed Allele Frequency")
    } else {
        fig <- ggplot(muts,aes(x=Population, y=Count, fill=Mutation)) +
            ylim(0,30) +
            ylab("Count")
    }
    fig <- fig +
        geom_bar(stat='identity') +
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
                  legend.text=element_text(size=12),
                  legend.position="bottom")
    } else {
        fig <- fig + guides(fill = "none")
    }
    
    return(fig)
}

make.mutation.class.df <- function(evolved.mutations.df) {
    evolved.mutations.df %>%
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
}

make.Fig3 <- function(evolved.mutations, weight.by.freq) {
    B30.evolved.mutations <- evolved.mutations %>%
        filter(Transposon == "B30")
    B20.evolved.mutations <- evolved.mutations %>%
        filter(Transposon == "B20")
    B30.Tet0.evolved.mutations <- B30.evolved.mutations %>%
        filter(Tet==0)
    B30.Tet50.evolved.mutations <- B30.evolved.mutations %>%
        filter(Tet==50)
    B20.Tet0.evolved.mutations <- B20.evolved.mutations %>%
        filter(Tet==0)
    B20.Tet50.evolved.mutations <- B20.evolved.mutations %>%
        filter(Tet==50)
    
    B30.Tet50.mutation.class.df <- make.mutation.class.df(B30.Tet50.evolved.mutations)
    B20.Tet50.mutation.class.df <- make.mutation.class.df(B20.Tet50.evolved.mutations)
    B30.Tet0.mutation.class.df <- make.mutation.class.df(B30.Tet0.evolved.mutations)
    B20.Tet0.mutation.class.df <- make.mutation.class.df(B20.Tet0.evolved.mutations)
    
    ## make sure colors are the same across plots by setting levels.
    fig3A <- plot.mutations.stackbar(B30.Tet50.mutation.class.df, my.title="B30 no plasmid\nTet 50", FALSE, weight.by.freq)
    fig3B <- plot.mutations.stackbar(B30.Tet50.mutation.class.df, "B30 p15A\nTet 50", FALSE, weight.by.freq) 
    fig3C <- plot.mutations.stackbar(B30.Tet50.mutation.class.df, "B30 pUC\nTet 50", TRUE, weight.by.freq) 
    
    ## pull the legend from fig3C.
    fig3legend <- cowplot::get_legend(fig3C)
    ## remove the legend from fig3C.
    fig3C <- fig3C + guides(fill = "none")
    
    fig3D <- plot.mutations.stackbar(B20.Tet50.mutation.class.df, "B20 no plasmid\nTet 50", FALSE, weight.by.freq) 
    fig3E <- plot.mutations.stackbar(B20.Tet50.mutation.class.df, "B20 p15A\nTet 50", FALSE, weight.by.freq) 
    fig3F <- plot.mutations.stackbar(B20.Tet50.mutation.class.df, "B20 pUC\nTet 50", FALSE, weight.by.freq) 
    
    ## add plots for Tet 0 control populations. 
    fig3G <- plot.mutations.stackbar(B30.Tet0.mutation.class.df, "B30 no plasmid\nTet 0", FALSE, weight.by.freq) 
    fig3H <- plot.mutations.stackbar(B30.Tet0.mutation.class.df, "B30 p15A\nTet 0", FALSE, weight.by.freq) 
    fig3I <- plot.mutations.stackbar(B30.Tet0.mutation.class.df, "B30 pUC\nTet 0", FALSE, weight.by.freq) 
    
    fig3J <- plot.mutations.stackbar(B20.Tet0.mutation.class.df, "B20 no plasmid\nTet 0", FALSE, weight.by.freq) 
    fig3K <- plot.mutations.stackbar(B20.Tet0.mutation.class.df, "B20 p15A\nTet 0", FALSE, weight.by.freq) 
    fig3L <- plot.mutations.stackbar(B20.Tet0.mutation.class.df, "B20 pUC\nTet 0", FALSE, weight.by.freq)
    
    fig3 <- plot_grid(
        plot_grid(
            fig3A, fig3B, fig3C,
            fig3D, fig3E, fig3F,
            fig3G, fig3H, fig3I,
            fig3J, fig3K, fig3L,
            ncol = 3,
            nrow=4),
        fig3legend, ncol = 1, rel_heights = c(1,0.1))
    return(fig3)
}

fig3.output <- "../results/draft-manuscript-1A/Fig3.pdf")
fig3 <- make.Fig3(evolved.mutations, FALSE)
ggsave(fig3, file=fig3.output,width=8,height=8)

## Repeat, but weight by allele frequency.
weighted.fig3.output <- "../results/draft-manuscript-1A/weighted-Fig3.pdf"
weightedfig3 <- make.Fig3(evolved.mutations, TRUE)
ggsave(weightedfig3, file=weighted.fig3.output,width=8,height=8)

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
        theme(legend.position="bottom") +
        ylim(0,40)
}


make.facet.mut.plot <- function(df) {
    make.summed.plot(df) + facet_wrap(Plasmid~Tet,scales="free",nrow=4) 
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
        facet_grid(Plasmid~Tet) + 
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
        facet_grid(Plasmid~Tet,scales="free") + 
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
ggsave("../results/draft-manuscript-1A/Fig4.pdf",
       Fig4, width=7, height=5)


point.mut.spectrum.df <- evolved.snps %>%
    group_by(Sample, Transposon, Plasmid, Tet, Population, Spectrum) %>%
    summarize(Count=n(), WeightedCount = sum(Frequency)) %>%
    ungroup() %>%
    data.frame()

Fig5A <- make.spectrum.plot(point.mut.spectrum.df) + COL_SCALE + ylim(0, 40)

point.mut.spectrum.class.df <- evolved.snps %>%
    group_by(Sample, Transposon, Plasmid, Tet, Population, Spectrum.Class) %>%
    summarize(Count=n(), WeightedCount = sum(Frequency)) %>%
    ungroup() %>%
    data.frame()

Fig5B <- make.spectrum.class.plot(point.mut.spectrum.class.df) + ylim(0, 40)

Fig5 <- plot_grid(Fig5A, Fig5B, labels = c('A','B'),nrow=1)

ggsave("../results/draft-manuscript-1A/Fig5.pdf",
       Fig5, width=10, height=10)

## repeat, but weigh by allele frequency.
weighted.Fig5A <- make.spectrum.plot(point.mut.spectrum.df,weight.by.freq=TRUE) +
    COL_SCALE + ylim(0,4)

weighted.Fig5B <- make.spectrum.class.plot(point.mut.spectrum.class.df,
                                           weight.by.freq=TRUE) + ylim(0,4)

weighted.Fig5 <- plot_grid(weighted.Fig5A, weighted.Fig5B,
                           labels = c('A','B'),nrow=1) 

ggsave("../results/draft-manuscript-1A/weighted-Fig5.pdf",
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

MakeMutCountMatrixFigure <- function(evolved.mutations, show.all=FALSE) {

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.mutations %>%
        ## unite the Transposon, Plasmid, Tet columns together.
        unite("Treatment", Transposon:Tet, sep="\n", remove = FALSE) %>%
        group_by(Gene, Sample, Transposon, Plasmid, Tet, Treatment) %>%
        summarize(mutation.count = n())
    
    total.muts <- matrix.data %>%
        group_by(Gene) %>%
        summarize(total.mutation.count = sum(mutation.count))
    
    matrix.data <- left_join(matrix.data, total.muts)

    if (!show.all) { ## then filter out genes that are only hit in one sample.
        matrix.data <- matrix.data %>% filter(total.mutation.count > 1)
    }
    
    ## sort genes by number of mutations in each row.
    ## also check out the alternate sorting method that follows.
    gene.hit.sort <- matrix.data %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(hits=sum(mutation.count)) %>%
        arrange(desc(hits))
    
    ## alternate sorting method: difference in hits between environments,
    ## AKA the (absolute value of the) difference in number of pops with hits
    ## between the Tet 50 and Tet 0 treatments.
    Tet50.hit.count.df <- filter(matrix.data,Tet==50) %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(Tet50.hit.count=n())
    
    Tet0.hit.count.df <- filter(matrix.data,Tet==0) %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(Tet0.hit.count=n())
    
    treatment.hit.sort <- full_join(Tet50.hit.count.df, Tet0.hit.count.df) %>%
        mutate(hit.diff = Tet50.hit.count - Tet0.hit.count) %>%
        arrange(desc(hit.diff))
    
    ## now sort genes.
    use.treatment.hit.sort <- TRUE
    if (use.treatment.hit.sort) {
        matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(treatment.hit.sort$Gene))
    } else {
        matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.hit.sort$Gene))
    }
    
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
            ggtitle(treatment) +
            theme_tufte(base_family='Helvetica') +
            theme(axis.ticks = element_blank(),
                  axis.text.x = element_text(size=10,angle=45,hjust=1),
                  axis.text.y = element_text(size=10,hjust=1,face="italic"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  ) +
            scale_y_discrete(drop=FALSE) + ## don't drop missing genes.
            scale_fill_manual(name="Mutations",
                              values = c("#ffdf00", "#bebada", "#fb8072", "#80b1d3", "#fdb462"))
        
        if (leg == FALSE) {
            fig <- fig + guides(fill= "none")
        }
        return(fig)
    }

    B20.noPlasmid.Tet50.matrix.panel <- make.matrix.panel(matrix.data, "B20\nNone\n50")
    ## Remove the gene labels to save space.
    B20.A31.Tet50.matrix.panel <- make.matrix.panel(matrix.data, "B20\np15A\n50") +
        theme(axis.text.y=element_blank())
    B20.A18.Tet50.matrix.panel <- make.matrix.panel(matrix.data, "B20\npUC\n50") +
        theme(axis.text.y=element_blank())

    B30.noPlasmid.Tet50.matrix.panel <- make.matrix.panel(matrix.data,"B30\nNone\n50") +
        theme(axis.text.y=element_blank())
    B30.A31.Tet50.matrix.panel <- make.matrix.panel(matrix.data, "B30\np15A\n50") +
        theme(axis.text.y=element_blank())
    B30.A18.Tet50.matrix.panel <- make.matrix.panel(matrix.data,"B30\npUC\n50") +
        theme(axis.text.y=element_blank())

    ## Tet 0 panels.
    B20.noPlasmid.Tet0.matrix.panel <- make.matrix.panel(matrix.data, "B20\nNone\n0") +
        theme(axis.text.y=element_blank())
    B20.A31.Tet0.matrix.panel <- make.matrix.panel(matrix.data, "B20\np15A\n0") +
        theme(axis.text.y=element_blank())
    B20.A18.Tet0.matrix.panel <- make.matrix.panel(matrix.data, "B20\npUC\n0") +
        theme(axis.text.y=element_blank())
    B30.noPlasmid.Tet0.matrix.panel <- make.matrix.panel(matrix.data,"B30\nNone\n0") +
        theme(axis.text.y=element_blank())
    B30.A31.Tet0.matrix.panel <- make.matrix.panel(matrix.data, "B30\np15A\n0") +
        theme(axis.text.y=element_blank())
    B30.A18.Tet0.matrix.panel <- make.matrix.panel(matrix.data,"B30\npUC\n0") +
        theme(axis.text.y=element_blank())

    ## Using the patchwork library for layout.
    matrix.figure <-
        B20.noPlasmid.Tet50.matrix.panel +
        B30.noPlasmid.Tet50.matrix.panel +
        B20.A31.Tet50.matrix.panel +
        B30.A31.Tet50.matrix.panel +
        B20.A18.Tet50.matrix.panel +
        B30.A18.Tet50.matrix.panel +
        B20.noPlasmid.Tet0.matrix.panel +
        B30.noPlasmid.Tet0.matrix.panel +
        B20.A31.Tet0.matrix.panel +
        B30.A31.Tet0.matrix.panel +
        B20.A18.Tet0.matrix.panel +
        B30.A18.Tet0.matrix.panel +
        plot_layout(nrow = 1)
    
    return(matrix.figure)
}

## Use summed allele frequency for the heatmap.
MakeSummedAlleleFrequencyMatrixFigure <- function(evolved.mutations,
                                                  allele.freq.threshold = 0.2,
                                                  show.all=FALSE) {

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.mutations %>%
        ## unite the Transposon, Plasmid, Tet columns together.
        unite("Treatment", Transposon:Tet, sep="\n", remove = FALSE) %>%
        group_by(Gene, Sample, Transposon, Plasmid, Tet, Treatment) %>%
        summarize(summed.Allele.Frequency = sum(Frequency))
    
    total.allele.freqs <- matrix.data %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(total.Allele.Frequency = sum(summed.Allele.Frequency))

    ## filter matrix.data for genes that pass the allele frequency threshold,
    ## based on total allele frequency summed across all pops.
    if (!show.all) {
        matrix.data <- left_join(matrix.data, total.allele.freqs) %>%
            filter(total.Allele.Frequency > allele.freq.threshold)
    }

    ## sort genes by the total allele frequency  in each row.
    ## also check out the alternate sorting method that follows.
    gene.freq.sort <- matrix.data %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(totalallelefreq = sum(summed.Allele.Frequency)) %>%
        arrange(desc(totalallelefreq))
    
    ## alternative sorting method:
    ## difference in allele frequency between the Tet 50 and Tet 0 treatments..
    Tet50.allele.freq.df <- filter(matrix.data, Tet==50) %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(Tet50.allele.frequency=sum(summed.Allele.Frequency))
    
    Tet0.allele.freq.df <- filter(matrix.data, Tet==0) %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(Tet0.allele.frequency = sum(summed.Allele.Frequency))
    
    treatment.freq.sort <- full_join(Tet50.allele.freq.df, Tet0.allele.freq.df) %>%
        replace_na(list(Tet50.allele.frequency = 0, Tet0.allele.frequency = 0)) %>%
        mutate(allele.diff = Tet50.allele.frequency - Tet0.allele.frequency) %>%
        arrange(desc(allele.diff))

    ## sort the genes.
    use.treatment.mut.sort <- TRUE
    if (use.treatment.mut.sort) {
        matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(treatment.freq.sort$Gene))
    } else {
        matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.freq.sort$Gene))
    }

    make.allele.freq.matrix.panel <- function(mdata, treatment, leg=FALSE) {
        fig <- ggplot(filter(mdata,Treatment==treatment),
                      aes(x=Sample,
                          y=Gene,
                          fill=summed.Allele.Frequency,
                          frame=Treatment)
                      ) +
            geom_tile(color="black",size=0.1) +
            ggtitle(treatment) +
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
    

    B20.noPlasmid.Tet50.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "B20\nNone\n50")
    ## Remove the gene labels for the additional matrices to save space.
    B20.A31.Tet50.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "B20\np15A\n50")  +
        theme(axis.text.y=element_blank())
    B20.A18.Tet50.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "B20\npUC\n50")  +
        theme(axis.text.y=element_blank())
    
    B30.noPlasmid.Tet50.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "B30\nNone\n50")  +
        theme(axis.text.y=element_blank())
    B30.A31.Tet50.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "B30\np15A\n50")  +
        theme(axis.text.y=element_blank())
    B30.A18.Tet50.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "B30\npUC\n50")  +
        theme(axis.text.y=element_blank())

    ## Tet 0 panels.
    B20.noPlasmid.Tet0.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "B20\nNone\n0") +
        theme(axis.text.y=element_blank())
    B20.A31.Tet0.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "B20\np15A\n0") +
        theme(axis.text.y=element_blank())
    B20.A18.Tet0.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "B20\npUC\n0") +
        theme(axis.text.y=element_blank())
    B30.noPlasmid.Tet0.matrix.panel <- make.allele.freq.matrix.panel(matrix.data,"B30\nNone\n0") +
        theme(axis.text.y=element_blank())
    B30.A31.Tet0.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "B30\np15A\n0") +
        theme(axis.text.y=element_blank())
    B30.A18.Tet0.matrix.panel <- make.allele.freq.matrix.panel(matrix.data,"B30\npUC\n0") +
        theme(axis.text.y=element_blank())

    ## Using the patchwork library for layout.
    matrix.figure <-
        B20.noPlasmid.Tet50.matrix.panel +
        B30.noPlasmid.Tet50.matrix.panel +
        B20.A31.Tet50.matrix.panel +
        B30.A31.Tet50.matrix.panel +
        B20.A18.Tet50.matrix.panel +
        B30.A18.Tet50.matrix.panel +
        B20.noPlasmid.Tet0.matrix.panel +
        B30.noPlasmid.Tet0.matrix.panel +
        B20.A31.Tet0.matrix.panel +
        B30.A31.Tet0.matrix.panel +
        B20.A18.Tet0.matrix.panel +
        B30.A18.Tet0.matrix.panel +
        plot_layout(nrow = 1)

    return(matrix.figure)

}

Fig6 <- MakeMutCountMatrixFigure(evolved.mutations)
matrix.outf <- "../results/draft-manuscript-1A/Fig6.pdf"
ggsave(matrix.outf, Fig6, height=8, width=12)

Fig6singles <- MakeMutCountMatrixFigure(evolved.mutations, show.all=T)
ggsave("../results/draft-manuscript-1A/Fig6-singles.pdf", Fig6singles, height=20, width=12)


Fig7 <- MakeSummedAlleleFrequencyMatrixFigure(evolved.mutations)
Fig7.matrix.outf <- "../results/draft-manuscript-1A/Fig7.pdf"
ggsave(Fig7.matrix.outf, Fig7, height=8, width=12)

Fig7singles <- MakeSummedAlleleFrequencyMatrixFigure(evolved.mutations, show.all=T)
ggsave("../results/draft-manuscript-1A/Fig7-singles.pdf", Fig7singles, height=20, width=12)

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
