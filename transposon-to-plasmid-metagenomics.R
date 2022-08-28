## transposon-to-plasmid-metagenomics.R by Rohan Maddamsetti.

## Here's a list of NEB reference genomes. I use the most recent version of NEB5-alpha as reference.
## https://international.neb.com/tools-and-resources/usage-guidelines/competent-e-coli-genome-sequences-tool

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggthemes)
library(viridis)
library(ggrepel)

## DH5a origin, based on aligning Jeff Barrick's manual annotation of the
## REL606 oriC sequence against NZ_CP017100 using NCBI BLAST.
NZ_CP017100_oriC_START = 3866291
NZ_CP017100_oriC_END = 3866522
NZ_CP017100_oriC_MID = (NZ_CP017100_oriC_START+NZ_CP017100_oriC_END)/2

## GC skew calculations using the webserver at:
## https://genskew.csb.univie.ac.at/webskew
## Also, see: https://skewdb.org/view/?seq=NZ_CP017100.1.
## GCskew_max <- 1494058
## GCskew_min <- 3858886
## GCskew_min is in atpA, so this is not exactly right.
## use Jeff Barrick's annotation.

## IDEA, perhaps for future work: verify GC skew and replication origin correlation,
## by comparing GC skew against actual replication of origin, based on looking
## at wave pattern in sequencing coverage in my Tet50 genome sequencing samples--
## some of these cultures were still in exponential phase when I sequenced them.
## what about for plasmids? see skewDB.

rotate.NEB5alpha.chr <- function(my.position) {
    #' function to rotate REL606 genome coordinates,
    #' setting oriC at the center of plots
    #' that examine mutation bias over the chromosome.
    ## we want to change coordinates so that c is the new origin.
    GENOME.LENGTH <- 4583637
    midpoint <- GENOME.LENGTH/2
    oriC <- 3886105
    
    if (oriC >= midpoint) {
        L <- oriC - midpoint
        ifelse(my.position > L, my.position - oriC, GENOME.LENGTH - oriC + my.position)
    } else { ## midpoint is greater than new.origin.
        L <- midpoint + oriC
        ifelse(my.position > L, my.position - GENOME.LENGTH - oriC, my.position - oriC)
    }
}


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

evolved.mutations <- read.csv(
    file.path(projdir,
              "results/genome-analysis/evolved_mutations.csv"),
    stringsAsFactors=FALSE) %>%
    ## rotate genome coordinates based on oriC.
    mutate(oriC.coordinate=rotate.NEB5alpha.chr(Position)) %>%
    mutate(Mbp.coordinate=oriC.coordinate/1000000)


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

###############################################
## Figure 2AB: Plot the distribution of measured allele frequencies in each population.

make.allele.freq.histogram <- function(evolved.mutations.df, my.title,annotate=FALSE) {
    p <- ggplot(evolved.mutations.df, aes(x=Frequency)) +
        geom_histogram(bins = 100) +
        theme_classic() +
        ylab("Count") +
        xlab("Allele Frequency") +
        scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1), limits = c(0,1.1)) +
        ggtitle(my.title) +
        facet_grid(Plasmid~.) +
    geom_vline(xintercept=0.10,color="red",linetype="dashed",size=0.2)

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


Tet0.evolved.mutations <- filter(evolved.mutations, Tet==0)
Tet50.evolved.mutations <- filter(evolved.mutations, Tet==50)
Fig2A <- make.allele.freq.histogram(Tet0.evolved.mutations, "Tet 0 populations",TRUE)
Fig2B <- make.allele.freq.histogram(Tet50.evolved.mutations, "Tet 50 populations",TRUE)

## This is a very important plot: what does this distribution say about
## the possibility of false positives? how can I interpret this?
## any theoretical basis in population genetics?

## Idea for an empirical control.
## 1) downsample reads from the treatment without plasmid, and re-run breseq
## to see if false positive mutation calls arise when coverage is ~40X rather than
## 300X.

Fig2AB <- Fig2A + Fig2B
fig2AB.output <- "../results/draft-manuscript-1A/Fig2AB.pdf"
ggsave(Fig2AB, file=fig2AB.output,width=10,height=4)

###############################################
## Figure 2C: make a stacked bar plot of the kinds of mutations in each treatment.
## Figure 2D: make a stacked bar plot of the kinds of mutations in each treatment, weighted by allele frequency.

## This function sums mutations per replicate population.
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
        group_by(Sample, Transposon, Plasmid, Tet, Population, Mutation) %>%
        summarize(Count=n(),WeightedCount = sum(Frequency)) %>%
        ungroup() %>%
        data.frame() %>%
        mutate(Mutation=as.factor(as.character(Mutation)))
}


plot.mutation.summary.stackbar <- function(mutation.class.df, leg=FALSE, weight.by.freq=FALSE) {

    if (weight.by.freq) {
        fig <- ggplot(mutation.class.df, aes(x=Plasmid, y=WeightedCount, fill=Mutation)) +
            ylab("Summed Allele Frequency")
    } else {
        fig <- ggplot(mutation.class.df, aes(x=Plasmid, y=Count, fill=Mutation)) +
            ylab("Count")
    }

    fig <- fig +
        ## show both tetracycline concentrations.
        facet_wrap(.~Tet) +
        geom_bar(stat='identity') +
        scale_fill_brewer(palette = "RdYlBu", direction=-1,drop=FALSE) +        
        theme_classic(base_family='Helvetica') +
        theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
              axis.text.y=element_text(size=12),
              panel.border=element_blank(),
              strip.background = element_blank(),
              panel.spacing.x=unit(1, "cm"),
              panel.spacing.y=unit(0.5, "cm"))

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

## Now make Figure 2CD.
mutation.class.df <- make.mutation.class.df(evolved.mutations) %>%
    mutate(Tet=recode(Tet,
                      `0` = "Tet 0",
                      `50` = "Tet 50"))

fig2C <- plot.mutation.summary.stackbar(mutation.class.df, FALSE, FALSE) 
## Repeat, but weight by allele frequency.
fig2D <- plot.mutation.summary.stackbar(mutation.class.df, TRUE, TRUE)
## pull the legend from fig2D
fig2CDlegend <- cowplot::get_legend(fig2D)
## and remove the legend from fig2D
fig2D <- fig2D + guides(fill = "none")

fig2CD <- fig2C + fig2D
full.fig2CD <- plot_grid(fig2CD, fig2CDlegend,nrow=2, rel_heights = c(1, 0.1))

## save figure 2CD.
fig2CD.output <- "../results/draft-manuscript-1A/Fig2CD.pdf"
ggsave(full.fig2CD, file=fig2CD.output,width=8,height=5)


#####################################################################################

## let's take a close look at the different kinds of evolved mutations.

## ALL of these MOB insertions are miniTn5-Tet insertions, either into the KanR gene
## on the plasmid, or into chromosomal genes in the no plasmid treatment.
evolved.MOB <- evolved.mutations %>% filter(Mutation == "MOB") 

## extremely strong parallel evolution in tetA-- seen only in no plasmid treatment.
## the parallel evolution is at the nucleotide level in the no plasmid treatment:
## (3100, 3102) in the tetA promoter, (3134 has two independent dS mutations, 3135
## has a nonsynonymous mutation,), and (4285, 4285, 4286 are independent nonsense mutations
## that truncate the very end of the protein.)
evolved.tetA <- evolved.mutations %>% filter(str_detect(Gene, "tetA")) %>%
    filter(Mutation != "MOB") %>%
    arrange(Sample, Position)

## parallel evolution of robA in no plasmid treatment.
evolved.nonsynonymous <- evolved.mutations %>% filter(Mutation == "nonsynonymous") %>%
    arrange(Gene, Position, Sample)

## most synonymous mutations in pops 3,4,5  of the pUC plasmid treatment.
## This will certainly be a significant association. perhaps some cryptic hypermutator
## clades here?
evolved.synonymous <- evolved.mutations %>% filter(Mutation == "synonymous") %>%
    arrange(Gene, Position, Sample)


evolved.INDEL <- evolved.mutations %>% filter(Mutation == "INS" | Mutation == "DEL")

evolved.intergenic <- evolved.mutations %>% filter(Mutation == "intergenic") %>%
    arrange(Gene, Position, Plasmid, Sample) %>%
    select(-Transposon, -Mutation, -Mutation_Category, -Population)



#####################################################################################
## examine DNA repair and DNA polymerase/replication genes for mutator and anti-mutator
## candidates.

DNA.repair.loci <- read.csv("../data/draft-manuscript-1A/DNA-repair-and-replication.csv",header=TRUE,as.is=TRUE)
## some mutations in DNA repair genes, but hard to conclude anything from this.
DNA.repair.muts <- filter(evolved.mutations, Gene %in% unique(DNA.repair.loci$Gene))

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

## check parallel evolution for synonymous mutations too.
parallel.dS.Table <- filter(evolved.mutations, Position %in% parallel.dS$Position) %>% arrange(Position)


##################################################################################
## analysis of parallel evolution at the gene level (including intergenic regions).

gene.level.parallel.mutations <- evolved.mutations %>% group_by(Gene) %>%
summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.genes <- gene.level.parallel.mutations %>%
    select(Gene, count, Plasmid, Transposon, Tet) %>%
    distinct() %>%
    arrange(desc(count))



################################################################################
### Figure 3: make a matrix plot of genes with mutations in two or more clones.
################################################################################
MakeMutCountMatrixFigure <- function(evolved.muts, show.all=FALSE, use.treatment.hit.sort=FALSE) {

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.muts %>%
        ## unite the Transposon, Plasmid, Tet columns together.
        unite("Treatment", Transposon:Tet, sep="\n", remove = FALSE) %>%
        group_by(Gene, Sample, Transposon, Plasmid, Tet, Treatment) %>%
        summarize(mutation.count = n()) %>%
        ## This is for sorting mutations.
        mutate(is.MOB = ifelse(str_detect(Gene,"tetA-Tn5"), TRUE, FALSE))
    
    total.muts <- matrix.data %>%
        group_by(Gene) %>%
        summarize(total.mutation.count = sum(mutation.count))
    
    matrix.data <- left_join(matrix.data, total.muts)
    
    if (!show.all) { ## then filter out genes that are only hit in one sample.
        matrix.data <- matrix.data %>% filter(total.mutation.count > 1)
    }
    
    ## sort genes by number of mutations in each row, but put all the transposon mutations together.
    ## put all the transposon mutations t
    ## also check out the alternate sorting method that follows.
    gene.hit.sort <- matrix.data %>%
        group_by(Gene, is.MOB, .drop = FALSE) %>%
        summarize(hits=sum(mutation.count)) %>%
        arrange(desc(is.MOB), desc(hits))
    
    ## now sort genes.
    if (use.treatment.hit.sort) {
        ## alternate sorting method: difference in hits between environments,
        ## AKA the (absolute value of the) difference in number of pops with hits
        ## between the None and pUC treatments.
        pUC.hit.count.df <- filter(matrix.data,Plasmid=="pUC") %>%
            group_by(Gene, is.MOB, .drop = FALSE) %>%
            summarize(pUC.hit.count=n())
        
        noPlasmid.hit.count.df <- filter(matrix.data,Plasmid=="None") %>%
            group_by(Gene, is.MOB, .drop = FALSE) %>%
            summarize(noPlasmid.hit.count=n())
        
        treatment.hit.sort <- full_join(pUC.hit.count.df, noPlasmid.hit.count.df) %>%
            mutate(hit.diff = pUC.hit.count - noPlasmid.hit.count) %>%
            arrange(desc(is.MOB), desc(hit.diff))

        ## cast Gene into a factor for plotting.
        matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(treatment.hit.sort$Gene))
    } else {
        matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.hit.sort$Gene))
    }
    
    ## cast mutation.count into a factor for plotting.
    matrix.data$mutation.count <- factor(matrix.data$mutation.count)

    make.matrix.panel <- function(mdata, treatment, leg=FALSE) {
        panel.data <- filter(mdata,Treatment==treatment)
        fig <- ggplot(panel.data,
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
                  axis.title.y = element_blank()) +
            scale_y_discrete(drop=FALSE) + ## don't drop missing genes.
            scale_fill_manual(name="Mutations",
                              values = c("#ffdf00", "#bebada", "#fb8072", "#80b1d3", "#fdb462"))
        
        if (leg == FALSE) {
            fig <- fig + guides(fill= "none")
        }
        return(fig)
    }

    
    ## make Tet50 panels.
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
    
    ## Using the patchwork library for layout.
    matrix.figure <-
        B20.noPlasmid.Tet50.matrix.panel +
        B30.noPlasmid.Tet50.matrix.panel +
        B20.A31.Tet50.matrix.panel +
        B30.A31.Tet50.matrix.panel +
        B20.A18.Tet50.matrix.panel +
        B30.A18.Tet50.matrix.panel +
        plot_layout(nrow = 1)
    return(matrix.figure)
}


## Use summed allele frequency for the heatmap.
MakeSummedAlleleFrequencyMatrixFigure <- function(evolved.muts,
                                                  allele.freq.threshold = 0.2,
                                                  show.all=FALSE) {

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.muts %>%
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

    ## sort genes by the total allele frequency in each row.
    ## also check out the alternate sorting method that follows.
    gene.freq.sort <- matrix.data %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(totalallelefreq = sum(summed.Allele.Frequency)) %>%
        arrange(desc(totalallelefreq))
    
    ## alternative sorting method:
    ## difference in allele frequency between the pUC and noPlasmid treatments..
    pUC.allele.freq.df <- filter(matrix.data, Plasmid=="pUC") %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(pUC.allele.frequency=sum(summed.Allele.Frequency))
    
    noPlasmid.allele.freq.df <- filter(matrix.data, Plasmid=="None") %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(noPlasmid.allele.frequency = sum(summed.Allele.Frequency))
    
    treatment.freq.sort <- full_join(pUC.allele.freq.df, noPlasmid.allele.freq.df) %>%
        replace_na(list(pUC.allele.frequency = 0, pUC.allele.frequency = 0)) %>%
        mutate(allele.diff = noPlasmid.allele.frequency - noPlasmid.allele.frequency) %>%
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
    
    ## Using the patchwork library for layout.
    matrix.figure <-
        B20.noPlasmid.Tet50.matrix.panel +
        B30.noPlasmid.Tet50.matrix.panel +
        B20.A31.Tet50.matrix.panel +
        B30.A31.Tet50.matrix.panel +
        B20.A18.Tet50.matrix.panel +
        B30.A18.Tet50.matrix.panel +
        plot_layout(nrow = 1)

    return(matrix.figure)
}


## examine genes that show parallelism across Tet 0 and Tet 50.
parallel.genes.across.Tet0.and.Tet50 <- parallel.genes %>%
    select(Gene, Tet) %>%
    distinct() %>%
    group_by(Gene) %>%
    summarize(num.tet.conc.found.in = n()) %>%
    filter(num.tet.conc.found.in == 2)

parallel.mutations.across.Tet0.and.Tet50 <- evolved.mutations %>%
    filter(Gene %in% parallel.genes.across.Tet0.and.Tet50$Gene)

## the intergenic mutations on the plasmid occur in the oriC sequences.
## IMPORTANT TODO: analyze these in a separate figure, to relate to the
## evolved changes in plasmid copy number.
plasmid.origin.muts <- parallel.mutations.across.Tet0.and.Tet50 %>%
    filter(Gene == "–/KanR")

## remove these from the other parallel mutations across the 2 treatments.
## the remaining mutations are not that interesting.
## 3 loci: mrcA, narU, rpsA, yeeJ, 2 parallel mutations in each
filtered.parallel.mutations.across.Tet0.and.Tet50 <- parallel.mutations.across.Tet0.and.Tet50 %>%
    filter(Gene != "–/KanR")


## genes that only show parallelism in Tet 0.
parallel.genes.in.Tet0 <- parallel.genes %>%
    select(Gene, Tet) %>%
    filter(Tet == 0) %>%
    distinct() %>%
    filter(!(Gene %in% parallel.genes.across.Tet0.and.Tet50$Gene))
## not a whole lot that is interesting. maybe the high frequency oxyR mutations.
## Off the top of my head, I believe that is a high level I-modulon regulator.
parallel.mutations.in.only.Tet0 <- evolved.mutations %>%
    filter((Gene %in% parallel.genes.in.Tet0$Gene))


## Figure 3.
## genes that only show parallelism in Tet 50, and all MOB (these are only in Tet50 treatment anyway).
parallel.genes.in.Tet50 <- parallel.genes %>%
    select(Gene, Tet) %>%
    filter(Tet == 50) %>%
    distinct() %>%
    filter(!(Gene %in% parallel.genes.across.Tet0.and.Tet50$Gene))

parallel.mutations.in.only.Tet50 <- evolved.mutations %>%
    filter(Gene %in% parallel.genes.in.Tet50$Gene)

Fig3.data <- full_join(evolved.MOB,
                       filter(parallel.mutations.in.only.Tet50, Allele != "MOB")) %>%
    filter(Frequency > 0.10)

Fig3A <- MakeMutCountMatrixFigure(Fig3.data,
                                 show.all=TRUE, ## This is needed to show the MOB insertions too.
                                 use.treatment.hit.sort=FALSE)

Fig3.outf <- "../results/draft-manuscript-1A/Fig3A.pdf"
ggsave(Fig3.outf, Fig3A, height=6, width=12)

S1Fig <- MakeMutCountMatrixFigure(evolved.data, show.all=TRUE, use.treatment.hit.sort=FALSE)
S1matrix.outf <- "../results/draft-manuscript-1A/S1Fig.pdf"
ggsave(S1matrix.outf, S1Fig, height=8, width=12)


Fig3B <- MakeSummedAlleleFrequencyMatrixFigure(Fig3.data, show.all=T)
Fig3B.matrix.outf <- "../results/draft-manuscript-1A/Fig3B.pdf"
ggsave(Fig3B.matrix.outf, Fig3B, height=8, width=12)

Fig7singles <- MakeSummedAlleleFrequencyMatrixFigure(evolved.mutations, show.all=T)
ggsave("../results/draft-manuscript-1A/Fig7-singles.pdf", Fig7singles, height=20, width=12)

