##copy-number-analysis.R by Rohan Maddamsetti.

## 1) use xml2 to get negative binomial fit and dispersion from
## breseq output summary.html. This is H0 null distribution of 1x coverage.

## 2) Find intervals longer than max.read.len that reject H0 coverage in genome.
##    at an uncorrected alpha = 0.05. This is to have generous predicted boundaries for amplifications.

## 3) Do a more rigorous test for each region. take positions in the region separated by more than max.read.len,
## and determine the probability that all are independently significant under the null, compared to
## a corrected bonferroni. The max.read.len ensures positions cannot be spanned by a single Illumina read.

## 4) Estimate copy number by dividing mean coverage in each region by the mean
##   of the H0 1x coverage distribution.

## 5) return copy number and boundaries for each significant amplification.

## (so far, there are no significant amplifications).

library(xml2)
library(assertthat)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

library(tidyverse)
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(DT)          # prettier data.frame output
library(data.table)  # faster fread()
library(dtplyr)      # dplyr works with data.table now.
library(cowplot)     # layout figures nicely.

#' parse the summary.html breseq output file, and return the mean and dispersion
#' of the negative binomial fit to the read coverage distribution, returned as a
#' data.frame with columns {mean, dispersion}.
#' NOTE: this code has only been tested on the summary file
#' output by breseq 0.30.0. It might fail on earlier or later versions.

coverage.nbinom.from.html <- function(breseq.output.dir) {
    summary.html.f <- file.path(breseq.output.dir, "output", "summary.html")
    tree <- read_html(summary.html.f)
    ## print text in the table 'Reference Sequence Information.
    query <- '//table[./tr/th[contains(text(),"fit dispersion")]]'
    table <- xml_find_first(tree,query)
    table.data <- xml_find_all(table,'./tr/td')
    avg <- as.numeric(xml_text(table.data[5]))
    dispersion <- as.numeric(xml_text(table.data[6]))
    print(paste('mean coverage is:',avg))
    print(paste('dispersion is:',dispersion))
    return(data.frame('mean'=avg,'dispersion'=dispersion,'variance'=avg*dispersion))
}

#' get the maximum length of a sequencing read from the summary.html breseq
#' output file.
max.readlen.from.html <- function(breseq.output.dir) {
    summary.html.f <- file.path(breseq.output.dir, "output", "summary.html")
    tree <- read_html(summary.html.f)
    ## print text in the table 'Read File Information.
    query <- '//table[./tr/th[contains(text(),"longest")]]'
    table <- xml_find_first(tree,query)
    table.data <- xml_find_all(table,'./tr/td')
    readlen.index <- length(table.data) - 1
    max.readlen <- xml_integer(xml_find_all(table.data[readlen.index],".//b//text()"))
    print(paste('max read length is:',max.readlen))
    return(max.readlen)
}

#' Find intervals longer than max.read.len that reject H0 coverage in genome.
#' at an uncorrected alpha = 0.05. This is to have generous predicted boundaries for amplifications.
#' Then do a more rigorous test for each region. take positions in the region separated by more than max.read.len,
#' and determine the probability that all are independently significant under the null, compared to
#' a corrected bonferroni. max.read.len ensures positions cannot be spanned by a single Illumina read.
#' Estimate copy number by dividing mean coverage in each region by the mean of the H0 1x coverage distribution.
#' return mean copy number, and boundaries for each region that passes the amplification test.
find.amplifications <- function(breseq.output.dir, gnome) { #gnome is not a misspelling.
    
    gnome <- as.character(gnome)
    print(gnome)
    ## Use xml2 to get negative binomial fit and dispersion from
    ## breseq output summary.html. This is H0 null distribution of 1x coverage.
    nbinom.fit <- coverage.nbinom.from.html(breseq.output.dir)
    
    ## Use xml2 to get max read length from summary.html.
    max.read.len <- max.readlen.from.html(breseq.output.dir)
    genome.length <- 4584653 ## length of NZ_CP053607 reference.
    
    alpha <- 0.05
    uncorrected.threshold <- qnbinom(p=alpha,mu=nbinom.fit$mean,size=nbinom.fit$dispersion,lower.tail=FALSE)
    
    genome.coverage.file <- file.path(breseq.output.dir,"08_mutation_identification", "NZ_CP053607.coverage.tab")
    
    ## use dtplyr for speed!
    genome.coverage <- lazy_dt(fread(genome.coverage.file)) %>%
        select(position,unique_top_cov,unique_bot_cov) %>% mutate(coverage=unique_top_cov+unique_bot_cov)
    
    ## find candidate amplifications that pass the uncorrected threshold.
    candidate.amplifications <- genome.coverage %>%
        filter(coverage > uncorrected.threshold) %>%
        ## now finally turn into a dataframe (as using lazy_dt)
        as.data.frame()
    
    ## calculate intervals of candidate amplifications.
    boundaries <- candidate.amplifications %>%
        mutate(left.diff=position - lag(position)) %>%
        mutate(right.diff=lead(position) - position) %>%
        ## corner case: check for the NA values at the endpoints and set them as boundaries.
        mutate(is.right.boundary=is.na(right.diff)|ifelse(right.diff>1,TRUE,FALSE)) %>%
        mutate(is.left.boundary=is.na(left.diff)|ifelse(left.diff>1,TRUE,FALSE)) %>%
        filter(is.left.boundary==TRUE | is.right.boundary==TRUE)
 
    left.boundaries <- filter(boundaries,is.left.boundary==TRUE) %>%
        arrange(position)
        
    right.boundaries <- filter(boundaries,is.right.boundary==TRUE) %>%
        arrange(position)
    
    assert_that(nrow(left.boundaries) == nrow(right.boundaries))
    
    ## helper higher-order function to get min, max, mean coverage of each segment.
    get.segment.coverage <- function(left.bound,right.bound,coverage.table,funcx) {
        seg <- coverage.table %>% filter(position>left.bound) %>% filter(position<right.bound)
        return(funcx(seg$coverage))
    }
    
    amplified.segments <- data.frame(left.boundary=left.boundaries$position,right.boundary=right.boundaries$position) %>%
        ## filter out intervals less than 2 * max.read.len.
        mutate(len=right.boundary-left.boundary) %>% filter(len>(2*max.read.len)) %>% mutate(amplication.index=row_number())

    ## return empty dataframe  if there are no significant amplified segments.
    if (nrow(amplified.segments) == 0) return(data.frame())

    amplified.segments <- amplified.segments %>%
        ## find min, max, and mean coverage of each amplified segment.
        group_by(left.boundary,right.boundary) %>%
        summarise(coverage.min=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,min),
                  coverage.max=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,max),
                  coverage.mean=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,mean)) %>%
        mutate(len=right.boundary-left.boundary) %>%
        mutate(copy.number.min=coverage.min/nbinom.fit$mean,copy.number.max=coverage.max/nbinom.fit$mean,
               copy.number.mean=coverage.mean/nbinom.fit$mean)
    
    ##print(data.frame(amplified.segments))
    ## divide alpha by the number of tests for the bonferroni correction.
    bonferroni.alpha <- alpha/(genome.length + sum(amplified.segments$len))
    corrected.threshold <- qnbinom(p=bonferroni.alpha,mu=nbinom.fit$mean,size=nbinom.fit$dispersion,lower.tail=FALSE)
    
    ## This is my test: take the probability of the minimum coverage under H0 to the power of the number of
    ## uncorrelated sites in the amplification (sites more than max.read.len apart). Then see if this is smaller than the
    ## bonferroni corrected p-value for significance..
    significant.amplifications <- amplified.segments %>%
        mutate(pval=(pnbinom(q=coverage.min,
                             mu=nbinom.fit$mean,
                             size=nbinom.fit$dispersion,
                             lower.tail=FALSE))^(len%/%max.read.len)) %>%
        mutate(is.significant=ifelse(pval<bonferroni.alpha,TRUE,FALSE)) %>%
        filter(is.significant==TRUE) %>% mutate(Genome=as.character(gnome)) %>%
        mutate(bonferroni.corrected.pval=pval*alpha/bonferroni.alpha)
    
    return(significant.amplifications)
}

## input: LCA.gbk: file.path of the reference genome,
##    amplifications: data.frame returned by find.amplifications.
annotate.amplifications <- function(amplifications,LCA.gff) {
    
    ## create the IRanges object.
    amp.ranges <- IRanges(amplifications$left.boundary,amplifications$right.boundary)
    ## Turn into a GRanges object in order to find overlaps with REL606 genes.
    g.amp.ranges <- GRanges("NZ_CP053607",ranges=amp.ranges)
    ## and add the data.frame of amplifications as metadata.
    mcols(g.amp.ranges) <- amplifications
    
    ## find the genes within the amplifications.
    pLCA <- import.gff(LCA.gff)
    LCA.Granges <- as(pLCA, "GRanges")
    
    LCA.genes <- LCA.Granges[LCA.Granges$type == 'gene']
    ## find overlaps between annotated genes and amplifications.
    hits <- findOverlaps(LCA.genes,g.amp.ranges,ignore.strand=FALSE)
    
    ## take the hits, the LCA annotation, and the amplifications,
    ## and produce a table of genes found in each amplication.
    
    hits.df <- data.frame(query.index=queryHits(hits),subject.index=subjectHits(hits))
    
    query.df <- data.frame(query.index=seq_len(length(LCA.genes)),
                           gene=LCA.genes$Name,locus_tag=LCA.genes$ID,
                           start=start(ranges(LCA.genes)),end=end(ranges(LCA.genes)))
    
    subject.df <- bind_cols(data.frame(subject.index=seq_len(length(g.amp.ranges))),data.frame(mcols(g.amp.ranges)))
    
    amplified.genes.df <- left_join(hits.df,query.df) %>% left_join(subject.df) %>%
        ## if gene is NA, replace with locus_tag. have to change factors to strings!
        mutate(gene = ifelse(is.na(gene),as.character(locus_tag),as.character(gene)))
    
    return(amplified.genes.df)
}

plot.amp.segments <- function(annotated.amps,clone.labels) {
    
    ## for annotated.amps and clone.labels to play nicely with each other.
    clone.labels$Name <- as.character(clone.labels$Name)
    
    labeled.annotated.amps <- left_join(annotated.amps,clone.labels,by=c("Genome" = 'Name')) %>%
        select(-query.index,-subject.index,-is.significant,-SampleType, -Population) %>%
        ## replace 'sfcA' with 'maeA' in the plot.
        mutate(gene = replace(gene, gene == 'sfcA', 'maeA')) %>%
        mutate(log.pval=log(bonferroni.corrected.pval)) %>%
        mutate(log2.copy.number.mean=log2(copy.number.mean)) %>%
        transform(Population = PopulationLabel) %>%
        filter(!(Genome==ParentClone)) %>%
        mutate(left.boundary.MB = left.boundary/1000000) %>%
        mutate(right.boundary.MB = right.boundary/1000000) %>%
        mutate(Genome.Class=recode(Environment,
                                   DM0 = "DM0-evolved genomes",
                                   DM25 = "DM25-evolved genomes"))
    
    ## order the genes by start to get axes correct on heatmap.
    labeled.annotated.amps$gene <- with(labeled.annotated.amps, reorder(gene, start))
    ## reverse the order of genomes to make axes consistent with stacked barplot.
    labeled.annotated.amps$Genome <- factor(labeled.annotated.amps$Genome)
    labeled.annotated.amps$Genome <- factor(labeled.annotated.amps$Genome,
                                            levels=rev(levels(labeled.annotated.amps$Genome)))
    
    segmentplot <- ggplot(
        labeled.annotated.amps,
        aes(x=left.boundary.MB,
            xend=right.boundary.MB,
            y=Genome,
            yend=Genome,
            color=log2.copy.number.mean,
            size=20,
            frame=Genome.Class)) +
        geom_segment() +
        ## draw vertical lines at maeA, dctA.
        geom_vline(size=0.2,
                   linetype='dashed',
                   xintercept = c(1534704/1000000,3542785/1000000)
                   ) +
        xlab("Genomic position (Mb)") +
        ylab("") +
        scale_color_viridis(name=bquote(log[2]~"(copy number)"),option="plasma") +
        facet_wrap(~Genome.Class,nrow=2, scales = "free_y") +
        theme_classic(base_family='Helvetica') +
        guides(size=FALSE) +
        theme(legend.position="bottom") +
        theme(axis.ticks=element_line(size=0.1))
    return(segmentplot)
}

## assert that we are in the src directory, such that
## proj.dir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("transposon-plasmid-evolution","src")))
projdir <- file.path("..")

outdir <- file.path(projdir, "results", "draft-manuscript-1A")
breseq.output.dir <- file.path(projdir, "results", "draft-manuscript-1A", "genome-analysis", "mixed-pops")
RM6.200.6.gff3 <- file.path(projdir,"results", "draft-manuscript-1A", "genome-analyss", "RM6-200-6.gff3")
all.genomes <- list.files(breseq.output.dir,pattern='^RM')

all.genome.paths <- sapply(all.genomes, function(x) file.path(breseq.output.dir,x))
genome.input.df <- data.frame(Genome=all.genomes,path=all.genome.paths)

amps <- map2_df(genome.input.df$path,
                genome.input.df$Genome,
                find.amplifications) %>%
    ungroup()


