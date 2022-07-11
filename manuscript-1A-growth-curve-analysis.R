## manuscript-1A-growth-curve-analysis.R by Rohan Maddamsetti.
## This script plots growth over time.

## If all p15 plasmids within the cell has the transposon,
## but 50% of the pUC plasmids have the transposon,
## then the change in the frequency of the allele on the plasmid
## could allow the strains with the high copy number plasmid to have
## phenotypic plasticity in the resistance phenotype, depending
## on the growth conditions before the plating.

## I could test this hypothesis by preconditioning the p15A and pUC clones in
## LB, LB+Tet20, and LB+Tet50,
## and then measuring copy number by qPCR,
## and measuring growth curves again to see if there is any change.


## CRITICAL TODO: for growth rate vs. lag measurements,
## fit a linear model to estimate lag_time as follows:
## time to hit OD threshold = lag_time + growth_rate * time

library(tidyverse)
library(cowplot)
library(broom)

### Make a plot of the Tet antibiotic treatment over time.
experiment.design.df <- data.frame(Day = c(1,2,3,4,5,6,7,8,9),
                                   Tet.conc = c(2,4,6,8,10,20,30,40,50))
treatment.plot <- ggplot(experiment.design.df,
                         aes(x = Day, y = Tet.conc)) +
    theme_classic() + geom_line() +
    ylab("Tetracycline (ug/mL)") +
    xlab("Day")


subtract.OD600.blanks <- function(OD600.df) {
    blanks.df <- filter(OD600.df, Treatment == "Blank")
    media.blank <- mean(blanks.df$RawOD600)
    subtracted.df <- OD600.df %>%
        filter(Treatment != "Blank") %>%
        mutate(OD600 = RawOD600 - media.blank)
    return(subtracted.df)
}

################################################################################
## plot OD600 over time in the nine day evolution experiment.
december.11.OD600.df <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-12-11-Full-9-day-OD600.csv") %>%
    ## subtract LB blanks to get the true OD measurement
    subtract.OD600.blanks() %>%
    mutate(BiologicalReplicate = as.factor(BiologicalReplicate)) %>%
    filter(!(Treatment %in% c("Blank", "Empty")))
    
december.11.OD600.plot <- ggplot(december.11.OD600.df,
                              aes(x = Day,
                                  y = OD600,
                                  ##shape = BiologicalReplicate,
                                  color = Treatment)) +
    facet_wrap(.~Transposon) +
    theme_classic() + geom_point() + geom_smooth()
ggsave("../results/draft-manuscript-1A/growth-results/Full-9-day-evolution-experiment-OD600-2021-12-11.pdf", december.11.OD600.plot)

################################################################################
## analyze data from the large-scale growth curve experiments.
## CRITICAL NOTE: The layout of the p15A plasmid and pUC plasmids are SWITCHED
## on the plate relative to the experiment conducted in November 2021!!!!
well.to.treatment.for.big.clone.OD600.experiment <- function(Well) {
    ## This helper function maps the well in
    ## the 96-well plate to the biological replicate.

    rowletter <- substring(Well,1,1)
    colnum <- as.numeric(substring(Well,2))
    
    ## The first and last columns are filled with saline.
    if (colnum == 1 || colnum == 12)
        return("Empty")
    
    if (rowletter == "A" || rowletter == "B")
        return("no_plasmid")

    if (rowletter == "C") {
        if (colnum <= 6)
            return("no_plasmid")
        else
            return("p15A_plasmid")
    }

    if (rowletter == "D" || rowletter == "E")
        return("p15A_plasmid")
    if (rowletter == "F" || rowletter == "G")
        return("pUC_plasmid")

    if (rowletter == "H") {
        if (colnum <= 6)
            return("pUC_plasmid")
        else
            return("Blank")
    }

    return(my.treatment)
}


tidy.tecan.clone.plate.data <- function(long.format.data, Tet.conc) {
    long.format.data %>%
        pivot_longer(cols=V1, values_to="Well") %>%
        select(-name) %>% ## drop this useless column
        pivot_longer(!Well,
                     names_to = "Time",
                     values_to = "RawOD600",
                     names_pattern = "V(.+)") %>%
        mutate(Time = as.numeric(Time) - 1) %>%
        mutate(minutes = 10*Time) %>%
        mutate(hours = minutes/60) %>%
        mutate(Tet = Tet.conc) %>%
        mutate(Treatment = sapply(Well, well.to.treatment.for.big.clone.OD600.experiment)) %>%
        mutate(Treatment = as.factor(Treatment)) %>%
        ## subtract LB blanks to get the true OD measurement.
        subtract.OD600.blanks() %>%
        filter(!(Treatment %in% c("Blank", "Empty")))
}

plot.OD600.over.time <- function(tidy.data) {
    tidy.data %>%
        ggplot(
            aes(x = hours,
                y = OD600,
                color = Treatment)) +
        geom_point(size=0.5) +
        geom_vline(xintercept=24,linetype="dashed",color="red")
}


plot.delta.log2.OD600.over.time <- function(tidy.data) {
    ## let's plot the derivatives of log(OD600) over time.
    delta.log2.OD600.data <- tidy.data %>%
        group_by(Well, Treatment) %>%
        mutate(log2.OD600 = log2(OD600)) %>%
        mutate(delta.log2.OD600 = log2.OD600 - lag(log2.OD600))
    
    delta.log2.OD600.plot <- delta.log2.OD600.data %>%
        group_by(Well) %>%
        ## ignore the first hour of data; remove this noise.
        filter(hours > 1) %>%
        ## ignore data after 30 hours; remove this noise
        ## (maybe meaningful, but skip for now).
        filter(hours < 30) %>%
        ggplot(
            aes(x = hours,
                y = delta.log2.OD600,
                color = Well)) + geom_line() + theme_classic() +
        facet_grid(Treatment~.) +
        ylab("Delta log2(OD600)") +
        xlab("Time (hours)") +
        guides(color = "none")
}

make.max.growth.rate.df <- function(tidy.data) {
    tidy.max.growth.rate.data <- tidy.data %>%
        filter(OD600 >= 0.3) %>%
        filter(OD600 <= 0.4) %>%
        group_by(Well, Treatment) %>%
        nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
        mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
        ## get the fit parameters.
        ## see tutorial at:
        ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
        mutate(tidied = map(fit, tidy)) %>%
        select(Well, Treatment, tidied) %>%
        unnest(tidied) %>%
        filter(term == 'hours') ## we only care about the slope parameter.
}


summarize.time.lag <- function(well.df, OD.threshold = 0.125) {
    ## to compare time lags, measure point when OD600 hits the threshold
    lag.min.index <- max(min(which(well.df$OD600 > OD.threshold)), 1)
    t.OD.hit.lag.min <- well.df[lag.min.index,]$hours

    ret.df <- mutate(well.df, time.lag = t.OD.hit.lag.min)
    return(ret.df)
}


make.time.lag.df <- function(tidy.data, lag.threshold = 0.125) {

    ## make a one-variable function given the lag.threshold.
    .summarize.time.lag <- partial(.f = summarize.time.lag, OD.threshold = lag.threshold)
    
    tidy.data %>%
        ## filter out the first 30 minutes of data.
        filter(hours > 0.5) %>%
        split(.$Well) %>%
        map_dfr(.f=.summarize.time.lag) %>%
        select(Well, Tet, Treatment, time.lag) %>%
        distinct() %>%
        mutate(inverse.lag = 1/time.lag)
}


plot.growth.lag <- function(lag.df) {
    lag.df %>%
        ggplot(
            aes(x = Treatment,
                y = time.lag,
                color = Treatment)) + geom_point() + theme_classic() +
        ylab("Time to reach OD600 threshold")
}


plot.max.growth.rate <- function(tidy.max.growth.rate.data) {   
    tidy.max.growth.rate.data %>%
        ggplot(
            aes(x = Treatment,
                y = estimate,
                color = Treatment)) + geom_point() + theme_classic() +
        ylab("growth rate between OD600 0.3-0.4")
}


plot.max.growth.rate.histogram <- function(tidy.max.growth.rate.data) {
    tidy.max.growth.rate.data %>%
        ggplot(
            aes(fill = Treatment,
                x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
        facet_grid(Treatment~.) +
        xlab("growth rate between OD600 0.3-0.4") +
        labs(fill = "Presence/absence of plasmid") +
        theme(legend.position="top")
}

plot.24h.OD600.histogram <- function(tidy.data) {
    tidy.data %>%
        filter(hours == 24) %>%
        ggplot(
            aes(fill = Treatment,
                x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
        facet_grid(Treatment~.) +
        ylab("OD600 at 24h timepoint") +
        labs(fill = "Presence/absence of plasmid") +
        theme(legend.position="top")
}


plot.growth.rate.vs.lag <- function(lag.df, max.growth.rate.df) {
    lag.vs.growth.rate.df <- full_join(lag.df, max.growth.rate.df)
    lag.vs.growth.rate.df %>%
        ggplot(aes(fill = Treatment, y = time.lag, x = estimate, color = Treatment)) +
        geom_point() + theme_classic() +
        ylab("Time to OD Threshold (h)") + xlab("growth rate between OD600 0.3-0.4")
}

########################################################
## B30 clones, Replicate 1.
########################################################
## Import and tidy the data.

Jan.30.LB.evolved.clone.long.format.data <- read.csv(
    "../data/draft-manuscript-1A/OD600-data/2022-01-30-LB-B30-clones.csv", header=FALSE)
Jan.30.LB.tidy.data <- Jan.30.LB.evolved.clone.long.format.data %>%
    tidy.tecan.clone.plate.data(0)

Jan.30.Tet20.evolved.clone.long.format.data <- read.csv(
    "../data/draft-manuscript-1A/OD600-data/2022-01-30-LB-Tet20-B30-clones.csv",
    header=FALSE)
Jan.30.Tet20.tidy.data <- Jan.30.Tet20.evolved.clone.long.format.data %>%
    tidy.tecan.clone.plate.data(20)

Jan.30.Tet50.evolved.clone.long.format.data <- read.csv(
    "../data/draft-manuscript-1A/OD600-data/2022-01-30-LB-Tet50-B30-clones.csv",
    header=FALSE)
Jan.30.Tet50.tidy.data <- Jan.30.Tet50.evolved.clone.long.format.data %>%
    tidy.tecan.clone.plate.data(50)

## write out data for Emrah.
write.csv(x=Jan.30.LB.tidy.data, file="../results/draft-manuscript-1A/growth-results/2022-01-30-LB-Tet0-B30-clones-tidy.csv", row.names=FALSE, quote = FALSE)
write.csv(x=Jan.30.Tet20.tidy.data, file="../results/draft-manuscript-1A/growth-results/2022-01-30-LB-Tet20-B30-clones-tidy.csv", row.names=FALSE, quote = FALSE)
write.csv(x=Jan.30.Tet50.tidy.data, file="../results/draft-manuscript-1A/growth-results/2022-01-30-LB-Tet50-B30-clones-tidy.csv", row.names=FALSE, quote = FALSE)

############################
## Plot growth over time.

Jan.30.LB.plot <- plot.OD600.over.time(Jan.30.LB.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-LB-plot.pdf",
       Jan.30.LB.plot)

Jan.30.delta.log2.LB.plot <- plot.delta.log2.OD600.over.time(Jan.30.LB.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-LB-delta-growth.pdf",
       Jan.30.delta.log2.LB.plot)

Jan.30.Tet20.plot <- plot.OD600.over.time(Jan.30.Tet20.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet20-plot.pdf",
       Jan.30.Tet20.plot)

Jan.30.delta.log2.Tet20.plot <- plot.delta.log2.OD600.over.time(Jan.30.Tet20.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet20-delta-growth.pdf",
       Jan.30.delta.log2.Tet20.plot)

Jan.30.Tet50.plot <- plot.OD600.over.time(Jan.30.Tet50.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet50-plot.pdf",
       Jan.30.Tet50.plot)

Jan.30.delta.log2.Tet50.plot <- plot.delta.log2.OD600.over.time(Jan.30.Tet50.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet50-delta-growth.pdf",
       Jan.30.delta.log2.Tet50.plot)

#######################################
## Plot lag distributions.

Jan.30.LB.lag.df <- make.time.lag.df(Jan.30.LB.tidy.data)
Jan.30.LB.lag.plot <- plot.growth.lag(Jan.30.LB.lag.df)

Jan.30.Tet20.lag.df <- make.time.lag.df(Jan.30.Tet20.tidy.data)
Jan.30.Tet20.lag.plot <- plot.growth.lag(Jan.30.Tet20.lag.df)

Jan.30.Tet50.lag.df <- make.time.lag.df(Jan.30.Tet50.tidy.data)
Jan.30.Tet50.lag.plot <- plot.growth.lag(Jan.30.Tet50.lag.df)

#######################################
## Plot max growth rate distributions.

Jan.30.LB.max.growth.rate.df <- make.max.growth.rate.df(Jan.30.LB.tidy.data)
Jan.30.LB.max.growth.rate.plot <- plot.max.growth.rate(Jan.30.LB.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-LB-growth-rate.pdf",
       Jan.30.LB.max.growth.rate.plot)
Jan.30.LB.max.growth.rate.histogram <- plot.max.growth.rate.histogram(Jan.30.LB.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-LB-growth-rate-histogram.pdf",
       Jan.30.LB.max.growth.rate.histogram)

Jan.30.Tet20.max.growth.rate.df <- make.max.growth.rate.df(Jan.30.Tet20.tidy.data)
Jan.30.Tet20.max.growth.rate.plot <- plot.max.growth.rate(Jan.30.Tet20.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet20-growth-rate.pdf",
       Jan.30.Tet20.max.growth.rate.plot)
Jan.30.Tet20.max.growth.rate.histogram <- plot.max.growth.rate.histogram(Jan.30.Tet20.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet20-growth-rate-histogram.pdf", Jan.30.Tet20.max.growth.rate.histogram)


Jan.30.Tet50.max.growth.rate.df <- make.max.growth.rate.df(Jan.30.Tet50.tidy.data)
Jan.30.Tet50.max.growth.rate.plot <- plot.max.growth.rate(Jan.30.Tet50.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet50-growth-rate.pdf",
       Jan.30.Tet20.max.growth.rate.plot)
Jan.30.Tet50.max.growth.rate.histogram <- plot.max.growth.rate.histogram(Jan.30.Tet50.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet50-growth-rate-histogram.pdf", Jan.30.Tet50.max.growth.rate.histogram)

########################################
## compare growth rates to lags.

Jan.30.LB.lag.vs.growth.rate.plot <- plot.growth.rate.vs.lag(
    Jan.30.LB.lag.df,
    Jan.30.LB.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-LB-growth-rate-vs-lag.pdf",
       Jan.30.LB.lag.vs.growth.rate.plot)

Jan.30.Tet20.lag.vs.growth.rate.plot <- plot.growth.rate.vs.lag(
    Jan.30.Tet20.lag.df,
    Jan.30.Tet20.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet20-growth-rate-vs-lag.pdf",
       Jan.30.Tet20.lag.vs.growth.rate.plot)

Jan.30.Tet50.lag.vs.growth.rate.plot <- plot.growth.rate.vs.lag(
    Jan.30.Tet50.lag.df,
    Jan.30.Tet50.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet50-growth-rate-vs-lag.pdf",
       Jan.30.Tet50.lag.vs.growth.rate.plot)


########################################
## plot OD600 at 24 hour distributions.

Jan.30.LB.24h.OD600.histogram <- plot.24h.OD600.histogram(Jan.30.LB.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-LB-24h-OD600-histogram.pdf", Jan.30.LB.24h.OD600.histogram)

Jan.30.Tet20.24h.OD600.histogram <- plot.24h.OD600.histogram(Jan.30.Tet20.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet20-24h-OD600-histogram.pdf", Jan.30.Tet20.24h.OD600.histogram)

Jan.30.Tet50.24h.OD600.histogram <- plot.24h.OD600.histogram(Jan.30.Tet50.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Jan-30-Tet50-24h-OD600-histogram.pdf", Jan.30.Tet50.24h.OD600.histogram)

########################################################
## B30 clones, Replicate 2.
########################################################
## Import and tidy the data.

Feb.01.LB.evolved.clone.long.format.data <- read.csv(
    "../data/draft-manuscript-1A/OD600-data/2022-02-01-LB-B30-clones.csv", header=FALSE)
Feb.01.LB.tidy.data <- Feb.01.LB.evolved.clone.long.format.data %>%
    tidy.tecan.clone.plate.data(0)

Feb.01.Tet20.evolved.clone.long.format.data <- read.csv(
    "../data/draft-manuscript-1A/OD600-data/2022-02-01-LB-Tet20-B30-clones.csv",
    header=FALSE)
Feb.01.Tet20.tidy.data <- Feb.01.Tet20.evolved.clone.long.format.data %>%
    tidy.tecan.clone.plate.data(20)

Feb.01.Tet50.evolved.clone.long.format.data <- read.csv(
    "../data/draft-manuscript-1A/OD600-data/2022-02-01-LB-Tet50-B30-clones.csv",
    header=FALSE)
Feb.01.Tet50.tidy.data <- Feb.01.Tet50.evolved.clone.long.format.data %>%
    tidy.tecan.clone.plate.data(50)

## write out data for Emrah.
write.csv(x=Feb.01.LB.tidy.data, file="../results/draft-manuscript-1A/growth-results/2022-02-01-LB-Tet0-B30-clones-tidy.csv", row.names=FALSE, quote = FALSE)
write.csv(x=Feb.01.Tet20.tidy.data, file="../results/draft-manuscript-1A/growth-results/2022-02-01-LB-Tet20-B30-clones-tidy.csv", row.names=FALSE, quote = FALSE)
write.csv(x=Feb.01.Tet50.tidy.data, file="../results/draft-manuscript-1A/growth-results/2022-02-01-LB-Tet50-B30-clones-tidy.csv", row.names=FALSE, quote = FALSE)

############################
## Plot growth over time.

Feb.01.LB.plot <- plot.OD600.over.time(Feb.01.LB.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-LB-plot.pdf",
       Feb.01.LB.plot)
Feb.01.delta.log2.LB.plot <- plot.delta.log2.OD600.over.time(Feb.01.LB.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-LB-delta-growth.pdf",
       Feb.01.delta.log2.LB.plot)

Feb.01.Tet20.plot <- plot.OD600.over.time(Feb.01.Tet20.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet20-plot.pdf",
       Feb.01.Tet20.plot)
Feb.01.delta.log2.Tet20.plot <- plot.delta.log2.OD600.over.time(Feb.01.Tet20.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet20-delta-growth.pdf",
       Feb.01.delta.log2.Tet20.plot)

Feb.01.Tet50.plot <- plot.OD600.over.time(Feb.01.Tet50.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet50-plot.pdf",
       Feb.01.Tet50.plot)
Feb.01.delta.log2.Tet50.plot <- plot.delta.log2.OD600.over.time(Feb.01.Tet50.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet50-delta-growth.pdf",
       Feb.01.delta.log2.Tet50.plot)


#######################################
## Plot lag distributions.

Feb.01.LB.lag.df <- make.time.lag.df(Feb.01.LB.tidy.data)
Feb.01.LB.lag.plot <- plot.growth.lag(Feb.01.LB.lag.df)

Feb.01.Tet20.lag.df <- make.time.lag.df(Feb.01.Tet20.tidy.data)
Feb.01.Tet20.lag.plot <- plot.growth.lag(Feb.01.Tet20.lag.df)

Feb.01.Tet50.lag.df <- make.time.lag.df(Feb.01.Tet50.tidy.data)
Feb.01.Tet50.lag.plot <- plot.growth.lag(Feb.01.Tet50.lag.df)


#######################################
## Plot max growth rate distributions.

Feb.01.LB.max.growth.rate.df <- make.max.growth.rate.df(Feb.01.LB.tidy.data)
Feb.01.LB.max.growth.rate.plot <- plot.max.growth.rate(Feb.01.LB.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-LB-growth-rate.pdf",
       Feb.01.LB.max.growth.rate.plot)
Feb.01.LB.max.growth.rate.histogram <- plot.max.growth.rate.histogram(Feb.01.LB.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-LB-growth-rate-histogram.pdf",
       Feb.01.LB.max.growth.rate.histogram)

Feb.01.Tet20.max.growth.rate.df <- make.max.growth.rate.df(Feb.01.Tet20.tidy.data)
Feb.01.Tet20.max.growth.rate.plot <- plot.max.growth.rate(Feb.01.Tet20.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet20-growth-rate.pdf",
       Feb.01.Tet20.max.growth.rate.plot)
Feb.01.Tet20.max.growth.rate.histogram <- plot.max.growth.rate.histogram(Feb.01.Tet20.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet20-growth-rate-histogram.pdf", Feb.01.Tet20.max.growth.rate.histogram)


Feb.01.Tet50.max.growth.rate.df <- make.max.growth.rate.df(Feb.01.Tet50.tidy.data)
Feb.01.Tet50.max.growth.rate.plot <- plot.max.growth.rate(Feb.01.Tet50.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet50-growth-rate.pdf",
       Feb.01.Tet20.max.growth.rate.plot)
Feb.01.Tet50.max.growth.rate.histogram <- plot.max.growth.rate.histogram(Feb.01.Tet50.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet50-growth-rate-histogram.pdf", Feb.01.Tet50.max.growth.rate.histogram)

########################################
## compare growth rates to lags.

Feb.01.LB.lag.vs.growth.rate.plot <- plot.growth.rate.vs.lag(
    Feb.01.LB.lag.df,
    Feb.01.LB.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-LB-growth-rate-vs-lag.pdf",
       Feb.01.LB.lag.vs.growth.rate.plot)

Feb.01.Tet20.lag.vs.growth.rate.plot <- plot.growth.rate.vs.lag(
    Feb.01.Tet20.lag.df,
    Feb.01.Tet20.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet20-growth-rate-vs-lag.pdf",
       Feb.01.Tet20.lag.vs.growth.rate.plot)

Feb.01.Tet50.lag.vs.growth.rate.plot <- plot.growth.rate.vs.lag(
    Feb.01.Tet50.lag.df,
    Feb.01.Tet50.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet50-growth-rate-vs-lag.pdf",
       Feb.01.Tet50.lag.vs.growth.rate.plot)

########################################
## plot OD600 at 24 hour distributions.

Feb.01.LB.24h.OD600.histogram <- plot.24h.OD600.histogram(Feb.01.LB.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-LB-24h-OD600-histogram.pdf", Feb.01.LB.24h.OD600.histogram)

Feb.01.Tet20.24h.OD600.histogram <- plot.24h.OD600.histogram(Feb.01.Tet20.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet20-24h-OD600-histogram.pdf", Feb.01.Tet20.24h.OD600.histogram)

Feb.01.Tet50.24h.OD600.histogram <- plot.24h.OD600.histogram(Feb.01.Tet50.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-01-Tet50-24h-OD600-histogram.pdf", Feb.01.Tet50.24h.OD600.histogram)

########################################################
## B30 clones, Replicate 3.
########################################################
## Import and tidy the data.

Feb.03.LB.evolved.clone.long.format.data <- read.csv(
    "../data/draft-manuscript-1A/OD600-data/2022-02-03-LB-B30-clones.csv", header=FALSE)
Feb.03.LB.tidy.data <- Feb.03.LB.evolved.clone.long.format.data %>%
    tidy.tecan.clone.plate.data(0)

Feb.03.Tet20.evolved.clone.long.format.data <- read.csv(
    "../data/draft-manuscript-1A/OD600-data/2022-02-03-LB-Tet20-B30-clones.csv",
    header=FALSE)
Feb.03.Tet20.tidy.data <- Feb.03.Tet20.evolved.clone.long.format.data %>%
    tidy.tecan.clone.plate.data(20)

Feb.03.Tet50.evolved.clone.long.format.data <- read.csv(
    "../data/draft-manuscript-1A/OD600-data/2022-02-03-LB-Tet50-B30-clones.csv",
    header=FALSE)
Feb.03.Tet50.tidy.data <- Feb.03.Tet50.evolved.clone.long.format.data %>%
    tidy.tecan.clone.plate.data(50)

write.csv(x=Feb.03.LB.tidy.data, file="../results/draft-manuscript-1A/growth-results/2022-02-03-LB-Tet0-B30-clones-tidy.csv", row.names=FALSE, quote = FALSE)
write.csv(x=Feb.03.Tet20.tidy.data, file="../results/draft-manuscript-1A/growth-results/2022-02-03-LB-Tet20-B30-clones-tidy.csv", row.names=FALSE, quote = FALSE)
write.csv(x=Feb.03.Tet50.tidy.data, file="../results/draft-manuscript-1A/growth-results/2022-02-03-LB-Tet50-B30-clones-tidy.csv", row.names=FALSE, quote = FALSE)


############################
## Plot growth over time.

Feb.03.LB.plot <- plot.OD600.over.time(Feb.03.LB.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-LB-plot.pdf",
       Feb.03.LB.plot)
Feb.03.delta.log2.LB.plot <- plot.delta.log2.OD600.over.time(Feb.03.LB.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-LB-delta-growth.pdf",
       Feb.03.delta.log2.LB.plot)

Feb.03.Tet20.plot <- plot.OD600.over.time(Feb.03.Tet20.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet20-plot.pdf",
       Feb.03.Tet20.plot)
Feb.03.delta.log2.Tet20.plot <- plot.delta.log2.OD600.over.time(Feb.03.Tet20.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet20-delta-growth.pdf",
       Feb.03.delta.log2.Tet20.plot)

Feb.03.Tet50.plot <- plot.OD600.over.time(Feb.03.Tet50.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet50-plot.pdf",
       Feb.03.Tet50.plot)
Feb.03.delta.log2.Tet50.plot <- plot.delta.log2.OD600.over.time(Feb.03.Tet50.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet50-delta-growth.pdf",
       Feb.03.delta.log2.Tet50.plot)

#######################################
## Plot lag distributions.

Feb.03.LB.lag.df <- make.time.lag.df(Feb.03.LB.tidy.data)
Feb.03.LB.lag.plot <- plot.growth.lag(Feb.03.LB.lag.df)

Feb.03.Tet20.lag.df <- make.time.lag.df(Feb.03.Tet20.tidy.data)
Feb.03.Tet20.lag.plot <- plot.growth.lag(Feb.03.Tet20.lag.df)

Feb.03.Tet50.lag.df <- make.time.lag.df(Feb.03.Tet50.tidy.data)
Feb.03.Tet50.lag.plot <- plot.growth.lag(Feb.03.Tet50.lag.df)

#######################################
## Plot max growth rate distributions.

Feb.03.LB.max.growth.rate.df <- make.max.growth.rate.df(Feb.03.LB.tidy.data)
Feb.03.LB.max.growth.rate.plot <- plot.max.growth.rate(Feb.03.LB.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-LB-growth-rate.pdf",
       Feb.03.LB.max.growth.rate.plot)
Feb.03.LB.max.growth.rate.histogram <- plot.max.growth.rate.histogram(Feb.03.LB.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-LB-growth-rate-histogram.pdf",
       Feb.03.LB.max.growth.rate.histogram)

Feb.03.Tet20.max.growth.rate.df <- make.max.growth.rate.df(Feb.03.Tet20.tidy.data)
Feb.03.Tet20.max.growth.rate.plot <- plot.max.growth.rate(Feb.03.Tet20.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet20-growth-rate.pdf",
       Feb.03.Tet20.max.growth.rate.plot)
Feb.03.Tet20.max.growth.rate.histogram <- plot.max.growth.rate.histogram(Feb.03.Tet20.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet20-growth-rate-histogram.pdf", Feb.03.Tet20.max.growth.rate.histogram)


Feb.03.Tet50.max.growth.rate.df <- make.max.growth.rate.df(Feb.03.Tet50.tidy.data)
Feb.03.Tet50.max.growth.rate.plot <- plot.max.growth.rate(Feb.03.Tet50.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet50-growth-rate.pdf",
       Feb.03.Tet20.max.growth.rate.plot)
Feb.03.Tet50.max.growth.rate.histogram <- plot.max.growth.rate.histogram(Feb.03.Tet50.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet50-growth-rate-histogram.pdf", Feb.03.Tet50.max.growth.rate.histogram)

########################################
## compare growth rates to lags.

Feb.03.LB.lag.vs.growth.rate.plot <- plot.growth.rate.vs.lag(
    Feb.03.LB.lag.df,
    Feb.03.LB.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-LB-growth-rate-vs-lag.pdf",
       Feb.03.LB.lag.vs.growth.rate.plot)

Feb.03.Tet20.lag.vs.growth.rate.plot <- plot.growth.rate.vs.lag(
    Feb.03.Tet20.lag.df,
    Feb.03.Tet20.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet20-growth-rate-vs-lag.pdf",
       Feb.03.Tet20.lag.vs.growth.rate.plot)

Feb.03.Tet50.lag.vs.growth.rate.plot <- plot.growth.rate.vs.lag(
    Feb.03.Tet50.lag.df,
    Feb.03.Tet50.max.growth.rate.df)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet50-growth-rate-vs-lag.pdf",
       Feb.03.Tet50.lag.vs.growth.rate.plot)

########################################
## plot OD600 at 24 hour distributions.

Feb.03.LB.24h.OD600.histogram <- plot.24h.OD600.histogram(Feb.03.LB.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-LB-24h-OD600-histogram.pdf", Feb.03.LB.24h.OD600.histogram)

Feb.03.Tet20.24h.OD600.histogram <- plot.24h.OD600.histogram(Feb.03.Tet20.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet20-24h-OD600-histogram.pdf", Feb.03.Tet20.24h.OD600.histogram)

Feb.03.Tet50.24h.OD600.histogram <- plot.24h.OD600.histogram(Feb.03.Tet50.tidy.data)
ggsave("../results/draft-manuscript-1A/growth-results/Feb-03-Tet50-24h-OD600-histogram.pdf", Feb.03.Tet50.24h.OD600.histogram)


