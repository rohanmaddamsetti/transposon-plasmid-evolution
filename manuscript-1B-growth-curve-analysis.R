## manuscript-1B-growth-curve-analysis.R by Rohan Maddamsetti.
## This script plots growth over time.

library(tidyverse)
library(cowplot)
library(broom)

## for OD600 measurement analysis of 9 day evolution experiments.
subtract.OD600.blanks <- function(OD600.df) {

    blanks.df <- filter(OD600.df, Treatment == "Blank")
    media.blank <- mean(blanks.df$RawOD600)
    subtracted.df <- OD600.df %>%
        filter(Treatment != "Blank") %>%
        mutate(OD600 = RawOD600 - media.blank)
    return(subtracted.df)
}

############################################################
## 7/11/2021. I made a plot from
## OD600 measurements from 9 day B51/B75 evolution experiment.

july.1.B51.B75.OD600.df <- read.csv("../data/draft-manuscript-1B/OD600-data/2021-B51-B75-9-day-OD600.csv") %>%
    ## subtract LB blanks to get the true OD measurement.
    subtract.OD600.blanks() %>%
    mutate(BiologicalReplicate = as.factor(BiologicalReplicate))
    
july.1.B51.B75.OD600.plot <- ggplot(july.1.B51.B75.OD600.df,
                             aes(x = Day,
                                 y = OD600,
                                 ##shape = BiologicalReplicate,
                                 color = Treatment)) +                           
    theme_classic() + geom_point() + geom_smooth()

ggsave("../results/draft-manuscript-1B/growth-results/OD600-2021-B51-B75-9-day-exp.pdf", july.1.B51.B75.OD600.plot)

##################################################################################
## for DH5a + B51 and B75 strains,
## under the following conditions:
## Tet 0, 2, 4, 6, 8, 10, 20, 30, 40, 50.

B51.B75.long.format.dose.response.data <- read.csv("../data/draft-manuscript-1B/OD600-data/2021-06-19-B51-B75-dose-response-OD600.csv", header=FALSE)

well.to.Tet.2021.06.19 <- function(well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the tetracycline concentration.
    colnum <- as.numeric(substring(well,2))
    tet.conc.map.vec <- c(NA,0,2,4,6,8,10,20,30,40,50)
    return(tet.conc.map.vec[colnum])
}

well.to.strain.2021.06.19 <- function(well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the technical replicate.
    rowletter <- substring(well,1,1)
    strain <- NA
    if (rowletter %in% c("A", "B", "C", "D")) strain <- "B51"
    if (rowletter %in% c("E", "F", "G", "H")) strain <- "B75"
    return(strain)
}

well.to.rep.2021.06.19 <- function(well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the technical replicate.
    rowletter <- substring(well,1,1)
    rep <- NA
    
    if (rowletter == "A" || rowletter == "E") rep <- 1
    if (rowletter == "B" || rowletter == "F") rep <- 2
    if (rowletter == "C" || rowletter == "G") rep <- 3
    if (rowletter == "D" || rowletter == "H") rep <- 4
    return(rep)
}

B51.B75.tidy.dose.response.data <- B51.B75.long.format.dose.response.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = as.factor(well.to.Tet.2021.06.19(Well))) %>%
    mutate(strain = sapply(Well, well.to.strain.2021.06.19)) %>%
    mutate(replicate = sapply(Well, well.to.rep.2021.06.19))
                 
B51.B75.dose.response.plot <- ggplot(B51.B75.tidy.dose.response.data,
                             aes(x = hours,
                                 y = OD600,
                                 color = strain)) +
    geom_point(size=0.5) +
    facet_wrap(.~Tet)

## surprising! The B51 strain without the mini-Tn5 transposon is more resistant
## than B75, which has the transposon!
## maybe the presence of the transposon suppresses TetR expression compared
## to the control plasmid?
ggsave("../results/draft-manuscript-1B/growth-results/tet-B51-B75-dose-response.pdf", B51.B75.dose.response.plot)

####################################################################################
## 7/15/2021. Compare Tet50 growth of B51 and B75 evolved populations.

well.to.rep.2021.07.15 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the clone replicate.
    rowletter <- substring(Well,1,1)
    return(rowletter)
}

well.to.pop.2021.07.15 <- function(Well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the biological replicate.
    colnum <- as.numeric(substring(Well,2))
    pop.map.vec <- rep(c(1,2,3,4,5),2)
    return(pop.map.vec[colnum])
}

well.to.plasmid.2021.07.15 <- function(Well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the plasmid.
    colnum <- as.numeric(substring(Well,2))
    map.vec <- c(rep("B51",5),rep("B75",5))
    return(map.vec[colnum])
}

B51.B75.evolved.long.format.data <- read.csv(
    "../data/draft-manuscript-1B/OD600-data/2021-07-15-B51-B75-9-day-pop-Tet50-growth.csv", header=FALSE)

B51.B75.evolved.tidy.data <- B51.B75.evolved.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 50) %>%
    mutate(pop = sapply(Well, well.to.pop.2021.07.15)) %>%
    mutate(pop = as.factor(pop)) %>%
    mutate(replicate = sapply(Well, well.to.rep.2021.07.15)) %>%
    mutate(plasmid = sapply(Well, well.to.plasmid.2021.07.15))

B51.B75.evolved.pop.plot1 <- ggplot(
    B51.B75.evolved.tidy.data,
    aes(x = hours,
        y = OD600,
        shape = pop,
        color = plasmid)) + geom_point(size=0.5)

ggsave("../results/draft-manuscript-1B/growth-results/B51-B75-Tet50-evolved-pops.pdf",
       B51.B75.evolved.pop.plot1)

## let's plot the distribution of OD at 24h.
B51.B75.evolved.pop.tidy.24h.data <- B51.B75.evolved.tidy.data %>%
    filter(hours == 24)

B51.B75.evolved.pop.24h.OD600.plot <- ggplot(
    B51.B75.evolved.pop.tidy.24h.data,
    aes(x = plasmid,
        y = OD600,
        shape = replicate,
        color = pop)) + geom_point() + theme_classic() +
    ylab("OD600 at 24h timepoint") +
    xlab("Presence/absence of transposase on plasmid")

ggsave("../results/draft-manuscript-1B/growth-results/B51-B75-Tet50-evolved-pops-24h-OD600.pdf",
       B51.B75.evolved.pop.24h.OD600.plot)

## now, calculate the max growth rate (look at points between OD600 = 0.3-0.4).
B51.B75.Tet50.evolved.pop.tidy.max.growth.rate.data <- B51.B75.evolved.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, pop, replicate, plasmid) %>%
    nest() %>% ## nest the data by well, then calculate the OD slope using lm().
    mutate(fit = map(data, ~lm(OD600 ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well,pop,replicate,plasmid, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

B51.B75.Tet50.evolved.pop.24h.max.growth.rate.plot <- ggplot(
    B51.B75.Tet50.evolved.pop.tidy.max.growth.rate.data,
    aes(x = plasmid,
        y = estimate,
        shape = replicate,
        color = pop)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of tranposase on plasmid")

ggsave("../results/draft-manuscript-1B/growth-results/B51-B75-Tet50-evolved-pops-growth-rate.pdf",
       B51.B75.Tet50.evolved.pop.24h.max.growth.rate.plot)

summarize.time.lag <- function(well.df) {
    ## to compare time lags, measure point when OD600 hits 0.25.
    lag.min.index <- max(min(which(well.df$OD600 > 0.25)), 1)
    t.OD.hit.lag.min <- well.df[lag.min.index,]$hours

    ret.df <- mutate(well.df, time.lag = t.OD.hit.lag.min)
    return(ret.df)
}

## now, calculate the time lag (time to reach OD600 = 0.25).
B51.B75.Tet50.evolved.pop.tidy.lagtime.data <- B51.B75.evolved.tidy.data %>%
    split(.$Well) %>%
    map_dfr(.f=summarize.time.lag) %>%
    group_by(pop, plasmid) %>%
    ## clunky hack to so that time.lag is not shown for every single time point!
    ## TODO: come up with a simpler solution.
    ## TODO: DOUBLE CHECK WHETHER I SHOULD KEEP REPLICATE OR NOT!
    summarize(mean.time.lag = mean(time.lag))
    
B51.B75.Tet50.evolved.pop.timelag.plot <- ggplot(
    B51.B75.Tet50.evolved.pop.tidy.lagtime.data,
    aes(x = plasmid,
        y = mean.time.lag,
        color = pop)) + geom_point() + theme_classic() +
    ylab("Time to reach OD600 = 0.25") +
    xlab("Presence/absence of transposase on plasmid")

ggsave("../results/draft-manuscript-1B/growth-results/B51-B75-Tet50-evolved-pops-time-lag.pdf",
       B51.B75.Tet50.evolved.pop.timelag.plot)

## let's actually calculate statistics to compare these growth parameters.
B51.timelag.data <- filter(B51.B75.Tet50.evolved.pop.tidy.lagtime.data,
                           plasmid == "B51")
B75.timelag.data <- filter(B51.B75.Tet50.evolved.pop.tidy.lagtime.data,
                           plasmid == "B75")

## no significant difference in growth performance, overall.
wilcox.test(B51.timelag.data$mean.time.lag,
       B75.timelag.data$mean.time.lag)

combined.B51.B75.plot <- plot_grid(
    B51.B75.Tet50.evolved.pop.timelag.plot,
    B51.B75.Tet50.evolved.pop.24h.max.growth.rate.plot,
    B51.B75.evolved.pop.24h.OD600.plot,
    labels=c('A','B','C'))

ggsave("../results/draft-manuscript-1B/growth-results/B51-B75-Tet50-evolved-pops-combined-growth-plot.pdf",
       combined.B51.B75.plot)

########################################################################
