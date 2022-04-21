## manuscript-1A-analysis-of-older-data.R by Rohan Maddamsetti.
## This script contains analyses of data that will probably
## never be published, but that are worth keeping on hand
## (both for the data, as well as for code recipes)
## for now.

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

### Make a plot of the Tet antibiotic treatment over time.
experiment.design.df <- data.frame(Day = c(1,2,3,4,5,6,7,8,9),
                                   Tet.conc = c(2,4,6,8,10,20,30,40,50))
treatment.plot <- ggplot(experiment.design.df,
                         aes(x = Day, y = Tet.conc)) +
    theme_classic() + geom_line() +
    ylab("Tetracycline (ug/mL)") +
    xlab("Day")

## make a plot of OD600 over time in the nine day experiment.
april.14.OD600.df <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-04-14-B30-OD600.csv") %>%
    ## subtract LB blanks to get the true OD measurement.
    subtract.OD600.blanks() %>%
    mutate(BiologicalReplicate = as.factor(BiologicalReplicate))
    
april.14.OD600.plot <- ggplot(april.14.OD600.df,
                              aes(x = Day,
                                  y = OD600,
                                  ##shape = BiologicalReplicate,
                                  color = Treatment)) +                                
    theme_classic() + geom_point() + geom_smooth()
ggsave("../results/draft-manuscript-1A/growth-results/OD600-2021-04-14.pdf", april.14.OD600.plot)

##########################
## 6/15/2021. I made a plot using
## OD600 measurements from 9 day A18-B30 evolution experiment.

june.15.OD600.df <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-A18-B30-9-day-OD600.csv") %>%
    ## subtract LB blanks to get the true OD measurement.
    subtract.OD600.blanks() %>%
    mutate(BiologicalReplicate = as.factor(BiologicalReplicate))
    
june.15.OD600.plot <- ggplot(june.15.OD600.df,
                             aes(x = Day,
                                 y = OD600,
                                 ##shape = BiologicalReplicate,
                                 color = Treatment)) +                           
    theme_classic() + geom_point() + geom_smooth()

ggsave("../results/draft-manuscript-1A/growth-results/OD600-2021-06-15.pdf", june.15.OD600.plot)

############################################################
## 7/11/2021. I made a plot from
## OD600 measurements from 9 day A31 evolution experiment.

july.1.A31.OD600.df <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-A31-9-day-OD600.csv") %>%
    ## subtract LB blanks to get the true OD measurement.
    subtract.OD600.blanks() %>%
    mutate(BiologicalReplicate = as.factor(BiologicalReplicate))
    
july.1.A31.OD600.plot <- ggplot(july.1.A31.OD600.df,
                             aes(x = Day,
                                 y = OD600,
                                 shape = BiologicalReplicate,
                                 color = Treatment)) +                           
    theme_classic() + geom_point() + geom_smooth()

ggsave("../results/draft-manuscript-1A/growth-results/OD600-2021-A31-9-day-exp.pdf", july.1.A31.OD600.plot)

############################################################
## for DH5a + B30 initial strain,
## under the following conditions:
## Tet 0, 2, 4, 6, 8, 10, 20, 30, 40, 50.

B30.long.format.dose.response.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-04-27_Rohan-DH5a-B30-dose-response-OD600.csv", header=FALSE)

well.to.Tet.2021.04.27 <- function(well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the tetracycline concentration.
    colnum <- as.numeric(substring(well,2))
    tet.conc.map.vec <- c(NA,0,2,4,6,8,10,20,30,40,50)
    return(tet.conc.map.vec[colnum])
}

well.to.rep.2021.04.27 <- function(well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the technical replicate.
    rowletter <- substring(well,1,1)
    if (rowletter == "B") return(1)
    if (rowletter == "C") return(2)
    if (rowletter == "D") return(3)
    if (rowletter == "E") return(4)
    return(NA)
}

B30.tidy.dose.response.data <- B30.long.format.dose.response.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = as.factor(well.to.Tet.2021.04.27(Well))) %>%
    mutate(replicate = sapply(Well, well.to.rep.2021.04.27))
                 
B30.dose.response.plot <- ggplot(B30.tidy.dose.response.data,
                             aes(x = hours,
                                 y = OD600,
                                 color = Tet)) +
    geom_point(size=0.5) +
    facet_wrap(.~Tet)

ggsave("../results/draft-manuscript-1A/growth-results/tet-B30-dose-response.pdf", B30.dose.response.plot)

#######################################################################################
## for DH5a + B30 miniTn5 transposon. Comparison between -/+ A18 pUC plasmid,
## under the following conditions:
## Tet 0, 2, 4, 6, 8, 10, 20, 30, 40, 50.

## Plate 1 from Tecan J. I have not analyzed Plate 2 from Tecan K-- looks too different
## to compare. I think experimental design should stick to one machine.

A18.plate.1.long.format.dose.response.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-06-20_dose-response-plate-1-Tecan-J.csv", header=FALSE)

well.to.Tet.2021.06.20.plate.1 <- function(well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the tetracycline concentration.
    ## In this experiment, measurements started from column 1.
    colnum <- as.numeric(substring(well,2))
    tet.conc.map.vec <- c(0,2,4,6,8,10,20,30,40,50, NA, NA)
    return(tet.conc.map.vec[colnum])
}

well.to.strain.2021.06.20.plate.1 <- function(well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the technical replicate.
    rowletter <- substring(well,1,1)
    strain <- NA
    if (rowletter %in% c("A", "C", "E")) strain <- "A18_pUC_plasmid-"
    if (rowletter %in% c("B", "D", "F")) strain <- "A18_pUC_plasmid+"
    return(strain)
}

well.to.rep.2021.06.20.plate.1 <- function(well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the technical replicate.
    rowletter <- substring(well,1,1)
    rep <- NA
    
    if (rowletter == "A" || rowletter == "B") rep <- 1
    if (rowletter == "C" || rowletter == "D") rep <- 2
    if (rowletter == "E" || rowletter == "F") rep <- 3
    return(rep)
}

A18.plate.1.tidy.dose.response.data <- A18.plate.1.long.format.dose.response.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = as.factor(well.to.Tet.2021.06.20.plate.1(Well))) %>%
    mutate(strain = sapply(Well, well.to.strain.2021.06.20.plate.1)) %>%
    mutate(replicate = sapply(Well, well.to.rep.2021.06.20.plate.1)) %>%
    filter(!is.na(Tet))
                 
A18.plate.1.dose.response.plot <- ggplot(A18.plate.1.tidy.dose.response.data,
                             aes(x = hours,
                                 y = OD600,
                                 color = strain)) +
    geom_point(size=0.5) +
    facet_grid(.~Tet)

## surprising! The B51 strain without the mini-Tn5 transposon is more resistant
## than B75, which has the transposon!
## maybe the presence of the transposon suppresses TetR expression compared
## to the control plasmid?
ggsave("../results/draft-manuscript-1A/growth-results/tet-A18-plate-1-dose-response.pdf", A18.plate.1.dose.response.plot)

#######################################################################################
## for DH5a + B30 miniTn5 transposon. Comparison of growth in Tet50 between pops evolved for ~90 generations
## in LB (relaxed selection), compared to ancestor.

well.to.strain.2021.06.29 <- function(Well) {
    ## This helper function maps the columns of the well in
    ## the 96-well plate to the strain.
    col <- as.numeric(substring(Well,2))
    strain <- NA
    if (col == 1) strain <- "ancestor A"
    if (col == 2) strain <- "ancestor B"
    if (col == 3) strain <- "evolved A1"
    if (col == 4) strain <- "evolved A2"
    if (col == 5) strain <- "evolved A3"
    if (col == 6) strain <- "evolved A4"
    if (col == 7) strain <- "evolved A5"
    if (col == 8) strain <- "evolved B1"
    if (col == 9) strain <- "evolved B2"
    if (col == 10) strain <- "evolved B3"
    if (col == 11) strain <- "evolved B4"
    if (col == 12) strain <- "evolved B5"
    
    return(strain)
}

well.to.rep.2021.06.29 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the technical replicate.
    rowletter <- substring(Well,1,1)
    rep <- NA
    if (rowletter == "A") rep <- 1
    if (rowletter == "B") rep <- 2
    if (rowletter == "C") rep <- 3
    if (rowletter == "D") rep <- 4
    if (rowletter == "E") rep <- 5
    return(rep)
}


no.plasmid.B30.relaxed.selection.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-06-29-B30-clone-relaxed-selection-Tet50-growth.csv", header=FALSE)

no.plasmid.B30.relaxed.selection.tidy.data <-
    no.plasmid.B30.relaxed.selection.long.format.data %>%
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
    mutate(strain = sapply(Well, well.to.strain.2021.06.29)) %>%
    mutate(replicate = sapply(Well, well.to.rep.2021.06.29)) %>%
    filter(!is.na(replicate)) %>%
    mutate(is.ancestor = ifelse(str_detect(strain, "ancestor"), "ancestor", "evolved")) %>%
    mutate(rowletter = substring(Well,1,1)) %>%
    mutate(col = as.numeric(substring(Well,2))) %>%
    mutate(clone = ifelse(str_detect(strain, "A"), "A", "B"))



no.plasmid.B30.relaxed.selection.plot1 <- ggplot(no.plasmid.B30.relaxed.selection.tidy.data,
                                                 aes(x = hours,
                                                     y = OD600,
                                                     color = strain,
                                                     shape=clone)) + geom_point(size=0.5)

no.plasmid.B30.relaxed.selection.plot2 <- no.plasmid.B30.relaxed.selection.plot1 +
    facet_wrap(is.ancestor~clone)

no.plasmid.B30.relaxed.selection.plot1 <- no.plasmid.B30.relaxed.selection.plot1 +
    facet_wrap(.~clone)

## surprising! After 9 days of relaxed selection, the evolved strains are even more
## resistant to Tet, despite no selection for Tet resistance!
ggsave("../results/draft-manuscript-1A/growth-results/no-plasmid-B30-relaxed-selection-1.pdf",
       no.plasmid.B30.relaxed.selection.plot1)

ggsave("../results/draft-manuscript-1A/growth-results/no-plasmid-B30-relaxed-selection-2.pdf",
       no.plasmid.B30.relaxed.selection.plot2)

## one possible mechanism: miniTn5 insertion into genes that are under selection for being
## knocked out would 1) be beneficial in LB 2) cause higher Tet resistance.

##########################################################################

A18.plasmid.B30.relaxed.selection.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-07-10-B30-A18-pop-relaxed-selection-Tet50-growth.csv", header=FALSE)

well.to.rep.2021.07.10 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the technical replicate.
    rowletter <- substring(Well,1,1)
    if (rowletter == "A") return(1)
    if (rowletter == "B") return(2)
    if (rowletter == "C") return(3)
    if (rowletter == "D") return(4)
    return(NA)
}

well.to.pop.2021.07.10 <- function(Well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the biological replicate.
    colnum <- as.numeric(substring(Well,2))
    pop.map.vec <- rep(c(1,2,3,4,5),2)
    return(pop.map.vec[colnum])
}

well.to.anc.or.evol.2021.07.10 <- function(Well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the biological replicate.
    colnum <- as.numeric(substring(Well,2))
    anc.evol.map.vec <- c(rep("Anc",5),rep("Evol",5))
    return(anc.evol.map.vec[colnum])
}


A18.plasmid.B30.relaxed.selection.tidy.data <-
    A18.plasmid.B30.relaxed.selection.long.format.data %>%
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
    mutate(pop = sapply(Well, well.to.pop.2021.07.10)) %>%
    mutate(replicate = sapply(Well, well.to.rep.2021.07.10)) %>%
    filter(!is.na(replicate)) %>%
    mutate(anc.or.evol = sapply(Well, well.to.anc.or.evol.2021.07.10)) %>%
    mutate(rowletter = substring(Well,1,1)) %>%
    mutate(col = as.numeric(substring(Well,2)))

A18.plasmid.B30.relaxed.selection.plot <- ggplot(
    A18.plasmid.B30.relaxed.selection.tidy.data,
    aes(x = hours,
        y = OD600,
        color = anc.or.evol)) + geom_point(size=0.5) +
    facet_grid(.~pop)

ggsave("../results/draft-manuscript-1A/growth-results/A18-plasmid-relaxed-selection-Tet50-OD600.pdf", A18.plasmid.B30.relaxed.selection.plot)

####################################################################################

## 7/14/2021. initial data for Tet50 growth of isolated clones.

well.to.clone.2021.07.13 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the clone replicate.
    rowletter <- substring(Well,1,1)
    return(rowletter)
}

well.to.pop.2021.07.13 <- function(Well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the biological replicate.
    colnum <- as.numeric(substring(Well,2))
    pop.map.vec <- rep(c(1,2,3,4,5),2)
    return(pop.map.vec[colnum])
}

well.to.minus.plus.plasmid.2021.07.13 <- function(Well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the biological replicate.
    colnum <- as.numeric(substring(Well,2))
    minus.plus.map.vec <- c(rep("A18_pUC_plasmid-",5),rep("A18_pUC_plasmid+",5))
    return(minus.plus.map.vec[colnum])
}

Tet50.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-07-13-evolved-clones-Tet50-dose-OD600-response.csv", header=FALSE)

Tet50.evolved.clone.tidy.data <- Tet50.evolved.clone.long.format.data %>%
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
    mutate(pop = sapply(Well, well.to.pop.2021.07.13)) %>%
    mutate(pop = as.factor(pop)) %>%
    mutate(clone = sapply(Well, well.to.clone.2021.07.13)) %>%
    mutate(plasmid.minus.or.plus = sapply(Well, well.to.minus.plus.plasmid.2021.07.13))

Tet50.evolved.clones.plot <- ggplot(
    Tet50.evolved.clone.tidy.data,
    aes(x = hours,
        y = OD600,
        shape = clone,
        color = plasmid.minus.or.plus)) + geom_point(size=0.5)

## really nice result! The best clones without the plasmid out-compete the clones
## with the plasmid!
ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-clones.pdf", Tet50.evolved.clones.plot)


## let's plot the distribution of OD at 24h..
Tet50.evolved.clone.tidy.24h.data <- Tet50.evolved.clone.tidy.data %>%
    filter(hours == 24)

Tet50.evolved.clones.24h.OD600.plot <- ggplot(
    Tet50.evolved.clone.tidy.24h.data,
    aes(x = plasmid.minus.or.plus,
        y = OD600,
        shape = clone,
        color = pop)) + geom_point() + theme_classic() +
    ylab("OD600 at 24h timepoint") +
    xlab("Presence/absence of high-copy pUC target plasmid")

ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-clones-24h-OD600.pdf",
       Tet50.evolved.clones.24h.OD600.plot)

## now, calculate the max growth rate (look at points between OD600 = 0.3-0.4).
Tet50.evolved.clone.tidy.max.growth.rate.data <- Tet50.evolved.clone.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, pop, clone, plasmid.minus.or.plus) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well,pop,clone,plasmid.minus.or.plus, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Tet50.evolved.clones.24h.max.growth.rate.plot <- ggplot(
    Tet50.evolved.clone.tidy.max.growth.rate.data,
    aes(x = plasmid.minus.or.plus,
        y = estimate,
        shape = clone,
        color = pop)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid")

ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-clones-growth-rate.pdf",
       Tet50.evolved.clones.24h.max.growth.rate.plot)

## let's calculate some statistics.
## the proper unit of replication is the population.
## it actually may make more sense to calculate these statistics
## on the mixed populations themselves.

Tet50.evolved.clone.24h.OD600.summary <- Tet50.evolved.clone.tidy.data %>%
    group_by(pop, plasmid.minus.or.plus) %>%
    summarize(mean_OD600 = mean(OD600), var_OD600 = var(OD600))

## significant difference in OD600 at 24h.
wilcox.test(x=filter(Tet50.evolved.clone.24h.OD600.summary,
                     plasmid.minus.or.plus == 'A18_pUC_plasmid-')$mean_OD600,
            y=filter(Tet50.evolved.clone.24h.OD600.summary,
                     plasmid.minus.or.plus == 'A18_pUC_plasmid+')$mean_OD600)

## no significance difference in variance in OD600 across populations across
## treatments.
wilcox.test(x=filter(Tet50.evolved.clone.24h.OD600.summary,
                     plasmid.minus.or.plus == 'A18_pUC_plasmid-')$var_OD600,
            y=filter(Tet50.evolved.clone.24h.OD600.summary,
                     plasmid.minus.or.plus == 'A18_pUC_plasmid+')$var_OD600)

Tet50.evolved.clone.max.growth.rate.summary <-
    Tet50.evolved.clone.tidy.max.growth.rate.data %>%
    group_by(pop, plasmid.minus.or.plus) %>%
    summarize(mean_growth_rate = mean(estimate), var_growth_rate = var(estimate))

## no significant difference in growth rates (BUT note that data is bi-modal).
wilcox.test(x=filter(Tet50.evolved.clone.max.growth.rate.summary,
                     plasmid.minus.or.plus == 'A18_pUC_plasmid-')$mean_growth_rate,
            y=filter(Tet50.evolved.clone.max.growth.rate.summary,
                     plasmid.minus.or.plus == 'A18_pUC_plasmid+')$mean_growth_rate)

## no significant difference in variance in growth rates within pops
wilcox.test(x=filter(Tet50.evolved.clone.max.growth.rate.summary,
                     plasmid.minus.or.plus == 'A18_pUC_plasmid-')$var_growth_rate,
            y=filter(Tet50.evolved.clone.max.growth.rate.summary,
                     plasmid.minus.or.plus == 'A18_pUC_plasmid+')$var_growth_rate)

## replot these figures as histograms, per Lingchong's request.

Tet50.evolved.clones.24h.OD600.histogram <- ggplot(
    Tet50.evolved.clone.tidy.24h.data,
    aes(fill = plasmid.minus.or.plus,
        x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(plasmid.minus.or.plus~.) +
    ylab("OD600 at 24h timepoint") +
    labs(fill = "Presence/absence of high-copy pUC target plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-clones-24h-OD600-histogram.pdf",
       Tet50.evolved.clones.24h.OD600.histogram)

Tet50.evolved.clones.24h.max.growth.rate.histogram <- ggplot(
    Tet50.evolved.clone.tidy.max.growth.rate.data,
    aes(fill = plasmid.minus.or.plus,
        x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(plasmid.minus.or.plus~.) +
    xlab("growth rate between OD600 0.3-0.4") +
    labs(fill = "Presence/absence of high-copy pUC target plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-clones-growth-rate-histogram.pdf",
       Tet50.evolved.clones.24h.max.growth.rate.histogram)

## let's plot and compare variances in these parameters across populations.

####################################################################################

## 10/21/2021. initial data for Tet50 growth of evolved B30 (no plasmid vs. pUC plasmid).

well.to.pop.replicate.2021.10.21 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the pop replicate.
    rowletter <- substring(Well,1,1)
    return(rowletter)
}

well.to.pop.2021.10.21 <- function(Well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the biological replicate.
    colnum <- as.numeric(substring(Well,2))
    pop.map.vec <- c(c("1","2","3","4","5"),"Blank","Blank",c("1","2","3","4","5"))
    return(pop.map.vec[colnum])
}

well.to.minus.plus.plasmid.2021.10.21 <- function(Well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the treatment.
    colnum <- as.numeric(substring(Well,2))
    minus.plus.map.vec <- c(rep("no_plasmid",5),"Blank", "Blank",rep("pUC_plasmid",5))
    return(minus.plus.map.vec[colnum])
}

Tet50.evolved.pops.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-10-21-evolved-pops-Tet50-dose-OD600-response.csv", header=FALSE)

Tet50.evolved.pops.tidy.data <- Tet50.evolved.pops.long.format.data %>%
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
    mutate(pop = sapply(Well, well.to.pop.2021.10.21)) %>%
    mutate(pop = as.factor(pop)) %>%
    mutate(replicate = sapply(Well, well.to.pop.replicate.2021.10.21)) %>%
    mutate(plasmid.minus.or.plus = sapply(Well, well.to.minus.plus.plasmid.2021.10.21)) %>%
    filter(pop != "Blank")

Tet50.evolved.pops.plot <- ggplot(
    Tet50.evolved.pops.tidy.data,
    aes(x = hours,
        y = OD600,
        shape = pop,
        color = plasmid.minus.or.plus)) + geom_point(size=0.5)

## really nice result, but seems VERY different from the growth curves with the clones??
ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-pops.pdf", Tet50.evolved.pops.plot)


## TODO: (change to 24hr once those data are collected).
## let's plot the distribution of OD at 20h.
Tet50.evolved.pops.tidy.24h.data <- Tet50.evolved.pops.tidy.data %>%
    filter(hours == 20)

Tet50.evolved.pops.24h.OD600.plot <- ggplot(
    Tet50.evolved.pops.tidy.24h.data,
    aes(x = plasmid.minus.or.plus,
        y = OD600,
        color = pop)) + geom_point() + theme_classic() +
    ylab("OD600 at 20h timepoint") +
##    ylab("OD600 at 24h timepoint") +
    xlab("Presence/absence of high-copy pUC target plasmid")

ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-pops-24h-OD600.pdf",
       Tet50.evolved.pops.24h.OD600.plot)

## now, calculate the max growth rate (look at points between OD600 = 0.3-0.4).
Tet50.evolved.pops.tidy.max.growth.rate.data <- Tet50.evolved.pops.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, pop, replicate, plasmid.minus.or.plus) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well,pop, replicate, plasmid.minus.or.plus, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Tet50.evolved.pops.24h.max.growth.rate.plot <- ggplot(
    Tet50.evolved.pops.tidy.max.growth.rate.data,
    aes(x = plasmid.minus.or.plus,
        y = estimate,
        color = pop)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid")

ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-pops-growth-rate.pdf",
       Tet50.evolved.pops.24h.max.growth.rate.plot)

## let's calculate some statistics.
## the proper unit of replication is the population.
## it actually may make more sense to calculate these statistics
## on the mixed populations themselves.

Tet50.evolved.pops.24h.OD600.summary <- Tet50.evolved.pops.tidy.data %>%
    group_by(pop, plasmid.minus.or.plus) %>%
    summarize(mean_OD600 = mean(OD600), var_OD600 = var(OD600))

## no difference in OD600 at 24h.
wilcox.test(x=filter(Tet50.evolved.pops.24h.OD600.summary,
                     plasmid.minus.or.plus == 'no_plasmid')$mean_OD600,
            y=filter(Tet50.evolved.pops.24h.OD600.summary,
                     plasmid.minus.or.plus == 'pUC_plasmid')$mean_OD600)

## marginally non-significant difference in variance in OD600 across populations across
## treatments.
wilcox.test(x=filter(Tet50.evolved.pops.24h.OD600.summary,
                     plasmid.minus.or.plus == 'no_plasmid')$var_OD600,
            y=filter(Tet50.evolved.pops.24h.OD600.summary,
                     plasmid.minus.or.plus == 'pUC_plasmid')$var_OD600)

Tet50.evolved.pops.max.growth.rate.summary <-
    Tet50.evolved.pops.tidy.max.growth.rate.data %>%
    group_by(pop, plasmid.minus.or.plus) %>%
    summarize(mean_growth_rate = mean(estimate), var_growth_rate = var(estimate))

## significant difference in growth rates between OD 0.3-0.4.
wilcox.test(x=filter(Tet50.evolved.pops.max.growth.rate.summary,
                     plasmid.minus.or.plus == 'no_plasmid')$mean_growth_rate,
            y=filter(Tet50.evolved.pops.max.growth.rate.summary,
                     plasmid.minus.or.plus == 'pUC_plasmid')$mean_growth_rate)

## no significant difference in variance in growth rates within pops
wilcox.test(x=filter(Tet50.evolved.pops.max.growth.rate.summary,
                     plasmid.minus.or.plus == 'no_plasmid')$var_growth_rate,
            y=filter(Tet50.evolved.pops.max.growth.rate.summary,
                     plasmid.minus.or.plus == 'pUC_plasmid')$var_growth_rate)

## replot these figures as histograms, per Lingchong's request.

Tet50.evolved.pops.24h.OD600.histogram <- ggplot(
    Tet50.evolved.pops.tidy.24h.data,
    aes(fill = plasmid.minus.or.plus,
        x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(plasmid.minus.or.plus~.) +
    ##    ylab("OD600 at 24h timepoint") +
    ylab("OD600 at 20h timepoint") +
    labs(fill = "Presence/absence of high-copy pUC target plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-pops-24h-OD600-histogram.pdf",
       Tet50.evolved.pops.24h.OD600.histogram)

Tet50.evolved.pops.24h.max.growth.rate.histogram <- ggplot(
    Tet50.evolved.pops.tidy.max.growth.rate.data,
    aes(fill = plasmid.minus.or.plus,
        x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(plasmid.minus.or.plus~.) +
    xlab("growth rate between OD600 0.3-0.4") +
    labs(fill = "Presence/absence of high-copy pUC target plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-pops-growth-rate-histogram.pdf",
       Tet50.evolved.pops.24h.max.growth.rate.histogram)

## TODO: plot differences in growth lag!

summarize.time.lag.OD0.2 <- function(well.df) {
    ## to compare time lags, measure point when OD600 hits 0.2
    lag.min.index <- max(min(which(well.df$OD600 > 0.2)), 1)
    t.OD.hit.lag.min <- well.df[lag.min.index,]$hours

    ret.df <- mutate(well.df, time.lag = t.OD.hit.lag.min)
    return(ret.df)
}

## now, calculate the time lag (time to reach OD600 = 0.25).
Tet50.evolved.pops.tidy.lagtime.data <- Tet50.evolved.pops.tidy.data %>%
    split(.$Well) %>%
    map_dfr(.f=summarize.time.lag.OD0.2) %>%
    group_by(pop, plasmid.minus.or.plus) %>%
    ## clunky hack to so that time.lag is not shown for every single time point!
    ## TODO: come up with a simpler solution.
    ## TODO: DOUBLE CHECK WHETHER I SHOULD KEEP REPLICATE OR NOT!
    summarize(mean.time.lag = mean(time.lag))
    
Tet50.evolved.pops.timelag.plot <- ggplot(
    Tet50.evolved.pops.tidy.lagtime.data,
    aes(x = plasmid.minus.or.plus,
        y = mean.time.lag,
        color = pop)) + geom_point() + theme_classic() +
    ylab("Time to reach OD600 = 0.2") +
    xlab("Presence/absence of high-copy pUC target plasmid")

ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-pops-time-lag.pdf",
       Tet50.evolved.pops.timelag.plot)

## let's actually calculate statistics to compare these growth parameters.
no.plasmid.timelag.data <- filter(Tet50.evolved.pops.tidy.lagtime.data,
                           plasmid.minus.or.plus == "no_plasmid")
pUC.plasmid.timelag.data <- filter(Tet50.evolved.pops.tidy.lagtime.data,
                           plasmid.minus.or.plus == "pUC_plasmid")

## no significant difference in growth performance, overall.
wilcox.test(no.plasmid.timelag.data$mean.time.lag,
       pUC.plasmid.timelag.data$mean.time.lag)

combined.Tet50.evolved.pops.growth.parameter.plot <- plot_grid(
    Tet50.evolved.pops.timelag.plot,
    Tet50.evolved.pops.24h.max.growth.rate.plot,
    Tet50.evolved.pops.24h.OD600.plot,
    labels=c('A','B','C'))

ggsave("../results/draft-manuscript-1A/growth-results/Tet50-evolved-pops-combined-growth-plot.pdf",
       combined.Tet50.evolved.pops.growth.parameter.plot)

#######################################################################################

## look at preliminary growth curves for B30+A31 clones, B20 clones, B20+A31 clones,
## and B20+A18 clones.

well.to.treatment.2021.10.26 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the treatment.
    rowletter <- substring(Well,1,1)
    if (rowletter == "A" || rowletter == "B")
        treatment <- "B20+p15A"
    else if (rowletter == "C" || rowletter == "D")
        treatment <- "B20+pUC"
    else if (rowletter == "E" || rowletter == "F")
        treatment <- "B20+no plasmid"
    else
        treatment <- "B30+p15A"
    return(treatment)
}

Oct.26.Tet50.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-10-26-tet50-day9-clones.csv", header=FALSE)


Oct.26.Tet50.evolved.clone.tidy.data <- Oct.26.Tet50.evolved.clone.long.format.data %>%
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
    ## assign different 4 pairs of rows to different (unknown for now) treatments.
    mutate(treatment = sapply(Well, well.to.treatment.2021.10.26))


Oct.26.Tet50.evolved.clone.plot <- ggplot(
    Oct.26.Tet50.evolved.clone.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) + geom_point(size=0.5) + facet_grid(.~treatment)

## These data have exactly the same format as the Oct 26 2021 data.
Oct.27.Tet50.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-10-27-tet50-day9-clones.csv", header=FALSE)

Oct.27.Tet50.evolved.clone.tidy.data <- Oct.27.Tet50.evolved.clone.long.format.data %>%
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
    ## assign different 4 pairs of rows to different (unknown for now) treatments.
    mutate(treatment = sapply(Well, well.to.treatment.2021.10.26))


Oct.27.Tet50.evolved.clone.plot <- ggplot(
    Oct.27.Tet50.evolved.clone.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) + geom_point(size=0.5) + facet_grid(.~treatment)

ggsave("../results/draft-manuscript-1A/growth-results/Oct-26-Tet50-evolved-clones-24h-growth-curves.pdf", Oct.26.Tet50.evolved.clone.plot)

ggsave("../results/draft-manuscript-1A/growth-results/Oct-27-Tet50-evolved-clones-24h-growth-curves.pdf", Oct.27.Tet50.evolved.clone.plot)
## Further analysis of the B20 transposon growth data from Oct 26 and Oct 27.


#######################################################################################
## look at polished growth curves for all B30 clones, from 11/4/2021.

well.to.well.number <- function(Well) {
    ## first map the row and column to a number from 1-96.
    rowletter <- substring(Well,1,1)
    myLetters <- letters[1:26]
    rownum <- match(tolower(rowletter), myLetters)
    colnum <- as.numeric(substring(Well,2))
    wellnum <- 12*(rownum - 1) + colnum
    return(wellnum)
}

well.to.treatment.2021.11.03 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the treatment.
    wellnum <- well.to.well.number(Well)

    if (wellnum <= 25)
        treatment <- "B30+no plasmid"
    else if (wellnum <= 50)
        treatment <- "B30+pUC plasmid"
    else if (wellnum <= 75)
        treatment <- "B30+p15A plasmid"
    else
        treatment <- "Blank"
    
    return(treatment)    
}

well.to.clone.2021.11.03 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the treatment.
    wellnum <- well.to.well.number(Well)
    clone <- (wellnum %% 5)
    if (clone == 0) clone <- 5 ## 1-5 rather than 1,2,3,4,0 encoding.
    return(clone)
}


well.to.pop.2021.11.03 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the population replicate.
    wellnum <- well.to.well.number(Well)
    well.mod.25 <- (wellnum %% 25)
    if (well.mod.25 == 0) ## 25th well.
        well.mod.25 <- 25 
    
    pop <- ((well.mod.25 - 1) %/% 5) + 1
    return(pop)
}


Nov.3.Tet50.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-03-LB-Tet50-Day9-all-B30-clones.csv", header=FALSE)

Nov.3.Tet50.evolved.clone.tidy.data <- Nov.3.Tet50.evolved.clone.long.format.data %>%
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
    ## assign different 4 pairs of rows to different (unknown for now) treatments.
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.03)) %>%
    filter(treatment != "Blank") %>%
    mutate(clone = sapply(Well, well.to.clone.2021.11.03)) %>%
    mutate(pop = sapply(Well, well.to.pop.2021.11.03)) %>%
    mutate(pop = as.factor(pop))


Nov.3.Tet50.evolved.clone.plot <- ggplot(
    Nov.3.Tet50.evolved.clone.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) + geom_point(size=0.5) + facet_grid(pop~treatment)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-3-Tet50-B30-evolved-clones-48h-growth-curves.pdf", Nov.3.Tet50.evolved.clone.plot)

## let's calculate the growth rate between 0.19 and 0.21.
Nov.3.tidy.max.growth.rate.data <- Nov.3.Tet50.evolved.clone.tidy.data %>%
    filter(OD600 >= 0.19) %>%
    filter(OD600 <= 0.21) %>%
    group_by(Well, pop, clone, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well,pop, clone, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.3.growth.rate.plot <- ggplot(
    Nov.3.tidy.max.growth.rate.data,
    aes(x = treatment,
        y = estimate,
        color = pop)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.19-0.21") +
    xlab("treatment")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-3-Tet50-all-B30-clones-growth-rate.pdf",
       Nov.3.growth.rate.plot)

#######################################################################################
## Dose response data for clone 1 from each of populations 1-4 for B30, B30+A18, B30+A31.
## 12 clones total.

well.to.Tet.dose.2021.11.10 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the Tet dosage.
    rowletter <- substring(Well,1,1)
    if (rowletter == "A")
        Tet <- 0
    else if (rowletter == "B")
        Tet <- 2
    else if (rowletter == "C")
        Tet <- 4
    else if (rowletter == "D")
        Tet <- 8
    else if (rowletter == "E")
        Tet <- 10
    else if (rowletter == "F")
        Tet <- 20
    else if (rowletter == "G")
        Tet <- 30
    else if (rowletter == "H")
        Tet <- 40
    return(Tet)
}

well.to.treatment.2021.11.10 <- function(Well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the biological replicate.
    colnum <- as.numeric(substring(Well,2))
    plas.map.vec <- c(rep("no_plasmid",4), rep("pUC_plasmid",4), rep("p15A_plasmid",4))
    return(plas.map.vec[colnum])
}

well.to.pop.2021.11.10 <- function(Well) {
    ## This helper function maps the column of the well in
    ## the 96-well plate to the biological replicate.
    colnum <- as.numeric(substring(Well,2))
    pop <- colnum %% 4
    if (pop == 0) pop <- 4
    return(pop)
}


Nov.10.dose.response.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-09-LB-Tet-dose-response-Day9-some-B30-clones.csv", header=FALSE)

Nov.10.evolved.clone.tidy.data <- Nov.10.dose.response.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(clone = 1) %>% ## always measuring clone 1 out of the 5 clones for each pop.
    mutate(Tet = sapply(Well, well.to.Tet.dose.2021.11.10)) %>%
    mutate(pop = sapply(Well, well.to.pop.2021.11.10)) %>%
    ## assign the columns to the different treatments.
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.10)) %>%
    mutate(pop = as.factor(pop)) %>%
    mutate(treatment = as.factor(treatment))

Nov.10.dose.response.evolved.clone.plot1 <- ggplot(
    Nov.10.evolved.clone.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment,
        shape = pop)) +
    geom_point(size=0.5) +
    facet_wrap(.~Tet) + geom_vline(xintercept=24,linetype="dashed",color="red")

Nov.10.dose.response.evolved.clone.plot2 <- ggplot(
    filter(Nov.10.evolved.clone.tidy.data,
          treatment != "p15A_plasmid"),
    aes(x = hours,
        y = OD600,
        color = treatment,
        shape = pop)) +
    geom_point(size=0.5) +
    facet_wrap(.~Tet) + geom_vline(xintercept=24,linetype="dashed",color="red")

Nov.10.dose.response.evolved.clone.plot3 <- ggplot(
    filter(Nov.10.evolved.clone.tidy.data,
          treatment != "pUC_plasmid"),
    aes(x = hours,
        y = OD600,
        color = treatment,
        shape = pop)) +
    geom_point(size=0.5) +
    facet_wrap(.~Tet) + geom_vline(xintercept=24,linetype="dashed",color="red")

Nov.10.dose.response.evolved.clone.plot4 <- ggplot(
    filter(Nov.10.evolved.clone.tidy.data,
          treatment != "no_plasmid"),
    aes(x = hours,
        y = OD600,
        color = treatment,
        shape = pop)) +
    geom_point(size=0.5) +
    facet_wrap(.~Tet) + geom_vline(xintercept=24,linetype="dashed",color="red")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-10-dose-response-some-B30-clones-plot1.pdf",
       Nov.10.dose.response.evolved.clone.plot1)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-10-dose-response-some-B30-clones-plot2.pdf",
       Nov.10.dose.response.evolved.clone.plot2)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-10-dose-response-some-B30-clones-plot3.pdf",
       Nov.10.dose.response.evolved.clone.plot3)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-10-dose-response-some-B30-clones-plot4.pdf",
       Nov.10.dose.response.evolved.clone.plot4)



#### let's plot OD600 24 hours and growth rate between OD600 = 0.3-0.4.
Nov.10.evolved.clone.tidy.24h.data <- Nov.10.evolved.clone.tidy.data %>%
    filter(hours == 24)

Nov.10.evolved.clones.24h.OD600.plot <- ggplot(
    Nov.10.evolved.clone.tidy.24h.data,
    aes(x = treatment,
        y = OD600,
        color = pop)) + geom_point() + theme_classic() +
    ylab("OD600 at 24h timepoint") +
    xlab("plasmid copy number") + facet_wrap(.~Tet)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-10-evolved-clones-24h-OD600.pdf", Nov.10.evolved.clones.24h.OD600.plot)

## now, calculate the max growth rate (look at points between OD600 = 0.3-0.4).
Nov.10.evolved.clone.tidy.max.growth.rate.data <- Nov.10.evolved.clone.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, pop, Tet, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well,pop, Tet, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.10.evolved.clones.24h.max.growth.rate.plot <- ggplot(
    Nov.10.evolved.clone.tidy.max.growth.rate.data,
    aes(x = treatment,
        y = estimate,
        color = pop)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid") +
    facet_wrap(.~Tet)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-10-evolved-clones-growth-rate.pdf",
       Nov.10.evolved.clones.24h.max.growth.rate.plot)


## now, calculate the max growth rate (look at points between OD600 = 0.1-0.2,
## for comparison to the Tet 50 data).
Nov.10.evolved.clone.tidy.max.growth.rate.data2 <- Nov.10.evolved.clone.tidy.data %>%
    filter(OD600 >= 0.15) %>%
    filter(OD600 <= 0.25) %>%
    group_by(Well, pop, Tet, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well,pop, Tet, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.10.evolved.clones.24h.max.growth.rate.plot2 <- ggplot(
    Nov.10.evolved.clone.tidy.max.growth.rate.data2,
    aes(x = treatment,
        y = estimate,
        color = pop)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.1-0.2") +
    xlab("Plasmid copy number") +
    facet_wrap(.~Tet)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-10-evolved-clones-growth-rate2.pdf",
       Nov.10.evolved.clones.24h.max.growth.rate.plot2)


################################################################################

## analyze data from the first big scale growth curve data, started 11/17/2021.

well.to.treatment.2021.11.17 <- function(Well) {
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
            return("pUC_plasmid")
    }

    if (rowletter == "D" || rowletter == "E")
        return("pUC_plasmid")
    if (rowletter == "F" || rowletter == "G")
        return("p15A_plasmid")

    if (rowletter == "H") {
        if (colnum <= 6)
            return("p15A_plasmid")
        else
            return("Blank")
    }

    return(my.treatment)
}

Nov.17.LB.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-17-LB-B30-clones.csv", header=FALSE)

Nov.17.LB.tidy.data <- Nov.17.LB.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 0) %>%
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.17)) %>%
    mutate(treatment = as.factor(treatment)) %>%
    filter(!(treatment %in% c("Blank", "Empty")))

Nov.17.LB.plot <- ggplot(
    Nov.17.LB.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) +
    geom_point(size=0.5) +
    geom_vline(xintercept=24,linetype="dashed",color="red")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-17-LB-plot.pdf",
       Nov.17.LB.plot)


Nov.17.Tet20.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-17-LB-Tet20-B30-clones.csv", header=FALSE)

Nov.17.Tet20.tidy.data <- Nov.17.Tet20.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 0) %>%
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.17)) %>%
    mutate(treatment = as.factor(treatment)) %>%
    filter(!(treatment %in% c("Blank", "Empty")))

Nov.17.Tet20.plot <- ggplot(
    Nov.17.Tet20.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) +
    geom_point(size=0.2) +
    geom_vline(xintercept=24,linetype="dashed",color="red")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-17-Tet20-plot.pdf",
       Nov.17.Tet20.plot)

## These readings in Tecan Loaner are all messed up!!
Nov.17.Tet50.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-17-LB-Tet50-B30-clones.csv", header=FALSE)

Nov.17.Tet50.tidy.data <- Nov.17.Tet50.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 0) %>%
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.17)) %>%
    filter(!(treatment %in% c("Blank", "Empty")))
##    mutate(treatment = as.factor(treatment)) ##%>%

Nov.17.Tet50.plot1 <- ggplot(
    Nov.17.Tet50.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) +
    geom_point(size=0.5) +
    geom_vline(xintercept=24,linetype="dashed",color="red") +
    facet_wrap(.~treatment)

## let's plot the derivatives of log(OD600) over time.
Nov.17.LB.delta.log2.OD600.data <- Nov.17.LB.tidy.data %>%
    group_by(Well, treatment) %>%
    mutate(log2.OD600 = log2(OD600)) %>%
    mutate(delta.log2.OD600 = log2.OD600 - lag(log2.OD600))

Nov.17.LB.delta.log2.OD600.plot <- Nov.17.LB.delta.log2.OD600.data %>%
    group_by(Well) %>%
    filter(hours > 1) %>%
    ggplot(
    aes(x = hours,
        y = delta.log2.OD600,
        color = Well)) + geom_line() + theme_classic() +
    facet_grid(treatment~.) +
    ylab("Delta log2(OD600)") +
    xlab("Time (hours)") +
    guides(color=FALSE)

## let's plot the derivative of log(OD600) over time.
Nov.17.Tet20.delta.log2.OD600.data <- Nov.17.Tet20.tidy.data %>%
    group_by(Well, treatment) %>%
    mutate(log2.OD600 = log2(OD600)) %>%
    mutate(delta.log2.OD600 = log2.OD600 - lag(log2.OD600))

ggsave("../results/draft-manuscript-1A/growth-results/Nov-17-LB-delta-growth.pdf",
       Nov.17.LB.delta.log2.OD600.plot)

Nov.17.Tet20.delta.log2.OD600.plot <- Nov.17.Tet20.delta.log2.OD600.data %>%
    group_by(Well) %>%
    filter(hours > 1) %>%
    ggplot(
    aes(x = hours,
        y = delta.log2.OD600,
        color = Well)) + geom_line() + theme_classic() +
    facet_grid(treatment~.) +
    ylab("Delta log2(OD600)") +
    xlab("Time (hours)") +
    guides(color=FALSE)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-17-Tet20-delta-growth.pdf",
       Nov.17.Tet20.delta.log2.OD600.plot)


## plot growth rate and OD600 at 24 hour distributions.

## for LB + Tet0.
Nov.17.LB.tidy.max.growth.rate.data <- Nov.17.LB.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.17.LB.max.growth.rate.plot <- ggplot(
    Nov.17.LB.tidy.max.growth.rate.data,
    aes(x = treatment,
        y = estimate,
        color = treatment)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid")

Nov.17.LB.max.growth.rate.histogram <- ggplot(
        Nov.17.LB.tidy.max.growth.rate.data,
    aes(fill = treatment,
        x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    xlab("growth rate between OD600 0.3-0.4") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-17-LB-growth-rate-histogram.pdf",
       Nov.17.LB.max.growth.rate.histogram)


Nov.17.LB.24h.OD600.histogram <- Nov.17.LB.tidy.data %>%
    filter(hours == 24) %>%
    ggplot(
        aes(fill = treatment,
        x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    ylab("OD600 at 24h timepoint") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-17-LB-24h-OD600-histogram.pdf",
       Nov.17.LB.24h.OD600.histogram)


#### Now repeat for Tet20.

Nov.17.Tet20.tidy.max.growth.rate.data <- Nov.17.Tet20.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.17.Tet20.max.growth.rate.plot <- ggplot(
    Nov.17.Tet20.tidy.max.growth.rate.data,
    aes(x = treatment,
        y = estimate,
        color = treatment)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid")

Nov.17.Tet20.max.growth.rate.histogram <- ggplot(
        Nov.17.Tet20.tidy.max.growth.rate.data,
    aes(fill = treatment,
        x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    xlab("growth rate between OD600 0.3-0.4") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-17-Tet20-growth-rate-histogram.pdf", Nov.17.Tet20.max.growth.rate.histogram)

Nov.17.Tet20.24h.OD600.histogram <- Nov.17.Tet20.tidy.data %>%
    filter(hours == 24) %>%
    ggplot(
        aes(fill = treatment,
        x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    ylab("OD600 at 24h timepoint") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-17-Tet20-24h-OD600-histogram.pdf", Nov.17.Tet20.24h.OD600.histogram)


## single readings of LB+Tet50 plate in Tecan J, since the readings in Tecan Loaner
## were all messed up.

Nov.19.Tet50.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-19-Tet50-single-reading-B30-clones.csv", header=FALSE)

Nov.19.Tet50.tidy.data <- Nov.19.Tet50.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 0) %>%
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.17)) %>%
    filter(!(treatment %in% c("Blank", "Empty")))

Nov.19.Tet50.plot <- ggplot(
    Nov.19.Tet50.tidy.data,
    aes(x = treatment,
        y = OD600,
        color = treatment)) +
    geom_point(size=0.5) +
    theme_classic() +
    geom_vline(xintercept=24,linetype="dashed",color="red")

Nov.19.Tet50.histogram <- ggplot(
    Nov.19.Tet50.tidy.data,
    aes(fill = treatment,
        x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-19-Tet50-36h-OD600-histogram.pdf", Nov.19.Tet50.histogram)


#############################################
## results from Replicate 2: started on 11/21/2021.
## TODO: refactor all this code to get rid of the copy-pasting.
## NOTE: same format as Replicate 1 dated 11/17/2021.

Nov.21.LB.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-21-LB-B30-clones.csv", header=FALSE)

Nov.21.LB.tidy.data <- Nov.21.LB.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 0) %>%
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.17)) %>%
    mutate(treatment = as.factor(treatment)) %>%
    filter(!(treatment %in% c("Blank", "Empty")))

Nov.21.LB.plot <- ggplot(
    Nov.21.LB.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) +
    geom_point(size=0.5) +
    geom_vline(xintercept=24,linetype="dashed",color="red")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-LB-plot.pdf",
       Nov.21.LB.plot)


Nov.21.Tet20.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-21-LB-Tet20-B30-clones.csv", header=FALSE)

Nov.21.Tet20.tidy.data <- Nov.21.Tet20.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 0) %>%
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.17)) %>%
    mutate(treatment = as.factor(treatment)) %>%
    filter(!(treatment %in% c("Blank", "Empty")))

Nov.21.Tet20.plot <- ggplot(
    Nov.21.Tet20.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) +
    geom_point(size=0.2) +
    geom_vline(xintercept=24,linetype="dashed",color="red")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-Tet20-plot.pdf",
       Nov.21.Tet20.plot)

Nov.21.Tet50.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-21-LB-Tet50-B30-clones.csv", header=FALSE)

Nov.21.Tet50.tidy.data <- Nov.21.Tet50.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 0) %>%
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.17)) %>%
    filter(!(treatment %in% c("Blank", "Empty")))


Nov.21.Tet50.plot <- ggplot(
    Nov.21.Tet50.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) +
    geom_point(size=0.2) +
    geom_vline(xintercept=24,linetype="dashed",color="red")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-Tet50-plot.pdf",
       Nov.21.Tet50.plot)


## let's plot the derivative of log(OD600) over time.

Nov.21.LB.delta.log2.OD600.data <- Nov.21.LB.tidy.data %>%
    group_by(Well, treatment) %>%
    mutate(log2.OD600 = log2(OD600)) %>%
    mutate(delta.log2.OD600 = log2.OD600 - lag(log2.OD600))

Nov.21.LB.delta.log2.OD600.plot <- Nov.21.LB.delta.log2.OD600.data %>%
    group_by(Well) %>%
    filter(hours > 1) %>%
    ggplot(
    aes(x = hours,
        y = delta.log2.OD600,
        color = Well)) + geom_line() + theme_classic() +
    facet_grid(treatment~.) +
    ylab("Delta log2(OD600)") +
    xlab("Time (hours)") +
    guides(color=FALSE)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-LB-delta-growth.pdf",
       Nov.21.LB.delta.log2.OD600.plot)


Nov.21.Tet20.delta.log2.OD600.data <- Nov.21.Tet20.tidy.data %>%
    group_by(Well, treatment) %>%
    mutate(log2.OD600 = log2(OD600)) %>%
    mutate(delta.log2.OD600 = log2.OD600 - lag(log2.OD600))

Nov.21.Tet20.delta.log2.OD600.plot <- Nov.21.Tet20.delta.log2.OD600.data %>%
    group_by(Well) %>%
    filter(hours > 1) %>%
    ggplot(
    aes(x = hours,
        y = delta.log2.OD600,
        color = Well)) + geom_line() + theme_classic() +
    facet_grid(treatment~.) +
    ylab("Delta log2(OD600)") +
    xlab("Time (hours)") +
    guides(color=FALSE)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-Tet20-delta-growth.pdf",
       Nov.21.Tet20.delta.log2.OD600.plot)


Nov.21.Tet50.delta.log2.OD600.data <- Nov.21.Tet50.tidy.data %>%
    group_by(Well, treatment) %>%
    mutate(log2.OD600 = log2(OD600)) %>%
    mutate(delta.log2.OD600 = log2.OD600 - lag(log2.OD600))

Nov.21.Tet50.delta.log2.OD600.plot <- Nov.21.Tet50.delta.log2.OD600.data %>%
    group_by(Well) %>%
    filter(hours > 1) %>%
    ggplot(
    aes(x = hours,
        y = delta.log2.OD600,
        color = Well)) + geom_line() + theme_classic() +
    facet_grid(treatment~.) +
    ylab("Delta log2(OD600)") +
    xlab("Time (hours)") +
    guides(color=FALSE)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-Tet50-delta-growth.pdf",
       Nov.21.Tet50.delta.log2.OD600.plot)


## plot growth rate and OD600 at 24 hour distributions.

## for LB + Tet0.
Nov.21.LB.tidy.max.growth.rate.data <- Nov.21.LB.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.21.LB.max.growth.rate.plot <- ggplot(
    Nov.21.LB.tidy.max.growth.rate.data,
    aes(x = treatment,
        y = estimate,
        color = treatment)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid")

Nov.21.LB.max.growth.rate.histogram <- ggplot(
    Nov.21.LB.tidy.max.growth.rate.data,
    aes(fill = treatment,
        x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    xlab("growth rate between OD600 0.3-0.4") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-LB-growth-rate-histogram.pdf",
       Nov.21.LB.max.growth.rate.histogram)


Nov.21.LB.24h.OD600.histogram <- Nov.21.LB.tidy.data %>%
    filter(hours == 24) %>%
    ggplot(
        aes(fill = treatment,
            x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    ylab("OD600 at 24h timepoint") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-LB-24h-OD600-histogram.pdf",
       Nov.21.LB.24h.OD600.histogram)


#### Now repeat for Tet20.

Nov.21.Tet20.tidy.max.growth.rate.data <- Nov.21.Tet20.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.21.Tet20.max.growth.rate.plot <- ggplot(
    Nov.21.Tet20.tidy.max.growth.rate.data,
    aes(x = treatment,
        y = estimate,
        color = treatment)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid")

Nov.21.Tet20.max.growth.rate.histogram <- ggplot(
        Nov.21.Tet20.tidy.max.growth.rate.data,
    aes(fill = treatment,
        x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    xlab("growth rate between OD600 0.3-0.4") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-Tet20-growth-rate-histogram.pdf", Nov.21.Tet20.max.growth.rate.histogram)

Nov.21.Tet20.24h.OD600.histogram <- Nov.21.Tet20.tidy.data %>%
    filter(hours == 24) %>%
    ggplot(
        aes(fill = treatment,
            x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    ylab("OD600 at 24h timepoint") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-Tet20-24h-OD600-histogram.pdf", Nov.21.Tet20.24h.OD600.histogram)

#### Now repeat for Tet50.

Nov.21.Tet50.tidy.max.growth.rate.data <- Nov.21.Tet50.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.21.Tet50.max.growth.rate.plot <- ggplot(
    Nov.21.Tet50.tidy.max.growth.rate.data,
    aes(x = treatment,
        y = estimate,
        color = treatment)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid")

Nov.21.Tet50.max.growth.rate.histogram <- ggplot(
        Nov.21.Tet50.tidy.max.growth.rate.data,
    aes(fill = treatment,
        x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    xlab("growth rate between OD600 0.3-0.4") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-Tet50-growth-rate-histogram.pdf", Nov.21.Tet50.max.growth.rate.histogram)

Nov.21.Tet50.24h.OD600.histogram <- Nov.21.Tet50.tidy.data %>%
    filter(hours == 24) %>%
    ggplot(
        aes(fill = treatment,
            x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    ylab("OD600 at 24h timepoint") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-21-Tet50-24h-OD600-histogram.pdf", Nov.21.Tet50.24h.OD600.histogram)

#############################################
## results from Replicate 3: started on 11/23/2021.
## TODO: refactor all this code to get rid of the copy-pasting.
## NOTE: same format as Replicate 1 dated 11/17/2021.

Nov.23.LB.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-23-LB-B30-clones.csv", header=FALSE)

Nov.23.LB.tidy.data <- Nov.23.LB.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 0) %>%
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.17)) %>%
    mutate(treatment = as.factor(treatment)) %>%
    filter(!(treatment %in% c("Blank", "Empty")))

Nov.23.LB.plot <- ggplot(
    Nov.23.LB.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) +
    geom_point(size=0.5) +
    geom_vline(xintercept=24,linetype="dashed",color="red")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-LB-plot.pdf",
       Nov.23.LB.plot)


Nov.23.Tet20.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-23-LB-Tet20-B30-clones.csv", header=FALSE)

Nov.23.Tet20.tidy.data <- Nov.23.Tet20.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 0) %>%
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.17)) %>%
    mutate(treatment = as.factor(treatment)) %>%
    filter(!(treatment %in% c("Blank", "Empty")))

Nov.23.Tet20.plot <- ggplot(
    Nov.23.Tet20.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) +
    geom_point(size=0.2) +
    geom_vline(xintercept=24,linetype="dashed",color="red")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-Tet20-plot.pdf",
       Nov.23.Tet20.plot)

Nov.23.Tet50.evolved.clone.long.format.data <- read.csv("../data/draft-manuscript-1A/OD600-data/2021-11-23-LB-Tet50-B30-clones.csv", header=FALSE)

Nov.23.Tet50.tidy.data <- Nov.23.Tet50.evolved.clone.long.format.data %>%
    pivot_longer(cols=V1, values_to="Well") %>%
    select(-name) %>% ## drop this useless column
    pivot_longer(!Well,
                 names_to = "Time",
                 values_to = "OD600",
                 names_pattern = "V(.+)") %>%
    mutate(Time = as.numeric(Time) - 1) %>%
    mutate(minutes = 10*Time) %>%
    mutate(hours = minutes/60) %>%
    mutate(Tet = 0) %>%
    mutate(treatment = sapply(Well, well.to.treatment.2021.11.17)) %>%
    filter(!(treatment %in% c("Blank", "Empty")))


Nov.23.Tet50.plot <- ggplot(
    Nov.23.Tet50.tidy.data,
    aes(x = hours,
        y = OD600,
        color = treatment)) +
    geom_point(size=0.2) +
    geom_vline(xintercept=24,linetype="dashed",color="red")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-Tet50-plot.pdf",
       Nov.23.Tet50.plot)


Nov.23.LB.delta.log2.OD600.data <- Nov.23.LB.tidy.data %>%
    group_by(Well, treatment) %>%
    mutate(log2.OD600 = log2(OD600)) %>%
    mutate(delta.log2.OD600 = log2.OD600 - lag(log2.OD600))

Nov.23.LB.delta.log2.OD600.plot <- Nov.23.LB.delta.log2.OD600.data %>%
    group_by(Well) %>%
    filter(hours > 1) %>%
    ggplot(
    aes(x = hours,
        y = delta.log2.OD600,
        color = Well)) + geom_line() + theme_classic() +
    facet_grid(treatment~.) +
    ylab("Delta log2(OD600)") +
    xlab("Time (hours)") +
    guides(color=FALSE)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-LB-delta-growth.pdf",
       Nov.23.LB.delta.log2.OD600.plot)


Nov.23.Tet20.delta.log2.OD600.data <- Nov.23.Tet20.tidy.data %>%
    group_by(Well, treatment) %>%
    mutate(log2.OD600 = log2(OD600)) %>%
    mutate(delta.log2.OD600 = log2.OD600 - lag(log2.OD600))

Nov.23.Tet20.delta.log2.OD600.plot <- Nov.23.Tet20.delta.log2.OD600.data %>%
    group_by(Well) %>%
    filter(hours > 1) %>%
    ggplot(
    aes(x = hours,
        y = delta.log2.OD600,
        color = Well)) + geom_line() + theme_classic() +
    facet_grid(treatment~.) +
    ylab("Delta log2(OD600)") +
    xlab("Time (hours)") +
    guides(color=FALSE)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-Tet20-delta-growth.pdf",
       Nov.23.Tet20.delta.log2.OD600.plot)


Nov.23.Tet50.delta.log2.OD600.data <- Nov.23.Tet50.tidy.data %>%
    group_by(Well, treatment) %>%
    mutate(log2.OD600 = log2(OD600)) %>%
    mutate(delta.log2.OD600 = log2.OD600 - lag(log2.OD600))

Nov.23.Tet50.delta.log2.OD600.plot <- Nov.23.Tet50.delta.log2.OD600.data %>%
    group_by(Well) %>%
    filter(hours > 1) %>%
    ggplot(
    aes(x = hours,
        y = delta.log2.OD600,
        color = Well)) + geom_line() + theme_classic() +
    facet_grid(treatment~.) +
    ylab("Delta log2(OD600)") +
    xlab("Time (hours)") +
    guides(color=FALSE)

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-Tet50-delta-growth.pdf",
       Nov.23.Tet50.delta.log2.OD600.plot)

## plot growth rate and OD600 at 24 hour distributions.

## for LB + Tet0.
Nov.23.LB.tidy.max.growth.rate.data <- Nov.23.LB.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.23.LB.max.growth.rate.plot <- ggplot(
    Nov.23.LB.tidy.max.growth.rate.data,
    aes(x = treatment,
        y = estimate,
        color = treatment)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid")

Nov.23.LB.max.growth.rate.histogram <- ggplot(
    Nov.23.LB.tidy.max.growth.rate.data,
    aes(fill = treatment,
        x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    xlab("growth rate between OD600 0.3-0.4") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-LB-growth-rate-histogram.pdf",
       Nov.23.LB.max.growth.rate.histogram)


Nov.23.LB.24h.OD600.histogram <- Nov.23.LB.tidy.data %>%
    filter(hours == 24) %>%
    ggplot(
        aes(fill = treatment,
            x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    ylab("OD600 at 24h timepoint") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-LB-24h-OD600-histogram.pdf",
       Nov.23.LB.24h.OD600.histogram)


#### Now repeat for Tet20.

Nov.23.Tet20.tidy.max.growth.rate.data <- Nov.23.Tet20.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.23.Tet20.max.growth.rate.plot <- ggplot(
    Nov.23.Tet20.tidy.max.growth.rate.data,
    aes(x = treatment,
        y = estimate,
        color = treatment)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid")

Nov.23.Tet20.max.growth.rate.histogram <- ggplot(
        Nov.23.Tet20.tidy.max.growth.rate.data,
    aes(fill = treatment,
        x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    xlab("growth rate between OD600 0.3-0.4") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-Tet20-growth-rate-histogram.pdf", Nov.23.Tet20.max.growth.rate.histogram)

Nov.23.Tet20.24h.OD600.histogram <- Nov.23.Tet20.tidy.data %>%
    filter(hours == 24) %>%
    ggplot(
        aes(fill = treatment,
            x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    ylab("OD600 at 24h timepoint") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-Tet20-24h-OD600-histogram.pdf", Nov.23.Tet20.24h.OD600.histogram)

#### Now repeat for Tet50.

Nov.23.Tet50.tidy.max.growth.rate.data <- Nov.23.Tet50.tidy.data %>%
    filter(OD600 >= 0.3) %>%
    filter(OD600 <= 0.4) %>%
    group_by(Well, treatment) %>%
    nest() %>% ## nest the data by well, then calculate the log2(OD600) slope using lm().
    mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
    ## get the fit parameters.
    ## see tutorial at:
    ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
    mutate(tidied = map(fit, tidy)) %>%
    select(Well, treatment, tidied) %>%
    unnest(tidied) %>%
    filter(term == 'hours') ## we only care about the slope parameter.

Nov.23.Tet50.max.growth.rate.plot <- ggplot(
    Nov.23.Tet50.tidy.max.growth.rate.data,
    aes(x = treatment,
        y = estimate,
        color = treatment)) + geom_point() + theme_classic() +
    ylab("growth rate between OD600 0.3-0.4") +
    xlab("Presence/absence of high-copy pUC target plasmid")

Nov.23.Tet50.max.growth.rate.histogram <- ggplot(
        Nov.23.Tet50.tidy.max.growth.rate.data,
    aes(fill = treatment,
        x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    xlab("growth rate between OD600 0.3-0.4") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-Tet50-growth-rate-histogram.pdf", Nov.23.Tet50.max.growth.rate.histogram)

Nov.23.Tet50.24h.OD600.histogram <- Nov.23.Tet50.tidy.data %>%
    filter(hours == 24) %>%
    ggplot(
        aes(fill = treatment,
            x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
    facet_grid(treatment~.) +
    ylab("OD600 at 24h timepoint") +
    labs(fill = "Presence/absence of plasmid") +
    theme(legend.position="top")

ggsave("../results/draft-manuscript-1A/growth-results/Nov-23-Tet50-24h-OD600-histogram.pdf", Nov.23.Tet50.24h.OD600.histogram)
