 ## qPCR-analysis.R by Rohan Maddamsetti.
## This script analyzes qPCR data using Yi's transposon
## and antibiotic resistance marker system.

library(tidyverse)
library(cowplot)

calc.probe.fold.differences <- function(well.df) {
    ## this is a helper function for calculating probe fold differences
    ## per well.

    ## data analysis using constants calculated from Yi's standard curve calibration.
    T.per.C.constant <- 0.39071847356712
    K.per.C.constant <- 0.58657387456313
    T.per.K.constant <- 0.666102754504912

    C <- filter(well.df, probe == 'C')$cycle_at_threshold
    T <- filter(well.df, probe == 'T')$cycle_at_threshold
    K <- filter(well.df, probe == 'K')$cycle_at_threshold

    T.per.C <- 2^(C - T)/T.per.C.constant
    K.per.C <- 2^(C - K)/K.per.C.constant
    T.per.K <- 2^(K - T)/T.per.K.constant

    ## This is what Yi does on his spreadsheet.
    ## He subtracts 1 in the denominator to account
    ## for the copy of the transposon on the chromosome.
    Yi.transposon.on.plasmid.fraction.calc <- 1 - (K.per.C - T.per.C)/(K.per.C - 1)

    return.df <- data.frame(Well = unique(well.df$Well),
                            Transposon = unique(well.df$Transposon),
                            Treatment = unique(well.df$Treatment),
                            Replicate = unique(well.df$Replicate),
                            transposons.per.chromosome = T.per.C,
                            plasmids.per.chromosome = K.per.C,
                            transposons.per.plasmid = T.per.K,
                            Yi.transposon.frac = Yi.transposon.on.plasmid.fraction.calc
                            )
    
    return(return.df)
}

calc.probe.fold.differences2 <- function(well.df) {
    ## this is a helper function for calculating probe fold differences
    ## per well.

    ## data analysis using constants calculated from Yi's standard curve calibration.
    T.per.C.constant <- 0.39071847356712
    C <- filter(well.df, probe == 'C')$cycle_at_threshold
    T <- filter(well.df, probe == 'T')$cycle_at_threshold
    T.per.C <- 2^(C - T)/T.per.C.constant

    return.df <- data.frame(Well = unique(well.df$Well),
                            Transposon = unique(well.df$Transposon),
                            Treatment = unique(well.df$Treatment),
                            Replicate = unique(well.df$Replicate),
                            transposons.per.chromosome = T.per.C)
    
    return(return.df)
}

calc.tet.cm.CNV.with.B30.control <- function(df.with.B30.control) {
    ## IMPORTANT: this function assumes perfect amplification efficiency
    ## (m = 1, or 2x amplification each cycle.)
    ## I can use a standard curve to determine m more precisely in the future.

    m.cm <- 1 ## TODO: calibrate with a standard curve.
    m.tet <- 1 ## TODO: calibrate with a standard curve.

    control.data <- filter(df.with.B30.control, Treatment == "B30")
    tet.control.data <- filter(control.data, probe == "T")
    cm.control.data <- filter(control.data, probe == "C")

    ## question: maybe I should be using geometric mean here instead?
    control.tet.Cq <- mean(tet.control.data$cycle_at_threshold)
    control.cm.Cq <- mean(cm.control.data$cycle_at_threshold)
    
    beta <- 2^(m.tet * control.tet.Cq - m.cm * control.cm.Cq)

    helper.func <- function(well.df) {
        ## note: this helper uses the beta variable defined above.
        C <- filter(well.df, probe == 'C')$cycle_at_threshold
        T <- filter(well.df, probe == 'T')$cycle_at_threshold
        T.per.C <- beta * 2^(m.cm * C - m.tet * T)

        return.df <- data.frame(Well = unique(well.df$Well),
                                Transposon = unique(well.df$Transposon),
                                Treatment = unique(well.df$Treatment),
                                Sample = unique(well.df$Sample),
                                transposons.per.chromosome = T.per.C)
        return(return.df)
    }

        
    results <- df.with.B30.control %>%
        split(.$Well) %>%
        map_dfr(helper.func)
    return(results)
}

calc.tet.cm.CNV.with.INIT.control <- function(df.with.INIT.control) {
    ## IMPORTANT: this function assumes perfect amplification efficiency
    ## (m = 1, or 2x amplification each cycle.)
    ## I can use a standard curve to determine m more precisely in the future.

    m.cm <- 2^1.055 ## calibrated with May 6 2021 standard curve.
    m.tet <- 2^0.904 ## calibrated with May 6 2021 standard curve.

    control.data <- filter(df.with.INIT.control, Treatment == "INIT")
    tet.control.data <- filter(control.data, probe == "T")
    cm.control.data <- filter(control.data, probe == "C")

    ## question: maybe I should be using geometric mean here instead?
    mean.control.tet.Cq <- mean(tet.control.data$cycle_at_threshold)
    mean.control.cm.Cq <- mean(cm.control.data$cycle_at_threshold)
    
    beta <- 2^(m.tet * mean.control.tet.Cq - m.cm * mean.control.cm.Cq)

    helper.func <- function(well.df) {
        ## note: this helper uses the beta variable defined above.
        C <- filter(well.df, probe == 'C')$cycle_at_threshold
        T <- filter(well.df, probe == 'T')$cycle_at_threshold
        T.per.C <- beta * 2^(m.cm * C - m.tet * T)

        return.df <- data.frame(Well = unique(well.df$Well),
                                Transposon = unique(well.df$Transposon),
                                Treatment = unique(well.df$Treatment),
                                Replicate = unique(well.df$Replicate),
                                transposons.per.chromosome = T.per.C)
        return(return.df)
    }

    
    results <- df.with.INIT.control %>%
        split(.$Well) %>%
        map_dfr(helper.func)
    return(results)
}

######################################################################

## standard curve analysis. Data collected May 6 2021.

calibration.data.may.06 <- read.csv(
    "../data/draft-manuscript-1A/qPCR-data/2021-05-06-Rohan-B30-calibration.csv") %>%
    ## data from LB culture dilutions craps out at the low dilution end.
    ## the DNA data looks quite good across the whole range.
    ## filter out the low end.
    filter(Log2DilutionFactor > -3)
    

## plot a linear regression per Treatment.
## INIT := 1:100 dilution of LB culture of DH5a + B30 cells in PCR H20.
## B30 := dilution of purified B30 miniTn5 plasmid DNA in PCR H20.

calibration.plot <- ggplot(calibration.data.may.06,
                           aes(x = Log2DilutionFactor,
                               y = cycle_at_threshold,
                               color = probe)) +
    geom_point() + geom_smooth(method="lm") +
    facet_wrap(.~Treatment)

ggsave("../results/draft-manuscript-1A/qPCR-results/2021-05-06-B30-calibration-std-curve.pdf",calibration.plot)

## get the estimates of the slopes of the linear regression.
## the slope m is a measure of the probe amplification efficiency.
## slope = 1 with perfect efficiency.

## the purified plasmid sample is a control.
## in terms of parameter estimates, we only care about the culture dilution.

plasmid.calibration.data <- calibration.data.may.06 %>%
    filter(Treatment == "B30")

cell.calibration.data <- calibration.data.may.06 %>%
    filter(Treatment == "INIT")

T.plasmid.calibration.data <- plasmid.calibration.data %>%
    filter(probe == "T")

C.plasmid.calibration.data <- plasmid.calibration.data %>%
    filter(probe == "C")

T.cell.calibration.data <- cell.calibration.data %>%
    filter(probe == "T")

C.cell.calibration.data <- cell.calibration.data %>%
    filter(probe == "C")

T.plasmid.regression <- lm(cycle_at_threshold~Log2DilutionFactor,
                           data=T.plasmid.calibration.data)

summary(T.plasmid.regression)

C.plasmid.regression <- lm(cycle_at_threshold~Log2DilutionFactor,
                           data=C.plasmid.calibration.data)

summary(C.plasmid.regression)

T.cell.regression <- lm(cycle_at_threshold~Log2DilutionFactor,
                           data=T.cell.calibration.data)

summary(T.cell.regression)

C.cell.regression <- lm(cycle_at_threshold~Log2DilutionFactor,
                           data=C.cell.calibration.data)

summary(C.cell.regression)


######################################################################

##### data from replicating Yi's experiment with target plasmid.

data.march.9 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-03-09-reformatted.csv")

figure.march.9 <- ggplot(data.march.9, aes(x = Treatment,
                                   y = cycle_at_threshold,
                                   color = probe)) +
    geom_point() +
    theme_classic()

ggsave("../results/draft-manuscript-1A/qPCR-results/raw-qPCR-2021-03-09.pdf",
       figure.march.9, height = 3, width = 3)


results.march.9 <- data.march.9 %>%
    mutate(Transposon = "B30") %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    filter(!(Treatment == -1))

Fig2A <- ggplot(results.march.9, aes(x = Treatment,
                                   y = transposons.per.chromosome,
                                   color = Replicate)) +
    geom_point() +
    ##geom_line() +
    theme_classic() +
    guides(color = FALSE) +
    ggtitle("Transposons per chromosome")

Fig2B <- ggplot(results.march.9, aes(x = Treatment,
                                   y = plasmids.per.chromosome,
                                   color = Replicate)) +
    geom_point() +
    ##geom_line() +
    theme_classic() +
    guides(color = FALSE) +
    ggtitle("Plasmids per chromosome")

Fig2C <- ggplot(results.march.9, aes(x = Treatment,
                                   y = transposons.per.plasmid,
                                   color = Replicate)) +
    geom_point() +
    ##geom_line() +
    theme_classic() +
    guides(color = FALSE) +
    ggtitle("Transposons per plasmid")

Fig2 <- plot_grid(Fig2A, Fig2B, Fig2C, labels=c('A','B','C'),nrow=1)

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-03-09.pdf", Fig2, height = 3, width = 9)



##### data from initial no target plasmid experiments.

data.march.24 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-03-24-reformatted.csv")

figure.march.24 <- ggplot(data.march.24, aes(x = Treatment,
                                   y = cycle_at_threshold,
                                   color = probe)) +
    geom_point() +
    theme_classic()


results.march.24 <- data.march.24 %>%
    mutate(Transposon = "B30") %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    filter(!(Treatment == -1)) %>%
    mutate(Day = "Day 1")

Fig.march.24 <- ggplot(results.march.24, aes(x = Treatment,
                                   y = transposons.per.chromosome,
                                   shape = Replicate)) +
    geom_point() +
    theme_classic() +
    ggtitle("Transposons per chromosome: Day 1")
ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-03-24.pdf", Fig.march.24, height = 5, width = 5)

##### data from week-long no target plasmid experiment.

data.march.29 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-03-29-reformatted.csv")

figure.march.29 <- ggplot(data.march.29, aes(x = Treatment,
                                   y = cycle_at_threshold,
                                   color = probe)) +
    geom_point() +
    theme_classic()

results.march.29 <- data.march.29 %>%
    mutate(Transposon = "B30") %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    filter(!(Treatment == -1)) %>%
    mutate(Day = "Day 6")


initial.results <- full_join(results.march.24, results.march.29)

initial.fig <- ggplot(initial.results,
                      aes(x = Treatment,
                          y = transposons.per.chromosome,
                          shape = Replicate,
                          color = Day)) +
    geom_point() +
    theme_classic() +
    xlab("Tet concentration (ug/mL)") +
    ylab("Transposons per chromosome") +
    ggtitle("Selection for increased tet resistance causes Tn5 duplications")
ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-03-29-initial-result.pdf", initial.fig, width= 6.5, height=3)

##### data after one week-- analysis of B30 miniTn5 (after 1 week of daily transfers) and
## B107 IS1A (after 3 days without transfer).

data.march.31 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-03-31-reformatted.csv")

results.march.31 <- data.march.31 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    filter(!(Treatment == -1)) %>%
    mutate(Day = ifelse(Transposon == "B30","Day 8","96h"))


big.results <- results.march.24 %>%
    full_join(results.march.29) %>%
    full_join(results.march.31)

big.fig <- ggplot(big.results,
                      aes(x = Treatment,
                          y = transposons.per.chromosome,
                          shape = Replicate,
                          color = Day)) +
    facet_wrap(.~Transposon) +
    geom_point() +
    theme_classic() +
    xlab("Tet concentration (ug/mL)") +
    ylab("Transposons per chromosome") +
    ggtitle("Selection for tet resistance causes transposition-mediated duplications")
ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-03-31.pdf", big.fig, width= 6.5, height=3)


###########################
## Day 5 result of B30 experiment with 5x biological replicates (4/11/21).
## This qPCR experiment did not have technical replicates, so I will conduct
## a new qPCR experiment with the same cultures/cells today (4/12/21).


data.april.11 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-04-11-reformatted.csv")

results.april.11 <- data.april.11 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    filter(!(Treatment == -1)) %>%
    mutate(Day = ifelse(Transposon == "B30","Day 5","NA"))

april.11.fig <- ggplot(results.april.11,
                       aes(x = Treatment,
                           y = transposons.per.chromosome,
                           shape = Replicate,
                           color = Day)) +
    geom_point() +
    theme_classic() +
    xlab("Tet concentration (ug/mL)") +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tet resistance causes transposition-mediated duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-04-11.pdf", april.11.fig, width= 6.5, height=3)

###########################
## Day 5 result of B30 experiment with 5x biological replicates,
## and 3x technical replicates (4/12/21).

data.april.12 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-04-12-reformatted.csv")

results.april.12 <- data.april.12 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    filter(!(Treatment == -1)) %>%
    mutate(Day = ifelse(Transposon == "B30","Day 5","NA"))

april.12.fig <- ggplot(results.april.12,
                       aes(x = Treatment,
                           y = transposons.per.chromosome,
                           shape = Replicate,
                           color = Day)) +
    geom_point() +
    theme_classic() +
    xlab("Tet concentration (ug/mL)") +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tet resistance causes transposition-mediated duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-04-12.pdf", april.12.fig, width= 6.5, height=3)

###########################
## Day 9 result of B30 experiment with 5x biological replicates,
## and 3x technical replicates (4/14/21).

data.april.14 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-04-14-reformatted.csv")

results.april.14 <- data.april.14 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    filter(!(Treatment == -1)) %>%
    mutate(Day = ifelse(Transposon == "B30","Day 9","NA"))

april.14.fig <- ggplot(results.april.14,
                       aes(x = Treatment,
                           y = transposons.per.chromosome,
                           shape = Replicate,
                           color = Day)) +
    geom_point() +
    theme_classic() +
    xlab("Tet concentration (ug/mL)") +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tet resistance causes transposition-mediated duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-04-14.pdf", april.14.fig, width= 6.5, height=3)

###########################
## Day 9 result of B30 experiment with 5x biological replicates,
## and 3x technical replicates (4/14/21).
## RE-RUN the analysis using different cycle thresholds:
## Cm probe = 0.2,
## Tet probe = 0.04.
## These are based on what Yi did for his standard curve.

data.april.14.thres2 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-04-14-threshold2-reformatted.csv")

results.april.14.thres2 <- data.april.14.thres2 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences2) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    filter(!(Treatment == -1)) %>%
    mutate(Day = ifelse(Transposon == "B30","Day 9","NA"))

april.14.thres2.fig <- ggplot(results.april.14.thres2,
                       aes(x = Treatment,
                           y = transposons.per.chromosome,
                           shape = Replicate)) +
    geom_point() +
    theme_classic() +
    xlab("Tet concentration (ug/mL)") +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tet resistance causes transposition-mediated duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-04-14-thres2.pdf", april.14.thres2.fig, width= 6.5, height=3)


###########################################################################################

## Day 9 result of B30 experiment with 5x biological replicates,
## and 3x technical replicates (4/14/21), grown on 4/20/21.
## experiment done on 4/21/21.
## I used the automatic machine Ct thresholding this time.
## IMPORTANT: wells E4, E5, E6 were COMPLETELY EMPTY,
## due to pipette error of my 37x master mix, using my newly "calibrated" P20.

data.april.21 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-04-21-reformatted.csv")

results.april.21 <- data.april.21 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences2) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    filter(!(Treatment == -1)) %>%
    mutate(Day = ifelse(Transposon == "B30", 9,"NA")) %>%
    ## remove the empty wells.
    filter(!(Well %in% c('E4','E5','E6'))) %>%
    ## set the internal controls to Day 9.
    mutate(Day = ifelse(Treatment %in% c('INIT', 'B30'), 0, 9))

april.21.fig <- ggplot(results.april.21,
                       aes(x = Day,
                           y = transposons.per.chromosome,
                           color = Treatment,
                           shape = Replicate)) +
    geom_point() +
    theme_classic() +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tet resistance causes transposition-mediated duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-04-21.pdf", april.21.fig, width= 6.5, height=3.5)

##########################

## let's recalculate copy number variation, assuming perfect amplification efficiency,
## and using the B30 plasmid DNA internal control.

results2.april.21 <- data.april.21 %>%
    calc.tet.cm.CNV.with.B30.control() %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(Day = ifelse(Transposon == "B30", 9,"NA")) %>%
    ## remove the empty wells.
    filter(!(Well %in% c('E4','E5','E6'))) %>%
    ## set the internal controls to Day 0.
    mutate(Day = ifelse(Treatment %in% c('INIT', 'B30'), 0, 9))

results3.april.21 <- data.april.21 %>%
    calc.tet.cm.CNV.with.INIT.control() %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(Day = ifelse(Transposon == "B30", 9,"NA")) %>%
    ## remove the empty wells.
    filter(!(Well %in% c('E4','E5','E6'))) %>%
    ## set the internal controls to Day 0.
    mutate(Day = ifelse(Treatment %in% c('INIT', 'B30'), 0, 9))

april.21.fig2 <- ggplot(results2.april.21,
                       aes(x = Day,
                           y = transposons.per.chromosome,
                           color = Treatment,
                           shape = Replicate)) +
    geom_point() +
    theme_classic() +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tet resistance causes transposition-mediated duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-04-21-fig2.pdf", april.21.fig2, width= 6.5, height=3.5)

april.21.fig3 <- ggplot(results3.april.21,
                       aes(x = Day,
                           y = transposons.per.chromosome,
                           color = Treatment,
                           shape = Replicate)) +
    geom_point() +
    theme_classic() +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tet resistance causes transposition-mediated duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-04-21-fig3.pdf", april.21.fig3, width= 6.5, height=3.5)

##########################################################

## analyze Day 5 and Day 9 data, revived from glycerol stock.
## Day 5 data, using 1:100 cell dilutions in PCR H20,
## qPCR experiment conducted on April 27 2021.

## Day 9 data, using 1:100 cell dilutions in PCR H20,
## qPCR experiment conducted on April 26 2021.

## Day 5 data.
data.april.27 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-04-27-reformatted.csv")

## Day 9 data.
data.april.26 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-04-26-reformatted.csv")

## Important: do the analysis for each day separately,
## in order to calibrate for batch effects.

results.april.27 <- data.april.27 %>%
    calc.tet.cm.CNV.with.INIT.control() %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    ## set Day to 5 for samples.
    mutate(Day = ifelse(Transposon == "B30", 5,"NA")) %>%
    ## set the internal controls to Day 0, and keep samples at day 5.
    mutate(Day = ifelse(Treatment %in% c('INIT', 'B30'), 0, 5)) %>%
    mutate(Date = "April 27")

results.april.26 <- data.april.26 %>%
    calc.tet.cm.CNV.with.INIT.control() %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    ## set Day to 9 for samples.
    mutate(Day = ifelse(Transposon == "B30", 9,"NA")) %>%
    ## set the internal controls to Day 0, and keep samples at day 9.
    mutate(Day = ifelse(Treatment %in% c('INIT', 'B30'), 0, 9)) %>%
    mutate(Date = "April 26")


results.april.26.and.27 <- full_join(results.april.27,results.april.26) %>%
    mutate(ExperimentReplicate = 1)


april.26.and.27.fig <- ggplot(results.april.26.and.27,
                       aes(x = Day,
                           y = transposons.per.chromosome,
                           color = Treatment,
                           shape = Replicate)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tet resistance causes transposition-mediated duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-04-28-initial-culture-calibration.pdf",
       april.26.and.27.fig,
       width= 6.5, height=4)

########################################################################

## analyze Day 5 and Day 9 data, revived from glycerol stock.
## Day 5 data, using 1:100 cell dilutions in PCR H20,
## qPCR experiment conducted on May 4 2021.

data.may.04 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-05-04-reformatted.csv")


## Important: do the analysis for each day separately,
## in order to calibrate for batch effects.

results.may.04 <- data.may.04 %>%
    calc.tet.cm.CNV.with.INIT.control() %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    ## hack to put the Day information back in.
    left_join(select(data.may.04,Well,Day)) %>%
    mutate(Date = "May 4") %>%
    mutate(ExperimentReplicate = 2)


may.04.fig <- ggplot(results.may.04,
                       aes(x = Day,
                           y = transposons.per.chromosome,
                           color = Treatment,
                           shape = Replicate)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tet resistance causes transposition-mediated duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-05-04-initial-culture-calibration.pdf",
       may.04.fig,
       width= 6.5, height=4)

########################################################################

## analyze Day 5 and Day 9 data, revived from glycerol stock.
## using 1:100 cell dilutions in PCR H20,
## qPCR experiment conducted on May 5 2021.
## This is a completely separate replication (fresh overnight cultures)
## from the May 4 2021 experiment.

data.may.05 <- read.csv("../data/draft-manuscript-1A/qPCR-data/Yi-transposon-qPCR-2021-05-05-reformatted.csv")


## Important: do the analysis for each day separately,
## in order to calibrate for batch effects.

results.may.05 <- data.may.05 %>%
    calc.tet.cm.CNV.with.INIT.control() %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    ## hack to put the Day information back in.
    left_join(select(data.may.05,Well,Day)) %>%
    mutate(Date = "May 5") %>%
    mutate(ExperimentReplicate = 3)


may.05.fig <- ggplot(results.may.05,
                       aes(x = Day,
                           y = transposons.per.chromosome,
                           color = Treatment,
                           shape = Replicate)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tet resistance causes transposition-mediated duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-05-05-initial-culture-calibration.pdf",
       may.05.fig,
       width= 6.5, height=4)


########################################################################
## Screen for isolates from Day 9 Tet 50 population 3 with miniTn5 duplications.

screen.data.may.09 <- read.csv("../data/draft-manuscript-1A/qPCR-data/2021-05-09-tet50-pop3-screen.csv")

results.may.09 <- screen.data.may.09 %>%
    calc.tet.cm.CNV.with.B30.control()

## colonies 2 and 4 show evidence of duplications.
may.09.fig <- ggplot(results.may.09,
                       aes(x = Treatment,
                           y = transposons.per.chromosome,
                           color = Sample)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    ylab("Transposons per chromosome") +
    ggtitle("May9 screen for B30 transposon duplications")

########################################################################
## make polished figure for 2021 Hartwell fellowship application.

hartwell.results <- rbind(results.april.26.and.27,
                          results.may.04,
                          results.may.05) %>%
    filter(!(Treatment %in% c("B30")))

hartwell.fig <- ggplot(hartwell.results,
                       aes(x = Day,
                           y = transposons.per.chromosome,
                           color = Treatment,
                           shape = Replicate)) +
    geom_point() +
    theme_classic() +
    guides(shape = FALSE) +
    scale_color_discrete(name = "tetracycline concentration\n(ug/mL)") +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    ylab("Transposons per chromosome") +
ggtitle("Selection for tetracycline resistance causes transposition-mediated duplications")

########################################################################
## Screen for isolates from Day 9 Tet 50 population 3 with miniTn5 duplications.

screen.data.may.25 <- read.csv("../data/draft-manuscript-1A/qPCR-data/2021-05-25_tet50_pop3_screen.csv")

results.may.25 <- screen.data.may.25 %>%
    calc.tet.cm.CNV.with.B30.control()


may.25.fig <- ggplot(results.may.25,
                       aes(x = Treatment,
                           y = transposons.per.chromosome,
                           color = Sample)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    ylab("Transposons per chromosome") +
    ggtitle("May25 screen for B30 transposon duplications")

########################################################################
##  Verify colony 5 and 6 isolates from Day 9 Tet 50 population 3 with miniTn5 duplications,
## by qPCR on purified gDNA.

data.june.5 <- read.csv("../data/draft-manuscript-1A/qPCR-data/2021-06-05_tet50_pop3_gDNA-verification.csv")

results.june.5 <- data.june.5 %>%
    calc.tet.cm.CNV.with.B30.control()

june.5.fig <- ggplot(results.june.5,
                       aes(x = Treatment,
                           y = transposons.per.chromosome,
                           color = Sample)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    ylab("Transposons per chromosome") +
    ggtitle("June 5 verification for B30 transposon duplications")



########################################################################
## analyze Day 10 relaxed selection pops, revived from glycerol stock.
## using 1:10 cell dilutions in PCR H20,
## qPCR experiment conducted on June 29 2021.

calc.tet.cm.CNV.with.INIT.control.june.29 <- function(df.with.INIT.control) {
    ## IMPORTANT: this function assumes perfect amplification efficiency
    ## (m = 1, or 2x amplification each cycle.)
    ## I can use a standard curve to determine m more precisely in the future.

    m.cm <- 2^1.055 ## calibrated with May 6 2021 standard curve.
    m.tet <- 2^0.904 ## calibrated with May 6 2021 standard curve.

    control.data <- filter(df.with.INIT.control, Treatment == "INIT")
    tet.control.data <- filter(control.data, probe == "T")
    cm.control.data <- filter(control.data, probe == "C")

    ## question: maybe I should be using geometric mean here instead?
    mean.control.tet.Cq <- mean(tet.control.data$cycle_at_threshold)
    mean.control.cm.Cq <- mean(cm.control.data$cycle_at_threshold)
    
    beta <- 2^(m.tet * mean.control.tet.Cq - m.cm * mean.control.cm.Cq)

    helper.func <- function(well.df) {
        ## note: this helper uses the beta variable defined above.
        C <- filter(well.df, probe == 'C')$cycle_at_threshold
        T <- filter(well.df, probe == 'T')$cycle_at_threshold
        T.per.C <- beta * 2^(m.cm * C - m.tet * T)

        return.df <- data.frame(Well = unique(well.df$Well),
                                Transposon = unique(well.df$Transposon),
                                Treatment = unique(well.df$Treatment),
                                Replicate = unique(well.df$Replicate),
                                Clone = unique(well.df$Clone),
                                transposons.per.chromosome = T.per.C)
        return(return.df)
    }

    
    results <- df.with.INIT.control %>%
        split(.$Well) %>%
        map_dfr(helper.func)
    return(results)
}


data.june.29 <- read.csv("../data/draft-manuscript-1A/qPCR-data/2021-06-29_B30-clone-relaxed-selection.csv")

results.june.29 <- data.june.29 %>%
    mutate(Clone = ifelse(Sample %in% c("RM6.139.1", "RM6.161.1","RM6.161.2","RM6.161.3",
                                        "RM6.161.4","RM6.161.5"), "A", "B")) %>%
    mutate(Clone = ifelse(Sample == "RM6.102.1", "INIT", Clone)) %>%
    calc.tet.cm.CNV.with.INIT.control.june.29() %>%
    mutate(Replicate = as.factor(Replicate))
                                

june.29.fig <- ggplot(results.june.29,
                       aes(x = Treatment,
                           y = transposons.per.chromosome,
                           color = Replicate,
                           shape = Clone)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    ylab("Transposons per chromosome") +
ggtitle("10 days of relaxed selection on clones with Tet duplications")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-2021-06-29-B30-duplication-clone-relaxed-selection.pdf",
       june.29.fig,
       width= 6.5, height=4)

############################################################
## analyze clones isolated from Tet50 B30 no plasmid and Tet50 B30 + A18 plasmid
## populations.
## using 1:10 cell dilutions in PCR H20,
## qPCR experiment conducted on July 13 2021.
## I FORGOT TO INCLUDE INIT CULTURE CONTROL SAMPLES!!!
## I will normalize based on the standard curve that Yi initially made.

well.to.clone.2021.07.13 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the clone replicate.
    rowletter <- substring(Well,1,1)
    return(rowletter)
}

data.july.13 <- read.csv("../data/draft-manuscript-1A/qPCR-data/2021-07-13_evolved-clones-qPCR.csv")



july.13.raw.figure <- ggplot(data.july.13, aes(x = cycle_at_threshold,
                                   y = cycle_at_threshold,
                                   color = probe)) +
    geom_point() +
    ylim(0,30) +
    theme_classic() +
    facet_wrap(Sample ~ .)

results.july.13 <- data.july.13 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Population = as.factor(Replicate)) %>%
    mutate(Clone = well.to.clone.2021.07.13(Well))

no.plasmid.results.july.13 <- results.july.13 %>%
    filter(Treatment == "A18_pUC_plasmid-")

has.plasmid.results.july.13 <- results.july.13 %>%
    filter(Treatment == "A18_pUC_plasmid+")

has.plasmid.fig1A <- ggplot(has.plasmid.results.july.13,
                              aes(x = Population,
                                  y = transposons.per.chromosome,
                                  color = Clone)) +
    geom_point() +
    guides(color=FALSE) +
    theme_classic() +
    ggtitle("pUC plasmid+: Transposons per chromosome")


has.plasmid.fig1B <- ggplot(has.plasmid.results.july.13,
                              aes(x = Population,
                                 y = transposons.per.plasmid,
                                 color = Clone)) +
    geom_point() +
    theme_classic() +
    ggtitle("pUC plasmid+: Transposons per plasmid") 

has.plasmid.fig1 <- plot_grid(has.plasmid.fig1A, has.plasmid.fig1B,
                              nrow=2, labels=c('A','B'))
ggsave("../results/draft-manuscript-1A/qPCR-results/has-plasmid-qPCR-07-13-2021.pdf", has.plasmid.fig1)

no.plasmid.figure1 <- ggplot(no.plasmid.results.july.13,
                              aes(x = Population,
                                  y = transposons.per.chromosome,
                                  color = Clone)) +
    geom_point() +
    theme_classic() + 
    ggtitle("no target plasmid: Transposons per chromosome")

ggsave("../results/draft-manuscript-1A/qPCR-results/no-plasmid-qPCR-07-13-2021.pdf", no.plasmid.figure1)

#################################################################################
## analyze Tet0 and Tet50 populations starting from A31 chromosomal integration.
## using 1:10 cell dilutions in PCR H20,
## qPCR experiment conducted on July 14 2021.
## for now, normalize based on the standard curve that Yi initially made.

A31.data.july.14 <- read.csv("../data/draft-manuscript-1A/qPCR-data/2021-07-14_A31_pops.csv")


A31.july.14.raw.figure <- ggplot(A31.data.july.14, aes(x = cycle_at_threshold,
                                   y = cycle_at_threshold,
                                   color = probe)) +
    geom_point() +
    ylim(0,30) +
    theme_classic() +
    facet_wrap(Sample ~ .)

A31.results.july.14 <- A31.data.july.14 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Population = as.factor(Replicate))


A31.figure1 <- ggplot(A31.results.july.14,
                              aes(x = Treatment,
                                  y = transposons.per.chromosome,
                                  color = Population,
                                  shape = Treatment)) +
    geom_point() +
    theme_classic() +
    ggtitle("Transposons per chromosome")


A31.figure2 <- ggplot(A31.results.july.14,
                              aes(x = Treatment,
                                 y = transposons.per.plasmid,
                                 color = Population,
                                 shape = Treatment)) +
    geom_point() +
    theme_classic() +
    ggtitle("Transposons per plasmid")

#################################################################################
## analyze all Tet50 pops, using undiluted gDNA and 1:10 culture dilutions in PCR H20.
## next time, dilute gDNA 1:100 before running!
## qPCR experiment conducted on July 18 2021.
## for now, normalize based on the standard curve that Yi initially made.

paired.gDNA.culture.data.july.28 <- read.csv("../data/draft-manuscript-1A/qPCR-data/2021-07-28_paired_gDNA-culture_qPCR.csv")

july.28.metadata <- paired.gDNA.culture.data.july.28 %>%
    select(Well, Treatment, Sample)

results.july.28 <- paired.gDNA.culture.data.july.28 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Population = as.factor(Replicate)) %>%
    mutate(transposons.per.plasmid = ifelse(Treatment == "no_plasmid", 0, transposons.per.plasmid)) %>%
    left_join(july.28.metadata) %>%
    mutate(is.gDNA = str_detect(Sample, "gDNA")) %>%
    mutate(template = ifelse(
               is.gDNA,"gDNA template","Culture template")) %>%
    filter(template == "Culture template") %>%
    mutate(Treatment = replace(
               Treatment, Treatment == "A18_plasmid", "pUC_plasmid")) %>%
    mutate(Treatment = replace(
               Treatment, Treatment == "A31_plasmid", "p15A_plasmid"))


figure1.july.28 <- ggplot(results.july.28,
                              aes(x = Treatment,
                                  y = transposons.per.chromosome,
                                  color = Population,
                                  shape = Treatment)) +
    geom_point() +
    theme_classic() +
    facet_wrap(.~template) +
    ggtitle("Transposons per chromosome")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-07-28-2021-fig1.pdf",figure1.july.28)

figure2.july.28 <- ggplot(results.july.28,
                              aes(x = Treatment,
                                 y = transposons.per.plasmid,
                                 color = Population,
                                 shape = Treatment)) +
    geom_point() +
    theme_classic() +
    facet_wrap(.~template) +
    ggtitle("Transposons per plasmid")

ggsave("../results/draft-manuscript-1A/qPCR-results/qPCR-07-28-2021-fig2.pdf",figure2.july.28)

################################################################
