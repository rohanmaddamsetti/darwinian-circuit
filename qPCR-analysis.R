 ## qPCR-analysis.R by Rohan Maddamsetti.
## This script analyzes qPCR data using Yi's transposon and antibiotic resistance marker system.


library(tidyverse)
library(cowplot)


##################################################################
## functions for analysis of 2021 data.
## These data are used for Figure 2D in the manuscript.

calc.probe.fold.differences.for.2021.data <- function(well.df) {
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


well.to.clone.2021.07.13 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the clone replicate.
    rowletter <- substring(Well,1,1)
    return(rowletter)
}


######################################################################
## analyze clones isolated from Tet50 B30 no plasmid and Tet50 B30 + A18 plasmid
## populations.
## using 1:10 cell dilutions in PCR H20,
## qPCR experiment conducted on July 13 2021.
## Normalize based on the standard curve that Yi made.
## For the manuscript, we only care about the Tet50 B30 + A18 pUC plasmid clones.

results.july.13 <- read.csv("../data/qPCR-data/2021-07-13_evolved-clones-qPCR.csv") %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences.for.2021.data) %>%
    mutate(Population = as.factor(Replicate)) %>%
    mutate(Clone = well.to.clone.2021.07.13(Well))

## For the manuscript, we only care about the Tet50 B30 + A18 pUC plasmid clones.
Fig2E.data <- results.july.13 %>%
    filter(Treatment == "A18_pUC_plasmid+") %>%
    mutate(Population = as.factor(Population))

Fig2E <- Fig2E.data %>%
    ## update the names of the Population factor for a prettier plot.
    mutate(Population = fct_recode(as.factor(Population),
                                   `Population 1` = "1",
                                   `Population 2` = "2",
                                   `Population 3` = "3",
                                   `Population 4` = "4",
                                   `Population 5` = "5")) %>%
    ggplot(
        aes(x=Clone, y=Yi.transposon.frac)) +
    geom_point() +
    facet_grid(.~Population) +
    theme_classic() +
    geom_hline(yintercept=1, linetype="dashed", color="gray") +
    ylab("transposons per pUC plasmid") +
    scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
    theme(strip.background=element_blank(),
          axis.text.x = element_text(size=10),
          strip.text.x = element_text(size = 14))

ggsave("../results/Fig2E.pdf", Fig2E, width=7, height=2.25)


######################################################################
## analysis of 2025 Darwin circuit experiment qPCR data collected with Grayson Hamrick.

##################################################################
## functions for analysis of 2025 data

calc.probe.fold.differences.for.2025.data <- function(well.df) {
    ## this is a helper function for calculating probe fold differences
    ## per well.

    ## Use the calculation method in Yuanchi's paper.
    ## data analysis using constants calculated from Yuanchi's standard curve calibration,
    ## as reported in her paper.
    CmR.amplification.factor <- 2.01
    KanR.amplification.factor <- 2.03
    TetR.amplification.factor <- 1.99
    
    C <- filter(well.df, probe_target == 'cmR')$cycle_at_threshold
    T <- filter(well.df, probe_target == 'tetA')$cycle_at_threshold
    K <- filter(well.df, probe_target == 'kanR')$cycle_at_threshold

    ## subtract one to account for the tetA copy on the chromosome.
    K.per.C <- (CmR.amplification.factor^C)/(KanR.amplification.factor^K) ##2^(C - K)/K.per.C.constant
    T.per.C <- (CmR.amplification.factor^C)/(TetR.amplification.factor^T - 1) ##2^(C - T)/T.per.C.constant
    T.per.K <- (KanR.amplification.factor^K)/(TetR.amplification.factor^T - 1)##2^(K - T)/T.per.K.constant


    return.df <- data.frame(
        Well = unique(well.df$Well),
        Transposon = unique(well.df$Transposon),
        Sample = unique(well.df$Sample),
        Plasmid = unique(well.df$Plasmid),
        Treatment = unique(well.df$Treatment),
        Replicate = unique(well.df$Replicate),
        Block = unique(well.df$Block),
        transposons.per.chromosome = T.per.C,
        plasmids.per.chromosome = K.per.C,
        transposons.per.plasmid = T.per.K
    )
    
    return(return.df)
}


##################################################################
## analyze Darwin circuit experiment data from work with Grayson in August 2025.
## normalize based on the standard curve that Yuanchi made for her paper (see bioRxiv preprint).

## Day 0 data.
Day0.data.August.2.2025 <- read.csv("../data/qPCR-data/2025-08-02_Rohan_Grayson_Darwin_Expt_Day0.csv")

Day0.results <- Day0.data.August.2.2025 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences.for.2025.data) %>%
    mutate(Day = 0)

## Day 1 data.
Day1.data.August.3.2025 <- read.csv("../data/qPCR-data/2025-08-03_Rohan_Grayson_Darwin_Expt_Day1.csv")

Day1.results <- Day1.data.August.3.2025 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences.for.2025.data) %>%
    mutate(Day = 1)

## Day 2 data.
Day2.data.August.4.2025 <- read.csv("../data/qPCR-data/2025-08-04_Rohan_Grayson_Darwin_Expt_Day2.csv")

Day2.results <- Day2.data.August.4.2025 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences.for.2025.data) %>%
    mutate(Day = 2)

## Combine data across days.
full.results <- rbind(Day0.results, Day1.results, Day2.results) %>%
    mutate(Day = as.factor(Day)) %>%
    mutate(Replicate=as.factor(Replicate)) %>%
    mutate(Plasmid=factor(Plasmid, levels = c("pUC", "CloDF13", "pBR322", "p15A")))

plasmids.per.chromosome.fig <- full.results %>%
    ggplot(
        aes(x=Day, y=log10(plasmids.per.chromosome), color=Replicate, shape=Block)) +
    facet_grid(Treatment~Plasmid) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=1, linetype="dashed", color="gray") +
    ylab("log10(plasmids per chromosome)") +
    ##    scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
    theme(strip.background=element_blank(),
          axis.text.x = element_text(size=10),
          strip.text.x = element_text(size = 14)) +
    ggtitle("plasmids per chromosome")

plasmids.per.chromosome.fig
ggsave("../results/August2025-plasmids-per-chromosome-fig.pdf", plasmids.per.chromosome.fig)


transposons.per.plasmid.fig <- full.results %>%
    ggplot(
        aes(x=Day, y=transposons.per.plasmid, color=Replicate, shape=Block)) +
    facet_grid(Treatment~Plasmid) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=1, linetype="dashed", color="gray") +
    ylab("transposons per plasmid") +
    ##    scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
    theme(strip.background=element_blank(),
          axis.text.x = element_text(size=10),
          strip.text.x = element_text(size = 14)) +
    ggtitle("transposons per plasmid")

transposons.per.plasmid.fig
ggsave("../results/August2025-transposons-per-plasmid-fig.pdf", transposons.per.plasmid.fig)


transposons.per.chromosome.fig <- full.results %>%
    ggplot(
        aes(x=Day, y=log10(transposons.per.chromosome), color=Replicate, shape=Block)) +
    facet_grid(Treatment~Plasmid) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=1, linetype="dashed", color="gray") +
    ylab("log10(transposons per chromosome)") +
    ##    scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
    theme(strip.background=element_blank(),
          axis.text.x = element_text(size=10),
          strip.text.x = element_text(size = 14)) +
    ggtitle("transposons per chromosome")

transposons.per.chromosome.fig
ggsave("../results/August2025-transposons-per-chromosome-fig.pdf", transposons.per.chromosome.fig)
