 ## qPCR-analysis.R by Rohan Maddamsetti.
## This script analyzes qPCR data using Yi's transposon and antibiotic resistance marker system.


library(tidyverse)
library(cowplot)


## Days when Tet 10 selection was applied in the August 2025 experiment.
AUG_2025_TET_SELECTION_DAYS <- c(
    0,
    ##3, ## omit day 3 since we did not collect data on that day.
    10)


## Plasmid color scheme:
## https://colorbrewer2.org/#type=sequential&scheme=PuRd&n=5
PLASMID_COLORSCALE <- c(
    "pUC" = "#7a0177",
    "CloDF13" = "#c51b8a",
    "pBR322" = "#f768a1",
    "p15A" = "#fbb4b9"
)

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

## Day 9 data.
Day9.data.August.11.2025 <- read.csv("../data/qPCR-data/2025-08-11_Rohan_Grayson_Darwin_Expt_Day9-Run4.csv")

Day9.results <- Day9.data.August.11.2025 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences.for.2025.data) %>%
    mutate(Day = 9)

## Day 10 data.
Day10.data.August.12.2025 <- read.csv("../data/qPCR-data/2025-08-12_Rohan_Grayson_Darwin_Expt_Day10.csv")

Day10.results <- Day10.data.August.12.2025 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences.for.2025.data) %>%
    mutate(Day = 10)

## Day 12 data.
## Note that the original data from the qPCR machine has the wrong day (should be Day 12 but Day 13 written in title).
Day12.data.August.14.2025 <- read.csv("../data/qPCR-data/2025-08-14_Rohan_Grayson_Darwin_Expt_Day12.csv")

Day12.results <- Day12.data.August.14.2025 %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences.for.2025.data) %>%
    mutate(Day = 12)


## Generate NA data for the days that were not sampled.
NA.data.template <- Day12.results %>%
    mutate(transposons.per.chromosome = NA) %>%
    mutate(plasmids.per.chromosome = NA) %>%
    mutate(transposons.per.plasmid = NA) %>%
    mutate(Day = NA)

Day3.results <- NA.data.template %>%
    mutate(Day = 3)

Day4.results <- NA.data.template %>%
    mutate(Day = 4)

Day5.results <- NA.data.template %>%
    mutate(Day = 5)

Day6.results <- NA.data.template %>%
    mutate(Day = 6)

Day7.results <- NA.data.template %>%
    mutate(Day = 7)

Day8.results <- NA.data.template %>%
    mutate(Day = 8)



## Combine data across days.
full.results <- rbind(
    ## actual data here
    Day0.results, Day1.results, Day2.results,
    ## NA missing data here
    ##Day3.results, Day4.results, Day5.results, Day6.results, Day7.results, Day8.results,
    ## actual data here
    Day9.results, Day10.results, Day12.results) %>%
    mutate(Replicate=as.factor(Replicate)) %>%
    mutate(Day=as.factor(Day)) %>%
    mutate(Plasmid=factor(Plasmid, levels = c("pUC", "CloDF13", "pBR322", "p15A"))) %>%
    ## it looks like the transposon did not jump onto the plasmid of the p15A clone.
    ## let's remove this from the analysis.
    filter(Plasmid != "p15A")


Grayson.plasmids.per.chromosome.fig <- full.results %>%
    filter(Block == "Grayson") %>%
    ggplot(
        aes(x=Day, y=log10(plasmids.per.chromosome), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    facet_grid(Treatment~Replicate) +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    ylab("log10(plasmids per chromosome)") +
    ##    scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
    theme(strip.background=element_blank(),
          axis.text.x = element_text(size=10),
          strip.text.x = element_text(size = 14)) +
    ggtitle("Grayson's Block: plasmids per chromosome") +
    theme(legend.position = "bottom")

Grayson.plasmids.per.chromosome.fig
ggsave("../results/August2025-Grayson-plasmids-per-chromosome-fig.pdf", Grayson.plasmids.per.chromosome.fig)


Grayson.transposons.per.plasmid.fig <- full.results %>%
    filter(Block == "Grayson") %>%
    ggplot(
        aes(x=Day, y=transposons.per.plasmid, color=Plasmid, group = interaction(Replicate, Plasmid))) +
    facet_grid(Treatment~Replicate) +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    ylab("transposons per plasmid") +
    ##    scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
    theme(strip.background=element_blank(),
          axis.text.x = element_text(size=10),
          strip.text.x = element_text(size = 14)) +
    ggtitle("Grayson's Block: transposons per plasmid") +
    theme(legend.position = "bottom")

Grayson.transposons.per.plasmid.fig
ggsave("../results/August2025-Grayson-transposons-per-plasmid-fig.pdf", Grayson.transposons.per.plasmid.fig)


Grayson.transposons.per.chromosome.fig <- full.results %>%
    filter(Block == "Grayson") %>%
    ggplot(
        aes(x=Day, y=log10(transposons.per.chromosome), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    facet_grid(Treatment~Replicate) +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    ylab("log10(transposons per chromosome)") +
    ##    scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
    theme(strip.background=element_blank(),
          axis.text.x = element_text(size=10),
          strip.text.x = element_text(size = 14)) +
    ggtitle("Grayson's Block: transposons per chromosome") +
    theme(legend.position = "bottom")

Grayson.transposons.per.chromosome.fig
ggsave("../results/August2025-Grayson-transposons-per-chromosome-fig.pdf", Grayson.transposons.per.chromosome.fig)




Rohan.plasmids.per.chromosome.fig <- full.results %>%
    filter(Block == "Rohan") %>%
    ggplot(
        aes(x=Day, y=log10(plasmids.per.chromosome), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    facet_grid(Treatment~Replicate) +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    ylab("log10(plasmids per chromosome)") +
    ##    scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
    theme(strip.background=element_blank(),
          axis.text.x = element_text(size=10),
          strip.text.x = element_text(size = 14)) +
    ggtitle("Rohan's Block: plasmids per chromosome") +
    theme(legend.position = "bottom")

Rohan.plasmids.per.chromosome.fig
ggsave("../results/August2025-Rohan-plasmids-per-chromosome-fig.pdf", Rohan.plasmids.per.chromosome.fig)


Rohan.transposons.per.plasmid.fig <- full.results %>%
    filter(Block == "Rohan") %>%
    ggplot(
        aes(x=Day, y=transposons.per.plasmid, color=Plasmid, group = interaction(Replicate, Plasmid))) +
    facet_grid(Treatment~Replicate) +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    ylab("transposons per plasmid") +
    ##    scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
    theme(strip.background=element_blank(),
          axis.text.x = element_text(size=10),
          strip.text.x = element_text(size = 14)) +
    ggtitle("Rohan's Block: transposons per plasmid") +
    theme(legend.position = "bottom")

Rohan.transposons.per.plasmid.fig
ggsave("../results/August2025-Rohan-transposons-per-plasmid-fig.pdf", Rohan.transposons.per.plasmid.fig)


Rohan.transposons.per.chromosome.fig <- full.results %>%
    filter(Block == "Rohan") %>%
    ggplot(
        aes(x=Day, y=log10(transposons.per.chromosome), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    facet_grid(Treatment~Replicate) +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    ylab("log10(transposons per chromosome)") +
    ##    scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
    theme(strip.background=element_blank(),
          axis.text.x = element_text(size=10),
          strip.text.x = element_text(size = 14)) +
    ggtitle("Rohan's Block: transposons per chromosome") +
    theme(legend.position = "bottom")

Rohan.transposons.per.chromosome.fig
ggsave("../results/August2025-Rohan-transposons-per-chromosome-fig.pdf", Rohan.transposons.per.chromosome.fig)
