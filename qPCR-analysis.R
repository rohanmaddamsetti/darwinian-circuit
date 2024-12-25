 ## qPCR-analysis.R by Rohan Maddamsetti.
## This script analyzes qPCR data using Yi's transposon
## and antibiotic resistance marker system.

## These data are used for Figure 2D in the manuscript.

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
    map_dfr(calc.probe.fold.differences) %>%
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


