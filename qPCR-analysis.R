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


######################################################################
## analyze clones isolated from Tet50 B30 no plasmid and Tet50 B30 + A18 plasmid
## populations.
## using 1:10 cell dilutions in PCR H20,
## qPCR experiment conducted on July 13 2021.
## I will normalize based on the standard curve that Yi initially made.

well.to.clone.2021.07.13 <- function(Well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the clone replicate.
    rowletter <- substring(Well,1,1)
    return(rowletter)
}


data.july.13 <- read.csv("../data/qPCR-data/2021-07-13_evolved-clones-qPCR.csv")

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
    guides(color = "none") +
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
ggsave("../results/qPCR-results/has-plasmid-qPCR-07-13-2021.pdf", has.plasmid.fig1)

no.plasmid.figure1 <- ggplot(no.plasmid.results.july.13,
                              aes(x = Population,
                                  y = transposons.per.chromosome,
                                  color = Clone)) +
    geom_point() +
    theme_classic() + 
    ggtitle("no target plasmid: Transposons per chromosome")

ggsave("../results/qPCR-results/no-plasmid-qPCR-07-13-2021.pdf", no.plasmid.figure1)

#################################################################################
## analyze all Tet50 pops, using undiluted gDNA and 1:10 culture dilutions in PCR H20.
## next time, dilute gDNA 1:100 before running!
## qPCR experiment conducted on July 18 2021.
## normalize based on the standard curve that Yi initially made.

paired.gDNA.culture.data.july.28 <- read.csv("../data/qPCR-data/2021-07-28_paired_gDNA-culture_qPCR.csv")

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

ggsave("../results/qPCR-results/qPCR-07-28-2021-fig1.pdf",figure1.july.28)

figure2.july.28 <- ggplot(results.july.28,
                              aes(x = Treatment,
                                 y = transposons.per.plasmid,
                                 color = Population,
                                 shape = Treatment)) +
    geom_point() +
    theme_classic() +
    facet_wrap(.~template) +
    ggtitle("Transposons per plasmid")

ggsave("../results/qPCR-results/qPCR-07-28-2021-fig2.pdf",figure2.july.28)

