## full-nine-day-experiment-OD600-GFP.R by Rohan Maddamsetti.

## This script plots growth and GFP over time for the 30 evolving populations,
## 10 with no plasmid, 10 with p15A plasmid, 10 with pUC plasmid.
## with the tetA-GFP-barcode transposons.

## This script also plots OD600 and GFP for the 30 clones that I isolated.

library(tidyverse)
library(cowplot)
library(broom)

### Make a plot of the Tet antibiotic treatment over time.
experiment.design.df <- data.frame(Day = c(1,2,3,4,5,6,7,8,9,10),
                                   Tet.conc = c(0,2,4,6,8,10,20,30,40,50))
treatment.plot <- ggplot(experiment.design.df,
                         aes(x = Day, y = Tet.conc)) +
    theme_classic() + geom_line() +
    ylab("Tetracycline (ug/mL)") +
    xlab("Day")


subtract.OD600.and.GFP.blanks <- function(OD600.df) {
    blanks.df <- filter(OD600.df, Treatment == "Blank")
    OD600.media.blank <- mean(blanks.df$RawOD600)
    GFP55.media.blank <- mean(blanks.df$RawGFP55)
    GFP100.media.blank <- mean(blanks.df$RawGFP100)
    subtracted.df <- OD600.df %>%
        filter(Treatment != "Blank") %>%
        mutate(OD600 = RawOD600 - OD600.media.blank) %>%
        mutate(GFP55 = RawGFP55 - GFP55.media.blank) %>%
        mutate(GFP100 = RawGFP100 - GFP100.media.blank)
    return(subtracted.df)
}

################################################################################
## plot OD600 and GFP over time in the nine day evolution experiment conducted February 2023.

Feb.2023.nine.day.tetA.GFP.OD600.df <- read.csv("../data/2023-B30-tetA-GFP-barcode-9-day-evolution-OD600-GFP-measurements.csv") %>%
    ## subtract LB blanks to get the true OD measurement
    subtract.OD600.and.GFP.blanks() %>%
    mutate(BiologicalReplicate = as.factor(BiologicalReplicate)) %>%
    filter(!(Treatment %in% c("Blank"))) %>%
    mutate(normalized.GFP55 = GFP55/OD600) %>%
    mutate(normalized.GFP100 = GFP100/OD600)
    
Feb.2023.OD600.plot <- ggplot(Feb.2023.nine.day.tetA.GFP.OD600.df,
                              aes(x = Day, y = OD600, color = Treatment)) +
    theme_classic() + geom_point() + geom_smooth(se=FALSE) +
    scale_x_continuous(breaks=seq(0, 10, 1)) + guides(color="none")

Feb.2023.normalized.GFP55.plot <- ggplot(Feb.2023.nine.day.tetA.GFP.OD600.df,
                              aes(x = Day, y = normalized.GFP55, color = Treatment)) +
    theme_classic() + geom_point() + geom_smooth(se=FALSE) +
    scale_x_continuous(breaks=seq(0, 10, 1)) +
    ylab("GFP/OD600") +
    theme(legend.position="bottom")

## This is basically the same as the GFP55 plot, so omit from the combined figure.
Feb.2023.normalized.GFP100.plot <- ggplot(Feb.2023.nine.day.tetA.GFP.OD600.df,
                              aes(x = Day, y = normalized.GFP100, color = Treatment)) +
    theme_classic() + geom_point() + geom_smooth(se=FALSE) +
    scale_x_continuous(breaks=seq(0, 10, 1)) +
    ylab("GFP/OD600")

Feb.2023.evolution.figure <- plot_grid(Feb.2023.OD600.plot, Feb.2023.normalized.GFP55.plot,
                                       labels=c('E','F'), ncol=1)
ggsave("../results/Full-9-day-evolution-experiment-OD600-February-2023.pdf",
       Feb.2023.evolution.figure, width=4.5, height=4.5)

################################################################################
## let's examine the trajectory of GFP expression in the no plasmid populations.
no.plasmid.pop.GFP.OD.data <- Feb.2023.nine.day.tetA.GFP.OD600.df %>%
    filter(Treatment == "no_plasmid")

p15A.pop.GFP.OD.data <- Feb.2023.nine.day.tetA.GFP.OD600.df %>%
    filter(Treatment == "p15A")

pUC.pop.GFP.OD.data <- Feb.2023.nine.day.tetA.GFP.OD600.df %>%
    filter(Treatment == "pUC")


make.my.test.plot <- function(my.data) {
    ggplot(my.data,
           aes(x = Day, y = normalized.GFP55, color = Treatment)) +
        theme_classic() + geom_point() + geom_smooth(se=FALSE) +
        scale_x_continuous(breaks=seq(0, 10, 1)) +
        ylab("GFP/OD600") +
        theme(legend.position="bottom")
}

no.plasmid.Feb.2023.normalized.GFP55.plot <- make.my.test.plot(no.plasmid.pop.GFP.OD.data)
p15A.Feb.2023.normalized.GFP55.plot <- make.my.test.plot(p15A.pop.GFP.OD.data)
pUC.Feb.2023.normalized.GFP55.plot <- make.my.test.plot(pUC.pop.GFP.OD.data)

no.plasmid.Feb.2023.normalized.GFP55.plot 

p15A.Feb.2023.normalized.GFP55.plot

pUC.Feb.2023.normalized.GFP55.plot 



## let's look at all of them now.
faceted.Feb.2023.normalized.GFP55.plot <- ggplot(Feb.2023.nine.day.tetA.GFP.OD600.df,
                              aes(x = Day, y = normalized.GFP55, color = Treatment)) +
    theme_classic() + geom_point() + geom_smooth(se=FALSE) + facet_wrap(.~Treatment) +
    scale_x_continuous(breaks=seq(0, 10, 1)) +
    ylab("GFP/OD600") +
    theme(legend.position="bottom")

faceted.Feb.2023.normalized.GFP55.plot

################################################################################
## plot OD600 and GFP for the single clones isolated from each of the 30 populations of the
## evolution experiment conducted February 2023.

clones.Feb.2023.tetA.GFP.OD600.df <- read.csv("../data/2023-B30-tetA-GFP-barcode-Day9-clones-OD600-GFP-measurements.csv") %>%
    ## subtract LB blanks to get the true OD measurement
    subtract.OD600.and.GFP.blanks() %>%
    mutate(BiologicalReplicate = as.factor(BiologicalReplicate)) %>%
    filter(!(Treatment %in% c("Blank"))) %>%
    mutate(normalized.GFP55 = GFP55/OD600) %>%
    mutate(normalized.GFP100 = GFP100/OD600)

clones.Feb.2023.OD600.plot <- ggplot(clones.Feb.2023.tetA.GFP.OD600.df,
                              aes(x = Treatment, y = OD600, color = Treatment)) +
    theme_classic() + geom_point() + guides(color="none")

clones.Feb.2023.normalized.GFP55.plot <- ggplot(clones.Feb.2023.tetA.GFP.OD600.df,
                              aes(x = Treatment, y = normalized.GFP55, color = Treatment)) +
    theme_classic() + geom_point() +
    ylab("GFP/OD600") +
    theme(legend.position="bottom")

Feb.2023.clone.figure <- plot_grid(clones.Feb.2023.OD600.plot, clones.Feb.2023.normalized.GFP55.plot,
                                       labels=c('G','H'), ncol=1)
ggsave("../results/9-day-clones-OD600-February-2023.pdf",
       Feb.2023.clone.figure,
       width=4.5, height=4.5)

## Let's see if there's any correlation between OD600 and GFP/OD600.
## nothing to write home about here.
clones.normalized.GFP.OD600.correlation.plot <- ggplot(clones.Feb.2023.tetA.GFP.OD600.df,
                                                       aes(x = OD600, y = normalized.GFP55, color = Treatment)) +
    theme_classic() + geom_point() +
    ylab("GFP/OD600") +
    theme(legend.position="bottom")
ggsave("../results/9-day-clones-GFP-OD600-correlation-February-2023.pdf",
       clones.normalized.GFP.OD600.correlation.plot,
       width=4.5, height=4.5)


