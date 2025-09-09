## August-2025-Darwin-expt-OD600-GFP-analysis.R
## by Rohan Maddamsetti.

## Note: measurements "OVER" the GFP100 measurement range are encoded as 63000
## as a ceiling for all real Tecan Infinite Pro GFP100 readings.

library(tidyverse)


## Days when Tet 10 selection was applied in the August 2025 experiment.
AUG_2025_TET_SELECTION_DAYS <- c(0,3,10)


## Plasmid color scheme.
## https://colorbrewer2.org/#type=sequential&scheme=PuRd&n=5
PLASMID_COLORSCALE <- c(
    "pUC" = "#7a0177",
    "CloDF13" = "#c51b8a",
    "pBR322" = "#f768a1",
    "p15A" = "#fbb4b9"
)


## at first, I'm not subtracting the blanks.
## TODO: SUBTRACT BLANK READINGS!!!
August.2025.darwin.data <- read.csv("../data/OD600-GFP-data/2025-August-Darwin-expt-OD600-GFP-data.csv") %>%
    mutate(Plasmid = factor(Plasmid, levels=c("pUC", "CloDF13", "pBR322", "p15A"))) %>%
    mutate(normalized.GFP55 = GFP55/OD600) %>%
    mutate(normalized.GFP100 = GFP100/OD600) %>%
    mutate(log.normalized.GFP55 = log10(GFP55/OD600)) %>%
    mutate(log.normalized.GFP100 = log10(GFP100/OD600)) %>%
    ## let's remove the blank measurements for now.
    filter(Sample != "Blank") %>%
    mutate(Tet = as.factor(Tet)) 
    

Rohan.GFP55.timecourse.plot <- August.2025.darwin.data %>%
    filter(Block == "Rohan") %>%
    ggplot(
        aes(x=Day, y=log10(normalized.GFP55), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Rohan's Block") +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    theme(legend.position = "bottom")

Rohan.GFP100.timecourse.plot <- August.2025.darwin.data %>%
    filter(Block == "Rohan") %>%
    ggplot(
        aes(x=Day, y=log10(normalized.GFP100), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Rohan's Block") +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    theme(legend.position = "bottom")

Grayson.GFP55.timecourse.plot <- August.2025.darwin.data %>%
    filter(Block == "Grayson") %>%
    ggplot(
        aes(x=Day, y=log10(normalized.GFP55), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Grayson's Block") +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    theme(legend.position = "bottom")

Grayson.GFP100.timecourse.plot <- August.2025.darwin.data %>%
    filter(Block == "Grayson") %>%
    ggplot(
        aes(x=Day, y=log10(normalized.GFP100), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Grayson's Block") +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    theme(legend.position = "bottom")


ggsave("../results/August-2025-Rohan-GFP55-plot.pdf",Rohan.GFP55.timecourse.plot, height=4, width=7)
ggsave("../results/August-2025-Rohan-GFP100-plot.pdf",Rohan.GFP100.timecourse.plot, height=4, width=7)

ggsave("../results/August-2025-Grayson-GFP55-plot.pdf",Grayson.GFP55.timecourse.plot, height=4,  width=7)
ggsave("../results/August-2025-Grayson-GFP100-plot.pdf",Grayson.GFP100.timecourse.plot, height=4, width=7)


#####################################################
## normalize GFP values by the p15A samples, which don't have tetA-GFP on plasmid.
## assume that there is one tetA-GFP copy on the chromosome.
## This can be an internal control for double-checking qPCR calculations.

p15A.darwin.data.for.normalization <- August.2025.darwin.data %>%
    filter(Plasmid == "p15A") %>%
    select(Replicate, Treatment, Block, uL_measurement, transfer_dilution, Tet, Day, Date,
           normalized.GFP55, normalized.GFP100, log.normalized.GFP55, log.normalized.GFP100) %>%
    rename(
        p15A.normalized.GFP55 = normalized.GFP55,
        p15A.normalized.GFP100 = normalized.GFP100,
        p15A.log.normalized.GFP55 = log.normalized.GFP55,
        p15A.log.normalized.GFP100 = log.normalized.GFP100)


p15A.normalized.Rohan.GFP55.timecourse.plot <- August.2025.darwin.data %>%
    ## add additional columns from p15A.darwin.data.for.normalization
    left_join(p15A.darwin.data.for.normalization) %>%
    mutate(p15A.normalized.GFP55 = 10^(log.normalized.GFP55 - p15A.log.normalized.GFP55)) %>%
    mutate(p15A.normalized.GFP100 = 10^(log.normalized.GFP100 - p15A.log.normalized.GFP100)) %>%   
    filter(Block == "Rohan") %>%
    ggplot(
        aes(x=Day, y=p15A.normalized.GFP55, color=Plasmid, group = interaction(Replicate, Plasmid))) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Rohan's Block") +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    theme(legend.position = "bottom")

p15A.normalized.Rohan.GFP100.timecourse.plot <- August.2025.darwin.data %>%
    ## add additional columns from p15A.darwin.data.for.normalization
    left_join(p15A.darwin.data.for.normalization) %>%
    mutate(p15A.normalized.GFP55 = 10^(log.normalized.GFP55 - p15A.log.normalized.GFP55)) %>%
    mutate(p15A.normalized.GFP100 = 10^(log.normalized.GFP100 - p15A.log.normalized.GFP100)) %>%
    filter(Block == "Rohan") %>%
    ggplot(
        aes(x=Day, y=p15A.normalized.GFP100, color=Plasmid, group = interaction(Replicate, Plasmid))) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Rohan's Block") +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    theme(legend.position = "bottom")

p15A.normalized.Grayson.GFP55.timecourse.plot <- August.2025.darwin.data %>%
    ## add additional columns from p15A.darwin.data.for.normalization
    left_join(p15A.darwin.data.for.normalization) %>%
    mutate(p15A.normalized.GFP55 = 10^(log.normalized.GFP55 - p15A.log.normalized.GFP55)) %>%
    mutate(p15A.normalized.GFP100 = 10^(log.normalized.GFP100 - p15A.log.normalized.GFP100)) %>%
    filter(Block == "Grayson") %>%
    ggplot(
        aes(x=Day, y=p15A.normalized.GFP55, color=Plasmid, group = interaction(Replicate, Plasmid))) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Grayson's Block") +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    theme(legend.position = "bottom")

p15A.normalized.Grayson.GFP100.timecourse.plot <- August.2025.darwin.data %>%
    left_join(p15A.darwin.data.for.normalization) %>%
    mutate(p15A.normalized.GFP55 = 10^(log.normalized.GFP55 - p15A.log.normalized.GFP55)) %>%
    mutate(p15A.normalized.GFP100 = 10^(log.normalized.GFP100 - p15A.log.normalized.GFP100)) %>%
    filter(Block == "Grayson") %>%
    ggplot(
        aes(x=Day, y=p15A.normalized.GFP100, color=Plasmid, group = interaction(Replicate, Plasmid))) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Grayson's Block") +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    theme(legend.position = "bottom")


ggsave("../results/August-2025-Rohan-p15A-normalized-GFP55-plot.pdf",p15A.normalized.Rohan.GFP55.timecourse.plot, height=4, width=7)
ggsave("../results/August-2025-Rohan-p15A-normalized-GFP100-plot.pdf",p15A.normalized.Rohan.GFP100.timecourse.plot, height=4, width=7)

ggsave("../results/August-2025-Grayson-p15A-normalized-GFP55-plot.pdf",p15A.normalized.Grayson.GFP55.timecourse.plot, height=4,  width=7)
ggsave("../results/August-2025-Grayson-p15A-normalized-GFP100-plot.pdf",p15A.normalized.Grayson.GFP100.timecourse.plot, height=4, width=7)

