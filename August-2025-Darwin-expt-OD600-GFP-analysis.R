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
    mutate(Tet = as.factor(Tet)) %>%
    ## rename Grayson's block to G-block and Rohan's block to R-block.
    mutate(Block = recode(Block, "Grayson" = "G-block", "Rohan" = "R-block"))


## Kanamycin- treatment goes into Figure 5F.
Fig5F <- August.2025.darwin.data %>%
    filter(Treatment == "Kan+") %>%
    ggplot(
        aes(x=Day, y=log10(normalized.GFP55), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    ylab("log10(GFP/OD600)") +
    ggtitle("August 2025 experiment, Kanamycin+ treatment") +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    facet_grid(Block~Replicate) +
    theme(legend.position = "top") +
    theme(axis.text.x = element_text(size = 6)) +
    scale_x_continuous(breaks = 0:12)


## Kanamycin+ treatment goes into Supplementary Figure S7.
S7Fig <- August.2025.darwin.data %>%
    filter(Treatment == "Kan-") %>%
    ggplot(
        aes(x=Day, y=log10(normalized.GFP55), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    geom_line() +
    geom_vline(xintercept = AUG_2025_TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme_classic() +
    ylab("log10(GFP/OD600)") +
    ggtitle("August 2025 experiment, Kanamycin- treatment") +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    facet_grid(Block~Replicate) +
    theme(legend.position = "top") +
    theme(axis.text.x = element_text(size = 6)) +
    scale_x_continuous(breaks = 0:12)

## save figures.
ggsave("../results/Fig5F.pdf", Fig5F, height=3.5, width=8)
ggsave("../results/S7Fig.pdf", S7Fig, height=3.5,width=8)


