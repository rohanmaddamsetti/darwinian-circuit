## August-2025-Darwin-expt-OD600-GFP-analysis.R
## by Rohan Maddamsetti.

library(tidyverse)

## let's analyze data collected with Grayson.
## Note: measurements "OVER" the GFP100 measurement range are encoded as 80000
## as a guessed maximum-ish value.

## at first, I'm not subtracting the blanks.
## TODO: SUBTRACT BLANK READINGS!!!
August.2025.darwin.data <- read.csv("../data/OD600-GFP-data/2025-August-Darwin-expt-OD600-GFP-data.csv") %>%
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
    aes(x=Day, y=log10(normalized.GFP55), color=Plasmid, shape=Block)) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Rohan's Block") +
    geom_line() + theme_classic()

Rohan.GFP100.timecourse.plot <- August.2025.darwin.data %>%
    filter(Block == "Rohan") %>%
    ggplot(
    aes(x=Day, y=log10(normalized.GFP100), color=Plasmid, shape=Block)) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Rohan's Block") +
    geom_line() + theme_classic()

Grayson.GFP55.timecourse.plot <- August.2025.darwin.data %>%
    filter(Block == "Grayson") %>%
    ggplot(
    aes(x=Day, y=log10(normalized.GFP55), color=Plasmid, shape=Block)) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Grayson's Block") +
    geom_line() + theme_classic()

Grayson.GFP100.timecourse.plot <- August.2025.darwin.data %>%
    filter(Block == "Grayson") %>%
    ggplot(
    aes(x=Day, y=log10(normalized.GFP100), color=Plasmid, shape=Block)) +
    facet_grid(Treatment~Replicate) +
    theme(legend.position = "bottom") +
    ggtitle("Grayson's Block") +
    geom_line() + theme_classic()



ggsave("../results/August-2025-Rohan-GFP55-plot.pdf",Rohan.GFP55.timecourse.plot, height=2.5, width=7)
ggsave("../results/August-2025-Rohan-GFP100-plot.pdf",Rohan.GFP100.timecourse.plot, height=2.5, width=7)


ggsave("../results/August-2025-Grayson-GFP55-plot.pdf",Grayson.GFP55.timecourse.plot, height=2.5,  width=7)
ggsave("../results/August-2025-Grayson-GFP100-plot.pdf",Grayson.GFP100.timecourse.plot, height=2.5, width=7)

