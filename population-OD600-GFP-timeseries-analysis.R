## population-OD600-GFP-timeseries-analysis.R
## by Rohan Maddamsetti.

library(tidyverse)

## let's analyze initial data collected with Yuanchi, finished 12/14/2022.
## Note: measurements "OVER" the GFP100 measurement range are encoded as 80000
## as a guessed maximum-ish value.

## at first, I'm not subtracting the blanks.
## TODO: SUBTRACT BLANK READINGS!!!
initial.darwin.data <- read.csv("../data/initial-test-experiment.csv") %>%
    mutate(normalized.GFP55 = GFP55/OD600) %>%
    mutate(normalized.GFP100 = GFP100/OD600) %>%
    mutate(log.normalized.GFP55 = log10(GFP55/OD600)) %>%
    mutate(log.normalized.GFP100 = log10(GFP100/OD600)) %>%
    ## let's remove the blank measurements for now.
    filter(Sample != "Blank") %>%
    mutate(Tet = as.factor(Tet))
    

initial.GFP55.timecourse.plot <- ggplot(data=initial.darwin.data,
                                        aes(x=Day, y=log10(normalized.GFP55), color=Tet)) +
    geom_point() + theme_classic()

initial.GFP100.timecourse.plot <- ggplot(data=initial.darwin.data,
                                        aes(x=Day, y=log10(normalized.GFP100), color=Tet)) +
    geom_point() + theme_classic()

ggsave("../results/initial-GFP55-plot.pdf",initial.GFP55.timecourse.plot, height=4, width=4)
ggsave("../results/initial-GFP100-plot.pdf",initial.GFP100.timecourse.plot, height=4, width=4)

