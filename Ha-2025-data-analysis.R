## Ha-2025-data-analysis.R by Rohan Maddamsetti.

## This script re-analyzes relevant data from Ha et al. (2025)
## "Transposon-plasmid nesting enables fast response to fluctuating environments"
## https://www.biorxiv.org/content/10.1101/2025.06.04.657954v1

## This analysis only considers the data from experiments with strains RM8.36.1-4,
## which contain active Tn5 integrated into the chromosome.

## Experimental details from Ha et al. (2025)
## We exposed the different strains to periodic removal and addition of tetracycline,
## in three iterations for a total duration of 21 days, with five replicates per strain.
## Each iteration consisted of serial passaging of 30,000-fold dilution in fresh media,
## which corresponds to log2(30000) = ~15 generations, without tetracycline for five days
## followed by two days with tetracycline (5 μg/mL), as the response to the tetracycline addition
## occurs more rapidly than that of following the removal of tetracycline.


library(tidyverse)
library(cowplot)


## Days when Tet 5 selection was applied.
TET_SELECTION_DAYS <- c(0,6,7,13,14,20,21)

## Plasmid color scheme:
## https://colorbrewer2.org/#type=sequential&scheme=PuRd&n=5
PLASMID_COLORSCALE <- c(
    "pUC" = "#7a0177",
    "CloDF13" = "#c51b8a",
    "pBR322" = "#f768a1"
)


## Get OD600 and GFP data.
Ha2025.OD.GFP.data <- read.csv("../data/OD600-GFP-data/Ha2025/withTransposase_GFP-OD.csv") %>%
    mutate(Plasmid = if_else(Plasmid == "CloDF", "CloDF13", Plasmid)) %>%
    mutate(Plasmid = factor(Plasmid, levels=c("pUC", "CloDF13", "pBR322"))) %>%
    rename(OD600 = OD) %>%
    mutate(Tet_selection = if_else(Day %in% TET_SELECTION_DAYS, TRUE, FALSE)) %>%
    ## select the "raw data" to recalculate summary statistics here.
    select(Well, Replicate, Kan, Plasmid, Day, OD600, GFP) %>%
    ## don't need the chromosome control here.
    filter(Plasmid != "Chromosome") %>%
    ## calculate normalized GFP directly from GFP and OD600 measurements.
    mutate(normalized_GFP = GFP/OD600)


## Get qPCR data.
Ha2025.qPCR.data <- read.csv("../data/qPCR-data/Ha2025/withTransposase_qPCR.csv") %>%
    rename(Plasmid = PCN) %>%
    mutate(Plasmid = ifelse(Plasmid == "CloDF", "CloDF13", Plasmid)) %>%
    mutate(Plasmid = factor(Plasmid, levels=c("pUC", "CloDF13", "pBR322"))) %>%
    rename(Replicate = ExpRep) %>%
    mutate(Tet_selection = if_else(Day %in% c(0,6,7,13,14,20,21), TRUE, FALSE)) %>%
    ## select the "raw data" to recalculate summary statistics here.
    select(Replicate, Kan, Plasmid, Day, transposon.copy, plasmid.copy) %>%
    ## don't need the chromosome control here.
    filter(Plasmid != "Chromosome") %>%
    ## calculate transposons per plasmid.
    mutate(transposons.per.plasmid = transposon.copy/plasmid.copy)


## split the data into the Kan+ and Kan- treatments.
Kan.OD.GFP.data <- filter(Ha2025.OD.GFP.data, Kan == TRUE)
noKan.OD.GFP.data <- filter(Ha2025.OD.GFP.data, Kan == FALSE)

Kan.qPCR.data <- filter(Ha2025.qPCR.data, Kan == TRUE)
noKan.qPCR.data <- filter(Ha2025.qPCR.data, Kan == FALSE)


## Plot the Kan+ data.
## Figure 5C: log10(GFP/OD600)
Fig5C <- Kan.OD.GFP.data %>%
    ggplot(aes(x = Day, y = log10(normalized_GFP), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    ylab("log10(GFP/OD600)") +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    geom_line() +
    theme_classic() +
    geom_vline(xintercept = TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme(legend.position = "bottom")

## get the figure legend
Fig5CDE_legend <- get_legend(Fig5C)
## and remove from the panel.
Fig5C <- Fig5C + theme(legend.position = "none")


## Figure 5D: transposons by qPCR.
Fig5D <- Kan.qPCR.data %>%
    ggplot(aes(x = Day, y = log10(transposon.copy), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    ylab(expression("log10(" * italic("tetA") *" copies)")) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    geom_line() +
    theme_classic() +
    geom_vline(xintercept = TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme(legend.position = "none")

## Figure 5E: transposons per plasmid by qPCR.
Fig5E <- Kan.qPCR.data %>%
    ggplot(aes(x = Day, y = transposons.per.plasmid, color=Plasmid, group = interaction(Replicate, Plasmid))) +
    ylab(expression(italic("tetA") * " per plasmid")) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    geom_line() +
    theme_classic() +
    geom_vline(xintercept = TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme(legend.position = "none")

## Make the title.
## see cowplot documentation here: https://wilkelab.org/cowplot/articles/plot_grid.html
Fig5CDE_title <- ggdraw() + 
    draw_label(
        "Ha et al. (2025)\nkanamycin+ treatment",
        fontface = 'bold',
        x = 0,
        hjust = 0
    ) +
    theme(
        ## add margin on the left of the drawing canvas,
        ## so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
    )


Fig5CDE <- plot_grid(Fig5CDE_title, Fig5CDE_legend, Fig5C, Fig5D, Fig5E, labels = c("", "","C", "D", "E"), ncol = 1, rel_heights=c(0.25, 0.15, 1,1,1))
## save the figure.
ggsave("../results/Fig5CDE.pdf",Fig5CDE, height=7.75, width=3)


## Now plot the Kan- data for supplement.
## Supplementary Figure S6 A: log10(GFP/OD600)
S6FigA <- noKan.OD.GFP.data %>%
    ggplot(aes(x = Day, y = log10(normalized_GFP), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    ylab("log10(GFP/OD600)") +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    geom_line() +
    theme_classic() +
    geom_vline(xintercept = TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme(legend.position = "bottom")


## get the figure legend
S6Fig_legend <- get_legend(S6FigA)
## and remove from the panel.
S6FigA <- S6FigA + theme(legend.position = "none")


## Supplementary Figure S6 B: transposons by qPCR.
S6FigB <- noKan.qPCR.data %>%
    ggplot(aes(x = Day, y = log10(transposon.copy), color=Plasmid, group = interaction(Replicate, Plasmid))) +
    ylab(expression("log10(" * italic("tetA") *" copies)")) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    geom_line() +
    theme_classic() +
    geom_vline(xintercept = TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme(legend.position = "none")
 

## Supplementary Figure S6 C: transposons per plasmid by qPCR.
S6FigC <- noKan.qPCR.data %>%
    ggplot(aes(x = Day, y = transposons.per.plasmid, color=Plasmid, group = interaction(Replicate, Plasmid))) +
    ylab(expression(italic("tetA") * " per plasmid")) +
    scale_color_manual(values = PLASMID_COLORSCALE) +
    geom_line() +
    theme_classic() +
    geom_vline(xintercept = TET_SELECTION_DAYS, linetype = "dashed", color = "gray") +
    theme(legend.position = "none")


## Make the title.
## see cowplot documentation here: https://wilkelab.org/cowplot/articles/plot_grid.html
S6Fig_title <- ggdraw() + 
    draw_label(
        "Ha et al. (2025)\nkanamycin- treatment",
        fontface = 'bold',
        x = 0,
        hjust = 0
    ) +
    theme(
        ## add margin on the left of the drawing canvas,
        ## so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
    )


S6Fig <- plot_grid(S6Fig_title, S6Fig_legend, S6FigA, S6FigB, S6FigC, labels = c("","","A", "B", "C"), ncol = 1, rel_heights=c(0.25,0.15, 1,1,1))
## save the figure.
ggsave("../results/S6Fig.pdf", S6Fig, height=7.75, width=3)

 

