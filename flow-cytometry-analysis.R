## flow-cytometry-analysis.R by Rohan Maddamsetti.
## This script analyzes the raw flow cytometry data in .fcs format
## that represents GFP per cell in the Darwin circuit clones,
## and a K-12 MG1655 negative control with no GFP,
## as measured with assistance from Bin Li on 12/12/2023.

## report the GFP distributions for the 5 clones + negative control in this experiment.

##This script uses packages from Bioconductor that are amenable for flow
##cytometry analysis.

## Load packages.
library(tidyverse)
library(flowCore)

##########################
### GFP EXPRESSION ANALYSIS
##########################

## We collected 100000 events before gating.

K12.MG1655.negative.control.clone.data <- read.FCS("../data/tetA-GFP-flow-cytometry-data/Rohan-Maddamsetti-4512012-GFP-staff-assist-2023-12-12/Specimen_001_Neg_Ctrl_k12.fcs")

pUC.clone.data <- read.FCS("../data/tetA-GFP-flow-cytometry-data/Rohan-Maddamsetti-4512012-GFP-staff-assist-2023-12-12/Specimen_001_Sample_A.fcs")

CloDF13.clone.data <- read.FCS("../data/tetA-GFP-flow-cytometry-data/Rohan-Maddamsetti-4512012-GFP-staff-assist-2023-12-12/Specimen_001_Sample_B.fcs")

ColE1.clone.data <- read.FCS("../data/tetA-GFP-flow-cytometry-data/Rohan-Maddamsetti-4512012-GFP-staff-assist-2023-12-12/Specimen_001_Sample_C.fcs")

p15A.clone.data <- read.FCS("../data/tetA-GFP-flow-cytometry-data/Rohan-Maddamsetti-4512012-GFP-staff-assist-2023-12-12/Specimen_001_Sample_D.fcs")

SC101.clone.data <- read.FCS("../data/tetA-GFP-flow-cytometry-data/Rohan-Maddamsetti-4512012-GFP-staff-assist-2023-12-12/Specimen_001_Sample_E.fcs")


## Set a GFP gate based on the K12 negative control and pUC.
gfp_min <- 400
gfp_max <- Inf
gfp_gate <- rectangleGate("GFP-A" = c(gfp_min, gfp_max))
gated_negative_control <- Subset(K12.MG1655.negative.control.clone.data, gfp_gate)


## Gate the data based on GFP.
gated.pUC.flowdata <- Subset(pUC.clone.data, gfp_gate)
gated.CloDF13.flowdata <- Subset(CloDF13.clone.data, gfp_gate)
gated.ColE1.flowdata <- Subset(ColE1.clone.data, gfp_gate)
gated.p15A.flowdata <- Subset(p15A.clone.data, gfp_gate)
gated.pSC101.flowdata <- Subset(SC101.clone.data, gfp_gate)


## convert to data frames.
gated.pUC.df <- as_tibble(exprs(gated.pUC.flowdata)) %>%
    mutate(Plasmid = "pUC")
gated.CloDF13.df <- as_tibble(exprs(gated.CloDF13.flowdata)) %>%
    mutate(Plasmid = "CloDF13")
gated.ColE1.df <- as_tibble(exprs(gated.ColE1.flowdata)) %>%
    mutate(Plasmid = "ColE1")
gated.p15A.df <- as_tibble(exprs(gated.p15A.flowdata)) %>%
    mutate(Plasmid = "p15A")
gated.pSC101.df <- as_tibble(exprs(gated.pSC101.flowdata)) %>%
    mutate(Plasmid = "pSC101")


## merge into a big dataframe.
Fig4D.df <- bind_rows(
    gated.pUC.df,
    gated.CloDF13.df,
    gated.ColE1.df,
    gated.p15A.df,
    gated.pSC101.df
) %>%
    ## order the plasmid factor
    mutate(Plasmid = factor(Plasmid, levels=c("pUC", "CloDF13","ColE1","p15A","pSC101")))


Fig4D <- ggplot(Fig4D.df, aes(x = log10(`GFP-A`), fill=Plasmid)) +
    geom_histogram(bins=500, alpha = 0.5) +
    theme_classic() +
    facet_grid(Plasmid~.) +
    guides(fill="none") +
    xlab("log10(GFP) (arbitrary units)") +
    ylab("cell counts")
    
ggsave("../results/flow-cytometry/Fig4D.pdf", Fig4D, width=3, height=4)
