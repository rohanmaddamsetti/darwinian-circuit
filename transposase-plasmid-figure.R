## transposase-plasmid-figure.R by Rohan Maddamsetti

## TODO: for the first panel, make two bars: one for the density of transposases
## per million bp on the chromosome, and the density of transposases per million bp
## on the plasmid.

## TODO: for the second panel / graph, do a binned average analysis.

library(tidyverse)
library(cowplot)

transposase.plasmid.summary <- read.csv("../results/transposase-counts.csv")

## 250,337 transposases on plasmids
transposase.on.plasmids.count <- sum(filter(
    transposase.plasmid.summary,
    SeqType == "plasmid")$transposase_count)

## 779,481 transposases on chromosomes
transposase.on.chromosomes.count <- sum(filter(
    transposase.plasmid.summary,
    SeqType == "chromosome")$transposase_count)

## 4,429,801 proteins on plasmids
CDS.on.plasmids.count <- sum(filter(
    transposase.plasmid.summary,
    SeqType == "plasmid")$CDS_count)

## 68,540,522 proteins on chromosome.
CDS.on.chromosomes.count <- sum(filter(
    transposase.plasmid.summary,
    SeqType == "chromosome")$CDS_count)

Fig1A.df <- data.frame(
    value=c(transposase.on.plasmids.count, transposase.on.chromosomes.count),
    SeqType=c("plasmid", "chromosome"))

Fig1B.df <- data.frame(
    value=c(CDS.on.plasmids.count, CDS.on.chromosomes.count),
    SeqType=c("plasmid", "chromosome"))

Fig1A <- ggplot(Fig1A.df, aes(x = "", y = value, fill = SeqType)) +
  geom_col(color = "white") +
  geom_text(aes(label = SeqType), size=5,
            position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle("transposases on chromosome and plasmids\n(1,049,818 transposases)")

Fig1B <- ggplot(Fig1B.df, aes(x = "", y = value, fill = SeqType)) +
  geom_col(color = "white") +
  geom_text(aes(label = SeqType), size=5,
            position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle("all genes on chromosomes and plasmids\n(72,970,323 genes)")

Fig1 <- plot_grid(Fig1A, Fig1B, ncol=1)
ggsave("../results/transposase-plasmid-piechart-Fig1.pdf", Fig1)

############################################################################
## revised figure design.
## for the first panel, make two bars: one for the density of transposases
## per million bp on the chromosome, and the density of transposases per million bp
## on the plasmid.
revised.Fig1A.df <- transposase.plasmid.summary %>%
    group_by(SeqType) %>%
    summarize(
        total_transposases = sum(transposase_count),
        total_length = sum(SeqLength)) %>%
    mutate(transposase_density = total_transposases / total_length) %>%
    mutate(transposase_density_per_Mbp = transposase_density * 1000000) %>%
    filter(SeqType != "unknown") ## just consider plasmids and chromosomes
    
revised.Fig1A <- revised.Fig1A.df %>%
    ggplot(aes(x=SeqType, y = transposase_density_per_Mbp)) +
    geom_bar(stat="identity") +
    theme_cowplot() +
    xlab("") +
    ylab("transposase density per Mbp") +
    coord_flip()

ggsave("../results/transposase-plasmid-barchart-Fig1A.pdf", revised.Fig1A)

############################################################################

PIRA.data <- read.csv("../data/Maddamsetti2025-S2Data-PIRA-PCN-estimates-with-normalization.csv")

PIRA.plasmid.transposase.data <- PIRA.data %>%
    inner_join(transposase.plasmid.summary) %>%
    filter(SeqType == "plasmid") %>%
    ## use 1 for TRUE and 0 for FALSE??
    mutate(contains_transposase = ifelse(transposase_count >=1, TRUE, FALSE)) %>%
    arrange(PIRACopyNumber) %>%
    select(AnnotationAccession, SeqID, SeqType, PIRACopyNumber, transposase_count, contains_transposase) %>%
    mutate(PCN_floor = floor(PIRACopyNumber))

write.csv(PIRA.plasmid.transposase.data,file="../results/PIRA-plasmid-transposase-summary.csv")

head(PIRA.plasmid.transposase.data)

Fig1C.data <- PIRA.plasmid.transposase.data %>%
    group_by(PCN_floor) %>%
    summarize(transposase_frequency= sum(contains_transposase)/n())

Fig1C <- Fig1C.data %>%
    ggplot(aes(x=PCN_floor, y = transposase_frequency)) +
    geom_point() +
    theme_classic()

Fig1C

log10.Fig1C <- Fig1C.data %>%
    ggplot(aes(x=log10(PCN_floor), y = transposase_frequency)) +
    geom_point() +
    #geom_bar(stat = "identity") +
    theme_classic()

log10.Fig1C

############################################################################
## revised figure design.

## make a figure showing the fraction of plasmids per bin
## TODO: set bins by powers of 2 in terms of copy number.
## repeat the same figures, using length.
## aim for 10 bins to cover the whole range.
## I can re-use my binomial confidence interval code to put ranges on the intervals.
