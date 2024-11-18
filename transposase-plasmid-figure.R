## transposase-plasmid-figure.R by Rohan Maddamsetti

## TODO: change the Figure 1B axis to show logarithmically spaced ticks, but using normal PCN units
## (not log-transformed).

## TODO: make the symbol size in Figure 1B to show the relative plasmid length for plasmids in each bin.
## to do this, set some scale transformation for plasmid size, and add these relevant data to the data frame
## for figure 1B.

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

############################################################################
## Figure 1A.
## for the first panel, make two bars: one for the density of transposases
## per million bp on the chromosome, and the density of transposases per million bp
## on the plasmid.
Fig1A.df <- transposase.plasmid.summary %>%
    group_by(SeqType) %>%
    summarize(
        total_transposases = sum(transposase_count),
        total_length = sum(SeqLength)) %>%
    mutate(transposase_density = total_transposases / total_length) %>%
    mutate(transposase_density_per_Mbp = transposase_density * 1000000) %>%
    filter(SeqType != "unknown") ## just consider plasmids and chromosomes
    
Fig1A <- Fig1A.df %>%
    ggplot(aes(x=SeqType, y = transposase_density_per_Mbp)) +
    geom_bar(stat="identity") +
    theme_cowplot() +
    xlab("") +
    ylab("number of transposases per Mbp") +
    coord_flip()

############################################################################
## Figure 1B.
## make a figure showing the fraction of plasmids per bin
## set bins by powers of 2 in terms of copy number.
## repeat the same figures, using length.
## aim for 10 bins to cover the whole range.
## I can re-use my binomial confidence interval code to put ranges on the intervals.

## See Wikipedia reference:
## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

## Make Z-distributed confidence intervals for the fraction of isolates with
## duplicated ARGs (panel A),
## the fraction of isolates with single-copy ARGs (panel B),
## the fraction of isolates with duplicated genes (panel C).


calc.plasmid.confints <- function(df) {
    df %>%
        ## use the normal approximation for binomial proportion conf.ints
        mutate(se = sqrt(p*(1-p)/total_plasmids)) %>%
        ## See Wikipedia reference:
        ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
        mutate(Left = p - 1.96*se) %>%
        mutate(Right = p + 1.96*se) %>%
        ## truncate confidence limits to interval [0,1].
        rowwise() %>% mutate(Left = max(0, Left)) %>%
        rowwise() %>% mutate(Right = min(1, Right))
}


make.plasmid.confint.figure.panel <- function(Table, title, no.category.label = FALSE) {
    
    Fig.panel <- Table %>%
        ## back-transform to get into regulator copy number space.
        mutate(exp_log2_PCN_bin = 2^(log2_PCN_bin)) %>%
        mutate(label_text = paste(total_contains_transposase, total_plasmids, sep = "/")) %>%
        ggplot(aes(y = log2_PCN_bin, x = p)) +
        geom_point(size=1) +
        ylab("PCN bin") +
        xlab("proportion of plasmids") +
        theme_classic() +
        ggtitle(title) +
        ## plot CIs.
        geom_errorbarh(aes(xmin=Left,xmax=Right), height=0.2, size=0.2) +
        geom_text(aes(label=label_text), nudge_x = 0.1)
    
    if (no.category.label)
        Fig.panel <- Fig.panel +
            theme(axis.text.y=element_blank())
    
    return(Fig.panel)
}


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


Fig1B.df <- PIRA.plasmid.transposase.data %>%
    mutate(log2_PCN = log2(PIRACopyNumber)) %>%
    mutate(log2_PCN_bin = floor(log2_PCN)) %>%
    group_by(log2_PCN_bin) %>%
    summarize(
        total_contains_transposase = sum(contains_transposase),
        total_plasmids = n(), ## total number of plasmids in this bin.
        ## p is the frequency of plasmids containing a transposase in this bin.
        p = total_contains_transposase / total_plasmids) %>%
    calc.plasmid.confints() %>%
    ## remove bins with fewer than 50 plasmids.
    filter(total_plasmids >= 50)


Fig1B <- Fig1B.df %>%
    make.plasmid.confint.figure.panel("fraction of plasmids containing transposases")

Fig1 <- plot_grid(Fig1A, Fig1B, ncol=1, rel_heights=c(0.4,1), labels=c('A','B'))
ggsave("../results/transposase-plasmid-barplot-Fig1.pdf", Fig1, height = 4, width=9)

## TODO: take the bins in the panel B, and calculate the density of transposases
## for each bin, and calculate the mean length of plasmids in each bin, etc.
## do this as a sanity check for the testfigure below.

## If the following figure is correct, then this analysis indicates that smaller plasmids,
## per capita, have more transposases per Mbp. Larger plasmids have a higher
## chance of having some transposase, but normalizing for size, smaller plasmids
## may be enriched with transposases-- CAREFULLY DOUBLE-CHECK THIS!!

## TODO: plot the number of transposases per Mbp per bin defined by plasmid length
## CRITICAL TODO: CHECK IF THIS CALCULATION IS CORRECT!!


test.df <- transposase.plasmid.summary %>%
    filter(SeqType == "plasmid") %>%
    mutate(SeqLengthBin = floor(log10(SeqLength))) %>%
    group_by(SeqLengthBin) %>%
    summarize(
        total_transposases = sum(transposase_count),
        total_length = sum(SeqLength)) %>%
    mutate(transposase_density = total_transposases / total_length) %>%
    mutate(transposase_density_per_Mbp = transposase_density * 1000000)
    
testfig <- test.df %>%
    ggplot(aes(x=SeqLengthBin, y = transposase_density_per_Mbp)) +
    geom_bar(stat="identity") +
    theme_cowplot() +
    ylab("number of transposases per Mbp") +
    coord_flip()
