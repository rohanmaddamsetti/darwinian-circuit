##DH5a-expt-copy-number-analysis.R by Rohan Maddamsetti.

## 1) use xml2 to get negative binomial fit from
## breseq output summary.html. This is H0 null distribution of 1x coverage.

## 2) Find intervals longer than max.read.len that reject H0 coverage in genome.
##    at an uncorrected alpha = 0.05. This is to have generous predicted boundaries for amplifications.

## 3) Do a more rigorous test for each region. take positions in the region separated by more than max.read.len,
## and determine the probability that all are independently significant under the null, compared to
## a corrected bonferroni. The max.read.len ensures positions cannot be spanned by a single Illumina read.

## 4) Estimate copy number by dividing mean coverage in each region by the mean
##   of the H0 1x coverage distribution.

## 5) return copy number and boundaries for each significant amplification.

library(tidyverse)
library(xml2)
library(cowplot)
library(data.table)
library(dtplyr)
library(assertthat)
library(viridis)

## Bioconductor dependencies
##library(IRanges)
##library(GenomicRanges)
##library(rtracklayer)


#' parse the summary.html breseq output file, and return the mean and relative variance
#' of the negative binomial fit to the read coverage distribution, returned as a
#' data.frame with columns {mean, relative.variance}.
#' NOTE: this code has only been tested on the summary file
#' output by breseq 0.35.0. It will fail on breseq 0.37 and later, which uses the term "relative variance".

coverage.nbinom.from.html <- function(breseq.output.dir, sample.has.plasmid=TRUE) {
    summary.html.f <- file.path(breseq.output.dir, "output", "summary.html")
    tree <- read_html(summary.html.f)
    ## print text in the table 'Reference Sequence Information.
    query <- '//table[./tr/th[contains(text(),"fit dispersion")]]'
    table <- xml_find_first(tree,query)
    table.data <- xml_find_all(table,'./tr/td')
    chromosome.avg <- as.numeric(xml_text(table.data[5]))
    chromosome.relative.variance <- as.numeric(xml_text(table.data[6]))
    ## all samples should have these data.
    coverage.df <- data.frame('Sample' = basename(breseq.output.dir),
                              'mean'=c(chromosome.avg),
                              'relative.variance'=c(chromosome.relative.variance),
                              'variance'=c(chromosome.avg * chromosome.relative.variance),
                              'replicon'=c("chromosome"))
    if (sample.has.plasmid) {
            plasmid.avg <- as.numeric(xml_text(table.data[21]))
            plasmid.relative.variance <- as.numeric(xml_text(table.data[22]))
            plasmid.coverage.df <- data.frame('Sample' = basename(breseq.output.dir),
                                              'mean' = plasmid.avg,
                                              'relative.variance' = plasmid.relative.variance,
                                              'variance' = plasmid.avg * plasmid.relative.variance,
                                              'replicon' = "plasmid")
            ## now join the plasmid coverage data.
            coverage.df <- rbind(coverage.df, plasmid.coverage.df)
    }
    return(coverage.df)
}

#' get the maximum length of a sequencing read from the summary.html breseq
#' output file.
max.readlen.from.html <- function(breseq.output.dir) {
    summary.html.f <- file.path(breseq.output.dir, "output", "summary.html")
    tree <- read_html(summary.html.f)
    ## print text in the table 'Read File Information.
    query <- '//table[./tr/th[contains(text(),"longest")]]'
    table <- xml_find_first(tree,query)
    table.data <- xml_find_all(table,'./tr/td')
    readlen.index <- length(table.data) - 1
    max.readlen <- xml_integer(xml_find_all(table.data[readlen.index],".//b//text()"))
    return(max.readlen)
}

recode.treatment.factor.names <- function(df) {
    ## update the names of the Transposon, Plasmid, and Tet factors
    ## for a prettier plot.
    df %>%
        ## update the names of the Transposon factor for a prettier plot.
        mutate(Transposon_factor = fct_recode(as.factor(Transposon),
                                              `TetA++ (strong promoter)` = "B30",
                                              `TetA+ (weak promoter)` = "B20",
                                              `Tn5-` = "B59")) %>%
        ## update the names of the Plasmid factor for a prettier plot.
        mutate(Plasmid_factor = fct_recode(as.factor(Plasmid),
                                           `No plasmid` = "None",
                                           `p15A plasmid` = "p15A",
                                           `pUC plasmid` = "pUC")) %>%
        ## update the names of the Tet factor for a prettier plot.
        mutate(Tet_factor = fct_recode(as.factor(Tet),
                                       `Tet 0` = "0",
                                       `Tet 50` = "50"))
}

#######################################
## Analysis time!

## assert that we are in the src directory, such that
## proj.dir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("darwinian-circuit","src")))
projdir <- file.path("..")

## make a dataframe of samples and paths to the breseq output for each sample.
mixedpop.results.dir <- file.path(projdir, "results", "DH5a-expt-genome-analysis", "mixed-pops")
all.mixedpops <- list.files(mixedpop.results.dir,pattern='^RM')
all.mixedpop.paths <- sapply(all.mixedpops, function(x) file.path(mixedpop.results.dir,x))
mixedpop.input.df <- data.frame(Sample=all.mixedpops, path=all.mixedpop.paths)

DH5a.expt.metadata <- read.csv(
    file.path(
        projdir,
        "data",
        "DH5a-genome-sequencing",
        "DH5a-evolved-populations-and-clones.csv"))

######################################################################
## Plot the plasmid/chromosome and transposon/chromosome ratio in each sample.

## Get the actual coverage for the transposons in each mixed population.
## This is calculated by get-DH5a-expt-transposon-coverage.py.
transposon.coverage.file <- file.path(projdir, "results", "DH5a-expt-genome-analysis", "transposon-coverage.csv")
transposon.coverage.df <- read.csv(transposon.coverage.file) %>%
    ## let's add metadata and rename columns for compatibility.
    dplyr::rename(mean = TransposonCoverage) %>%
    dplyr::mutate(replicon = "transposon") 


evolved.replicon.coverage.df <- map_dfr(.x = mixedpop.input.df$path, .f = coverage.nbinom.from.html) %>%
    ## I am not examining dispersion or variance at this point.
    select(Sample, mean, replicon) %>%
    ## add transposon coverage data
    bind_rows(transposon.coverage.df) %>%
    ## set NA coverage values to zero.
    mutate(mean=sapply(mean, function(x) ifelse(is.na(x), 0,x))) %>%
    ## and add metadata.
    inner_join(DH5a.expt.metadata)


evolved.replicon.coverage.ratio.df <- evolved.replicon.coverage.df %>%
    pivot_wider(names_from = replicon, values_from = mean, names_prefix = "mean_") %>%
    group_by(Sample, Plasmid, Population) %>%
    summarise(transposons.per.chromosome = (mean_transposon/mean_chromosome),
              plasmids.per.chromosome = (mean_plasmid/mean_chromosome),
              transposons.per.plasmid = (mean_transposon/mean_plasmid)) %>%
    pivot_longer(cols = c(transposons.per.chromosome,plasmids.per.chromosome,transposons.per.plasmid),
                 names_to = "ratio_type", values_to = "ratio") %>%
    ungroup() %>%
    ## add metadata again.
    inner_join(DH5a.expt.metadata)
    
copy.number.plot.df <- evolved.replicon.coverage.ratio.df %>%
    ## we don't need transposons per plasmid, since we can get
    ## that from the other two ratios.
    filter(ratio_type != "transposons.per.plasmid") %>%
    mutate(ratio_type = fct_recode(as.factor(ratio_type),
                                       `Transposon copy number` = "transposons.per.chromosome",
                                   `Plasmid copy number` = "plasmids.per.chromosome")) %>%
    mutate(`Copy number` = ratio) %>%
    mutate(Population = as.factor(Population)) %>%
    ## log-transform copy number.
    mutate(`log(copy number)` = log2(ratio)) %>%
    ## update the names of the Transposon, Plasmid, and Tet factors
    ## for a prettier plot.
    recode.treatment.factor.names()

Fig2C.df <- copy.number.plot.df %>%
    filter(Transposon != "B59") %>%
    filter(Tet == 50) %>%
    filter(Plasmid == "pUC")

Fig2C <- Fig2C.df %>%
    ggplot(
        aes(x=Population, y=`Copy number`, fill=ratio_type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_classic() +
    facet_wrap(Transposon_factor~Plasmid_factor, scales="free") +
    theme(
        legend.title=element_blank(),
        legend.position="bottom",
        strip.background = element_blank())
ggsave("../results/Fig2C.pdf", Fig2C, width=4, height=2.5)

## let's write out the table.
 write.csv(evolved.replicon.coverage.ratio.df, "../results/DH5a-expt-genome-analysis/DH5a-expt-plasmid-transposon-coverage-ratios.csv", quote=F, row.names=FALSE)

########################################################
## 1) plot y-axis in log-scale.

Fig1I.variation1 <- ggplot(data=Fig1I.df, aes(x=Population, y=`log(copy number)`, fill=ratio_type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_classic() +
    facet_wrap(.~Plasmid, scales="free") +
    theme(legend.title=element_blank(), legend.position="bottom")
ggsave("../results/Fig1I-log-scale.pdf", Fig1I.variation1, width=8, height=3.5)

## 2) Rank order isolates in each group by transposon copy number.
transposon.copy.ranks <- Fig1I.df %>%
    filter(ratio_type == "Transposon copy number") %>%
    arrange(Plasmid, `Copy number`) %>%
    group_by(Plasmid) %>% 
    mutate(Rank = rank(`Copy number`, ties.method = "first")) %>%
    mutate(Rank = as.factor(Rank)) %>%
    ## have to drop columns to get the merge to work right.
    select(Sample, Plasmid, Rank)

rank.ordered.Fig1I.df <- Fig1I.df %>%
    full_join(transposon.copy.ranks)

rank.ordered.Fig1I <- ggplot(data=rank.ordered.Fig1I.df, aes(x=Rank, y=`Copy number`, fill=ratio_type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_classic() +
    facet_wrap(.~Plasmid, scales="free") +
    xlab("Ranked populations") +
    theme(legend.title=element_blank(), legend.position="bottom")
ggsave("../results/Fig1I-rank-ordered.pdf", rank.ordered.Fig1I, width=8, height=3.5)

## 3) Plot plasmid-copy-number against transposon-copy-number.
## let's calculate lines of best fit for p15A and pUC separately.
p15A.transposon.plasmid.copy.df <- evolved.replicon.coverage.ratio.df %>%
    filter(Plasmid == "p15A")

p15A.transposon.plasmid.correlation <- lm(
    transposons.per.chromosome~plasmids.per.chromosome,
    data=p15A.transposon.plasmid.copy.df)
summary(p15A.transposon.plasmid.correlation)

p15A.transposon.plasmid.correlation.plot <- ggplot(
    data=p15A.transposon.plasmid.copy.df,
    aes(x=plasmids.per.chromosome,y=transposons.per.chromosome,color=Plasmid)) +
    geom_point(color="blue") +
    ylim(0,180) +
    geom_abline(slope=1,
                intercept=0,
                color="gray",linetype="dashed",size=0.1) +
    geom_abline(slope=p15A.transposon.plasmid.correlation$coefficients[2],
                intercept=p15A.transposon.plasmid.correlation$coefficients[1],
                color="red",linetype="dashed",size=0.1) +
    facet_wrap(.~Plasmid) +
    theme_classic() +
    guides(color="none")


pUC.transposon.plasmid.copy.df <- evolved.replicon.coverage.ratio.df %>%
    filter(Plasmid == "pUC")

pUC.transposon.plasmid.correlation <- lm(
    transposons.per.chromosome~plasmids.per.chromosome,
    data=pUC.transposon.plasmid.copy.df)
summary(pUC.transposon.plasmid.correlation)

pUC.transposon.plasmid.correlation.plot <- ggplot(
    data=pUC.transposon.plasmid.copy.df,
    aes(x=plasmids.per.chromosome,y=transposons.per.chromosome,color=Plasmid)) +
    geom_point(color="orange") +
    ylim(0,2000) +
    geom_abline(slope=1,
                intercept=0,
                color="gray",linetype="dashed",size=0.1) +
    geom_abline(slope=pUC.transposon.plasmid.correlation$coefficients[2],
                intercept=pUC.transposon.plasmid.correlation$coefficients[1],
                color="red",linetype="dashed",size=0.1) +
    facet_wrap(.~Plasmid) +
    theme_classic() +
    guides(color="none")

transposon.plasmid.correlation.plot <- plot_grid(
    p15A.transposon.plasmid.correlation.plot,
    pUC.transposon.plasmid.correlation.plot,nrow=1)
ggsave("../results/evolved-transposon-plasmid-correlation.pdf", transposon.plasmid.correlation.plot, width=6, height=3.5)


########################################################
## let's examine copy number in the ancestral clones.
## these were directly streaked from glycerol stock onto LB plates, and then
## the plates were sent for sequencing. So no selection was applied to maintain
## the p15A or pUC plasmids in the ancestral clones on the plate.

ancestral.transposon.coverage.df <- filter(transposon.coverage.df, Sample %in% ancestral.clone.input.df$Sample)

ancestral.replicon.coverage.df <- map_dfr(.x = ancestral.clone.input.df$path, .f = coverage.nbinom.from.html) %>%
    ## I am not examining dispersion or variance at this point.
    select(Sample, mean, replicon) %>%
    ## add transposon coverage data
    bind_rows(ancestral.transposon.coverage.df) %>%
    ## set NA coverage values to zero.
    mutate(mean=sapply(mean, function(x) ifelse(is.na(x), 0,x))) %>%
    ## and add metadata.
    full_join(ancestral.clone.metadata)


ancestral.replicon.coverage.ratio.df <- ancestral.replicon.coverage.df %>%
    pivot_wider(names_from = replicon, values_from = mean, names_prefix = "mean_") %>%
    group_by(Sample, Plasmid, Population) %>%
    summarise(transposons.per.chromosome = (mean_transposon/mean_chromosome),
              plasmids.per.chromosome = (mean_plasmid/mean_chromosome),
              transposons.per.plasmid = (mean_transposon/mean_plasmid)) %>%
    pivot_longer(cols = c(transposons.per.chromosome,plasmids.per.chromosome,transposons.per.plasmid),
                 names_to = "ratio_type", values_to = "ratio")


ancestralFig1I.df <- ancestral.replicon.coverage.ratio.df %>%
    ## we don't need transposons per plasmid, since we can get
    ## that from the other two ratios.
    filter(ratio_type != "transposons.per.plasmid") %>%
    mutate(ratio_type = fct_recode(as.factor(ratio_type),
                                       `Transposon copy number` = "transposons.per.chromosome",
                                   `Plasmid copy number` = "plasmids.per.chromosome")) %>%
    mutate(`Copy number` = ratio) %>%
    mutate(Population = as.factor(Population)) %>%
    ## log-transform copy number.
    mutate(`log(copy number)` = log2(ratio))

ancestralFig1I <- ggplot(data=ancestralFig1I.df, aes(x=Population, y=`Copy number`, fill=ratio_type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_classic() +
    facet_wrap(.~Plasmid, scales="free") +
    theme(legend.title=element_blank(), legend.position="bottom")
ggsave("../results/ancestralFig1I.pdf", ancestralFig1I, width=8, height=3.5)

########################################################
## Make the full transposon-plasmid correlation plot that LC asked for,
## including DH5a mixed population data.

## add columns to merge datasets.
K12.ancestral.clone.plasmid.transposon.ratio.df <- ancestral.replicon.coverage.ratio.df %>%
    mutate(Strain = "K12") %>%
    mutate(SampleType = "Clone") %>%
    mutate(Tet=0) %>%
    mutate(Transposon="B30")

K12.evolved.clone.plasmid.transposon.ratio.df <- evolved.replicon.coverage.ratio.df %>%
    mutate(Strain = "K12") %>%
    mutate(SampleType = "Clone") %>%
    mutate(Tet=50) %>%
    mutate(Transposon="B30")

DH5a.evolved.mixed.pop.replicon.coverage.ratio.df <- read.csv(
    "../../transposon-plasmid-evolution/results/draft-manuscript-1A/plasmid-transposon-coverage-ratios.csv") %>%
    mutate(SampleType = "MixedPop") %>%
    mutate(Strain = "DH5a")

full.replicon.coverage.ratio.df <-
    K12.ancestral.clone.plasmid.transposon.ratio.df %>%
    bind_rows(K12.evolved.clone.plasmid.transposon.ratio.df) %>%
    bind_rows(DH5a.evolved.mixed.pop.replicon.coverage.ratio.df) %>%
    mutate(Population = as.factor(Population)) %>%
    mutate(Tet = as.factor(Tet)) %>%
    ## this is to fix discrepancies between labeling of DH5a and K-12 samples.
    mutate(Plasmid = replace(Plasmid, Plasmid == "None", "No plasmid"))


newFig1.df <- full.replicon.coverage.ratio.df %>%
    ## we don't need transposons per plasmid, since we can get
    ## that from the other two ratios.
    filter(ratio_type != "transposons.per.plasmid") %>%
    mutate(ratio_type = fct_recode(as.factor(ratio_type),
                                   `Transposon copy number` = "transposons.per.chromosome",
                                   `Plasmid copy number` = "plasmids.per.chromosome")) %>%
    pivot_wider(names_from = ratio_type, values_from = ratio) %>%
    ## set no plasmid treatment copy number to 0.
    mutate(`Plasmid copy number` = replace_na(`Plasmid copy number`, 0)) %>%
    ## unite the Strain, SampleType, Tet, Transposon columns together
        unite("Treatment", Strain:Transposon, sep="-", remove = FALSE)


newFig1 <- ggplot(data=newFig1.df,
                  aes(x=`Plasmid copy number`, y=`Transposon copy number`, color=Treatment, shape=Plasmid)) +
    geom_point() +
    theme_classic() +
#    facet_wrap(SampleType~Tet) +
    geom_abline(slope=1,
                intercept=0,
                color="gray",linetype="dashed",size=0.5)
ggsave("../results/newFig1.pdf", newFig1, height=3.5)

log.newFig1 <- newFig1 +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10')
ggsave("../results/log-newFig1.pdf", log.newFig1, height=3.5)

              
