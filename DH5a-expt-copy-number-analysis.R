##DH5a-expt-copy-number-analysis.R by Rohan Maddamsetti.

library(tidyverse)
## use xml2 to get negative binomial fit from
## breseq output summary.html. 
library(xml2)
library(assertthat)


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


recode.DH5a.expt.treatment.factor.names <- function(df) {
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
    ## update the names of the Transposon, Plasmid, and Tet factors
    ## for a prettier plot.
    recode.DH5a.expt.treatment.factor.names()

## For Figure 2C, just show pUC in the Tet50 treatment with the working Tn5 transposase.
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
        legend.position="top",
        strip.background = element_blank())
ggsave("../results/Fig2C.pdf", Fig2C, width=4, height=2.5)

## let's write out the table.
 write.csv(evolved.replicon.coverage.ratio.df, "../results/DH5a-expt-genome-analysis/DH5a-expt-plasmid-transposon-coverage-ratios.csv", quote=F, row.names=FALSE)

########################################################
## Supplementary Figure S1. Plot transposon copy number against plasmid copy number.

S1Fig.df <- evolved.replicon.coverage.ratio.df %>%
    pivot_wider(names_from = ratio_type, values_from = ratio) %>%
    ## update the names of the Transposon, Plasmid, and Tet factors
    ## for a prettier plot.
    recode.DH5a.expt.treatment.factor.names() %>%
    select(-Tet, -Plasmid, -Transposon) %>%
    rename(Tet = "Tet_factor") %>%
    rename(Plasmid = "Plasmid_factor") %>%
    rename(Transposon = "Transposon_factor")

S1Fig <- S1Fig.df %>%
    ggplot(
        aes(x=plasmids.per.chromosome, y=transposons.per.chromosome, color=Transposon, shape=Plasmid)) +
    geom_point() +
    theme_classic() +
    facet_wrap(~Tet) +
    xlab("plasmid copy number") +
    ylab("transposon copy number") +
    geom_abline(slope=1,
                intercept=0,
                color="gray",
                linetype="dashed",
                size=0.5)
ggsave("../results/S1Fig.pdf", S1Fig, height=3)

