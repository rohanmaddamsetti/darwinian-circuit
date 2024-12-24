##one-day-expt-copy-number-analysis.R by Rohan Maddamsetti.

library(tidyverse)
library(xml2)
## use xml2 to get negative binomial fit from
## breseq output summary.html.
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


recode.one.day.expt.treatment.factor.names <- function(df) {
    ## update the names of the Transposon, Plasmid, and Tet factors
    ## for a prettier plot.
    df %>%
        ## organize these columns in the dataframe first.
        relocate(Sample, Population, Transposon, Tet, Plasmid) %>%
        ## update the names of the Transposon factor for a prettier plot.
        mutate(Transposon_factor = fct_recode(as.factor(Transposon),
                                              `TetA++ (strong promoter)` = "B30",
                                              `Tn5-` = "B59")) %>%
        ## update the names of the Plasmid factor for a prettier plot.
        mutate(Plasmid_factor = fct_recode(as.factor(Plasmid),
                                           `No plasmid` = "None",
                                           `p15A plasmid` = "p15A",
                                           `pUC plasmid` = "pUC")) %>%
        ## update the names of the Tet factor for a prettier plot.
        mutate(Tet_factor = fct_recode(as.factor(Tet),
                                       `Tet 0` = "0",
                                       `Tet 5` = "5")) ##%>%
    ## unite the Transposon_factor, Tet_factor columns together.
    ##unite("Treatment", Transposon_factor:Tet_factor, sep="\n", remove = FALSE)
}


#######################################
## Analysis time!

## assert that we are in the src directory, such that
## proj.dir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("darwinian-circuit","src")))
projdir <- file.path("..")

## get metadata for all the evolved population metagenomes.
metagenome.metadata <- read.csv("../data/one-day-expt-evolved-sample-metadata.csv")

mixedpop.output.dir <- file.path(projdir, "results", "one-day-expt-genome-analysis", "one-day-expt-evolved-pops")
all.mixedpops <- list.files(mixedpop.output.dir,pattern='^RM')
all.mixedpop.paths <- sapply(all.mixedpops, function(x) file.path(mixedpop.output.dir,x))
mixedpop.input.df <- data.frame(Sample=all.mixedpops, path=all.mixedpop.paths) %>%
    ## skip the two clone samples for now.
    inner_join(metagenome.metadata)

## get metadata for the ancestral clones
ancestralclone.metadata <- read.csv("../data/one-day-expt-ancestral-sample-metadata.csv")

## get corresponding inputs for the ancestral clones.
ancestralclone.output.dir <- file.path(projdir, "results", "one-day-expt-genome-analysis")
all.ancestralclones <- list.files(ancestralclone.output.dir, pattern='^remapped-RM')
## we have to cut the "remapped-" prefix from the sample names.
remapped.prefix.len <- 10
all.ancestralclone.samples <- substring(all.ancestralclones, remapped.prefix.len)

all.ancestralclone.paths <- sapply(all.ancestralclones, function(x) file.path(ancestralclone.output.dir,x))

ancestralclone.input.df <- data.frame(
    Sample=all.ancestralclone.samples,
    path=all.ancestralclone.paths) %>%
    inner_join(ancestralclone.metadata)

ancestral.clones.df <- ancestralclone.metadata %>%
    ## these clones are the ancestors, so add a Ancestor column to sample ID.
    ## IMPORTANT: ancestral.clones.df is needed downstream is some places,
    ## in other places ancestralclone.metadata is what is needed.
   dplyr::rename(Ancestor = Sample) %>%
    select(-SampleType) %>%
    mutate(gff_name = paste0(Ancestor, ".gff3")) %>%
    mutate(gff_path = file.path(projdir, "results", "one-day-expt-genome-analysis", gff_name))

######################################################################
## Plot the plasmid/chromosome and transposon/chromosome ratio in each sample.

## Get the actual coverage for the B30 and B59 transposons.
## This is calculated by get-one-day-expt-transposon-coverage.py.
transposon.coverage.file <- file.path(projdir, "results", "one-day-expt-genome-analysis", "transposon-coverage.csv")
transposon.coverage.df <- read.csv(transposon.coverage.file) %>%
    ## let's add metadata and rename columns for compatibility.
    dplyr::rename(mean = TransposonCoverage) %>%
    dplyr::mutate(replicon = "transposon")

## data for the ancestral samples.
ancestral.transposon.coverage.df <- ancestralclone.metadata %>%
    left_join(transposon.coverage.df)

ancestral.replicon.coverage.df <- map_dfr(.x = ancestralclone.input.df$path, .f = coverage.nbinom.from.html) %>%
    ## fix the names of the samples.
    mutate(Sample = substring(Sample, remapped.prefix.len)) %>%
    ## add metadata.
    full_join(ancestralclone.metadata) %>%
    ## select the columns that matter in this analysis.
    select(Sample, SampleType, Transposon, Plasmid, replicon, mean) %>%
    ## add rows for the transposon coverage data.
    full_join(ancestral.transposon.coverage.df) %>%
    ## set sensible default values for these columns
    ## for comparison to the evolved populations.
    mutate(Population = 1) %>%
    mutate(Tet = 0)

ancestral.replicon.coverage.ratio.df <- ancestral.replicon.coverage.df %>%
    pivot_wider(names_from = replicon, values_from = mean, names_prefix = "mean_") %>%
    group_by(Sample, Transposon, Plasmid, Population, Tet) %>%
    summarise(transposons.per.chromosome = (mean_transposon/mean_chromosome),
              plasmids.per.chromosome = (mean_plasmid/mean_chromosome),
              transposons.per.plasmid = (mean_transposon/mean_plasmid)) %>%
    pivot_longer(cols = c(transposons.per.chromosome,plasmids.per.chromosome,transposons.per.plasmid),
                 names_to = "ratio_type", values_to = "ratio") %>%
    ## add a Day column to compare to evolved samples.
    mutate(Day = 0)

## data for the evolved samples.
evolved.transposon.coverage.df <- metagenome.metadata %>%
    left_join(transposon.coverage.df) %>%
    ## we don't need this when joining to evolved.replicon.coverage.df.
    select(-SampleType)

evolved.replicon.coverage.df <- map_dfr(.x = mixedpop.input.df$path, .f = coverage.nbinom.from.html) %>%
    inner_join(metagenome.metadata) %>%
    ## I am not examining dispersion or variance at this point.
    select(Sample, mean, replicon, Transposon, Plasmid, Population, Tet) %>%
    ## add rows for the transposon coverage data.
    full_join(evolved.transposon.coverage.df)    

evolved.replicon.coverage.ratio.df <- evolved.replicon.coverage.df %>%
    pivot_wider(names_from = replicon, values_from = mean, names_prefix = "mean_") %>%
    group_by(Sample, Transposon, Plasmid, Population, Tet) %>%
    summarise(transposons.per.chromosome = (mean_transposon/mean_chromosome),
              plasmids.per.chromosome = (mean_plasmid/mean_chromosome),
              transposons.per.plasmid = (mean_transposon/mean_plasmid)) %>%
    ungroup() %>%
    pivot_longer(cols = c(transposons.per.chromosome,plasmids.per.chromosome,transposons.per.plasmid),
                 names_to = "ratio_type", values_to = "ratio") %>%
    ## add a Day column to compare to ancestral samples.
    mutate(Day = 1)

Tet5.ratio.plot <- evolved.replicon.coverage.ratio.df %>%
    filter(Tet == 5) %>%
    ggplot(aes(y = ratio, x = Plasmid, color = ratio_type, shape = Transposon)) +
    geom_point() +
    theme_classic() +
    facet_wrap(ratio_type~Transposon, scales = "free") +
    theme(strip.background = element_blank()) +
    ggtitle("5 ug/mL tetracycline, Day 1") +
    guides(color= "none", shape = "none")

Tet0.ratio.plot <- evolved.replicon.coverage.ratio.df %>%
    filter(Tet == 0) %>%
    ggplot(aes(y = ratio, x = Plasmid, color = ratio_type, shape = Transposon)) +
    geom_point() +
    theme_classic() +
    facet_wrap(ratio_type~Transposon, scales = "free") +
    theme(strip.background = element_blank()) +
    ggtitle("0 ug/mL tetracycline, Day 1") +
    guides(color = "none", shape = "none")

library(cowplot)

ratio.figure <- plot_grid(Tet5.ratio.plot, Tet0.ratio.plot, labels=c('A','B'),nrow=2)
## The warning messages in saving the plot correspond to the missing ratios
## for the 'No plasmid' strains for which transposon.per.plasmid
## and plasmid.per.chromosome are undefined (NA).



ggsave("../results/one-day-expt-coverage-ratios.pdf", ratio.figure)

## let's write out the table too.
write.csv(evolved.replicon.coverage.ratio.df, "../results/one-day-expt-plasmid-transposon-coverage-ratios.csv",
          quote=F, row.names=FALSE)

################################################################################
## Figure 2A.

## This is the big data frame for making Figure 2A.
ancestral.and.evolved.replicon.coverage.ratio.df <- full_join(
    ancestral.replicon.coverage.ratio.df, evolved.replicon.coverage.ratio.df) %>%
    ungroup()

Fig2A.df <- ancestral.and.evolved.replicon.coverage.ratio.df %>%
    filter(ratio_type == "transposons.per.chromosome") %>%
    ## update the names of the Plasmid factor for a prettier plot.
    mutate(Plasmid = fct_recode(as.factor(Plasmid),
                                `No plasmid` = "None",
                                p15A = "p15A")) %>%
    ## update the names of the Transposon factor for a prettier plot.
    mutate(Transposon = fct_recode(as.factor(Transposon),
                                   Inactive = "B59",
                                   Active = "B30")) %>%
    ## turn Tet, Population, Day into discrete factors for plotting.
    mutate(Tet = as.factor(Tet)) %>%
    mutate(Population = as.factor(Population)) %>%
    mutate(Day = as.factor(Day))

## Write Figure 2A Source Data.
write.csv(Fig2A.df, "../results/Source-Data/Fig2A-Source-Data.csv", row.names=FALSE, quote=FALSE)

## Make Figure 2A.
Fig2A <- ggplot(Fig2A.df, aes(x = Day, y = ratio, color = Transposon, shape = Tet)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    theme(strip.background = element_blank()) +
    ylab("tetA-transposons per chromosome")
ggsave("../results/Fig2A.pdf", Fig2A, width=4.5, height=3.25)


