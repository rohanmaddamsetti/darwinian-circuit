##copy-number-analysis.R by Rohan Maddamsetti.

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
library(IRanges)
library(GenomicRanges)
library(rtracklayer)


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


#' Find intervals longer than max.read.len that reject H0 coverage in genome.
#' at an uncorrected alpha = 0.05. This is to have generous predicted boundaries for amplifications.
#' Then do a more rigorous test for each region. take positions in the region separated by more than max.read.len,
#' and determine the probability that all are independently significant under the null, compared to
#' a corrected bonferroni. max.read.len ensures positions cannot be spanned by a single Illumina read.
#' Estimate copy number by dividing mean coverage in each region by the mean of the H0 1x coverage distribution.
#' return mean copy number, and boundaries for each region that passes the amplification test.

find.K12.candidate.amplifications <- function(breseq.output.dir, gnome) { #gnome is not a misspelling.

    gnome <- as.character(gnome)
    print(gnome)
    ## Use xml2 to get negative binomial fit and relative variance from
    ## breseq output summary.html. This is H0 null distribution of 1x coverage.
    nbinom.fit <- coverage.nbinom.from.html(breseq.output.dir) %>%
        filter(replicon=="chromosome")
    
    ## Use xml2 to get max read length from summary.html.
    max.read.len <- max.readlen.from.html(breseq.output.dir)
    genome.length <- 4641652 ## length of K-12 MG1655 reference.
    my.size.parameter <- nbinom.fit$mean^2/(nbinom.fit$variance - nbinom.fit$mean)
    alpha <- 0.05
   
    uncorrected.threshold <- qnbinom(p=alpha, mu=nbinom.fit$mean, size=my.size.parameter, lower.tail=FALSE)
    
    genome.coverage.file <- file.path(breseq.output.dir,"08_mutation_identification", "NC_000913.coverage.tab")
    
    ## use dtplyr for speed!
    genome.coverage <- lazy_dt(fread(genome.coverage.file)) %>%
        select(position,unique_top_cov,unique_bot_cov) %>% mutate(coverage=unique_top_cov+unique_bot_cov)
    
    ## find candidate amplifications that pass the uncorrected threshold.
    candidate.amplifications <- genome.coverage %>%
        filter(coverage > uncorrected.threshold) %>%
        ## now finally turn into a tibble.
        as_tibble()
    
    ## calculate intervals of candidate amplifications.
    boundaries <- candidate.amplifications %>%
        mutate(left.diff=position - lag(position)) %>%
        mutate(right.diff=lead(position) - position) %>%
        ## corner case: check for the NA values at the endpoints and set them as boundaries.
        mutate(is.right.boundary=is.na(right.diff)|ifelse(right.diff>1,TRUE,FALSE)) %>%
        mutate(is.left.boundary=is.na(left.diff)|ifelse(left.diff>1,TRUE,FALSE)) %>%
        filter(is.left.boundary==TRUE | is.right.boundary==TRUE)
 
    left.boundaries <- filter(boundaries,is.left.boundary==TRUE) %>%
        arrange(position)
        
    right.boundaries <- filter(boundaries,is.right.boundary==TRUE) %>%
        arrange(position)
    
    assert_that(nrow(left.boundaries) == nrow(right.boundaries))
    
    ## helper higher-order function to get min, max, mean coverage of each segment.
    get.segment.coverage <- function(left.bound,right.bound,coverage.table,funcx) {
        seg <- coverage.table %>% filter(position>left.bound) %>% filter(position<right.bound)
        return(funcx(seg$coverage))
    }
    
    amplified.segments <- data.frame(left.boundary=left.boundaries$position,right.boundary=right.boundaries$position) %>%
        ## filter out intervals less than 2 * max.read.len.
        mutate(len=right.boundary-left.boundary) %>%
    filter(len>(2*max.read.len)) %>% mutate(amplication.index=row_number())

    ## return empty dataframe  if there are no significant amplified segments.
    if (nrow(amplified.segments) == 0) return(data.frame())

    amplified.segments <- amplified.segments %>%
        ## find min, max, and mean coverage of each amplified segment.
        group_by(left.boundary,right.boundary) %>%
        summarise(coverage.min=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,min),
                  coverage.max=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,max),
                  coverage.mean=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,mean)) %>%
        mutate(len=right.boundary-left.boundary) %>%
        mutate(copy.number.min=coverage.min/nbinom.fit$mean,copy.number.max=coverage.max/nbinom.fit$mean,
               copy.number.mean=coverage.mean/nbinom.fit$mean) %>%
        ## annotate with the sample name.
        mutate(Sample=as.character(gnome))

    return(amplified.segments)
}


find.K12.chromosomal.amplifications <- function(breseq.output.dir, gnome) { #gnome is not a misspelling.
  
    amplified.segments <- find.K12.candidate.amplifications(breseq.output.dir, gnome)
    ## handle the case that there are no amplified.segments (empty dataframe).
    if (nrow(amplified.segments) == 0) return(amplified.segments)
    
    ## Use xml2 to get negative binomial fit and relative variance from
    ## breseq output summary.html. This is H0 null distribution of 1x coverage.
    nbinom.fit <- coverage.nbinom.from.html(breseq.output.dir) %>%
        filter(replicon=="chromosome")  
    ## Use xml2 to get max read length from summary.html.
    max.read.len <- max.readlen.from.html(breseq.output.dir)
    genome.length <- 4641652 ## length of K-12 MG1655 reference.
    my.size.parameter <- nbinom.fit$mean^2/(nbinom.fit$variance - nbinom.fit$mean)
    alpha <- 0.05

    ## divide alpha by the number of tests for the bonferroni correction.
    bonferroni.alpha <- alpha/(genome.length + sum(amplified.segments$len))
    corrected.threshold <- qnbinom(p = bonferroni.alpha, mu = nbinom.fit$mean, size = my.size.parameter, lower.tail=FALSE)
    
    ## This is my test: take the probability of the minimum coverage under H0 to the power of the number of
    ## uncorrelated sites in the amplification (sites more than max.read.len apart). Then see if this is smaller than the
    ## bonferroni corrected p-value for significance..
    significant.amplifications <- amplified.segments %>%
        mutate(pval=(pnbinom(q = coverage.min,
                             mu = nbinom.fit$mean,
                             size = my.size.parameter,
                             lower.tail=FALSE))^(len%/%max.read.len)) %>%
        mutate(is.significant = ifelse(pval < bonferroni.alpha, TRUE, FALSE)) %>%
        filter(is.significant==TRUE) %>%
        mutate(bonferroni.corrected.pval=pval*alpha/bonferroni.alpha)
    
    return(significant.amplifications)
}


annotate.sample.amplifications <- function(sample.amplifications) {

    ancestor.gff <- unique(sample.amplifications$gff_path)
    
    ## create the IRanges object.
    amp.ranges <- IRanges(sample.amplifications$left.boundary,
                          sample.amplifications$right.boundary)
    ## Turn into a GRanges object in order to find overlaps with K12 genes.
    g.amp.ranges <- GRanges("NC_000913", ranges=amp.ranges)
    ## and add the data.frame of sample.amplifications as metadata.
    mcols(g.amp.ranges) <- sample.amplifications
    
    ## find the genes within the amplifications.
    ancestor.gff.data <- import.gff(ancestor.gff)
    ancestor.Granges <- as(ancestor.gff.data, "GRanges")
    
    ancestor.genes <- ancestor.Granges[ancestor.Granges$type == 'gene']
    ## find overlaps between annotated genes and amplifications.
    hits <- findOverlaps(ancestor.genes,g.amp.ranges,ignore.strand=FALSE)
    
    ## take the hits, the ancestor annotation, and the amplifications,
    ## and produce a table of genes found in each amplication.
    
    hits.df <- data.frame(query.index=queryHits(hits),subject.index=subjectHits(hits))
    
    query.df <- data.frame(query.index=seq_len(length(ancestor.genes)),
                           gene=ancestor.genes$Name,locus_tag=ancestor.genes$ID,
                           start=start(ranges(ancestor.genes)),end=end(ranges(ancestor.genes)))
    
    subject.df <- bind_cols(data.frame(subject.index=seq_len(length(g.amp.ranges))),data.frame(mcols(g.amp.ranges)))
    
    amplified.genes.df <- left_join(hits.df,query.df) %>% left_join(subject.df) %>%
        ## if gene is NA, replace with locus_tag. have to change factors to strings!
        mutate(gene = ifelse(is.na(gene),as.character(locus_tag),as.character(gene)))
    
    return(amplified.genes.df)
}


annotate.amplifications <- function(amps.with.ancestors) {
    amps.with.ancestors %>% split(.$Sample) %>%
        map_dfr(.f = annotate.sample.amplifications)    
}


plot.amp.segments <- function(annotated.amps) {
    
    labeled.annotated.amps <- annotated.amps %>%
        mutate(log2.copy.number.mean=log2(copy.number.mean)) %>%
        mutate(left.boundary.MB = left.boundary/1000000) %>%
        mutate(right.boundary.MB = right.boundary/1000000)
    
    ## order the genes by start to get axes correct on heatmap.
    labeled.annotated.amps$gene <- with(labeled.annotated.amps, reorder(gene, start))
    ## reverse the order of genomes to make axes consistent with stacked barplot.
    labeled.annotated.amps$Sample <- factor(labeled.annotated.amps$Sample)
    labeled.annotated.amps$Sample <- factor(labeled.annotated.amps$Sample,
                                            levels=rev(levels(labeled.annotated.amps$Sample)))
    
    segmentplot <- ggplot(
        labeled.annotated.amps,
        aes(x=left.boundary.MB,
            xend=right.boundary.MB,
            y=Sample,
            yend=Sample,
            color=copy.number.mean,
            size=20,
            frame=Plasmid)) +
        geom_segment() +
        ## draw vertical lines at acrBAR starts.
        geom_vline(size=0.2,
                   color = 'red',
                   linetype = 'dashed',
                   xintercept = c(480478/1000000, 484843/1000000, 484985/1000000)) +
        xlab("K-12 chromosomal position (Mb)") +
        ylab("") +
        scale_color_viridis(name="copy number",option="plasma") +
        facet_wrap(.~Plasmid,ncol=1, scales = "free_y") +
        theme_classic(base_family='Helvetica') +
        guides(size= "none") +
        theme(legend.position="bottom") +
        theme(axis.ticks=element_line(size=0.1))
    return(segmentplot)
}

#######################################
## Analysis time!

## assert that we are in the src directory, such that
## proj.dir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("darwinian-circuit","src")))
projdir <- file.path("..")


## helper function for adding metadata for RM7.140.31-60.
Evolved.Sample.ID.to.Plasmid <- function(sample.id) {
    ## get rid of the prefix and turn the rest into a number.
    sample.num <- strtoi(str_replace(sample.id, "^RM7-140-", "" ))
    if (sample.num <= 40)
        plasmid.type = "No plasmid"
    else if (sample.num <= 50)
        plasmid.type = "p15A"
    else if (sample.num <= 60)
        plasmid.type = "pUC"
    else
        stop("sample from nine-day tetA-GFP barcode experiment out of range")
    return(plasmid.type)
}


## helper function for adding metadata for RM7.140.31-60.
Evolved.Sample.ID.to.Population <- function(sample.id) {
    ## get rid of the prefix and turn the rest into a number.
    sample.num <- strtoi(str_replace(sample.id, "^RM7-140-", "" ))
    population.num <- ((sample.num - 30) %% 10)
    if (population.num == 0)
        population.num <- 10
    return(population.num)
}


clone.data.dir <- file.path(projdir, "results", "nine-day-GFP-barcode-expt-genome-analysis")
## don't match the GFF files.
all.clones <- Filter(function(x) return(!str_detect(x, "gff")), list.files(clone.data.dir,pattern='^RM7'))
all.clone.paths <- sapply(all.clones, function(x) file.path(clone.data.dir,x))
clone.input.df <- data.frame(Sample=all.clones, path=all.clone.paths)

ancestral.clone.input.df <- filter(clone.input.df, !str_detect(Sample, "^RM7-140-"))
evolved.clone.input.df <- filter(clone.input.df, str_detect(Sample, "^RM7-140-"))

## generate evolved clone metadata from the Sample column.
evolved.clone.metadata <- evolved.clone.input.df %>%
    select(Sample) %>%
    ## get the plasmid type from the sample ID.
    mutate(Plasmid = sapply(Sample, Evolved.Sample.ID.to.Plasmid)) %>%
    ## get the population number from the sample ID.
    mutate(Population = sapply(Sample, Evolved.Sample.ID.to.Population)) %>%
    ## need a GFF file for annotating chromosomal amplifications.
    mutate(gff_path = "../results/nine-day-GFP-barcode-expt-genome-analysis/LCA.gff3")

## The no plasmid ancestor clones are RM7.107.3,4,5,6,7,8,9,10,11,13.
## The p15A ancestor clones are RM7.106.3,5,6,7,8,10,11,12,13,15.
## The pUC ancestor clones are RM7.107.38,39,41,42,44,45,46,48,49,56.
ancestral.clone.metadata <- data.frame(
    Sample = c(
        "RM7-107-3","RM7-107-4","RM7-107-5","RM7-107-6","RM7-107-7",
        "RM7-107-8","RM7-107-9","RM7-107-10","RM7-107-11","RM7-107-13",
        "RM7-106-3","RM7-106-5","RM7-106-6","RM7-106-7","RM7-106-8",
        "RM7-106-10","RM7-106-11","RM7-106-12","RM7-106-13","RM7-106-15",
        "RM7-107-38","RM7-107-39","RM7-107-41","RM7-107-42","RM7-107-44",
        "RM7-107-45","RM7-107-46","RM7-107-48","RM7-107-49","RM7-107-56"),
    Plasmid = c(
        rep("No plasmid",10), rep("p15A",10), rep("pUC",10)),
    Population = c(
        seq(1:10), seq(1:10), seq(1:10)),
    ## need a GFF file for annotating chromosomal amplifications.
    gff_path = rep("../results/nine-day-GFP-barcode-expt-genome-analysis/LCA.gff3", 30))


######################################################################
## Plot the plasmid/chromosome and transposon/chromosome ratio in each sample.

## Get the actual coverage for the B31-N20 transposons in each clone.
## This is calculated by get-nine-day-GFP-barcode-expt-transposon-coverage.py.
transposon.coverage.file <- file.path(projdir, "results", "nine-day-GFP-barcode-expt-genome-analysis", "transposon-coverage.csv")
transposon.coverage.df <- read.csv(transposon.coverage.file) %>%
    ## let's add metadata and rename columns for compatibility.
    dplyr::rename(mean = TransposonCoverage) %>%
    dplyr::mutate(replicon = "transposon") 

ancestral.transposon.coverage.df <- filter(transposon.coverage.df, Sample %in% ancestral.clone.input.df$Sample)
evolved.transposon.coverage.df <- filter(transposon.coverage.df, Sample %in% evolved.clone.input.df$Sample)

evolved.replicon.coverage.df <- map_dfr(.x = evolved.clone.input.df$path, .f = coverage.nbinom.from.html) %>%
    ## I am not examining dispersion or variance at this point.
    select(Sample, mean, replicon) %>%
    ## add transposon coverage data
    bind_rows(evolved.transposon.coverage.df) %>%
    ## set NA coverage values to zero.
    mutate(mean=sapply(mean, function(x) ifelse(is.na(x), 0,x))) %>%
    ## and add metadata.
    full_join(evolved.clone.metadata)


evolved.replicon.coverage.ratio.df <- evolved.replicon.coverage.df %>%
    pivot_wider(names_from = replicon, values_from = mean, names_prefix = "mean_") %>%
    group_by(Sample, Plasmid, Population) %>%
    summarise(transposons.per.chromosome = (mean_transposon/mean_chromosome),
              plasmids.per.chromosome = (mean_plasmid/mean_chromosome),
              transposons.per.plasmid = (mean_transposon/mean_plasmid)) %>%
    pivot_longer(cols = c(transposons.per.chromosome,plasmids.per.chromosome,transposons.per.plasmid),
                 names_to = "ratio_type", values_to = "ratio")
    

Fig1I.df <- evolved.replicon.coverage.ratio.df %>%
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

Fig1I <- ggplot(data=Fig1I.df, aes(x=Population, y=`Copy number`, fill=ratio_type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_classic() +
    facet_wrap(.~Plasmid, scales="free") +
    theme(legend.title=element_blank(), legend.position="bottom")
ggsave("../results/Fig1I.pdf", Fig1I, width=8, height=3.5)

## let's write out the table.
write.csv(evolved.replicon.coverage.ratio.df, "../results/nine-day-GFP-barcode-expt-genome-analysis/nine-day-GFP-barcode-expt-plasmid-transposon-coverage-ratios.csv", quote=F, row.names=FALSE)

########################################################
## Lingchong also asked for the following variations of Figure 1I.

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
## Supplementary Figure 1C. Heatmap figure for examining amplifications.

## Find chromosomal amplifications in the samples, and annotate with their reference gff file.
evolved.amps <- map2_df(evolved.clone.input.df$path,
                        evolved.clone.input.df$Sample,                              
                        ##find.K12.candidate.amplifications) %>% ## for uncorrected p-values
                        find.K12.chromosomal.amplifications) %>% ## for corrected p-values.
    ungroup() %>%
    ## add GFF column for annotate.amplifications function.
    mutate(gff_path = "../results/nine-day-GFP-barcode-expt-genome-analysis/LCA.gff3")

annotated.amps <- annotate.amplifications(evolved.amps) %>%
    ## get the plasmid type from the sample ID.
    ## we need this for the function plot.amp.segments.
    mutate(Plasmid = sapply(Sample, Evolved.Sample.ID.to.Plasmid))

## make a figure of the significant amplifications found with the bonferroni method.
amp.segment.plot <- plot.amp.segments(annotated.amps)
## show the plot.
ggsave("../results/S1FigC.pdf",
       amp.segment.plot, height=6.5, width=4)

parallel.amplified.genes <- annotated.amps %>%
    group_by(gene, locus_tag, start, end) %>%
    summarize(parallel.amplifications = n()) %>%
    arrange(desc(parallel.amplifications)) %>%
    filter(parallel.amplifications > 2)

acrABR.amps <- annotated.amps %>%
    filter(str_detect(gene, "acr"))

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

              
