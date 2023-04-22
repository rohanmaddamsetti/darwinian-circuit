## nine-day-expt-tetA-GFP-barcode-genomics.R by Rohan Maddamsetti.

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggthemes)
library(viridis)
library(ggrepel)

## K-12 MG1655 oriC replication origin annotation
## annotated as rep_origin in the genbank file.
## Also see: https://biocyc.org/ECOLI/NEW-IMAGE?type=EXTRAGENIC-SITE&object=G0-10506
## 3,925,744 -> 3,925,975
K12_oriC_START = 3925744
K12_oriC_END = 3925975
K12_oriC_MID = (K12_oriC_START+K12_oriC_END)/2


rotate.K12.chr <- function(my.position) {
    #' function to rotate genome coordinates,
    #' setting oriC at the center of plots
    #' that examine mutation bias over the chromosome.
    ## we want to change coordinates so that c is the new origin.
    GENOME.LENGTH <- 4641652
    midpoint <- GENOME.LENGTH/2
    oriC <- 3925860
    
    if (oriC >= midpoint) {
        L <- oriC - midpoint
        ifelse(my.position > L, my.position - oriC, GENOME.LENGTH - oriC + my.position)
    } else { ## midpoint is greater than new.origin.
        L <- midpoint + oriC
        ifelse(my.position > L, my.position - GENOME.LENGTH - oriC, my.position - oriC)
    }
}


## colorblind-friendly palette.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## assert that we are in the src directory, such that
## projdir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("darwinian-circuit","src")))
projdir <- file.path("..")

## This is the key data file for the analysis.
evolved.mutations <- read.csv(
    file.path(projdir,
              "results/nine-day-GFP-barcode-expt-genome-analysis/evolved_mutations.csv"),
    stringsAsFactors=FALSE) %>%
    mutate(Mbp.coordinate=Position/1000000) %>%
    ## update the names of the Plasmid factor for a prettier plot.
    mutate(Plasmid_factor = fct_recode(as.factor(Plasmid),
                                       `No plasmid` = "None",
                                       p15A = "p15A",
                                       pUC = "pUC"))


evolved.MOB <- filter(evolved.mutations, Mutation=='MOB')
###############################################
## Supplementary Figure S1A: make a stacked bar plot of the kinds of mutations in each treatment.

## This function sums mutations per replicate population.
make.mutation.class.df <- function(evolved.mutations.df) {
    evolved.mutations.df %>%
        ## give nicer names for mutation classes.
        mutate(Mutation=recode(Mutation,
                               MOB = "Mobile element transposition",
                               DEL = "Indel",
                               INS = "Indel",
                               SUB = "Multiple-base substitution",
                               nonsynonymous = "Nonsynonymous",
                               synonymous = "Synonymous",
                               nonsense = "Nonsense",
                               pseudogene = "Pseudogene",
                               intergenic = "Intergenic",
                               )) %>%
        group_by(Sample, Plasmid, Population, Mutation, Plasmid_factor) %>%
        summarize(Count=n()) %>%
        ungroup() %>%
        data.frame() %>%
        mutate(Mutation=as.factor(as.character(Mutation)))
}


plot.mutation.summary.stackbar <- function(mutation.class.df, leg=FALSE) {
    fig <- ggplot(mutation.class.df, aes(x=Population, y=Count, fill=Mutation)) +
        ylab("Count") +
        facet_wrap(. ~ Plasmid) +
        geom_bar(stat='identity') +
        scale_x_continuous(breaks=seq(0, 10, 1)) +
        scale_fill_brewer(palette = "RdYlBu", direction=-1,drop=FALSE) +        
        theme_classic(base_family='Helvetica') +
        theme(
            axis.text.y=element_text(size=12),
            panel.border=element_blank(),
            strip.background = element_blank(),
            panel.spacing.x=unit(1, "cm"),
            panel.spacing.y=unit(0.5, "cm"))
    
    if (leg == TRUE) {
        fig <- fig +
            theme(legend.title=element_text(size=8, face="bold"),
                  legend.title.align=0.5,
                  legend.text=element_text(size=8),
                  legend.position="bottom")
    } else {
        fig <- fig + guides(fill = "none")
    }
    
    return(fig)
}

## Now make Figure S1A.
mutation.class.df <- make.mutation.class.df(evolved.mutations)

S1FigA <- plot.mutation.summary.stackbar(mutation.class.df, TRUE)
S1figA.output <- "../results/S1FigA.pdf"
ggsave(S1FigA, file=S1figA.output,width=6,height=5)


#################################################################################
## analysis of parallel evolution at the same nucleotide.
## discuss numbers and finding in the text (no figure.).
## This could be a Supplementary Table.

 bp.parallel.mutations <- evolved.mutations %>% group_by(Position) %>%
summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.MOB <- filter(bp.parallel.mutations,Mutation=='MOB')
## no parallel DEL muts at bp level.
parallel.INS <- filter(bp.parallel.mutations,Mutation=='INS')
parallel.dN <- filter(bp.parallel.mutations,Mutation=='nonsynonymous')
parallel.dS <- filter(bp.parallel.mutations,Mutation=='synonymous')

## examine parallel evolution at amino acid level
parallel.AA.dN <- evolved.mutations %>% filter(Mutation=='nonsynonymous') %>% group_by(Position) %>% summarize(count=n()) %>% filter(count > 1)
parallel.dN.Table <- filter(evolved.mutations, Position %in% parallel.AA.dN$Position) %>% arrange(Position)

## check parallel evolution for synonymous mutations too.
parallel.dS.Table <- filter(evolved.mutations, Position %in% parallel.dS$Position) %>% arrange(Position)


#################################################################################
## analysis of parallel evolution at the gene level (including intergenic regions).

gene.level.parallel.mutations <- evolved.mutations %>% group_by(Gene) %>%
summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.genes <- gene.level.parallel.mutations %>%
    select(Gene, count, Plasmid) %>%
    distinct() %>%
    arrange(desc(count))

#################################################################################
### S1 Supplementary Figure S1B : make a matrix plot of genes with mutations in two or more clones.
################################################################################
MakeMutCountMatrixFigure <- function(evolved.muts, show.all=FALSE) {

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.muts %>%
        group_by(Gene, Sample, Plasmid, Plasmid_factor) %>%
        summarize(mutation.count = n()) %>%
        ## This is for sorting mutations.
        mutate(is.MOB = ifelse(str_detect(Gene,"tetA-Tn5"), TRUE, FALSE))
    
    total.muts <- matrix.data %>%
        group_by(Gene) %>%
        summarize(total.mutation.count = sum(mutation.count))
    
    matrix.data <- left_join(matrix.data, total.muts)
    
    if (!show.all) { ## then filter out genes that are only hit in one sample.
        matrix.data <- matrix.data %>%
            filter(total.mutation.count > 1)
    }
    
    ## sort genes by number of mutations in each row, but put all the transposon mutations together.
    ## also check out the alternate sorting method that follows.
    gene.hit.sort <- matrix.data %>%
        group_by(Gene, is.MOB, .drop = FALSE) %>%
        summarize(hits=sum(mutation.count)) %>%
        arrange(desc(is.MOB), desc(hits))
    ## now sort genes.
    matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.hit.sort$Gene))
    
    ## cast mutation.count into a factor for plotting.
    matrix.data$mutation.count <- factor(matrix.data$mutation.count)

    make.matrix.panel <- function(mdata, plasmid, leg=FALSE) {
        panel.data <- filter(mdata, Plasmid_factor == plasmid)
        fig <- ggplot(panel.data,
                      aes(x=Sample,
                          y=Gene,
                          fill=mutation.count,
                          frame= Plasmid_factor)
                      ) +
            geom_tile(color="black",size=0.1) +
            ggtitle(plasmid) +
            theme_tufte(base_family='Helvetica') +
            theme(axis.ticks = element_blank(),
                  axis.text.x = element_text(size=10,angle=45,hjust=1),
                  axis.text.y = element_text(size=10,hjust=1,face="italic"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank()) +
            scale_y_discrete(drop=FALSE) + ## don't drop missing genes.
            scale_fill_manual(name="Mutations",
                              values = c("#ffdf00", "#bebada", "#fb8072", "#80b1d3", "#fdb462"))
        
        if (leg == FALSE) {
            fig <- fig + guides(fill = "none")
        }
        return(fig)
    }

    
    ## make panels.
    noPlasmid.matrix.panel <- make.matrix.panel(matrix.data, "No plasmid")
    ## Remove the gene labels to save space.
    p15A.matrix.panel <- make.matrix.panel(matrix.data,"p15A") +
        theme(axis.text.y=element_blank())
    pUC.matrix.panel <- make.matrix.panel(matrix.data, "pUC") +
        theme(axis.text.y=element_blank())
    
    ## Using the patchwork library for layout.
    matrix.figure <-
        noPlasmid.matrix.panel +
        p15A.matrix.panel +
        pUC.matrix.panel +
        plot_layout(nrow = 1)
    return(matrix.figure)
}


## S1 Supplementary Figure Panel B.
S1FigB <- MakeMutCountMatrixFigure(evolved.mutations, show.all=TRUE)
S1FigB.outf <- "../results/S1FigB.pdf"
ggsave(S1FigB.outf, S1FigB, height=8, width=7)

S2Fig <- MakeMutCountMatrixFigure(evolved.mutations, show.all=FALSE)
S2Fig.outf <- "../results/S2Fig.pdf"
ggsave(S2Fig.outf, S2Fig, height=5, width=7)

        
