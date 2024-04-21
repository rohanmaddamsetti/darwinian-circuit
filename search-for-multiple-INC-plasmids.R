## search-for-multiple-INC-plasmids.R by Rohan Maddamsetti.
## use Hye-in's MobSuite results to look for multiple INC plasmids in the same host genome.
## see documentation here about columns for replicon typing: https://github.com/phac-nml/mob-suite


library(tidyverse)

mob.suite.results <- read.csv("../data/Maddamsetti2024_FileS5-MOBTyper-plasmid-annotations.csv")

simple.mob.suite.results <- mob.suite.results %>%
    select(sample_id, rep_type.s., rep_type_accession.s.) %>%
    filter(rep_type.s. != "-") %>%
    filter(rep_type_accession.s. != "-") %>%
    separate(sample_id, into = c("GenomeID", "SeqID"), sep = "_genomic_", remove = FALSE)

write.csv(simple.mob.suite.results, "../results/simple-MOBTyper-plasmid-annotations.csv", row.names=FALSE)

multiple.inc.group.genomes.df <- simple.mob.suite.results %>%
    group_by(GenomeID, rep_type.s., rep_type_accession.s.) %>%
    summarize(plasmid_count = n()) %>%
    filter(plasmid_count > 1) %>%
    arrange(desc(plasmid_count))
write.csv(multiple.inc.group.genomes.df, "../results/candidate_multiple_INC_group_genomes.csv", row.names=FALSE)

relevant.mob.suite.results <- simple.mob.suite.results %>%
    filter(GenomeID %in% multiple.inc.group.genomes.df$GenomeID)
write.csv(relevant.mob.suite.results, "../results/candidate_multiple_INC_group_plasmids.csv", row.names=FALSE)


## TODO: examine PCN copy number -- are any of these plasmids high copy number?
