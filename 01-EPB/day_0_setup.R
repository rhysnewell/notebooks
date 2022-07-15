library(data.table)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Read in checkm, coverm, and GTDB

com <- fread("/srv/projects/abisko/rhys/01-projects/04-bioreactor_slamM/00-d0_assembly/data/scaffolds1/checkm_table.txt")
rel <- fread("/srv/projects/abisko/rhys/01-projects/04-bioreactor_slamM/00-d0_assembly/data/scaffolds1/coverm_rel_abundances.tsv")
arc <- fread("/srv/projects/abisko/rhys/01-projects/04-bioreactor_slamM/00-d0_assembly/data/scaffolds1/gtdbtk/gtdbtk.ar122.summary.tsv")
bac <- fread("/srv/projects/abisko/rhys/01-projects/04-bioreactor_slamM/00-d0_assembly/data/scaffolds1/gtdbtk/gtdbtk.bac120.summary.tsv")

# Combine bac and rc gtdbtk
pro <- rbind(bac, arc)

# Join to completeness
com <- data.table(inner_join(inner_join(com, pro[, .(user_genome, classification, fastani_ani)], by=c("bin_id"="user_genome")), rel, by=c("bin_id"="Genome")))

# Join to relative abundances
rel <- data.table(full_join(rel, com[, .(bin_id, marker_lineage, completeness, contamination, strain_heterogeneity, classification)], by=c("Genome"="bin_id")))
