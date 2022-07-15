library(data.table)

comm_two_ani <- fread("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/02-community_two/02-parsnp_out/00-ani/ANIm_percentage_identity.tab")
comm_four_ani <- fread("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/04-community_four/02-parsnp_out/00-ani/ANIm_percentage_identity.tab")
comm_five_ani <- fread("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/05-community_five/02-parsnp_out/00-ani/ANIm_percentage_identity.tab")

cami_genome4 <- fread("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/13-camisim/00-small_metagenomes/source_genomes/00-ani/ani/ANIm_percentage_identity.tab")

comm_two_abd <- fread("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/02-community_two/abundances.tsv")
comm_four_abd <- fread("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/04-community_four/abundances.tsv")
comm_five_abd <- fread("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/05-community_five/abundances.tsv")