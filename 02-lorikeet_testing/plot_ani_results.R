library(ape)
library(ggtree)
library(data.table)
library(ggplot2)
library(heatmaply)
library(ggpubr)
library(grid)
library(gridExtra)

viral_ani <- read.table("/srv/projects2/debruijn_profiling/02-mocks/09-lorikeet_tests/01.1-viral_only_one_sample/00-ani/ANIm_percentage_identity.tab")
comm_two_ani <- read.table("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/02-community_two/02-parsnp_out/00-ani/ANIm_percentage_identity.tab")
comm_four_ani <- read.table("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/04-community_four/02-parsnp_out/00-ani/ANIm_percentage_identity.tab")
comm_five_ani <- read.table("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/05-community_five/02-parsnp_out/00-ani/ANIm_percentage_identity.tab")

comm_two_abd <- fread("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/02-community_two/abundances.tsv")
comm_four_abd <- fread("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/04-community_four/abundances.tsv")
comm_five_abd <- fread("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/05-community_five/abundances.tsv")

#### COMMUNITY TWO ####
comm_two_plt <- ggheatmap(as.matrix(comm_two_ani), colors = inferno(n=256)) 
ggsave("comm_two_plt.png", comm_two_plt, device="png", width = 10, height = 6)


heatmaply(comm_two_ani, 
          file = "/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/02-community_two/comm_two_heatmap.html",
          k_row = 2, k_col = 2)

comm_two_tab <- ggtexttable(t(comm_two_abd),
                             cols = c("Sample 1", "Sample 2", 
                                      "Sample 3", "Sample 4", 
                                      "Sample 5", "Sample 6", 
                                      "Sample 7", "Sample 8",
                                      "Sample 9", "Sample 10"),
                             theme=ttheme(
                               colnames.style = colnames_style(size = 8),
                               rownames.style = rownames_style(size = 8)))
ggsave("comm_two_tab.png", comm_two_tab, device="png", width = 10, height = 6)


comm_two_grid <- ggarrange(comm_two_plt, comm_two_tab, ncol = 1, nrow = 2)
comm_two_grid
ggsave("comm_two_grid.png", comm_two_grid, device="png", width = 12, height = 6)


#### COMMUNITY FOUR ####
comm_four_plt <- ggheatmap(as.matrix(comm_four_ani), colors = inferno(n=256)) 
ggsave("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/04-community_four/comm_four_plt.png", 
       comm_four_plt, device="png", width = 10, height = 6)


heatmaply(scale(comm_four_ani), 
          file = "/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/04-community_four/comm_four_heatmap.html",
          k_row = 2, k_col = 2)

comm_four_tab <- ggtexttable(t(comm_four_abd),
                             cols = c("Sample 1", "Sample 2", 
                                      "Sample 3", "Sample 4"),
                             theme=ttheme(
                               colnames.style = colnames_style(size = 8),
                               rownames.style = rownames_style(size = 8)))
ggsave("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/04-community_four/comm_four_tab.png", 
       comm_four_tab, device="png", width = 10, height = 6)


comm_four_grid <- ggarrange(comm_four_plt, comm_four_tab)
comm_four_grid
ggsave("comm_four_grid.png", comm_four_grid, device="png", width = 12, height = 6)

#### COMMUNITY FIVE #####
hm_5 <- heatmaply(scale(comm_five_ani), 
          file = "/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/05-community_five/comm_five_heatmap.html",
          k_row = 2, k_col = 2, plot_method="ggplot")

comm_five_plt <- ggheatmap(as.matrix(comm_five_ani), colors = inferno(n=256)) 
ggsave("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/05-community_five/comm_five_plt.png", 
       comm_five_plt, device="png", width = 10, height = 6)

comm_five_tab <- ggtexttable(t(comm_five_abd),
                             cols = c("Sample 1", "Sample 2", 
                                      "Sample 3", "Sample 4"),
                             theme=ttheme(
                               colnames.style = colnames_style(size = 8),
                               rownames.style = rownames_style(size = 8)))
ggsave("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/05-community_five/comm_five_tab.png", 
       comm_five_tab, device="png", width = 10, height = 6)

comm_five_grid <- ggarrange(comm_five_plt, comm_five_tab)
comm_five_grid
ggsave("/srv/projects/abisko/rhys/01-projects/03-lorikeet_testing/05-strains/05-community_five/comm_five_grid.png", 
       comm_four_grid, device="png", width = 12, height = 6)

## viral
viral_plt <- ggheatmap(as.matrix(viral_ani), colors = inferno(n=256)) 
viral_abd <- data.table()
viral_abd[, GCF_000919375.1:=c(10, 20, 0)]
viral_abd[, GCF_001505895.1:=c(10, 0, 20)]
viral_tab <- ggtexttable(t(viral_abd),
                             cols = c("Sample 1", "Sample 2", 
                                      "Sample 3"),
                             theme=ttheme(
                               colnames.style = colnames_style(size = 8),
                               rownames.style = rownames_style(size = 8)))

ggsave("/srv/projects2/debruijn_profiling/02-mocks/09-lorikeet_tests/01.1-viral_only_one_sample/viral_plt.png", 
       viral_plt, device="png", width = 10, height = 6)
ggsave("/srv/projects2/debruijn_profiling/02-mocks/09-lorikeet_tests/01.1-viral_only_one_sample/viral_tab.png", 
       viral_tab, device="png", width = 10, height = 6)
