
library(data.table)
source('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/abisko-stuff/abisko.R')
library(reshape2)
library(pheatmap)
library(ggplot2)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# from https://github.com/joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
gm_mean = function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}
# from https://github.com/joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
clr = function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}

read_uid_and_bin_names = function(){
  uid_and_bin_names = fread('/srv/projects/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/gtdb/tree_bacterial_set/export_U_numbers_and_bin_names.csv',header=F)
  setnames(uid_and_bin_names, names(uid_and_bin_names), c('genome_uid','fasta'))
  uid_and_bin_names[, genome := gsub('.fa$',x=fasta,'',perl=T)]
  uid_and_bin_names[, genome := gsub('.fasta.metabat-bins-','',genome)]
  return(uid_and_bin_names[, .(genome_uid, genome)])
}

read_phylum_wise_taxonomy = function(uid_and_bin_names){
  phylum_wise_taxonomy = fread('/srv/projects/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/gtdb/tree_bacterial_set/overall_poster_tree08092016.tree.taxonomy',header=F)
  setnames(phylum_wise_taxonomy, names(phylum_wise_taxonomy), c('genome_uid', 'taxonomy'))
  phylum_wise_taxonomy = split_taxonomy(phylum_wise_taxonomy)

  m = merge(phylum_wise_taxonomy, uid_and_bin_names, by='genome_uid')
  setnames(m, 'domain', 'manual_phylum')
  return(m[, .(genome, manual_phylum)])
}

read_graftm_s5_on_aterrible12 = function(){
  aterrible_graftm = fread('/srv/projects/abisko/shotgun_abundance/103_contig_tree_relative_abundances/9_s5.new.gpkg_vs_aterrible12_proteins.prokk3.fixed/prokka3_kingdom_specific.fixed/prokka3_kingdom_specific.fixed_read_tax.tsv',header=F)
  setnames(aterrible_graftm, c('V1','V2'), c('protein','s5_taxonomy'))
  aterrible_graftm[, genome := gsub('(.*)_.*','\\1',perl=T,x=protein)]
  aterrible_graftm = aterrible_graftm[genome %in% aterrible_graftm[, .N, by=genome][N==1]$genome]
  aterrible_graftm[, small_family := split_taxonomy(.SD, taxonomy_field='s5_taxonomy')[, substrRight(trim(paste(domain, phylum, class_name, order_name, family)), 40)]]
  return(aterrible_graftm[,.(genome, small_family, s5_taxonomy)])
}

## link to dereplicated genome list
read_dereplicated_genomes_clusters = function(all_genomes){
  cluster_table = fread('/srv/projects/abisko/aterrible_bins/14_dereplicated647/cluster_table_groups.tsv', header=T)
  cluster_table[order(-quality), representative := .SD$genome[1], by=group]
  cluster_table[, genome := gsub('.fasta.metabat-bins-','',genome)]
  g2 = merge(cluster_table, all_genomes, by='genome', all.y=T)
  g2[is.na(representative), representative := genome] #singletons aren't in the cluster file
  g2[, representative := gsub('.fasta.metabat-bins-','',representative)]
  return(g2[, .(genome, representative)])
}

## Link KO table to genome
read_targeted_gene_and_ko = function(){
  gene_and_ko = fread('/srv/projects/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/prokka3_kingdom_specific.ko_via_uniprot.csv',header=F)
  setnames(gene_and_ko, names(gene_and_ko), c('gene','uniprot_id','ko'))
  target_ko_list = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/key_genes20160908.csv',header=T)
  setnames(target_ko_list, names(target_ko_list)[ncol(target_ko_list)], 'crud')
  target_ko_list[, crud := NULL]
  target_ko_list[, spreadsheet_order := 1:.N]
  setnames(target_ko_list, 'KO', 'ko')
  gene_and_ko_targeted = merge(gene_and_ko, target_ko_list, by='ko', all.y=T)
  gene_and_ko_targeted[, genome := gsub('.unique_contig_names.*','',x=gene)]
  gene_and_ko_targeted[,genome := gsub('(.*)_.*','\\1',genome,perl=T)]
  #setnames(gene_and_ko_targeted, grep('Note:', names(gene_and_ko_targeted), value=T), 'note')
  #gene_and_ko_targeted[, note := NULL]
  gene_and_ko_targeted[, gpkg_complete := NULL]
  gene_and_ko_targeted[, gpkg_name := NULL]
  return(gene_and_ko_targeted)
}


generate_small_family = function(df){
  return(df[, substrRight(paste(domain,phylum,class_name,order_name,family),40)])
}

## Read in auto-families
read_genome_tree_autotax = function(uid_and_bin_names){
  autotax = fread('/srv/projects/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/gtdb/tree_bacterial_set/overall_poster_tree08092016.tree.autotaxonomy',header=F)
  setnames(autotax, names(autotax), c('genome_uid','autotax'))
  autotax[, genome_uid := gsub(' ','_',genome_uid)]
  autotax = merge(autotax, uid_and_bin_names, by='genome_uid')
  setnames(autotax, 'autotax', 'taxonomy')
  autotax = split_taxonomy(autotax, has_root=F)
  autotax[,small_family := generate_small_family(autotax)]
  return(autotax)
}

read_a_gh = function(hmmtblout_retabbed){
  ## Use cutoffs as described at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt
  min_evalue = 1e-18
  min_hmm_coverage = 0.35
  dbcan_hits = fread(hmmtblout_retabbed,header=F)
  newnames = c('gene','target_accession','tlen','query','query_accession','qlen','sequence_Evalue','sequence_score','sequence_bias','num','of','cEvalue','iEvalue','domain_score','domain_bias','hmm_from','hmm_to','ali_from','ali_to','env_from','env_to','acc','description')
  setnames(dbcan_hits, names(dbcan_hits)[1:length(newnames)], newnames)
  dbcan_hits2 = dbcan_hits[iEvalue <= 1e-18 & (hmm_to - hmm_from)/qlen >= min_hmm_coverage]
  dbcan_hits2[, query := gsub('.hmm','',query)]
  dbcan_hits2[, genome := gsub('_\\d{5}$','',gene)]
  dbcan_hits2 = dbcan_hits2[grep('^GH',query,perl=T)]
  firsts = dbcan_hits2[order(iEvalue)][, .SD[1,], by='gene']
  setnames(firsts, 'query', 'gh')
  return(firsts[, .(genome, gene, gh, sequence_Evalue, description)])
}

## Read in GH families
read_gh = function(){
  return(read_a_gh('/srv/projects/abisko/annotation/16_dbCAN_aterrible12/dbCAN-fam-HMMs.txt.v5_vs_aterrible12_prokka3.hmmsearch.domtblout.retabbed.csv'))
}

read_backgroud_ghs = function(){
  return(read_a_gh('/srv/projects/abisko/annotation/25_dbcan_domtblout/dbCAN-fam-HMMs.txt.v5_vs_random1000.prokka.hmmsearch.domtblout.retabbed.csv'))
}

## Read in abundance info of clustered genomes
read_abundances = function(){
  abundances = fread('/srv/projects/abisko/mapping/47_aterrible12_dereplicated647/derep_genome_tpmean_filtered_coverages.csv',header=T)
  ## read genome sizes and get an average
  genome_sizes = fread('/srv/projects/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/genome_sizes.txt', header=F)
  average_genome_size = mean(genome_sizes$V1)
  ## read number of reads in each sample
  read_counts = fread('/srv/projects/abisko/data/flat20150213/flat20150213_read_counts.csv',header=T)
  setnames(read_counts, 'Number of reads R1 (line count/4)', 'num_reads')
  setnames(read_counts, 'Sample','sample')
  ## read number of mapped reads to the dereplicated set
  mapped_reads = fread('/srv/projects/abisko/mapping/47_aterrible12_dereplicated647/filtered0.75_0.95.read1_counts.csv',header=F)
  setnames(mapped_reads, names(mapped_reads), c('sample','num_mapped_reads'))
  ## relative abundance is (coverage divided by total coverage) * num_mapped_reads / total_reads
  mr = merge(mapped_reads, read_counts, by='sample')
  ab2 = merge(abundances, mr, by='sample')
  ab2[, relative_abundance := coverage/sum(coverage)*num_mapped_reads/num_reads, by='sample']
  #ab2[, genome:= gsub('_(\\d+)','.\\1',genome,perl=T)]
  #ab2[, genome:= gsub('^','73.',genome,perl=T)]
  ab2[, relative_abundance_of_recovered := coverage/sum(coverage), by='sample']
  return(ab2[, .(genome, sample, coverage, relative_abundance, relative_abundance_of_recovered)])
}

read_gpkg_assignments = function(){
  gpkg_in = fread('/srv/projects/abisko/annotation/21_gpkg_aterrible12_sweep/run1.read_tax.tsv',header=F)
  setnames(gpkg_in, names(gpkg_in), c('gpkg','protein','taxonomy'))
  gpkg_in = split_taxonomy(gpkg_in,has_root=T)
  gpkg_in[, unsplit_protein := gsub('_split.*','',x=protein)]
  gpkg_in[, genome := gsub('(.*)_.*','\\1',perl=T,x=protein)]
  gpkg_in[, genome := gsub('(.*)_(.....)_split.*','\\1',perl=T,x=genome)]
  setnames(gpkg_in,
           c('taxonomy','domain','phylum','class_name','order_name','family','genus','species'),
           c('gpkg.taxonomy','gpkg.domain','gpkg.phylum','gpkg.class_name','gpkg.order_name','gpkg.family','gpkg.genus','gpkg.species'))
  return(gpkg_in)
}

construct_genome_types = function(graftm_s5_on_aterrible12){
  g2 = split_taxonomy(graftm_s5_on_aterrible12, taxonomy_field='s5_taxonomy')
  g2[, small_family := generate_small_family(g2)]
  acidoflorens_genomes = g2[phylum=='P__gCandidatus_Koribacter.sCandidatus_Koribacter_versatilis.2']$genome
  methanoflorens_genomes = g2[small_family == 'ales_1.fMethanosaetaceae.gMethanosaeta.1']$genome
  ad3_target_family_genomes = g2[small_family=='ia_4.1 O__dBacteria_4.4 F__dBacteria_4.3']$genome
  ad3_family2_genomes = g2[family == 'F__dBacteria_4.8']$genome
  ad3_family3_genomes = g2[family == 'F__dBacteria_4.4']$genome
  ad3_family4_genomes = g2[domain == 'K__dBacteria_4.30']$genome
  genome_types = data.table(genome=methanoflorens_genomes, taxon='Methanoflorentaceae')
  genome_types = rbind(genome_types, data.table(genome=acidoflorens_genomes, taxon='Acidofloren'))
  genome_types = rbind(genome_types, data.table(genome=ad3_target_family_genomes, taxon='ad3_target_family'))
  genome_types = rbind(genome_types, data.table(genome=ad3_family2_genomes, taxon='ad3_family2'))
  genome_types = rbind(genome_types, data.table(genome=ad3_family3_genomes, taxon='ad3_family3'))
  genome_types = rbind(genome_types, data.table(genome=ad3_family4_genomes, taxon='ad3_family4'))
  return(genome_types)
}

read_prokka_annotations = function(){
  d = fread('/srv/projects/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/prokka3_kingdom_specific/prokka3.EC.csv',header=F)
  setnames(d, names(d), c('protein','ec','annotation'))
  d[ec == '-', ec := NA]
  d[, genome := gsub('(.*)_.*','\\1',perl=T,x=protein)]
  return(d)
}

## read_gh_to_ec_mappings = function() {
##   d = fread('/srv/projects/abisko/annotation/22_gh_cazy_munging/gh_to_ec_mapping.csv',header=F)
##   setnames(d, names(d), c('gh','ec'))
##   return(d)
## }

read_a_gh_best_hit_and_ec = function(gh_file){
  d = fread(gh_file,header=F)
  setnames(d, c('V1','V2'), c('protein','genbank_id'))
  genbank_and_gc = fread('/srv/projects/abisko/annotation/22_gh_cazy_munging/genbank_and_ec.manual.csv',header=F)
  setnames(genbank_and_gc, names(genbank_and_gc), c('ec','genbank_id'))
  m = merge(d[, .(genbank_id, protein)], genbank_and_gc, all.x=T, by='genbank_id', allow.cartesian=T)
  m[,genome := gsub('(.*)_.*','\\1',protein,perl=T)]
  return(data.table(m))
}

read_gh_best_hit_and_ecs = function(){
  return(read_a_gh_best_hit_and_ec('/srv/projects/abisko/annotation/22_gh_cazy_munging/gh_proteins_to_blast.faaVcharacterised.with_ec'))
}

read_background_gh_best_hit_and_ecs = function(){
  return(read_a_gh_best_hit_and_ec('/srv/projects/abisko/annotation/24_gtdb_background_genome_set/gh_proteins_to_blast.faaVcharacterised'))
}

read_ec_definitions = function(){
  d = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/ec_definitions.csv',header=T)
  d = d[, .SD[1], by=ec]
  return(d)
}

read_s5_read_abundances = function(){
  d = fread('/srv/projects/abisko/shotgun_abundance/103_contig_tree_relative_abundances/3_again/s5.abisko.csv',
            header=T)
  setnames(d,'#OTU_ID','otu_id')
  d[, otu_id := NULL]
  setnames(d, names(d), gsub(x=names(d), '\\.1$', '',perl=T))
  setnames(d, 'ConsensusLineage', 'taxonomy')
  d = split_taxonomy(d)
  d[, full_family := paste(domain, phylum, class_name, order_name, family)]
  d[, small_family := trim(substrRight(as.character(full_family), 40))]
  d[, taxonomy := NULL]
  d[, domain := NULL]
  d[, phylum:= NULL]
  d[, class_name := NULL]
  d[, order_name := NULL]
  d[, family := NULL]
  d[, genus := NULL]
  d[, species := NULL]
  d[, full_family := NULL]
  m = melt(d, id.vars='small_family', value.name='count', variable.name='sample')[count > 0]
  m[, relative_abundance := count/sum(count), by=sample]
  return(m[, .(sample, count, relative_abundance, small_family)])
}

read_complete_pathways = function(){
  d = fread('/srv/projects/abisko/annotation/23_kegg_via_hmms/pathways3.complete.csv',header=F)
  setnames(d, names(d), c('genome','pathway'))
  return(d)
}

read_background_complete_pathways = function(){
  d = fread('/srv/projects/abisko/annotation/24_gtdb_background_genome_set/ko_assignment/pathways3.complete.csv')
  setnames(d, names(d), c('genome','pathway'))
  return(d)
}

                                        # Commented out because it is behind the times.
## get_pathway_metabolites = function(){
##   d = fread('/srv/projects/abisko/annotation/23_kegg_via_hmms/pathways1.molecules.csv',header=T)
##   m = melt(d, id.vars='pathway', value.name='chemical')[,.(pathway,chemical)][chemical != '']
##   return(m)
## }

get_ec_pathway_metabolites = function(genome_and_ec){
  d = data.table(rbind(
    genome_and_ec[ec == '3.2.1.4', 'cellulose', by=genome], # endo-cellulase
    genome_and_ec[ec == '3.2.1.4', 'cellobiose', by=genome],
    genome_and_ec[ec == '3.2.1.8', 'xylan', by=genome], # endo-xylanase
    genome_and_ec[ec == '3.2.1.8', 'xylose', by=genome],
    genome_and_ec[ec == '3.2.1.20', 'cellulose', by=genome], # cellobiose phosphorylase
    genome_and_ec[ec == '3.2.1.20', 'glucose', by=genome],
    genome_and_ec[ec == '3.2.1.21', 'cellobiose', by=genome], # cellobiose phosphorylase
    genome_and_ec[ec == '3.2.1.21', 'glucose', by=genome],
    genome_and_ec[ec == '3.2.1.91', 'cellulose', by=genome], # exo-cellulase
    genome_and_ec[ec == '3.2.1.91', 'glucose', by=genome]
    ))
  setnames(d, names(d), c('genome','chemical'))
  return(d)
}

read_direct_taxonomy = function(uid_and_bin_names){
  d = fread('/srv/projects/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/gtdb/tree_bacterial_set/decorated.rerooted.decorated.tree.my_sisters.csv', header=F, sep="\t")
  setnames(d, names(d), c('genome_uid','taxonomy','sister_taxonomies'))
  d[, genome_uid := gsub(' ','_',genome_uid)]
  d[, sister_taxonomies := NULL]
  d2 = merge(d, uid_and_bin_names, by='genome_uid')
  d2[, genome_uid := NULL]
}

get_methanogenesis_metabolites = function(direct_taxonomy){
  d3 = split_taxonomy(direct_taxonomy, has_root = F)[phylum == 'p__Euryarchaeota']
  to_return = data.table(rbind(
    d3[, .(genome,chemical='methane')],
    d3[family == 'f__Methanosaetaceae', .(genome, chemical='acetate')],
    d3[family != 'f__Methanosaetaceae', .(genome, chemical='CO2')],
    d3[family != 'f__Methanosaetaceae', .(genome, chemical='H2')]
  ))
  return(to_return)
}

get_bulk_geochem = function(){
  geochem = fread('/srv/projects/abisko/shotgun_abundance/103_contig_tree_relative_abundances/DATA SUMMARY.csv',header=T, na.strings=c('na',''))
  caitlin = fread('/srv/projects/abisko/shotgun_abundance/103_contig_tree_relative_abundances/Caitlin_pmoABC_singleM_raw_20151027.csv',header=T, na.strings=c('','na'))
  setnames(caitlin, 'UID - genomics', 'sample')
  cait = caitlin[!is.na(Name_f) & !is.na(sample)]
  cait[, d13C_CH4 := as.numeric(d13C_CH4)]
  cait[, d13C_DIC := as.numeric(d13C_DIC)]
  cait[, DIC.mM := as.numeric(DIC.mM)]
  cait[, CH4.mM := as.numeric(CH4.mM)]
  cait[, alpha := (d13C_DIC+1000)/(d13C_CH4+1000)]
  cait[, pH_pw := as.numeric(pH_pw)]
  return(cait[,.(sample,d13C_CH4,d13C_DIC,DIC.mM,CH4.mM,alpha,pH_pw,WTD.cm_neg,sulfate.uM,nitrate.uM,acetate.uM,T_3cm.degC,T_13cm.degC)])
}

get_hyddb_classifications = function(){
  hyddb = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/hydrogenases.hyddb.csv',header=F)
  setnames(hyddb, names(hyddb), c('protein','hyddb_classification'))
  hyddb[, genome := gsub('_.....$','',protein)]
  return(hyddb)
}

read_acidoflorens_gene_list = function(){
  d = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/Acidoflorens_gene_revised_list_191216.csv',header=T)
  d[, original_order := 1:nrow(d)]
  return(d)
}

read_acidoflorens_diamond_for_paul = function(){
  d = fread('/srv/projects/abisko/annotation/20_Acidobacteria_clustering/paul_diamond_shortcut/diamond_results3.csv',header=F)
  setnames(d, c('V1','V2','V11'),c('query','subject','evalue'))
  return(d)
}

read_ad3_diamond_for_paul = function(){
  d = fread('/srv/projects/abisko/annotation/28_ad3/diamond_results.csv',header=F)
  setnames(d, c('V1','V2','V11'),c('query','subject','evalue'))
  return(d)
}

read_ad3_gene_list = function(){
  d = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/AD3_gene_list_121216.csv',header=T)
  setnames(d,'Locus Tag','protein')
  setnames(d,'Name','paul_gene_name')
  d[, protein := gsub('73.20111000_P2D.25', '73.20111000_P2D.5', protein)]
  d[, original_order := 1:nrow(d)]
  return(d)
}

read_joel_pathways = function(){
  d = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/network_analyser_output_annotations.csv',header=T)
  setnames(d, 'Genome_name', 'genome')
  setnames(d, 'Module_id', 'module_id')
  setnames(d, 'Module_name', 'module_name')
  setnames(d, 'Steps_needed', 'num_steps')
  return(d[,.(genome,module_id,module_name,num_steps)])
}

read_ko_annotation = function(){
  d = fread('zcat /srv/projects/abisko/annotation/23_kegg_via_hmms/run1/3columns.best_bitscore.csv.gz',header=T)
  d[, genome := gsub('_.....$','',protein)]
  return(d)
}

read_ko_annotation_background = function(){
  d = fread('zcat /srv/projects/abisko/annotation/24_gtdb_background_genome_set/ko_assignment/run1/3columns.best_bitscore.csv.gz',header=T)
  d[, genome := gsub('_.....$','',protein)]
  return(d)
}

read_dereplicated647 = function(){
  d = fread('/srv/projects/abisko/aterrible_bins/14_dereplicated647/647_genomes_derep_97_70_for_paral.txt',header=F)
  setnames(d, 'V1','genome')
  return(d)
}

read_xylose_oxidoreductase_genomes = function(){
  d = fread('/srv/projects/abisko/Joel/99_phd/01_Projects/41-Population_genome_paper/experiments/03-Xylose_in_popn_genomes/oxidoreductase_pathway.txt',header=F)
  setnames(d,'V1','genome')
  return(d)
}

read_genome_taxonomy = function(){
  d = fread('/srv/projects/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/gtdb/tree_bacterial_set/bin_numbers_and_conservative_taxonomy.non_uniq_taxons.taxonomy',header=F)
  setnames(d, names(d), c('fasta','taxonomy'))
  d[, genome := gsub('.fa$',x=fasta,'',perl=T)]
  d[, genome := gsub('.fasta.metabat-bins-','',genome)]
  d2 = split_taxonomy(d, has_root=F)
  return(d2[,.(genome,domain,phylum,class_name,order_name,family,genus,species)])
}

list_acetoclastic_methanogens = function(genome_taxonomy){
  ## this list is from Joel, then just returning those that are in Methanosarcinales
  possible_acetoclasts = data.table(genome= c('73.20110700_S2D.54',
                                              '73.20120600_E3D.137',
                                              '73.20120600_E3D.31',
                                              '73.20110800_E1M.18',
                                              '73.20110600_S3M.49',
                                              '73.20120800_E3X.20',
                                              '73.20110800_S2D.59',
                                              '73.20100900_E3D.38',
                                              '73.20120700_S1D.117'))
  return(merge(possible_acetoclasts, genome_taxonomy, by='genome')[order_name=='o__Methanosarcinales'][,.(genome,pathway='acetoclastic_methanogenesis')])
}

list_methanotrophs = function(){
  methanotrophs = data.table(fasta = c('73.20120600_E3D.fasta.metabat-bins-.145.fa',
                                       '73.20120600_S3S.fasta.metabat-bins-.5.fa',
                                       '73.20120700_E2X.fasta.metabat-bins-.256.fa',
                                       '73.20120800_S3M.fasta.metabat-bins-.30.fa',
                                       '73.20120700_E3M.fasta.metabat-bins-.239.fa',
                                       '73.20110600_S1M.fasta.metabat-bins-.35.fa',
                                       '73.20120600_S2M.fasta.metabat-bins-.29.fa',
                                       '73.20120800_S3D.fasta.metabat-bins-.23.fa',
                                       '73.20111000_S3D.fasta.metabat-bins-.29.fa',
                                       '73.20120600_S3D.fasta.metabat-bins-.22.fa',
                                       '73.20111000_S2M.fasta.metabat-bins-.28.fa',
                                       '73.20120700_P3D.fasta.metabat-bins-.113.fa'))
  methanotrophs[, genome := gsub('.fa$',x=fasta,'',perl=T)]
  methanotrophs[, genome := gsub('.fasta.metabat-bins-','',genome)]
  return(methanotrophs[,.(genome,pathway='methanotrophy')])
}

read_sugar_names = function(){
  sugar_names = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/sugar_breakdown.csv',header=F)
  setnames(sugar_names, names(sugar_names), c('sugar','gene_name','ko'))
}

read_more_modules = function(){
  dir = '/srv/projects/abisko/Joel/99_phd/01_Projects/41-Population_genome_paper/experiments/02-KO_enrichment/more_modules'
  d = fread(paste(dir,'allM.txt',sep='/'),header=F)
  setnames(d, c('V1','V3'),c('genome','pathway'))
  return(d[,.(genome,pathway)])
}

read_more_modules_background = function(){
  dir = '/srv/projects/abisko/Joel/99_phd/01_Projects/41-Population_genome_paper/experiments/02-KO_enrichment/more_modules'
  d = fread(paste(dir,'allM.background.txt',sep='/'),header=F)
  setnames(d, c('V1','V3'),c('genome','pathway'))
  return(d[,.(genome,pathway)])
}

read_acidoflorens_families = function(){
  d = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/acidoflorens_families.csv',header=T)
  setnames(d, 'acidoflorens_family','acidoflorens_clade')
  return(d)
}
read_methanoflorens_families = function(){
  return(fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/methanoflorens_clades.csv',header=T))
}

read_phylogeny2 = function(){
  d = fread('/srv/projects/abisko/Caitlin/006_population_genomes/03_final_trees/1529_pop_genomes_tax_Ben_final_proteo_split_correct.csv',
            header=F, sep='\t')
  setnames(d, c('V1','V2','V3'), c('uid','phylum','genome'))
  d[,genome := gsub('.fa$','',genome)]
  d[, phylum_no_p := gsub('p__','',phylum)]
  d[phylum_no_p == 'Deltaproteobacteria', phylum_no_p := 'Proteobacteria']
  return(d[,.(genome,phylum=phylum_no_p)])
}

read_tax2tree_phylogeny2 = function(uid_and_bin_names){
  d = fread('/srv/projects/abisko/Caitlin/006_population_genomes/03_final_trees/abisko_bacteria_tree/gtdb.decorated.tree-consensus-strings',header=F,sep='\t')
  setnames(d, names(d), c('genome_uid','taxonomy'))
  translations = fread('/srv/projects/abisko/Caitlin/006_population_genomes/03_final_trees/abisko_bacteria_tree/01_genomes_to_tax/bacterial_genomes_taxonomy.tsv', header=F)
  setnames(translations, names(translations), c('genome_uid','phil_phylum','caitlin_phylum','genome'))
  translations[caitlin_phylum == 'Deltaproteobacteria', caitlin_phylum := 'Proteobacteria']
  translations[, genome := NULL]
  bacteria = merge(d, uid_and_bin_names, by='genome_uid')
  b2 = merge(bacteria, translations, by='genome_uid')
  d = fread('/srv/projects/abisko/Caitlin/006_population_genomes/03_final_trees/abisko_archaea_tree/gtdb_archaea_decorated.tree-consensus-strings', header=F,sep='\t')
  setnames(d, names(d), c('genome_uid','taxonomy'))
  translations = fread('/srv/projects/abisko/Caitlin/006_population_genomes/03_final_trees/abisko_archaea_tree/01_genomes_to_tax/archaea_95_taxonomy.tsv', header=F, sep='\t')
  translations[, V5 := NULL]
  setnames(translations, names(translations), c('genome_uid','phil_phylum','caitlin_phylum','genome'))
  translations[, genome := NULL]
  archaea = merge(d, uid_and_bin_names, by='genome_uid')
  a2 = merge(archaea, translations, by='genome_uid')
  return(rbind(b2,a2)[,.(genome,phil_phylum,caitlin_phylum,taxonomy)])
}

read_graftm_16S = function(){
  d = fread('/srv/projects/abisko/shotgun_abundance/100_graftm_gg97_flat20150213/graftm.gg97.otu_table.csv',header=T)
  setnames(d, names(d), gsub(x=names(d), '.1$', ''))
  d2 = melt(d, id.vars=c('#OTU_ID','ConsensusLineage'), variable.name='sample', value.name='count')[count > 0]
  setnames(d2, 'ConsensusLineage', 'taxonomy')
  d3 = split_taxonomy(d2)
  d3[, taxonomy := NULL]
  return(d3)
}

read_singlem_rarefield = function() {
  d = fread('/srv/projects/abisko/shotgun_abundance/94_beta_diversity_ok_singlem_on_abisko/pplacer1.rarefied100.otu_table.csv', header = T)
  ## exclude samples that do not have all genes with 100 seqs
  d2 = d[sample %in% d[, sum(num_hits), by='sample'][V1==1500]$sample]
  d2[, sample := gsub('.1$','',perl=T,sample)]
  return(d2)
}

read_networkanalyzer_pathways = function(){
  d = fread('/srv/projects/abisko/annotation/30_networkanalyzer/my_annotations.tsv',header=F)
  setnames(d, c('V1','V2','V6'), c('genome','pathway','completeness'))
  return(d[completeness==100, .(genome, pathway)])
}

read_all_networkanalyzer_pathways = function(){
  d = fread('/srv/projects/abisko/annotation/30_networkanalyzer/networkanalyzer_out_abisko_check4_annotations.tsv',header=T)
  setnames(d, 'Module_id', 'module')
  return(d[Percent_Steps_found==100, .(genome=Genome_name, module, pathway=Module_name)])
}

read_all_networkanalyzer_pathways_background = function(){
  d = fread('/srv/projects/abisko/annotation/30_networkanalyzer/networkanalyzer_out4_annotations.tsv',header=T)
  setnames(d, 'Module_id', 'module')
  return(d[Percent_Steps_found==100, .(genome=Genome_name, module, pathway=Module_name)])
}

read_actino_depth_profile_aai = function(){
  d = fread('/srv/projects/abisko/annotation/32_actino_aai/comparem_out/aai/aai_summary.tsv',header=T)
  setnames(d, c('Mean AAI','Orthologous fraction (OF)'), c('aai','orthologous_fraction'))
  setnames(d, c('Genome A','Genome B'), c('genome1','genome2'))
  return(d)
}

read_brite_lists = function(){
  t1 = data.table(group=gsub('.txt','',dir('/srv/projects/abisko/annotation/32_kegg_pathway_lists/lists')))
  total = t1[, .(module=fread(paste('/srv/projects/abisko/annotation/32_kegg_pathway_lists/lists/',group,".txt",sep=''),header=F)$V1), by=group]
  return(total)
}

read_proteobacteria_specific_taxonomy = function(){
  d = fread('/srv/projects/abisko/Caitlin/006_population_genomes/03_final_trees/1529_pop_genomes_tax_Ben_final_proteo_split.tsv')
  setnames(d, names(d), c('uid','phylum','genome'))
  d[, genome := gsub('.fa$','',genome)]
  return(d)
}

read_acidoflorens_aai = function(){
  d = fread('/srv/projects/abisko/annotation/26_acidoflorens/comparem_aai_output/aai/aai_summary.tsv',header=T)
  setnames(d, names(d), c('genome1','genes_in_1','genome2','genes_in_2','num_orthologous','aai','aai_sd','orthologous_fraction'))
  d[, genome_pair_index := paste(sort(c(genome1,genome2)),collapse=' '),by='genome1,genome2']
  return(d)
}

read_mash = function(){
  d = fread('/srv/projects/abisko/aterrible_bins/12_assembly73_individuals_flat20160115/mash/aterrible12_vs_aterrible12.700_or_more.csv',header=F)
  setnames(d, names(d), c('genome1','genome2','pvalue','thousand'))
  d[, genome1:= gsub('.fna','',genome1)]
  d[, genome2:= gsub('.fna','',genome2)]
  d[, ani := 1-pvalue]
  return(d[,.(genome1,genome2,ani)])
}

list_methanomassiliicoccales_genomes = function(){
  # read manually from the Caitlin's arb tree
  return(data.table(genome=c('73.20110800_E3D.115','73.20120700_E3M.288')))
}

read_methanoflorens_aai = function(){
  d = fread('/srv/projects/abisko/annotation/35_methanoflorens_aai/comparem_out/aai/aai_summary.tsv',header=T)
  setnames(d, names(d), c('genome1','genes_in_1','genome2','genes_in_2','num_orthologous','aai','aai_sd','orthologous_fraction'))
  d[, genome_pair_index := paste(sort(c(genome1,genome2)),collapse=' '),by='genome1,genome2']
  return(d)
}

read_virus_host_linkages = function(){
  d = fread('/srv/projects/abisko/annotation/37_virus_linkages/Supplementary_Tables_AVMM_v12.csv',header=T)
  setnames(d, names(d), c('viral_contig','viral_contig_length','genome','host_taxonomy','evidence','additional_evidence'))
  d[, genome := gsub('.fasta.metabat-bins-','',genome)]
  d2 = d[,.(viral_contig, genome)][,unlist(strsplit(gsub('Ties-','',genome),';')),by='genome,viral_contig']
  d2[, genome := NULL]
  setnames(d2, 'V1','genome')
  return(d2)
}

read_genome_publication_names = function(){
  d = fread('/srv/projects/abisko/Caitlin/006_population_genomes/03_final_trees/02_newicks_for_publication/working/publication_genome_names.tsv',header=F)
  setnames(d, names(d), c('uid','genome','genome_publication_name'))
  d[, genome := gsub('.fa$','',genome)]
  return(d[,.(genome,genome_publication_name)])
}

generate_metabolism_names_and_sources = function(){
  funcs = rbind(
    data.table(type = 'networkanalyzer_pathways_of_interest',
               name = c(
                 'galacturonic_acid_deg_bacteria',
                 'galacturonic_acid_deg_fungi',
                 'galactose_deg_leloir_pathway',
                 'lactose_deg',
                 'sucrose_deg',
                 'trehalose_deg',
                 'fructose_degradation_via_1P_to_fructose-16P',
                 'fructose_degradation_via_1P_to_glyceraldehyde',
                 'canon_fructose_deg',
                 'mannose_deg',
                 'fucose_deg',
                 'glycolysis'
               )),
    data.table(type = 'all_networkanalyzer_pathways_of_interest',
               name = c('Pectin degradation',
                        'D-Galacturonate degradation (fungi), D-galacturonate => glycerol',
                        'D-Galacturonate degradation (bacteria), D-galacturonate => pyruvate + D-glyceraldehyde 3P',
                        'Methanogenesis, methanol => methane')),
    #data.table(type = 'all_networkanalyzer_pathways_of_interest',
    #           name = 'Methanogenesis, methanol => methane'),
    data.table(type = 'all_networkanalyzer_pathways_of_interest_module_names',
               name = c('pyruvate_to_ethanol',
                        'pyruvate_to_acetate',
                        'pyruvate_to_butanoate',
                        'pyruvate_to_1butanol',
                        'pyruvate_to_acetone',
                        'pyruvate_to_2propanol',
                        'pyruvate_to_lactate',
                        'pyruvate_to_acrylate',
                        'pyruvate_to_propanoate',
                        'pyruvate_to_co2',
                        'reduced_ferredoxin_to_h2',
                        'ferredoxin_from_nadph',
                        'ferredoxin_and_nad_reduction',
                        #'acetoclastic_methanogenesis',
                        #'hydrogenotrophic_methanogenesis',
                        #'methanol_methanogenesis',
                        'amine_methanogenesis',
                        'canonical_xylose_degradation',
                        'fungal_xylose_to_xylulose_5phosphate',
                        'fungal_xylose_to_xylulose_5phosphate_xr_xylulose_reductase',
                        'xylonate_hydratase_pathway',
                        'xylonate_hydratase_pathway_weimberg')),
    data.table(type='complete_pathways_by_ko',
               name='glycolysis'),
    data.table(type='paul_annotation',
               name=c('acetoclastic_methanogenesis',
                      'hydrogenotrophic_methanogenesis',
                      'methylmethanogenesis',
                      'methanol_methanogenesis'))
  )

  return(funcs)
}

generate_metabolic_pathways = function(genome_types,
                                       protein_and_gh_ec, sugar_names, ko_annotation, more_pathways, networkanalyzer_pathways,
                                       all_networkanalyzer_pathways, methanotroph_functions,
                                       xylan_ecs, cellulose_ecs, beta_glucosidase_ecs, complete_pathways_by_ko,
                                       pauls_methanogen_functions,
                                       metabolism_names_and_sources){
  ##gp = complete_pathways_by_ko[pathway %in% pathways_of_interest]
  ## manually Add Methanoflorens genomes
  ##gp2 = rbind(gp, data.table(genome=genome_types[taxon=='Methanoflorentaceae',genome],pathway='hydrogenotrophic_methanogenesis'))
  gp2 = pauls_methanogen_functions
  ## gp2 = data.table(genome=genome_types[taxon=='Methanoflorentaceae',genome],pathway='hydrogenotrophic_methanogenesis')
  ## Add cellulose and xylose degradation
  gp2a = rbind(gp2,
               rbind(data.table(genome=protein_and_gh_ec[ec %in% xylan_ecs, unique(genome)], pathway='xylan_degradation'),
                     data.table(genome=protein_and_gh_ec[ec %in% cellulose_ecs, unique(genome)], pathway='cellulose_degradation'),
                     data.table(genome=protein_and_gh_ec[ec %in% beta_glucosidase_ecs, unique(genome)], pathway='beta_glucosidase')
                     ))
  ##gp2b = rbind(gp2a, acetoclastic_methanogen_functions)
  gp2c = rbind(gp2a, merge(sugar_names, ko_annotation, by='ko')[,.(pathway=paste(sugar,gene_name)),by='genome,gene_name'][,.SD[1],by='genome,pathway'][,.(genome,pathway)])
  gp2d = rbind(gp2c, more_pathways)
  gp2e = rbind(gp2d, networkanalyzer_pathways[pathway %in% metabolism_names_and_sources[type == 'networkanalyzer_pathways_of_interest', name]])
  gp2f = rbind(gp2e, all_networkanalyzer_pathways[pathway %in% metabolism_names_and_sources[type == 'all_networkanalyzer_pathways_of_interest', name], .(genome,pathway)])
  gp2g = rbind(gp2f, all_networkanalyzer_pathways[module %in% metabolism_names_and_sources[type == 'all_networkanalyzer_pathways_of_interest_module_names', name], .(genome,pathway=module)])
  gp2h = rbind(gp2g, complete_pathways_by_ko[pathway %in% metabolism_names_and_sources[type == 'complete_pathways_by_ko'], .(genome,pathway)])
  gp3 = rbind(gp2h, methanotroph_functions)
  return(gp3)
}

generate_background_pathway_prevalence = function(protein_and_gh_background_ec,
                                                  xylan_ecs, cellulose_ecs,
                                                  beta_glucosidase_ecs,
                                                  complete_pathways_by_ko_background,
                                                  metabolism_names_and_sources){
  backgrounds_from_ec = rbind(
    data.table(genome=protein_and_gh_background_ec[ec %in% xylan_ecs, unique(genome)], pathway='xylan_degradation'),
    data.table(genome=protein_and_gh_background_ec[ec %in% cellulose_ecs, unique(genome)], pathway='cellulose_degradation'),
    data.table(genome=protein_and_gh_background_ec[ec %in% beta_glucosidase_ecs, unique(genome)], pathway='beta_glucosidase'))
  background_pathway_prevalence = rbind(
    ##ko_annotation_background[ko==xylonate_dehydratase_ko,.(genome,pathway='xylonate_dehydratase_pathway')],
    ##complete_pathways_by_ko_background[pathway %in% pathways_of_interest],
    all_networkanalyzer_pathways_background[module %in% metabolism_names_and_sources[type == 'all_networkanalyzer_pathways_of_interest_module_names', name], .(genome,pathway=module)],
    all_networkanalyzer_pathways_background[pathway %in% metabolism_names_and_sources[type == 'all_networkanalyzer_pathways_of_interest', name], .(genome,pathway)],
    complete_pathways_by_ko_background[pathway %in% metabolism_names_and_sources[type == 'complete_pathways_by_ko'], .(genome,pathway)],
    backgrounds_from_ec
  )
  return(background_pathway_prevalence)
}

generate_used_pathways_and_better_names = function(){
  return(data.table(
    pathway = c("pyruvate_to_acetate",
                "canonical_xylose_degradation",
                "pyruvate_to_lactate",
                "fungal_xylose_to_xylulose_5phosphate",
                "pyruvate_to_ethanol",
                "xylonate_hydratase_pathway",
                "hydrogenotrophic_methanogenesis",
                "acetoclastic_methanogenesis",
                "pyruvate_to_2propanol",
                "fungal_xylose_to_xylulose_5phosphate_xr_xylulose_reductase",
                "pyruvate_to_propanoate",
                "methanol_methanogenesis",
                "methylmethanogenesis",
                "xylan_degradation",
                "cellulose_degradation",
                "methanotrophy",
                'canon_fructose_deg',
                'fucose_deg',
                'galactose_deg_leloir_pathway',
                'galacturonic_acid_deg_bacteria',
                'glycolysis',
                'lactose_deg',
                'mannose_deg'
                ),
    pathway_publication_name = c("acetate fermentation",
                                 "xylose degradation (isomerase pathway)",
                                 "lactate fermentation",
                                 "xylose degradation (previously-fungal pathway)",
                                 "ethanol fermentation",
                                 "xylose degradation (xylanate hydratase pathway)",
                                 "hydrogenotrophic methanogenesis",
                                 "acetoclastic methanogenesis",
                                 "2propanol fermentation",
                                 "xylose degradation (xylose reductase and xylulose reductase)",
                                 "propanoate fermentation",
                                 "methanol methanogenesis",
                                 "methylmethanogenesis",
                                 "xylan degradation",
                                 "cellulose degradation",
                                 "methanotrophy",
                                 'fructose degradation',
                                 'fucose degradation',
                                 'galactose degradation',
                                 'galacturonic acid degradation',
                                 'glycolysis',
                                 'lactose degradation',
                                 'mannose degradation'
                                 )))
}

read_dereplicated_set_aai = function(){
  d = fread('/srv/projects/abisko/annotation/40_aai_of_mapped_genomes/aai.at_least_50.tsv')
  setnames(d, names(d), c('genome1','genome2','aai'))
  d$relationship = 'higher'
  d[aai > 79, relationship := 'genus']
  d[aai > 95, relationship := 'species']
  return(d)
}

read_cton = function(){
  d = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/cton.tsv',na.strings = 'na')[!is.na(CtoN_wt)]
  return(d)
}

read_species_of_interest_list = function(){
  d = fread('/srv/projects/abisko/Caitlin/006_population_genomes/03_final_trees/03_chang_acid_mcrillii/chang_acido_mcrill.tsv', header=F)
  setnames(d, names(d), c('sublineage','genome'))
  d[, genome := gsub('.fa','',genome)]
  d$lineage = 'fail'
  d[grep('lineage',sublineage), lineage := 'Acidoflorens']
  d[grep('lineage',sublineage), lineage := 'Acidoflorens']
  d$lineage
  d[grep('Changshen',sublineage), lineage := 'Changshengia']
  d[grep('Mstor',sublineage), lineage := 'Methanoflorens']
  d[grep('Mcrillii',sublineage), lineage := 'Methanoflorens']
  d[grep('Acidobacteriales',sublineage), lineage := 'Koribacteriaceae?']
  d[lineage == 'fail']
  return(d)
}

read_ncbi_biosamples = function(){
  d = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/ncbi_biosamples.tsv.csv',header=T)
  d = d[,c(1:15),with=F]
  setnames(d, 'sample_name','sample')
  return(d)
}

read_sample_base_counts = function(){
  d = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/sample_base_counts.csv')
  return(d)
}

read_transcript_tpm = function(){
  d = fread('/srv/projects/abisko/Joel/99_phd/01_Projects/41-Population_genome_paper/experiments/25-TPM_directionality/10-Redoing_TPM_no_RNA/tpm_measurements.one_sided.tsv', header=T)
  setnames(d, 'TPM','tpm')
  setnames(d,'gene_id','gene')
  transcriptome_samples = unique(d[!(
    sample %in% c(
                  '20120700_S1S',# - too shallow
                  '20120800_E2S',# - too shallow
                  '20110600_E1D',# - dna contam to high
                  '20110600_E1M'# - dna contam too high
                  ))]$sample) #= 26 samples
  d2 = d[sample %in% transcriptome_samples]
  return(d2)
}

read_genome_mcras = function(){
  d = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/prokka3_mcrA_graftm/prokka3_kingdom_specific.fixed/prokka3_kingdom_specific.fixed_read_tax.tsv',header=F)
  setnames(d, 'V1', 'gene')
  d[, V2 := NULL] # These annotations aren't very specific
  d[, genome := gsub('_.....$','',gene)]
  return(d)
}

read_all_pathway_kos = function(){
  d = fread('/srv/projects/abisko/annotation/43_methanogenesis_comparison_tpm/module_to_description.bare.tsv',header=F)
  setnames(d, names(d), c('module','kos'))
  d2 = d[, strsplit(kos,' '), by=module]
  setnames(d2, 'V1','ko')
  return(d2)
}

read_new_archaea_taxonomy = function(uid_and_bin_names){
  d = fread('/srv/projects/abisko/Caitlin/006_population_genomes/03_final_trees/00_tree_redo_20170911/02_abisko_archaea_tree/gtdb_archaea_decorated.tree-consensus-strings',header=F,sep="\t")
  setnames(d, names(d), c('genome_uid','taxonomy'))
  d2 = merge(d, uid_and_bin_names, by='genome_uid')
  return(d2[,.(genome,taxonomy)])
}

read_pathway_expression = function(){
  d = fread('/srv/projects/abisko/Joel/99_phd/01_Projects/41-Population_genome_paper/experiments/25-TPM_directionality/08-Redrawing_plots/pathway_expression.no_rna.tsv',header=T)
  return(d[,.(sample,tpm=TPM,pathway)])
}

read_pathway_expression_per_genome = function(){
  d = fread('/srv/projects/abisko/Joel/99_phd/01_Projects/41-Population_genome_paper/experiments/25-TPM_directionality/08-Redrawing_plots/pathway_expression_per_genome.no_rna.tsv',header=T)
  setnames(d,names(d),c('sample','genome','tpm','pathway'))
  d[, sample:=gsub('\\+AF8-','_',sample)]
  d[, genome:=gsub('\\+AF8-','_',genome)]
  d[, tpm:=as.numeric(gsub('\\+AF8-','_',tpm))]
  d[, pathway:=gsub('\\+AF8-','_',pathway)]
  return(d)
}

read_pauls_methanogen_functions = function(genome_dereplication){
  d = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/paul_methanogen_classification.csv',header=T)
  d2 = melt(d,id.vars = 'genome', measure.vars = c('function1','function2','function3'))
  d3 = merge(d2, genome_dereplication[,.(un=genome,representative)], by.x='genome', by.y='representative') # Paul only annotated the representatives
  return(d3[value != '',.(genome=un,pathway=value)])
}

read_proteomics11 = function(){
  d= fread('/srv/projects/abisko/proteomics/11_genomes_paper_representatives_subset/expressed_list.csv',header=T)
  setnames(d, 'Genome', 'protein')
}

read_proteomics_variants = function(){
  d = fread('tail -n+2 /srv/projects/abisko/proteomics/10_genomes_paper/Variants_mapped_to_the_representative.csv |awk \'$3>1\' |sed \'s/;.*\t/\t/; s/;.*$//\' |cut -f4-5', header = F)
  setnames(d, c('protein','proteomic_representative'))
  return(d)
}

read_all_networkanalyzer_pathways_kos = function(){
  d = fread('/srv/projects/abisko/annotation/30_networkanalyzer/pathways2.tsv',header=F)
  d2 = d[, .(ko=unique(grep('.',unlist(tstrsplit(d$V2[1], '[\\(\\), ]')),value=T))), by=V1]
  setnames(d2,'V1','pathway')
  return(d2)
}

### Actual code
## Read in raw data
uid_and_bin_names = read_uid_and_bin_names()
phylum_wise_taxonomy = read_phylum_wise_taxonomy(uid_and_bin_names)
graftm_s5_on_aterrible12 = read_graftm_s5_on_aterrible12()
genome_dereplication = read_dereplicated_genomes_clusters(data.table(genome=phylum_wise_taxonomy[, genome]))
targeted_gene_and_ko = read_targeted_gene_and_ko()
#genome_tree_autotax = read_genome_tree_autotax(uid_and_bin_names)
#gh = read_gh()
#background_gh = read_backgroud_ghs()
mapping_abundances = read_abundances()
gpkg_assignments = read_gpkg_assignments()
genome_types = construct_genome_types(graftm_s5_on_aterrible12) #need to fix family names in function.
gpkg_functions = fread('/srv/projects/abisko/annotation/18_aterrible12_prokka_to_KO/gpkg_functions.csv',header=T)
protein_annotations = read_prokka_annotations()
##gh_to_ec = read_gh_to_ec_mappings()
protein_and_gh_ec = read_gh_best_hit_and_ecs()
protein_and_gh_background_ec = read_background_gh_best_hit_and_ecs()
ec_definitions = read_ec_definitions()
s5_read_abundances = read_s5_read_abundances()
complete_pathways_by_ko = read_complete_pathways()
complete_pathways_by_ko_background = read_background_complete_pathways()
## pathway_metabolites = get_pathway_metabolites()[pathway != 'pyruvate_to_co2'] #don't include this at the moment since there's too much
## pathway_metabolites = pathway_metabolites[chemical != 'CO2'] # exclude CO2 from fermentation
metabolites_from_ec = get_ec_pathway_metabolites(protein_and_gh_ec)
direct_taxonomy = read_direct_taxonomy(uid_and_bin_names)
methanogenesis_metabolites = get_methanogenesis_metabolites(direct_taxonomy)
geochem = get_bulk_geochem()
hyddb = get_hyddb_classifications()
acidoflorens_gene_list = read_acidoflorens_gene_list()
acidoflorens_diamond = read_acidoflorens_diamond_for_paul()
ad3_diamond = read_ad3_diamond_for_paul()
ad3_gene_list = read_ad3_gene_list()
joel_pathways = read_joel_pathways()
ko_annotation = read_ko_annotation()
ko_annotation_background = read_ko_annotation_background()
representatives647 = read_dereplicated647()
xylose_oxidoreductase_genomes = read_xylose_oxidoreductase_genomes()
genome_taxonomy = read_genome_taxonomy()
##acetoclastic_methanogen_functions = list_acetoclastic_methanogens(genome_taxonomy) # Paul's effort is superior
methanotroph_functions = list_methanotrophs()
site_colours = c('#703C1B', # palsa
                 '#058000', # sphag
                 '#0001FF')  # erio
sugar_names = read_sugar_names()
more_pathways = read_more_modules()
more_pathways_background = read_more_modules()
acidoflorens_families = read_acidoflorens_families()
methanoflorens_families = read_methanoflorens_families()
phylogeny2 = read_phylogeny2()
phylogeny2_proteos_together = phylogeny2
phylogeny2_proteos_together[grep('roteobacteria',phylum), phylum:='Proteobacteria']
graftm16S = read_graftm_16S()
singlem_rarefied = read_singlem_rarefield()
networkanalyzer_pathways = read_networkanalyzer_pathways()
phylogeny2_full = read_tax2tree_phylogeny2(uid_and_bin_names)
actino_depth_profile_aai = read_actino_depth_profile_aai()
brite_lists = read_brite_lists()
all_networkanalyzer_pathways = read_all_networkanalyzer_pathways()
all_networkanalyzer_pathways_background = read_all_networkanalyzer_pathways_background()
proteobacteria_taxonomy = read_proteobacteria_specific_taxonomy()
acidoflorens_aai = read_acidoflorens_aai()
mash = read_mash()
mmass_genomes = list_methanomassiliicoccales_genomes()
methanoflorens_aai = read_methanoflorens_aai()
viral_host_linkages = read_virus_host_linkages()
genome_publication_names = read_genome_publication_names()
xylan_ecs = c('3.2.1.8','3.2.1.151') #'3.2.1.20','3.2.1.21',
cellulose_ecs = c('3.2.1.4','3.2.1.74')
beta_glucosidase_ecs = c('3.2.1.21') # directly breaking glycosidic bond
pauls_methanogen_functions = read_pauls_methanogen_functions(genome_dereplication)
metabolism_names_and_sources = generate_metabolism_names_and_sources()

metabolism_pathway_annotations = generate_metabolic_pathways(
  genome_types, #acetoclastic_methanogen_functions,
  protein_and_gh_ec, sugar_names, ko_annotation, more_pathways, networkanalyzer_pathways,
  all_networkanalyzer_pathways, methanotroph_functions,
  xylan_ecs, cellulose_ecs, beta_glucosidase_ecs, complete_pathways_by_ko,
  pauls_methanogen_functions,
  metabolism_names_and_sources)
background_pathway_prevalences = generate_background_pathway_prevalence(
  protein_and_gh_background_ec,
  xylan_ecs, cellulose_ecs,
  beta_glucosidase_ecs,
  complete_pathways_by_ko_background,
  metabolism_names_and_sources)
used_paths_and_publication_names = generate_used_pathways_and_better_names()
used_paths = used_paths_and_publication_names$pathway
dereplicated_set_aai = read_dereplicated_set_aai()
cton = read_cton()
interesting_species_list = read_species_of_interest_list()
ncbi_biosamples = read_ncbi_biosamples()
sample_base_counts = read_sample_base_counts()
mcras = read_genome_mcras()
all_pathway_kos = read_all_pathway_kos()
new_archaea_taxonomy = read_new_archaea_taxonomy(uid_and_bin_names)
## Transcriptomics
tpm = read_transcript_tpm()
pathway_expression_per_genome = read_pathway_expression_per_genome()
pathway_expression = read_pathway_expression()
proteomics = read_proteomics11()
proteomic_variants = read_proteomics_variants()
all_networkanalyzer_pathways_kos = read_all_networkanalyzer_pathways_kos()
