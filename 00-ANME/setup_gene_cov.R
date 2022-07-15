library(data.table)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# All samples at once
filt_gene_cov <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/dirseq_out/filtrate_gene_cov.tsv")
retn_gene_cov <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/dirseq_out/retentate_gene_cov.tsv")
filt_gene_cov[, length:=end-start]
retn_gene_cov[, length:=end-start]

# Fitler out low coverage misassembly region
filt_gene_cov <- filt_gene_cov[, mean_cov:=(forward_average_coverage+reverse_average_coverage)/2][mean_cov>15,]
retn_gene_cov[, mean_cov:=(forward_average_coverage+reverse_average_coverage)/2]

# Normalize by the mean coverage. Can't be bothered to do CLR here, doesn't seem necessary yet
filt_gene_cov[, norm_cov:=(mean_cov/mean(filt_gene_cov$mean_cov))]
retn_gene_cov[, norm_cov:=(mean_cov/mean(retn_gene_cov$mean_cov))]

filt_gene_cov[, sample:="filtrate"]
retn_gene_cov[, sample:="retentate"]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Per sample gene coverages. Only associated fractions samples with their associated fraction MAG
#All samples are available in folder
filt_gene_cov_1 <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/dirseq_out/filtrate_anme_removed_chimeras.fna.SB8774_S41_R1_001.fastq.gz_gene_cov.tsv")[, mean_cov:=(forward_average_coverage+reverse_average_coverage)/2]
filt_gene_cov_2 <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/dirseq_out/filtrate_anme_removed_chimeras.fna.SC5420_S10_R1_001.fastq.gz_gene_cov.tsv")[, mean_cov:=(forward_average_coverage+reverse_average_coverage)/2]
filt_gene_cov_3 <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/dirseq_out/filtrate_anme_removed_chimeras.fna.SC5422_S12_R1_001.fastq.gz_gene_cov.tsv")[, mean_cov:=(forward_average_coverage+reverse_average_coverage)/2]
#Normalize coverages by mean coverage
filt_gene_cov_1[, norm_cov:=(mean_cov/mean(filt_gene_cov$mean_cov))]
filt_gene_cov_2[, norm_cov:=(mean_cov/mean(filt_gene_cov$mean_cov))]
filt_gene_cov_3[, norm_cov:=(mean_cov/mean(filt_gene_cov$mean_cov))]
# Add CDS length
filt_gene_cov_1[, length:=end-start]
filt_gene_cov_2[, length:=end-start]
filt_gene_cov_3[, length:=end-start]
# Melt tables to be more manageable
filt_melt_1 <- melt(filt_gene_cov_1[annotation!="hypothetical protein"], id.vars=c("annotation", "length"), measure.vars = c("mean_cov"), value.name="filt_sample_1")
filt_melt_2 <- melt(filt_gene_cov_2[annotation!="hypothetical protein"], id.vars=c("annotation", "length"), measure.vars = c("mean_cov"), value.name="filt_sample_2")
filt_melt_3 <- melt(filt_gene_cov_3[annotation!="hypothetical protein"], id.vars=c("annotation", "length"), measure.vars = c("mean_cov"), value.name="filt_sample_3")
# Join tables
filt_join <- data.table(inner_join(filt_melt_1, filt_melt_2, by=c("annotation", "length", "variable")))
filt_join <- data.table(inner_join(filt_join, filt_melt_3, by=c("annotation", "length", "variable")))
  
  
retn_gene_cov_1 <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/dirseq_out/retentate_anme_nanopore_redo.fna.SB8775_S42_R1_001.fastq.gz_gene_cov.tsv")[, mean_cov:=(forward_average_coverage+reverse_average_coverage)/2]
retn_gene_cov_2 <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/dirseq_out/retentate_anme_nanopore_redo.fna.SC5421_S11_R1_001.fastq.gz_gene_cov.tsv")[, mean_cov:=(forward_average_coverage+reverse_average_coverage)/2]
retn_gene_cov_3 <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/dirseq_out/retentate_anme_nanopore_redo.fna.SC5423_S13_R1_001.fastq.gz_gene_cov.tsv")[, mean_cov:=(forward_average_coverage+reverse_average_coverage)/2]
retn_gene_cov_4 <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/dirseq_out/retentate_anme_nanopore_redo.fna.SC5423_S48_R1_001.fastq.gz_gene_cov.tsv")[, mean_cov:=(forward_average_coverage+reverse_average_coverage)/2]
retn_gene_cov_1[, norm_cov:=(mean_cov/mean(retn_gene_cov$mean_cov))]
retn_gene_cov_2[, norm_cov:=(mean_cov/mean(retn_gene_cov$mean_cov))]
retn_gene_cov_3[, norm_cov:=(mean_cov/mean(retn_gene_cov$mean_cov))]
retn_gene_cov_4[, norm_cov:=(mean_cov/mean(retn_gene_cov$mean_cov))]
retn_gene_cov_1[, length:=end-start]
retn_gene_cov_2[, length:=end-start]
retn_gene_cov_3[, length:=end-start]
retn_gene_cov_4[, length:=end-start]
retn_melt_1 <- melt(retn_gene_cov_1[annotation!="hypothetical protein"], id.vars=c("annotation", "length"), measure.vars = c("mean_cov"), value.name="retn_sample_1")
retn_melt_2 <- melt(retn_gene_cov_2[annotation!="hypothetical protein"], id.vars=c("annotation", "length"), measure.vars = c("mean_cov"), value.name="retn_sample_2")
retn_melt_3 <- melt(retn_gene_cov_3[annotation!="hypothetical protein"], id.vars=c("annotation", "length"), measure.vars = c("mean_cov"), value.name="retn_sample_3")
retn_melt_4 <- melt(retn_gene_cov_4[annotation!="hypothetical protein"], id.vars=c("annotation", "length"), measure.vars = c("mean_cov"), value.name="retn_sample_4")
retn_join <- data.table(inner_join(retn_melt_1, retn_melt_2, by=c("annotation", "length", "variable")))
retn_join <- data.table(inner_join(retn_join, retn_melt_3, by=c("annotation", "length", "variable")))
retn_join <- data.table(inner_join(retn_join, retn_melt_4, by=c("annotation", "length", "variable")))
retn_join

# ~~~~~~~~~~~~~~ VCF Files ~~~~~~~~~~~~~~~~~~

retn_vcf <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/00-lorikeet_out/01-fractions/retentate_anme_nanopore_redo.vcf")
filt_vcf <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/00-lorikeet_out/01-fractions/filtrate_anme_removed_chimeras.vcf")

retn_vcf_rev <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/00-lorikeet_out/03-reversed/retentate_anme_nanopore_redo.vcf")
filt_vcf_rev <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/00-lorikeet_out/03-reversed/filtrate_anme_removed_chimeras.vcf")

# ~~~~~~~~~~~~~~ Gene info files ~~~~~~~~~~~~

retentate_nano <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/prokka_out/retentate_anme_nanopore_redo.tsv", header=T, sep = "\t", fill=TRUE)
filtrate_nano <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/prokka_out/filtrate_anme_removed_chimeras.tsv", fill=T)

retentate_illumina <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/prokka_out/ANME_Retentate_fraction.reassembled.tsv", fill=T)
filtrate_illumina <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/prokka_out/ANME_Filter_fraction.reassembled.tsv", fill=T)
fausi_illumina <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/prokka_out/GCA_000685155.1_ANME2D_V10_genomic.tsv", fill=T)

retentate_nano_gff <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/prokka_out/retentate_anme_nanopore_redo.gff", sep2=";")
filtrate_nano_gff <- fread("/work/microbiome/abisko/rhys/01-projects/02-anme/03-genomes/prokka_out/filtrate_anme_removed_chimeras.gff", sep2=";")