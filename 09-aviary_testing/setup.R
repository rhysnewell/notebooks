library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw())

set_binner <- function(bin_stats, bin) {
    bin_stats[bin_id %like% bin, binner:=bin]
}

set_binner_das <- function(bin_stats, bin) {
    bin_stats[, binner:=bin]
}

set_binner_for_all <- function(bin_stats) {
    set_binner(bin_stats, 'rosella')
    set_binner(bin_stats, 'vamb')
    set_binner(bin_stats, 'metabat2')
    set_binner(bin_stats, 'metabat_sens')
    set_binner(bin_stats, 'metabat_ssens')
    set_binner(bin_stats, 'metabat_spec')
    set_binner(bin_stats, 'metabat_sspec')
    set_binner(bin_stats, 'maxbin')
    set_binner(bin_stats, 'concoct')    
}

generate_ranks <- function(checkm_out, bin) {
    checkm_out[`Bin Id` %like% bin, binner:=bin]
    checkm_out[`Bin Id` %like% bin & Completeness >= 50 & Contamination <= 10, com_rank := frank(checkm_out[`Bin Id` %like% bin & Completeness >= 50 & Contamination <= 10], -Completeness, Contamination)]
    checkm_out[`Bin Id` %like% bin & Completeness >= 50 & Contamination <= 10, con_rank := frank(checkm_out[`Bin Id` %like% bin & Completeness >= 50 & Contamination <= 10], Contamination, -Completeness)]
}

generate_ranks_das <- function(checkm_out, bin) {
    checkm_out[, binner := bin]
    checkm_out[Completeness >= 50 & Contamination <= 10, com_rank := frank(checkm_out[Completeness >= 50 & Contamination <= 10], -Completeness, Contamination)]
    checkm_out[Completeness >= 50 & Contamination <= 10, con_rank := frank(checkm_out[Completeness >= 50 & Contamination <= 10], Contamination, -Completeness)]    
}

change_names <- function(checkm_out) {
    checkm_out[binner=='metabat2', binner:='MetaBAT2']
    checkm_out[binner=='metabat_sens', binner:='MetaBAT Sens.']
    checkm_out[binner=='metabat_ssens', binner:='MetaBAT SuperSens.']
    checkm_out[binner=='metabat_spec', binner:='MetaBAT Spec.']
    checkm_out[binner=='metabat_sspec', binner:='MetaBAT SuperSpec.']
    checkm_out[binner=='maxbin', binner:='MaxBin2']
    checkm_out[binner=='concoct', binner:='CONCOCT']
    checkm_out[binner=='vamb', binner:='VAMB']
    checkm_out[binner=='rosella', binner:='Rosella']
}

generate_ranks_refined <- function(checkm_out) {
    checkm_out[Completeness >= 50 & Contamination <= 10, com_rank := frank(-Completeness), by=binner]
    checkm_out[Completeness >= 50 & Contamination <= 10, con_rank := frank(Contamination), by=binner]
}

get_thresholds <- function(checkm_out, check_vamb=FALSE) {
#     output <- checkm_out[, .N, by=binner][order(binner)]
    output <- checkm_out[Completeness>=95 & Contamination<=5, .N, by=binner][order(binner)]
    output[, t1:=`N`]
    output[, N:=NULL]
    output <- full_join(output, checkm_out[Completeness>=90 & Contamination<=10, .N, by=binner][order(binner)])
    output[, t2:=`N`]
    output[, N:=NULL]
    output <- full_join(output, checkm_out[Completeness>=80 & Contamination<=10, .N, by=binner][order(binner)])
    output[, t3:=`N`]
    output[, N:=NULL]
    output <- full_join(output, checkm_out[Completeness>=70 & Contamination<=10, .N, by=binner][order(binner)])
    output[, t4:=`N`]
    output[, N:=NULL]
    output <- full_join(output, checkm_out[Completeness>=50 & Contamination<=10, .N, by=binner][order(binner)])
    output[, t5:=`N`]
    output[, N:=NULL]
    if (check_vamb) {
        if ('VAMB' %in% output$binner) {
            return(output[order(binner)])
        } else {
            output <- rbind(output, list('VAMB', 0, 0, 0, 0, 0))
        }
    } else {
        return(output)
    }
}

get_near_complete_counts <- function(checkm_out, contig_count=1, check_vamb=FALSE) {
#     output <- checkm_out[, .N, by=binner][order(binner)]
    output <- checkm_out[completeness>=95 & contamination<=10 & n_contigs == contig_count, .N, by=binner][order(binner)]
    output[, t1:=`N`]
    output[, N:=NULL]
    output <- full_join(output, checkm_out[completeness>=90 & contamination<=10 & n_contigs == contig_count, .N, by=binner][order(binner)])
    output[, t2:=`N`]
    output[, N:=NULL]
    if (check_vamb) {
        if ('VAMB' %in% output$binner) {
            return(output[order(binner)])
        } else {
            output <- rbind(output, list('VAMB', 0, 0))
        }
    } else {
        return(output)
    }
}

# Color generating function found here:
# http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

generate_plots <- function(aviary_path, metawrap_path) {
    aviary_checkm <- paste0(aviary_path, '/bins/checkm.out')
    atlas_checkm <- paste0(aviary_path, 'data/atlas_dastool/checkm.out') 
    das_checkm <- paste0(aviary_path, '/data/das_tool_bins_no_refine/checkm.out')
    metawrap_checkm <- paste0(metawrap_path, '/checkm.out')

    aviary <- fread(aviary_checkm)
    atlas <- fread(atlas_checkm)
    das <- fread(das_checkm)
    metawrap <- fread(metawrap_checkm)

    generate_ranks_das(das, 'DASTool')
    generate_ranks_das(aviary, 'Aviary')
    generate_ranks_das(atlas, 'ATLAS')
    generate_ranks_das(metawrap, 'MetaWRAP')
    
    all <- rbind(das, aviary, atlas, metawrap)
    adj_names = sort(unique(setdiff(c('DASTool', 'Aviary', 'ATLAS', 'MetaWRAP'), 'Aviary')))
    values = gg_color_hue(length(adj_names))
    names(values) = adj_names
    values = c(values, c(Aviary="#000000"))
    
    com <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=com_rank, y=Completeness, color=binner)) +
        ylim(50, 100) +
        scale_colour_manual("", values=values) +
        labs(x="Com. rank", color="", y="") +
        guides(color = guide_legend(override.aes = list(color = values))) +
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="bottom" ,
              legend.direction="horizontal", 
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    con <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=con_rank, y=Contamination, color=binner)) +   
        ylim(0, 10) +
        scale_colour_manual("", values=values) +
        labs(x="Con. rank", color="", y="") +
        guides(color = guide_legend(override.aes = list(color = values))) +
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="bottom" ,
              legend.direction="horizontal", 
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())


    t_all <- get_thresholds(all, check_vamb=FALSE)
    
    return(list(com, con, t_all, all))
}

get_total_plots <- function(...) {
    all <- lapply(..., function(x){x[[4]]}) # get results for all binners
    
    all <- do.call(rbind, all)
    

    generate_ranks_refined(all)
#     change_names(all)
    
    adj_names = sort(unique(setdiff(c('DASTool', 'Aviary', 'ATLAS', 'MetaWRAP'), 'Aviary')))
    values = gg_color_hue(length(adj_names))
    names(values) = adj_names
    values = c(values, c(Aviary="#000000"))
    
    com <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=com_rank, y=Completeness, color=binner)) +
        ylim(50, 100) +
        scale_colour_manual("", values=values) +
        labs(x="Com. rank", color="", y="") +
        guides(color = guide_legend(override.aes = list(color = values))) +
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="bottom" ,
              legend.direction="horizontal", 
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    con <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=con_rank, y=Contamination, color=binner)) +   
        ylim(0, 10) +
        scale_colour_manual("", values=values) +
        labs(x="Con. rank", color="", y="") +
        guides(color = guide_legend(override.aes = list(color = values))) +
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="bottom" ,
              legend.direction="horizontal", 
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())


    t_all <- get_thresholds(all, check_vamb=FALSE)

    return(list(com, con, t_all, all))
}

'%!in%' <- function(x,y)!('%in%'(x,y))



process_stats <- function(results, file_structure= '/www/metaquast/summary/TSV/Genome_fraction.tsv', genomes=c("B. subtilis","C. neoformans","En. faecalis","E. coli","L. fermentum","L. monocytogenes","P. aeruginosa","S. cerevisiae","Sa. enterica","St. aureus"), euk=c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), use_genomes=TRUE) {
    stats <- lapply(results, function(x) {
        result <- fread(paste0(x, file_structure))
        result$sample <- gsub('\\/', '', tstrsplit(x, 'final_assemblies', keep=2L))
        if (use_genomes) {
            result$Assemblies <- genomes
            result$euk <- euk
        }

        return(result)
    })
    stats <- do.call(rbind, stats)
    if (use_genomes) {
        stats <- melt(stats, id.vars=c("Assemblies", "sample", "euk"))
        stats <- transform(stats, value = as.numeric(value))
        stats[, value:=ifelse(is.na(value), 0.0, value)]
        return(stats)
    } else {
        stats <- melt(stats, id.vars=c("Assemblies", "sample"))
        stats <- transform(stats, value = as.numeric(value))
        stats[, value:=ifelse(is.na(value), 0.0, value)]
        return(stats)
    }
}

get_quast_stats <- function(directory = "/mnt/hpccs01/scratch/microbiome/n10853499/03-aviary_testing/00-zymo_hybrid_assembly_benchmark/final_assemblies/", use_genomes=FALSE) {
    results <- list.dirs(directory, recursive=FALSE)
#     print(results)
    genomes <- c("B. subtilis","C. neoformans","En. faecalis","E. coli","L. fermentum","L. monocytogenes","P. aeruginosa","S. cerevisiae","Sa. enterica","St. aureus")
    genomes_alt <- c("B. subtilis","C. neoformans","En. faecalis","E. coli","L. fermentum","L. monocytogenes","P. aeruginosa","S. cerevisiae","Sa. enterica","St. aureus", "Not Aligned")
    euk <- c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)
    euk_alt <- c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE)
    fractions <- process_stats(results)
    largest_alignment <- process_stats(results,'/www/metaquast/summary/TSV/Largest_alignment.tsv')
    largest_contig <- process_stats(results, '/www/metaquast/summary/TSV/Largest_contig.tsv', genomes=genomes_alt, euk=euk_alt)
    contig_count <- process_stats(results, '/www/metaquast/summary/TSV/num_contigs.tsv', genomes=genomes_alt, euk=euk_alt)
    num_indels <- process_stats(results, '/www/metaquast/summary/TSV/num_indels_per_100_kbp.tsv')
    num_mismatches <- process_stats(results, '/www/metaquast/summary/TSV/num_mismatches_per_100_kbp.tsv')
    num_misassemblies <- process_stats(results, '/www/metaquast/summary/TSV/num_misassemblies.tsv')
    
    
    return(list("fraction_recovered"=fractions,"largest_alignment"=largest_alignment,"largest_contig"=largest_contig,"contig_count"=contig_count, "indels"=num_indels, "mismatches"=num_mismatches, "misassemblies"=num_misassemblies))
}

get_cami_stats <- function(directory = "/mnt/hpccs01/microbiome/n10853499/00-rosella_testing/01-CAMI_II/CAMI_Airways/binning/assemblies/results/") {
    results <- list.dirs(directory, recursive=FALSE)
    fractions <- process_cami_stats(results)
    largest_alignment <- process_cami_stats(results,'/www/metaquast/summary/TSV/Largest_alignment.tsv')
    largest_contig <- process_cami_stats(results, '/www/metaquast/summary/TSV/Largest_contig.tsv')
    contig_count <- process_cami_stats(results, '/www/metaquast/summary/TSV/num_contigs.tsv')
#     num_indels <- process_cami_stats(results, '/www/metaquast/summary/TSV/num_indels_per_100_kbp.tsv')
#     num_mismatches <- process_cami_stats(results, '/www/metaquast/summary/TSV/num_mismatches_per_100_kbp.tsv')
#     num_misassemblies <- process_cami_stats(results, '/www/metaquast/summary/TSV/num_misassemblies.tsv')
    return(list("fraction_recovered"=fractions, "largest_alignment"=largest_alignment, "largest_contig"=largest_contig,"contig_count"=contig_count))
}

process_cami_stats <- function(results, file_structure= '/www/metaquast/summary/TSV/Genome_fraction.tsv') {
    stats <- lapply(results, function(x) {
#         print(paste0(x, file_structure))
        result <- fread(paste0(x, file_structure))
        result$sample <- gsub('sample_', 'Sample ', gsub('\\/', '', tstrsplit(x, 'results', keep=2L)))
        return(result)
    })
    stats <- do.call(rbind, stats)
    
    stats <- melt(stats, id.vars=c("Assemblies", "sample"))
    stats <- transform(stats, value = as.numeric(value))
    stats[, value:=ifelse(is.na(value), 0.0, value)]
    stats[, variable:=gsub("_", "-", variable)]
    return(stats)
  
}

retrieve_single_stats <- function(file = "/mnt/hpccs01/scratch/microbiome/n10853499/03-aviary_testing/00-zymo_hybrid_assembly_benchmark/stats/SRR10084338/stats_aviary.txt", tool_name = "Aviary") {
    test_aviary <- fread(file, skip='-')
#     print(readLines(file))
    max_scaffold <- as.integer(gsub("KB", "", gsub(" MB", "000", gsub("\\.", "", tstrsplit(readLines(file)[12], '\t', keep=2L)[[1]]))))
    colnames(test_aviary) <- c("min_scaffold_length", "number_of_scaffolds", "number_of_contigs", "total_scaffold_length", "total_contig_length", "scaffold_contig_coverage")
    test_aviary[, assembler:=tool_name]
    test_aviary[, assembly_size:=test_aviary[min_scaffold_length=="All"]$total_contig_length]
    test_aviary <- test_aviary[min_scaffold_length %!in% c("1 KB", "500", "250", "100", "50"), ]
    test_aviary[, min_scaffold_length:=ifelse(min_scaffold_length == "All", "1 KB", min_scaffold_length)]
    entire_file <- readLines(file)
    test_aviary[, N50:=str_trim(rev(str_split(entire_file[grep("Main genome contig N/L50", entire_file)], '/', simplify=TRUE))[1])]
    test_aviary[, max_contig_size:=str_trim(rev(str_split(entire_file[grep("Max contig length", entire_file)], ':', simplify=TRUE))[1])]
    test_aviary[, assembly_size_short:=str_trim(rev(str_split(entire_file[grep("Main genome scaffold sequence total", entire_file)], ':', simplify=TRUE))[1])]
    test_aviary[, max_scaffold_size:=max_scaffold]
    return(test_aviary)
}

retrieve_all_stats <- function(directory = "/mnt/hpccs01/scratch/microbiome/n10853499/03-aviary_testing/00-zymo_hybrid_assembly_benchmark/stats/", sample_name = "SRR10084338", with_spades=TRUE) {
    test_aviary <- retrieve_single_stats(paste0(directory, sample_name, "/stats_aviary.txt"), "Aviary")
#     test_metaspades_contigs <- retrieve_single_stats(paste0(directory, sample_name, "/stats_metaspades_contigs.txt"), "metaSPAdes contigs")
    if (with_spades) {
        test_metaspades_scaffolds <- retrieve_single_stats(paste0(directory, sample_name, "/stats_metaspades_scaffolds.txt"), "metaSPAdes")
        test_flye <- retrieve_single_stats(paste0(directory, sample_name, "/stats_flye.txt"), "metaFlye")
        test_operams <- retrieve_single_stats(paste0(directory, sample_name, "/stats_operams.txt"), "OPERA-MS")

        test_stats <- rbind(test_aviary, test_flye, test_metaspades_scaffolds, test_operams)
        test_stats[, number_of_scaffolds:=as.numeric(gsub(",", "", number_of_scaffolds))]
        test_stats[, number_of_contigs:=as.numeric(gsub(",", "", number_of_contigs))]
        test_stats[, total_scaffold_length:=as.numeric(gsub(",", "", total_scaffold_length))]
        test_stats[, total_contig_length:=as.numeric(gsub(",", "", total_contig_length))]
        test_stats[, assembly_size:=as.numeric(gsub(",", "", assembly_size))]
        test_stats <- test_stats[order(nrow(test_stats):1)]
        test_stats[, gained_size := total_scaffold_length - shift(total_scaffold_length), by = assembler]
        test_stats[, gained_size := ifelse(is.na(gained_size), total_scaffold_length, gained_size)]
        test_stats[, min_scaffold_length := factor(min_scaffold_length, levels = c("1 KB", "2.5 KB", "5 KB", "10 KB", "25 KB", "50 KB", "100 KB", "250 KB", "500 KB", "1 MB", "2.5 MB", "5 MB", "10 MB"))]
        test_stats[, sample_name := sample_name]

        return(test_stats)
    } else {
        test_flye <- retrieve_single_stats(paste0(directory, sample_name, "/stats_flye.txt"), "metaFlye")
        test_operams <- retrieve_single_stats(paste0(directory, sample_name, "/stats_operams.txt"), "OPERA-MS")

        test_stats <- rbind(test_aviary, test_flye, test_operams)
        test_stats[, number_of_scaffolds:=as.numeric(gsub(",", "", number_of_scaffolds))]
        test_stats[, number_of_contigs:=as.numeric(gsub(",", "", number_of_contigs))]
        test_stats[, total_scaffold_length:=as.numeric(gsub(",", "", total_scaffold_length))]
        test_stats[, total_contig_length:=as.numeric(gsub(",", "", total_contig_length))]
        test_stats[, assembly_size:=as.numeric(gsub(",", "", assembly_size))]
        test_stats <- test_stats[order(nrow(test_stats):1)]
        test_stats[, gained_size := total_scaffold_length - shift(total_scaffold_length), by = assembler]
        test_stats[, gained_size := ifelse(is.na(gained_size), total_scaffold_length, gained_size)]
        test_stats[, min_scaffold_length := factor(min_scaffold_length, levels = c("1 KB", "2.5 KB", "5 KB", "10 KB", "25 KB", "50 KB", "100 KB", "250 KB", "500 KB", "1 MB", "2.5 MB", "5 MB", "10 MB"))]
        test_stats[, sample_name := sample_name]

        return(test_stats)    
    
    }

}

concatenate_all_results <- function(directory = "/mnt/hpccs01/scratch/microbiome/n10853499/03-aviary_testing/00-zymo_hybrid_assembly_benchmark/stats/", with_spades=TRUE) {
    sample_names <- lapply(list.dirs(directory, recursive=FALSE), function(x) rev(str_split(x, '/', simplify=TRUE))[1])
#     sample_names <- c("SRR10084338", "SRR10084340")
    bind_rows(lapply(sample_names, function(x) retrieve_all_stats(directory, sample_name=x, with_spades=with_spades)))
}                     