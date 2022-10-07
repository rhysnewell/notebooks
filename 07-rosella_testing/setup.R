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

generate_ranks_refined <- function(checkm_out, contamination_threshold=10) {
    checkm_out[Completeness >= 50 & Contamination <= contamination_threshold, com_rank := frank(-Completeness), by=binner]
    checkm_out[Completeness >= 50 & Contamination <= contamination_threshold, con_rank := frank(Contamination), by=binner]
}

get_thresholds <- function(checkm_out, check_vamb=FALSE, contamination_threshold=10) {
#     output <- checkm_out[, .N, by=binner][order(binner)]
    output <- checkm_out[Completeness>=95 & Contamination<=contamination_threshold, .N, by=binner][order(binner)]
    output[, t1:=`N`]
    output[, N:=NULL]
    output <- full_join(output, checkm_out[Completeness>=90 & Contamination<=contamination_threshold, .N, by=binner][order(binner)])
    output[, t2:=`N`]
    output[, N:=NULL]
    output <- full_join(output, checkm_out[Completeness>=80 & Contamination<=contamination_threshold, .N, by=binner][order(binner)])
    output[, t3:=`N`]
    output[, N:=NULL]
    output <- full_join(output, checkm_out[Completeness>=70 & Contamination<=contamination_threshold, .N, by=binner][order(binner)])
    output[, t4:=`N`]
    output[, N:=NULL]
    output <- full_join(output, checkm_out[Completeness>=50 & Contamination<=contamination_threshold, .N, by=binner][order(binner)])
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


ranks_for_all <- function(checkm_out) {
    generate_ranks(checkm_out, 'rosella')
    generate_ranks(checkm_out, 'vamb')
    generate_ranks(checkm_out, 'metabat2')
    generate_ranks(checkm_out, 'metabat_sens')
    generate_ranks(checkm_out, 'metabat_ssens')
    generate_ranks(checkm_out, 'metabat_spec')
    generate_ranks(checkm_out, 'metabat_sspec')
    generate_ranks(checkm_out, 'maxbin')
    generate_ranks(checkm_out, 'concoct')
}


# Color generating function found here:
# http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

add_refine <- function(out_path, previous_plots) {
    rosella_refined <- fread(paste0(out_path, '/data/rosella_refine_rosella/checkm.out'))
    metabat_refined <- fread(paste0(out_path, '/data/rosella_refine_metabat2/checkm.out'))

    tryCatch(generate_ranks_das(rosella_refined, 'Rosella Refined'), error = function(e) NULL)
    tryCatch(generate_ranks_das(metabat_refined, 'MetaBAT2 Refined'), error = function(e) NULL)
   
    
    rosella_all <- previous_plots[[6]][binner == "Rosella"]
    metabat_all <- previous_plots[[6]][binner == "MetaBAT2"]
    
    
    rosella_all <- tryCatch(rbind(rosella_all, rosella_refined), error = function(e) rosella_all)
    metabat_all <- tryCatch(rbind(metabat_all, metabat_refined), error = function(e) metabat_all)
    
    
    generate_ranks_das(rosella_all, 'Rosella Refined')
    generate_ranks_das(metabat_all, 'MetaBAT2 Refined')

    dastool_refined <- tryCatch(fread(paste0(out_path, '/data/rosella_refine_das_tool/checkm.out')), error = function(e) NULL)
    if (!is.null(dastool_refined) && nrow(dastool_refined) > 0) {
        generate_ranks_das(dastool_refined, 'DASTool Refined')
        dastool_all <- previous_plots[[7]]
        dastool_all <- rbind(dastool_all, dastool_refined)
        generate_ranks_das(dastool_all, 'DASTool Refined')
    } else {
        dastool_all <- data.table(previous_plots[[7]])
        generate_ranks_das(dastool_all, 'DASTool Refined')
    }
    
    all <- rbind(previous_plots[[6]], rosella_all, metabat_all, previous_plots[[7]], previous_plots[[8]], dastool_all)
    

    if ("VAMB" %in% all$binner) {
        values = c(CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/o Rosella`="#A3A500", `DASTool Refined`="#D89000", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", `MetaBAT2`="#E76BF3", `MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000")
#     values = c(values, c(`MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000"))
        print(values)
    
        all$binner = factor(all$binner, levels=c('CONCOCT', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'DASTool Refined', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.',  'MetaBAT2', 'MetaBAT2 Refined', 'VAMB', 'Rosella', 'Rosella Refined'))
        linetypes = c(3, 2, 2, 1, 3, 3, 3, 3, 3, 3, 1, 3, 3, 1)
    } else {
        values = c(CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/o Rosella`="#A3A500", `DASTool Refined`="#D89000", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", `MetaBAT2`="#E76BF3", `MetaBAT2 Refined`="#E76BF3", Rosella="#000000", `Rosella Refined`="#000000")
#     values = c(values, c(`MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000"))
        print(values)
    
        all$binner = factor(all$binner, levels=c('CONCOCT', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'DASTool Refined', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.',  'MetaBAT2', 'MetaBAT2 Refined', 'Rosella', 'Rosella Refined'))
        linetypes = c(3, 2, 2, 1, 3, 3, 3, 3, 3, 3, 1, 3, 1)
    }
    com <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=com_rank, y=Completeness, color=binner, linetype=binner)) +
        ylim(50, 100) +
        scale_colour_manual("Binner", values=values) +
        scale_linetype_manual("Binner", values = linetypes, guide="none") + 
        labs(x="Com. rank", color="Binner", y="", linetype="Binner") +
        guides(color = guide_legend(override.aes = list(color = values, linetype = linetypes))) +
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
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    con <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=con_rank, y=Contamination, color=binner, linetype=binner)) +    
        ylim(0, 10) +
        scale_colour_manual("", values=values) +
        scale_linetype_manual("", values = linetypes, guide="none") + 
        labs(x="Con. rank", color="Binner", y="", linetype="") +
        guides(color = guide_legend(override.aes = list(color = values, linetype = linetypes))) +
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


    t_all <- get_thresholds(all, check_vamb=TRUE)
    
    return(list(com, con, t_all, all))
}
                                

generate_plots <- function(out_path) {
    all_path <- paste0(out_path, '/data/all_bins/checkm.out')
    das_path <- paste0(out_path, '/data/das_tool_bins/checkm.out')
    backup_path <- paste0(out_path, '/data/checkm.out')
    das_no_rosella_path <- paste0(out_path, '/data/checkm_without_rosella.out')
    
    all <- fread(all_path)
#     tryCatch(generate_ranks_das(rosella_refined, 'Rosella Refined'), error = function(e) NULL)
    das_all <- tryCatch(fread(das_path), error = function(e) fread(backup_path))
    das_wor <- fread(das_no_rosella_path)
    
    ranks_for_all(all)
    generate_ranks_das(das_all, 'DASTool w/ Rosella')
    generate_ranks_das(das_wor, 'DASTool w/o Rosella')
    change_names(all)
    
    # Dynamically generate default color values, but have Rosella="black".
    # adj_names = sort(unique(setdiff(c(all$binner, 'DASTool w/ Rosella', 'DASTool w/o Rosella'), 'Rosella')))

    # Had to change this to be static because VAMB doesn't show up sometimes which 
    # breaks the manual linetype values I set in the plots
    adj_names = sort(unique(setdiff(c('CONCOCT', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.',  'MetaBAT2', 'VAMB', 'Rosella', 'DASTool w/ Rosella', 'DASTool w/o Rosella'), 'Rosella')))
    values = gg_color_hue(length(adj_names))
    names(values) = adj_names
    values = c(values, c(Rosella="#000000"))
    print(values)
    
    com <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=com_rank, y=Completeness, color=binner, linetype=binner)) +
        geom_line(data=das_all[Completeness >= 50 & Contamination <= 10,], aes(x=com_rank, y=Completeness, color=binner, linetype=binner)) +
        geom_line(data=das_wor[Completeness >= 50 & Contamination <= 10,], aes(x=com_rank, y=Completeness, color=binner, linetype=binner)) +
        ylim(50, 100) +
        scale_colour_manual("Binner", values=values) +
        scale_linetype_manual("Binner", values = c(1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1), guide="none") + 
        labs(x="Com. rank", color="Binner", y="", linetype="Binner") +
        guides(color = guide_legend(override.aes = list(color = values, linetype = c(1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1)))) +
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
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    con <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=con_rank, y=Contamination, color=binner, linetype=binner)) +
        geom_line(data=das_all[Completeness >= 50 & Contamination <= 10,], aes(x=con_rank, y=Contamination, color=binner, linetype=binner)) +
        geom_line(data=das_wor[Completeness >= 50 & Contamination <= 10,], aes(x=con_rank, y=Contamination, color=binner, linetype=binner)) +       
        ylim(0, 10) +
        scale_colour_manual("", values=values) +
        scale_linetype_manual("", values = c(1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1), guide="none") + 
        labs(x="Con. rank", color="Binner", y="", linetype="") +
        guides(color = guide_legend(override.aes = list(color = values, linetype = c(1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1)))) +
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


    t_all <- get_thresholds(all, check_vamb=TRUE)
    t_das <- get_thresholds(das_all)
    t_wor <- get_thresholds(das_wor)
    
    return(list(com, con, t_all, t_das, t_wor, all, das_all, das_wor))
}


get_total_plots <- function(...) {
    all <- lapply(..., function(x){x[[6]]}) # get results for all binners
    das_all <- lapply(..., function(x){x[[7]]})
    das_wor <- lapply(..., function(x){x[[8]]})
    
    all <- do.call(rbind, all)
    das_all <- do.call(rbind, das_all)
    das_wor <- do.call(rbind, das_wor)
    
    ranks_for_all(all)
    generate_ranks_das(das_all, 'DASTool w/ Rosella')
    generate_ranks_das(das_wor, 'DASTool w/o Rosella')
    change_names(all)
    
    values = c(CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/o Rosella`="#A3A500", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", `MetaBAT2`="#E76BF3", VAMB="#FF62BC", Rosella="#000000")
#     values = c(values, c(`MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000"))
    print(values)
    
    all$binner = factor(all$binner, levels=c('CONCOCT', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.',  'MetaBAT2', 'VAMB', 'Rosella'))
    

    com <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=com_rank, y=Completeness, color=binner, linetype=binner)) +
        geom_line(data=das_all[Completeness >= 50 & Contamination <= 10,], aes(x=com_rank, y=Completeness, color=binner, linetype=binner)) +
        geom_line(data=das_wor[Completeness >= 50 & Contamination <= 10,], aes(x=com_rank, y=Completeness, color=binner, linetype=binner)) +
        labs(x="Com. rank", color="Binner", y="", linetype="Binner") +
        ylim(50, 100) +
        scale_linetype_manual("Binner", values = c(1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1), guide="none") + 
        scale_colour_manual("Binner", values=values) +
        guides(color = guide_legend(override.aes = list(color = values, linetype = c(1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1)))) +
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
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    con <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=con_rank, y=Contamination, color=binner, linetype=binner)) +
        geom_line(data=das_all[Completeness >= 50 & Contamination <= 10,], aes(x=con_rank, y=Contamination, color=binner, linetype=binner)) +
        geom_line(data=das_wor[Completeness >= 50 & Contamination <= 10,], aes(x=con_rank, y=Contamination, color=binner, linetype=binner)) +
        labs(x="Con. rank", color="Binner", y="", linetype="Binner") +
        ylim(0, 10) +
        scale_linetype_manual("Binner", values = c(1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1), guide="none") + 
        scale_colour_manual("Binner", values=values) +
        guides(color = guide_legend(override.aes = list(color = values, linetype = c(1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1)))) +
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
    t_all <- get_thresholds(all, check_vamb=TRUE)
    t_das <- get_thresholds(das_all)
    t_wor <- get_thresholds(das_wor)
    return(list(com, con, t_all, t_das, t_wor, all, das_all, das_wor))
}

get_total_plots_refined <- function(...) {
    all <- lapply(..., function(x){x[[4]]}) # get results for all binners
    
    all <- do.call(rbind, all)
    

    generate_ranks_refined(all)
#     change_names(all)
    
    if ("VAMB" %in% all$binner) {
        values = c(CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/o Rosella`="#A3A500", `DASTool Refined`="#D89000", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", `MetaBAT2`="#E76BF3", `MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000")
#     values = c(values, c(`MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000"))
        print(values)
    
        all$binner = factor(all$binner, levels=c('CONCOCT', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'DASTool Refined', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.',  'MetaBAT2', 'MetaBAT2 Refined', 'VAMB', 'Rosella', 'Rosella Refined'))
        linetypes = c(3, 2, 2, 1, 3, 3, 3, 3, 3, 3, 1, 3, 3, 1)
    } else {
        values = c(CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/o Rosella`="#A3A500", `DASTool Refined`="#D89000", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", `MetaBAT2`="#E76BF3", `MetaBAT2 Refined`="#E76BF3", Rosella="#000000", `Rosella Refined`="#000000")
#     values = c(values, c(`MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000"))
        print(values)
    
        all$binner = factor(all$binner, levels=c('CONCOCT', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'DASTool Refined', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.',  'MetaBAT2', 'MetaBAT2 Refined', 'Rosella', 'Rosella Refined'))
        linetypes = c(3, 2, 2, 1, 3, 3, 3, 3, 3, 3, 1, 3, 1)
    }
    
    com <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=com_rank, y=Completeness, color=binner, linetype=binner)) +
        labs(x="Com. rank", color="Binner", y="", linetype="Binner") +
        ylim(50, 100) +
        scale_linetype_manual("Binner", values = linetypes, guide="none") + 
        scale_colour_manual("Binner", values=values) +
        guides(color = guide_legend(override.aes = list(color = values, linetype = linetypes))) +
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
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    con <- ggplot() + 
        geom_line(data=all[Completeness >= 50 & Contamination <= 10,], aes(x=con_rank, y=Contamination, color=binner, linetype=binner)) +
        labs(x="Con. rank", color="Binner", y="", linetype="Binner") +
        ylim(0, 10) +
        scale_linetype_manual("Binner", values = linetypes, guide="none") + 
        scale_colour_manual("Binner", values=values) +
        guides(color = guide_legend(override.aes = list(color = values, linetype = linetypes))) +
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
    t_all <- get_thresholds(all, check_vamb=TRUE)

    return(list(com, con, t_all, all))
}

generate_bin_stats <- function(out_path, sample_name) {
    all_path <- paste0(out_path, '/data/bin_stats/all_bin_stats.tsv')
    das_path <- paste0(out_path, '/data/bin_stats/das_tool_wr_stats.tsv')
    das_no_rosella_path <- paste0(out_path, '/data/bin_stats/das_tool_nr_stats.tsv')
    m2_refined <- paste0(out_path, '/data/bin_stats/m2_refined_bin_stats.tsv')
    ro_refined <- paste0(out_path, '/data/bin_stats/ro_refined_bin_stats.tsv')
    dt_refined <- paste0(out_path, '/data/bin_stats/dt_refined_bin_stats.tsv')
    
    all <- fread(all_path)
    das_all <- fread(das_path)
    das_wor <- fread(das_no_rosella_path)
    m2_refined <- fread(m2_refined)
    ro_refined <- fread(ro_refined)
    dt_refined <- fread(dt_refined)
    
    set_binner_for_all(all)
    
    rosella_all <- all[binner == "rosella"]
    print(nrow(rosella_all))
    set_binner_das(ro_refined, "Rosella Refined")
    rosella_all <- rbind(ro_refined, rosella_all)
    print(nrow(rosella_all))
    set_binner_das(rosella_all, "Rosella Refined")
    print(nrow(rosella_all))
    
    m2_all <- all[binner == "metabat2"]
    set_binner_das(m2_refined, "MetaBAT2 Refined")
    m2_all <- rbind(m2_refined, m2_all)
    set_binner_das(m2_all, "MetaBAT2 Refined")
    
    
    dt_all <- rbind(dt_refined, das_all)
    set_binner_das(dt_all, "DASTool Refined")
    
    set_binner_das(das_all, 'DASTool w/ Rosella')
    set_binner_das(das_wor, 'DASTool w/o Rosella')
    change_names(all)
    all <- rbind(all, m2_all, rosella_all, das_all, das_wor, dt_all)
    all[, sample_name:=sample_name]

    # Dynamically generate default color values, but have Rosella="black".
    # adj_names = sort(unique(setdiff(c(all$binner, 'DASTool w/ Rosella', 'DASTool w/o Rosella'), 'Rosella')))

    # Had to change this to be static because VAMB doesn't show up sometimes which 
    # breaks the manual linetype values I set in the plots
    values = c(CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/o Rosella`="#A3A500", `DASTool Refined`="#D89000", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", `MetaBAT2`="#E76BF3", `MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000")
#     values = c(values, c(`MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000"))
    print(values)
    
    all$binner = factor(all$binner, levels=c('CONCOCT', 'VAMB', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.',  'MetaBAT2', 'MetaBAT2 Refined', 'Rosella', 'Rosella Refined', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'DASTool Refined'))
    
    n_contigs <- ggplot(all[completeness >= 50 & contamination <= 10], aes(x=binner, y=n_contigs, color=binner)) +
        geom_violin() +
        geom_boxplot(width=0.1) +
        scale_colour_manual("Binner", values=values) +
        scale_y_continuous(trans='log10') +
        labs(x="Binner", color="Binner", y="No. of contigs (Log10)", linetype="Binner") +
        theme(axis.text=element_text(size=10),
              axis.text.x=element_text(angle=30, vjust=1, hjust=1),
              axis.title=element_text(size=10),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="bottom" ,
              legend.direction="horizontal", 
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    size <- ggplot(all[completeness >= 50 & contamination <= 10], aes(x=binner, y=size, color=binner)) +
        geom_violin() +
        geom_boxplot(width=0.1) +
        scale_colour_manual("Binner", values=values) +
        labs(x="Binner", color="Binner", y="Size (base pairs)", linetype="Binner") +
        theme(axis.text=element_text(size=10),
              axis.text.x=element_text(angle=30, vjust=1, hjust=1),
              axis.title=element_text(size=10),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="bottom" ,
              legend.direction="horizontal", 
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())
    
    fragmentation <- ggplot(all[completeness >= 50 & contamination <= 10], aes(x=completeness, y=n_contigs)) +
#         stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
        geom_hex(bins = 20) +
        scale_fill_continuous(type = "viridis") +
#         scale_colour_manual("Binner", values=values) +
#         labs(x="Binner", color="Binner", y="Size (base pairs)", linetype="Binner") +
        theme(axis.text=element_text(size=10),
              axis.text.x=element_text(angle=30, vjust=1, hjust=1),
              axis.title=element_text(size=10),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="bottom" ,
              legend.direction="horizontal", 
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank()) +
       facet_wrap(~binner)
    
    return(list(n_contigs, size, all, fragmentation))
}

get_total_bin_stats <- function(...) {
    all <- lapply(..., function(x){x[[3]]}) # get results for all binners
    
    all <- do.call(rbind, all)
    change_names(all)
    
    thresholds <- get_near_complete_counts(all, 1, TRUE)
    # Dynamically generate default color values, but have Rosella="black".
    # adj_names = sort(unique(setdiff(c(all$binner, 'DASTool w/ Rosella', 'DASTool w/o Rosella'), 'Rosella')))

    # Had to change this to be static because VAMB doesn't show up sometimes which 
    # breaks the manual linetype values I set in the plots
    values = c(CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/o Rosella`="#A3A500", `DASTool Refined`="#D89000", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", `MetaBAT2`="#E76BF3", `MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000")
#     values = c(values, c(`MetaBAT2 Refined`="#E76BF3", VAMB="#FF62BC", Rosella="#000000", `Rosella Refined`="#000000"))
    print(values)
    
    all$binner = factor(all$binner, levels=c('CONCOCT', 'VAMB', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.',  'MetaBAT2', 'MetaBAT2 Refined', 'Rosella', 'Rosella Refined', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'DASTool Refined'))
    
    n_contigs <- ggplot(all[completeness >= 50 & contamination <= 10], aes(x=binner, y=n_contigs, color=binner)) +
        geom_violin() +
        geom_boxplot(width=0.1) +
        scale_colour_manual("Binner", values=values) +
        scale_y_continuous(trans='log10') +
        labs(x="Binner", color="Binner", y="No. of contigs (Log10)", linetype="Binner") +
        theme(axis.text=element_text(size=10),
              axis.text.x=element_text(angle=30, vjust=1, hjust=1),
              axis.title=element_text(size=10),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="bottom" ,
              legend.direction="horizontal", 
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    size <- ggplot(all[completeness >= 50 & contamination <= 10], aes(x=binner, y=size, color=binner)) +
        geom_violin() +
        geom_boxplot(width=0.1) +
        scale_colour_manual("Binner", values=values) +
        labs(x="Binner", color="Binner", y="Size (base pairs)", linetype="Binner") +
        theme(axis.text=element_text(size=10),
              axis.text.x=element_text(angle=30, vjust=1, hjust=1),
              axis.title=element_text(size=10),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="bottom" ,
              legend.direction="horizontal", 
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())
    
    fragmentation <- ggplot(all[completeness >= 50 & contamination <= 10], aes(x=completeness, y=n_contigs)) +
#         stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
        geom_hex(bins = 20) +
        scale_fill_continuous(type = "viridis") +
#         scale_colour_manual("Binner", values=values) +
#         labs(x="Binner", color="Binner", y="Size (base pairs)", linetype="Binner") +
        theme(axis.text=element_text(size=10),
              axis.text.x=element_text(angle=30, vjust=1, hjust=1),
              axis.title=element_text(size=10),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="bottom" ,
              legend.direction="horizontal", 
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank()) +
       facet_wrap(~binner)
    
    return(list(n_contigs, size, all, fragmentation, thresholds))
}


get_amber_files <- function(path = "/mnt/hpccs01/scratch/microbiome/n10853499/00-rosella_testing/00-CAMI_I/00-low_complexity/binning/data/amber_out/benchmarks/genome/") {
    amber <- do.call(
        rbind,
        lapply(list.dirs(path, recursive=FALSE), function(x) {
            fread(paste0(x, '/metrics_per_bin.tsv'))
        })
    )

    amber[, binner:=ifelse(
        grepl("concoct", `Bin ID`),
        "CONCOCT",
        ifelse(
            grepl("_sens", `Bin ID`),
            "MetaBAT Sens.",
            ifelse(
                grepl("_spec", `Bin ID`),
                "MetaBAT Spec.",
                ifelse(
                    grepl("_ssens", `Bin ID`),
                    "MetaBAT SuperSens.",
                    ifelse(
                        grepl("_sspec", `Bin ID`),
                        "MetaBAT SuperSpec.",
                        ifelse(
                            grepl("metabat_bins_2", `Bin ID`),
                            "MetaBAT2",
                            ifelse(
                                grepl("maxbin2", `Bin ID`),
                                "MaxBin2",
                                ifelse(
                                    grepl("rosella_bins", `Bin ID`),
                                    "Rosella",
                                    ifelse(
                                        grepl("semibin_bins", `Bin ID`),
                                        "SemiBin",
                                        ifelse(
                                            grepl("vamb", `Bin ID`),
                                            "VAMB",
                                            ifelse(
                                                grepl("refine_rosella", `Bin ID`),
                                                "Rosella Refined",
                                                ifelse(
                                                    grepl("refine_metabat2", `Bin ID`),
                                                    "MetaBAT2 Refined",
                                                    ifelse(
                                                        grepl("refine_semibin", `Bin ID`),
                                                        "SemiBin Refined",
                                                        "GSA"
                                                    )
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )]

    amber[, 
        binner:=ifelse(
            grepl("das_tool_without", `Bin ID`),
            "DASTool w/o Rosella",
            ifelse(
                grepl("das_tool_bins_no_refine", `Bin ID`),
                "DASTool w/ Rosella",
                ifelse(
                    grepl("das_tool_bins_with_refine", `Bin ID`),
                    "DASTool w/ Refine",
                    binner
                )
            )
    )]

    amber[, `Contamination`:=round((1 - `Purity (bp)`) * 100, 2)]
    amber[, `Completeness`:=`Completeness (bp)` * 100]
    amber <- amber[binner!="GSA"]
    amber[, group:=ifelse(binner %like% "DASTool", "DASTool", ifelse(binner %like% "Refined", "Refined", "Independent"))]
    
    return(amber)
}

generate_plots_amber <- function(
    path = "/mnt/hpccs01/scratch/microbiome/n10853499/00-rosella_testing/00-CAMI_I/00-low_complexity/binning/data/amber_out/benchmarks/genome/",
    legend = TRUE,
    y_labels = TRUE,
    x_labels = TRUE,
    title = "CAMI I - Low complexity"
    ) {
    amber <- get_amber_files(path)
    
    generate_ranks_refined(amber, 5)
    # Had to change this to be static because VAMB doesn't show up sometimes which 
    # breaks the manual linetype values I set in the plots
    if ("VAMB" %in% amber$binner) {
        values = c(
            CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/ Refine`="#D99000", `DASTool w/o Rosella`="#A3A500", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", 
            `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", VAMB="#FF62BC", `MetaBAT2`="#E76BF3", `MetaBAT2 Refined`="#E76BF3", SemiBin="#880808", `SemiBin Refined`="#880808", 
            Rosella="#000000", `Rosella Refined`="#000000"
        )
        print(values)
    
        amber$binner = factor(amber$binner, levels=c('CONCOCT', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'DASTool w/ Refine', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.', 'VAMB', 'MetaBAT2', 'MetaBAT2 Refined', 'SemiBin', 'SemiBin Refined', 'Rosella', 'Rosella Refined'))
        linetypes = c(3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 1, 3, 1, 3, 1)
    } else {
        values = c(CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/ Refine`="#D99000", `DASTool w/o Rosella`="#A3A500", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", 
                   `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", `MetaBAT2`="#E76BF3", `MetaBAT2 Refined`="#E76BF3", SemiBin="#880808", `SemiBin Refined`="#880808",
                   Rosella="#000000", `Rosella Refined`="#000000")
        print(values)
    
        amber$binner = factor(amber$binner, levels=c('CONCOCT', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'DASTool w/ Refine', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.', 'MetaBAT2', 'MetaBAT2 Refined', 'SemiBin', 'SemiBin Refined', 'Rosella', 'Rosella Refined'))
        linetypes = c(3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 1, 3, 1, 3, 1)
    }
    
    amber$group = factor(amber$group, levels=c("Independent", "Refined", "DASTool"))
    
    com <- ggplot() + 
        geom_line(data=amber[Completeness >= 50 & Contamination <= 5,], aes(x=com_rank, y=Completeness, color=binner, linetype=binner)) +
        ylim(50, 100) +
        scale_colour_manual("Binner", values=values) +
        scale_linetype_manual("Binner", values = linetypes, guide="none") + 
        labs(x="Com. rank", color="Binner", y="", linetype="Binner") +
        guides(color = guide_legend(override.aes = list(color = values, linetype=linetypes))) +
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
            #   axis.text.y=ifelse(y_labels, element_text(), element_blank()),
            #   axis.text.x=ifelse(x_labels, element_text(), element_blank()),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position=ifelse(legend, "bottom", "none"),
              legend.direction="horizontal", 
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    con <- ggplot() + 
        geom_line(data=amber[Completeness >= 50 & Contamination <= 5,], aes(x=con_rank, y=Contamination, color=binner, linetype=binner)) +   
        ylim(0, 5) +
        scale_colour_manual("", values=values) +
        scale_linetype_manual("", values = linetypes, guide="none") + 
        labs(x="Con. rank", color="Binner", y="", linetype="") +
        guides(color = guide_legend(override.aes = list(color = values, linetype=linetypes))) +
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
            #   axis.text.y=ifelse(y_labels, element_text(), element_blank()),
            #   axis.text.x=ifelse(x_labels, element_text(), element_blank()),
              axis.line = element_line(size=0.25),
              axis.ticks=element_line(size=0.25),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position=ifelse(legend, "bottom", "none"),
              legend.direction="horizontal", 
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    amber[, mag_group:=ifelse(Completeness >= 90 & Contamination <= 5, "90%", 
                              ifelse(Completeness >= 80 & Contamination <= 5, "80%", 
                                     ifelse(Completeness >= 70 & Contamination <= 5, "70%", 
                                            ifelse(Completeness >= 60 & Contamination <= 5, "60%", 
                                                   ifelse(Completeness >= 50 & Contamination <= 5, "50%", "<50%")))))]

    # amber_low[, .N, by=c("group", "binner", "mag_group")]

    bar_chart <- ggplot(data=amber[mag_group != "<50%", .N, by=c("group", "binner", "mag_group")], aes(fill=mag_group, y=N, x=binner)) +
        geom_bar(position="stack", stat="identity") + 
        labs(title=title, fill="Completeness", x = "", y = "MAGs") +
    #     scale_fill_viridis(discrete=TRUE, option="mako")
        scale_fill_brewer() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position=ifelse(legend, "right", "none"),
            legend.direction="vertical", 
            legend.text=element_text(size=10),
            axis.text.y=element_text(size=ifelse(y_labels, 10, 0)),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1),
            plot.title=element_text(hjust=0.5)
            ) +
        facet_wrap(~group, nrow=3, ncol=1, scales="free_y") +
        coord_flip()
    
    gt = ggplot_gtable(ggplot_build(bar_chart))
    gt$heights[8] = 2*gt$heights[8]
    gt$heights[13] = 0.5*gt$heights[13]
    gt$heights[18] = 0.5*gt$heights[18]

    t_all <- get_thresholds(amber, check_vamb=TRUE)
    
    return(list(com, con, gt, t_all, amber))
}

combine_amber_results <- function(
    path="/mnt/hpccs01/scratch/microbiome/n10853499/00-rosella_testing/01-CAMI_II/CAMI_uro/binning/", 
    regexp="amber_sample_*", 
    inner_directory="/data/amber_out/benchmarks/genome/",
    legend = TRUE,
    y_labels = TRUE,
    x_labels = TRUE,
    title = "CAMI II - Urogenital"
    ) {
    results <- Sys.glob(paste0(path, regexp))

    amber <- do.call(rbind, lapply(results, function(x) {
        result <- get_amber_files(paste0(x, "/", inner_directory))
        return(result)
    }))

    generate_ranks_refined(amber, 5)

    if ("VAMB" %in% amber$binner) {
        values = c(
            CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/ Refine`="#D99000", `DASTool w/o Rosella`="#A3A500", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", 
            `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", VAMB="#FF62BC", `MetaBAT2`="#E76BF3", `MetaBAT2 Refined`="#E76BF3", SemiBin="#880808", `SemiBin Refined`="#880808", 
            Rosella="#000000", `Rosella Refined`="#000000"
        )
        print(values)
    
        amber$binner = factor(amber$binner, levels=c('CONCOCT', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'DASTool w/ Refine', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.', 'VAMB', 'MetaBAT2', 'MetaBAT2 Refined', 'SemiBin', 'SemiBin Refined', 'Rosella', 'Rosella Refined'))
        linetypes = c(3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 1, 3, 1, 3, 1)
    } else {
        values = c(CONCOCT="#F8766D",  `DASTool w/ Rosella`="#D89000", `DASTool w/ Refine`="#D99000", `DASTool w/o Rosella`="#A3A500", MaxBin2="#39B600", `MetaBAT Sens.`="#00BF7D", `MetaBAT Spec.`="#00BFC4", 
                   `MetaBAT SuperSens.`="#00B0F6", `MetaBAT SuperSpec.`="#9590FF", `MetaBAT2`="#E76BF3", `MetaBAT2 Refined`="#E76BF3", SemiBin="#880808", `SemiBin Refined`="#880808",
                   Rosella="#000000", `Rosella Refined`="#000000")
        print(values)
    
        amber$binner = factor(amber$binner, levels=c('CONCOCT', 'DASTool w/o Rosella', 'DASTool w/ Rosella', 'DASTool w/ Refine', 'MaxBin2', 'MetaBAT Sens.', 'MetaBAT Spec.', 'MetaBAT SuperSens.', 'MetaBAT SuperSpec.', 'MetaBAT2', 'MetaBAT2 Refined', 'SemiBin', 'SemiBin Refined', 'Rosella', 'Rosella Refined'))
        linetypes = c(3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 1, 3, 1, 3, 1)
    }
    
    amber$group = factor(amber$group, levels=c("Independent", "Refined", "DASTool"))
    
    com <- ggplot() + 
        geom_line(data=amber[Completeness >= 50 & Contamination <= 5,], aes(x=com_rank, y=Completeness, color=binner, linetype=binner)) +
        ylim(50, 100) +
        scale_colour_manual("Binner", values=values) +
        scale_linetype_manual("Binner", values = linetypes, guide="none") + 
        labs(x="Com. rank", color="Binner", y="", linetype="Binner") +
        guides(color = guide_legend(override.aes = list(color = values, linetype=linetypes))) +
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
              legend.title=element_blank(),
              legend.text=element_text(size=10), 
              legend.background=element_blank(), 
              legend.key=element_blank())

    con <- ggplot() + 
        geom_line(data=amber[Completeness >= 50 & Contamination <= 5,], aes(x=con_rank, y=Contamination, color=binner, linetype=binner)) +   
        ylim(0, 5) +
        scale_colour_manual("", values=values) +
        scale_linetype_manual("", values = linetypes, guide="none") + 
        labs(x="Con. rank", color="Binner", y="", linetype="") +
        guides(color = guide_legend(override.aes = list(color = values, linetype=linetypes))) +
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

    amber[, mag_group:=ifelse(Completeness >= 90 & Contamination <= 5, "90%", 
                              ifelse(Completeness >= 80 & Contamination <= 5, "80%", 
                                     ifelse(Completeness >= 70 & Contamination <= 5, "70%", 
                                            ifelse(Completeness >= 60 & Contamination <= 5, "60%", 
                                                   ifelse(Completeness >= 50 & Contamination <= 5, "50%", "<50%")))))]

    # amber_low[, .N, by=c("group", "binner", "mag_group")]

    bar_chart <- ggplot(data=amber[mag_group != "<50%", .N, by=c("group", "binner", "mag_group")], aes(fill=mag_group, y=N, x=binner)) +
        geom_bar(position="stack", stat="identity") + 
        labs(title=title, fill="Completeness", x = "", y = "MAGs") +
    #     scale_fill_viridis(discrete=TRUE, option="mako")
        scale_fill_brewer() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position=ifelse(legend, "right", "none"),
              legend.direction="vertical", 
              legend.text=element_text(size=10),
              axis.text.y=element_text(size=ifelse(y_labels, 10, 0)),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1),
              plot.title=element_text(hjust=0.5)
            ) +
        facet_wrap(~group, nrow=3, ncol=1, scales="free_y") +
        coord_flip()
    
    gt = ggplot_gtable(ggplot_build(bar_chart))
    gt$heights[8] = 2*gt$heights[8]
    gt$heights[13] = 0.5*gt$heights[13]
    gt$heights[18] = 0.5*gt$heights[18]

    t_all <- get_thresholds(amber, check_vamb=TRUE)
    
    return(list(com, con, gt, t_all, amber))
}

get_benchmarks <- function(benchmark_folder) {
    rosella_b <- fread(paste0(benchmark_folder, "/rosella.benchmark.txt"))
    rosella_b$binner <- "Rosella"
    semibin_b <- fread(paste0(benchmark_folder, "/semibin.benchmark.txt"))
    semibin_b$binner <- "SemiBin"
    metabat2_b <- fread(paste0(benchmark_folder, "/metabat_2.benchmark.txt"))
    metabat2_b$binner <- "MetaBAT2"
    maxbin2_b <- fread(paste0(benchmark_folder, "/maxbin2.benchmark.txt"))
    maxbin2_b$binner <- "MaxBin2"
    metabat_sens_b <- fread(paste0(benchmark_folder, "/metabat_sens.benchmark.txt"))
    metabat_sens_b$binner <- "MetaBAT Sens."
    metabat_ssens_b <- fread(paste0(benchmark_folder, "/metabat_ssens.benchmark.txt"))
    metabat_ssens_b$binner <- "MetaBAT SuperSens."
    metabat_spec_b <- fread(paste0(benchmark_folder, "/metabat_spec.benchmark.txt"))
    metabat_spec_b$binner <- "MetaBAT Spec."
    metabat_sspec_b <- fread(paste0(benchmark_folder, "/metabat_sspec.benchmark.txt"))
    metabat_sspec_b$binner <- "MetaBAT SuperSpec."
    concoct_b <- fread(paste0(benchmark_folder, "/concoct.benchmark.txt"))
    concoct_b$binner <- "CONCOCT"
    vamb_b <- fread(paste0(benchmark_folder, "/vamb.benchmark.txt"))
    vamb_b$binner <- "VAMB"

    # refine_rosella <- fread(paste0(benchmark_folder, "/rosella_refine_rosella.benchmark.txt"))
    # refine_rosella$binner <- "Rosella Refined"
    # refine_semibin <- fread(paste0(benchmark_folder, "/rosella_refine_semibin.benchmark.txt"))
    # refine_semibin$binner <- "SemiBin Refined"
    # refine_metabat2 <- fread(paste0(benchmark_folder, "/rosella_refine_metabat2.benchmark.txt"))
    # refine_metabat2$binner <- "MetaBAT2 Refined"
    # refine_dastool <- fread(paste0(benchmark_folder, "/rosella_refine_dastool.benchmark.txt"))
    # refine_dastool$binner <- "DASTool Refined"
    
    das_tool_with_refine <- fread(paste0(benchmark_folder, '/das_tool_refine.benchmark.txt'))
    das_tool_with_refine$binner <- "DASTool w/ Refine"
    das_tool_no_rosella <- fread(paste0(benchmark_folder, '/das_tool_no_rosella.benchmark.txt'))
    das_tool_no_rosella$binner <- "DASTool w/o Rosella"
    das_tool_no_refine <- fread(paste0(benchmark_folder, '/das_tool_no_refine.benchmark.txt'))
    das_tool_no_refine$binner <- "DASTool w/ Rosella"

    bound <- rbind(rosella_b, semibin_b, metabat2_b, maxbin2_b, metabat_sens_b, metabat_ssens_b, metabat_spec_b, 
    metabat_sspec_b, concoct_b, vamb_b, das_tool_with_refine, das_tool_no_rosella, das_tool_no_refine)

    bound[, m:=s/60]
    bound[, h:=m/60]
    bound[, max_rss:=max_rss/1024]
    bound[, io_out:=io_out/1024]
    bound[, io_in:=io_in/1024]

    return(bound)
}