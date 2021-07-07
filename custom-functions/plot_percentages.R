#' Functions to plot the percentage of cell populations
#' and display differences based on the previous analyses.

### Create frequency file
create_freq_file <- function(abs_props = T, settings){
    if (settings$norm) {
        normalised_folder <- "normalised/"
    } else {
        normalised_folder <- "original/"
    }
    if (settings$subsetting) {
        subset_folder <- "subset/"
    } else {
        subset_folder <- "whole/"
    }
    if (settings$zoom) {
        zoom_folder <- paste0(settings$zoom_pop, "/")
        if (abs_props == T) {
            name = "_total_prop_table.csv"
        } else {
            name = "_prop_table.csv"
        }
    } else {
        zoom_folder <- NULL
        name = "_prop_table.csv"
    }
    tissue <- settings$tissue
    
    freqs <- fread(paste0("./", tissue, "/", 
                          subset_folder, normalised_folder, zoom_folder, tissue, name))
    
    # Add condition
    freqs[, condition := ifelse(grepl("control", sample_id), "control",
                                ifelse(grepl("ms", sample_id), "ms", 
                                       ifelse(grepl("ad", sample_id), "ad", NA)))]
    
    freqs[, perc := prop*100]
    
    # Summarise medians
    med_sum <- freqs[,
                .(fq_median = median(prop)),
                by = c("pop", "condition")]
    
    # Keep only significant (before correction)
    ms <- fread(paste0("./", tissue, "/", 
                       subset_folder, normalised_folder, zoom_folder, "/", tissue,
                       "_topTable_da_res_ms.csv"))
    ad <- fread(paste0("./", tissue, "/", 
                       subset_folder, normalised_folder, zoom_folder, "/", tissue,
                       "_topTable_da_res_ad.csv"))
    adms <- fread(paste0("./", tissue, "/", 
                         subset_folder, normalised_folder, zoom_folder, "/", tissue,
                         "_topTable_da_res_advsms.csv"))
    
    p_vals_table <- rbindlist(list(
        ms[, comparison := "ms_con"],
        ad[, comparison := "ad_con"],
        adms[, comparison := "ad_ms"]))
    
    return(list(freqs = freqs, med_sum = med_sum, p_vals_table = p_vals_table))
}

### Organise freqs for plots
organise_freqs <- function(freqs, settings) {
    # Re-order from more to less abundant, for a cleaner plot
    freqs$pop <- reorder(freqs$pop, -freqs$perc)

    # Split for visibility
    if (settings$tissue == "cp") {
        abundant <- freqs[pop %in% c("Myeloid", "Granulocytes", "Naive B cells",
                                     "NK cells", "B cell lineage",  "CD4+ T cells", "CD8+ T cells",  
                                     "CD8 c3", "DN", "CD8 c4", "CD8 c1",  
                                     "NK CD56dim")]
        mid_ab <- freqs[pop %in% c("Monocytes", 
                                    "CD8 effector memory", "CD8 naive", "CD8 Temra", "CD8 Temra resident",
                                   "CD8 c4", "CD8 c3")]
        low_ab <- freqs[pop %in% c("B cell lineage", "gd T cells",
                                   "ASCs", "B cells", "NK CD56dim", "NK CD56bright",
                                   "CD8 transitional RA/RO", "CD4 effector memory PD1-",
                                   "DN", "CD4 naive", "gd T cells",
                                   "CD8 c5", "DN", "CD4 c3", "gd T"
                                   )] # excluding "unknown'
        non_ab <- freqs[!pop %in% c("Myeloid", "Granulocytes", "Naive B cells",
                                    "NK cells", "B cell lineage",  "CD4+ T cells", "CD8+ T cells",  
                                    "CD8 c3", "DN", "CD8 c4", "CD8 c1",  
                                    "NK CD56dim")]
    }
    if (settings$tissue == "septum") {
        abundant <- freqs[pop %in% c("Microglia", "Myeloid",
                                     "CD8 c3", "CD4 c2", "CD4 c1", "CD8 c1", "CD8 c2")]
        mid_ab <-  freqs[pop %in% c("CD8+ T cells", "Granulocytes", "CD4+ T cells", 
                                    "CD4 effector memory", "CD8 Temra",
                                    "CD8 effector memory",
                                    "CD8 c3", "CD8 c4", "CD4 c2"
                                    )]
        low_ab <- freqs[!(pop %in% c("Microglia", "Myeloid", "Undefined", "CD45-",
                                     "CD8 Trm", "CD4 Trm", "CD8 Trm CD103+", 
                                     "CD4 effector memory", "CD8 Temra", "CD8 effector memory",
                                     "CD8+ T cells", "Granulocytes", "CD4+ T cells",
                                     "CD8 c1", "CD8 c2", "CD4 c1", "CD8 c3", "CD8 c4", "CD4 c2"))]
        non_ab <- freqs[!pop %in% c("Microglia", "Myeloid",
                                    "CD8 c3", "CD4 c2", "CD4 c1", "CD8 c1", "CD8 c2")]
    }
    if (settings$tissue == "blood") {
        abundant <- freqs[pop %in% c("Granulocytes", "NK cells")]
        mid_ab <- freqs[pop %in% c("CD4+ T cells", "CD8+ T cells", "Monocytes")]
        low_ab <- freqs[pop %in% c("mDCs", "pDCs", "B cells", "TCRgd T cells")]
    }
    
    freqs$condition <- factor(freqs$condition, levels = c("control", "ms", "ad")) 
    
    # Factor so the plot recognises the position per group of abundance:
    abundant$pop <- factor(abundant$pop)
    if (settings$tissue != "cp") {
        mid_ab$pop <- factor(mid_ab$pop)
    }
    low_ab$pop <- factor(low_ab$pop)
    
    
    return(list(freqs = freqs, abundant = abundant, 
                mid_ab = mid_ab, low_ab = low_ab, non_ab = non_ab))
}


### Create function for plots
plot_perc <- function(populations, p_vals_table, x = pop, y = perc, fill = condition, 
                      settings, abs_prop = abs_prop) {
    if (settings$zoom) {
        if (abs_prop == T) {
            main_pop <- "CD45+"
        } else {
            main_pop <- settings$zoom_pop   
        }
    } else {
        main_pop <- "CD45+"
    }
    # Plot colors
    colors = c("#043741", "#189cb3", "#e79d24")
    linecol = "gray30"
    black = c("black", "black", "black")
    #Create plot
    populations$condition <- factor(populations$condition, levels = c("control", "ad", "ms"))
    p <- ggplot(populations, aes(x = pop, y = perc, fill = condition)) +
        geom_boxplot(width = 0.75, color = linecol, lwd = 0.25, 
                     outlier.shape = NA, alpha = 0.5) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                      jitter.height = 0,
                                      dodge.width = 0.75), 
                   # position = position_dodge(width=0.5),
                   # position = position_jitter(w = 0.15, h = 0),
                   size = 1, shape = 21, stroke = 0.1,
                   # color = "black",
                   aes(color = condition, fill = condition)
                   # fill = "white"
                   # color = "black"
                   ) +
        xlab(NULL) +
        ylab(paste0("Percentage of ", main_pop ," cells (%)")) +
        expand_limits(y=0) +
        # scale_y_continuous(expand = c(0,0)) +
        scale_fill_manual(values=colors, labels = c('Control','Dementia','MS'), name = NULL) +
        scale_color_manual(values=colors, labels = c('Control','Dementia','MS'), name = NULL) +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 1)
              , legend.position = "none"
        )
    max_height <- 0
    for (cluster in populations[, unique(pop)]) {
        
        height <- populations[pop == cluster, max(perc, na.rm = T)]
        for (comp in c('ad_con', 'ad_ms', 'ms_con')) {
            cond_a <- unlist(str_split(comp, "_"))[1]
            cond_b <- unlist(str_split(comp, "_"))[2]
            if (cond_a == "con") {cond_a <- "control"}
            if (cond_b == "con") {cond_b <- "control"}
            
            cond_pos <- data.table(cond = c("control", "ad", "ms"), pos = c(-0.33, 0, +0.33))
            pos_a <- cond_pos[cond == cond_a, pos]
            pos_b <- cond_pos[cond == cond_b, pos]
            
            p_label <- data.table(ref_p = c(0.1, 0.05, 0.01), label = c("*", "**", "***"))
            my_p_value <- p_vals_table[comparison == comp & cluster_id == cluster, p_adj]
            label <- p_label[my_p_value <= ref_p][.N][, label]
            if (length(label) == 0) {
                print(paste("Skipping", comp, cluster))
                next
            } else {
                height <- height * 1.1
            }
            
            p <- p + geom_segment(
                x=as.numeric(factor(cluster, levels = levels(as.factor(populations$pop)))) + pos_a,
                y=height * 1.1,
                xend=as.numeric(factor(cluster, levels = levels(as.factor(populations$pop)))) + pos_b,
                yend=height * 1.1
            ) + geom_text(
                x = (as.numeric(factor(cluster, levels = levels(as.factor(populations$pop)))) + pos_a + 
                         as.numeric(factor(cluster, levels = levels(as.factor(populations$pop)))) + pos_b) / 2,
                y = height * 1.1 + 0.02,
                label = label
            )
        }
        max_height <- max(max_height, height)
    }
    p <- p + scale_y_continuous(expand = c(0,0), limits = c(0, max_height*1.2))
    
    print(p)
}


### Create function for separate plots per pop
plot_perc_sep <- function(freqs, p_vals_table, x = pop, y = perc, fill = condition, 
                      settings, abs_prop = abs, accuracy = 0.01) {
    if (settings$zoom) {
        if (abs_prop == T) {
            main_pop <- "CD45+"
        } else {
            main_pop <- settings$zoom_pop   
        }
    } else {
        main_pop <- "CD45+"
    }
    # Plot colors
    colors = c("#043741", "#189cb3", "#e79d24")
    linecol = "gray30"
    black = c("black", "black", "black")
    #Create plot
    populations <- unique(freqs$pop)
    freqs$condition <- factor(freqs$condition, levels = c("control", "ad", "ms"))
    for (population in populations) {
        p <- ggplot(freqs[pop == population], aes(x = pop, y = perc, fill = condition)) +
            geom_boxplot(width = 0.75, color = linecol, lwd = 0.25, 
                         outlier.shape = NA, alpha = 0.5) +
            geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                                       jitter.height = 0,
                                                       dodge.width = 0.75), 
                       # position = position_dodge(width=0.5),
                       # position = position_jitter(w = 0.15, h = 0),
                       size = 1, shape = 21, stroke = 0.1,
                       # color = "black",
                       aes(color = condition, fill = condition)
                       # fill = "white"
                       # color = "black"
            ) +
            xlab(NULL) +
            ylab(paste0("Percentage of ", main_pop ," cells (%)")) +
            expand_limits(y=0) +
            # scale_y_continuous(expand = c(0,0)) +
            scale_fill_manual(values=colors, labels = c('Control','Dementia','MS'), name = NULL) +
            scale_color_manual(values=colors, labels = c('Control','Dementia','MS'), name = NULL) +
            theme_classic() + 
            theme(#axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none"
            )
    
    max_height <- 0
    height <- freqs[pop == population, max(perc, na.rm = T)]
    for (comp in c('ad_con', 'ad_ms', 'ms_con')) {
        cond_a <- unlist(str_split(comp, "_"))[1]
        cond_b <- unlist(str_split(comp, "_"))[2]
        if (cond_a == "con") {cond_a <- "control"}
        if (cond_b == "con") {cond_b <- "control"}
            
        cond_pos <- data.table(cond = c("control", "ad", "ms"), pos = c(0.75, 1, 1.25))
        pos_a <- cond_pos[cond == cond_a, pos]
        pos_b <- cond_pos[cond == cond_b, pos]
            
        p_label <- data.table(ref_p = c(0.1, 0.05, 0.01), label = c("*", "**", "***"))
        my_p_value <- p_vals_table[comparison == comp & cluster_id == population, p_adj]
        label <- p_label[my_p_value <= ref_p][.N][, label]
        if (length(label) == 0) {
            print(paste("Skipping", comp, population))
            next
        } else {
            height <- height * 1.1
        }
            
        p <- p + geom_segment(
            x= pos_a,
            y=height * 1.1,
            xend=pos_b,
            yend=height * 1.1
        ) + geom_text(
            x =  as.numeric(pos_a + pos_b) / 2,
            y = height * 1.1 + 0.02,
            label = label
        )
        }
        max_height <- max(max_height, height)
        if (population %in% abundant) {
            p <- p + scale_y_continuous(expand = c(0,0), limits = c(0, max_height*1.2), labels = scales::number_format(accuracy = accuracy))
        } else {
            p <- p + scale_y_continuous(expand = c(0,0), limits = c(0, max_height*1.2), labels = scales::number_format(accuracy = accuracy))

        }
        
        set_dev(plot_path(settings, paste0("perc_", population)), height=2, width = 1.5)
        print(p)
        unset_dev()
        
    }
    
}



