#' Functions to plot the marker expressions per cell population
#' and display differences based on the previous analyses.

### Create frequency file with only combinations of marker-cell pop that we want to plot
create_expr_file <- function(settings){
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
        zoom_folder <- settings$zoom_pop
    } else {
        zoom_folder <- NULL
    }
    tissue <- settings$tissue
    
    ms <- fread(paste0(settings$experiment_dir, "/", tissue, "/", 
                        subset_folder, normalised_folder, zoom_folder, "/", tissue,
                        "_ds_res_ms_results.csv"))
    ad <- fread(paste0(settings$experiment_dir, "/", tissue, "/", 
                        subset_folder, normalised_folder, zoom_folder, "/", tissue,
                        "_ds_res_ad_results.csv"))
    adms <- fread(paste0(settings$experiment_dir, "/", tissue, "/", 
                          subset_folder, normalised_folder, zoom_folder, "/", tissue,
                          "_ds_res_advsms_results.csv"))
    
    expression <- rbind(ms, ad, adms) # Pool the comparisons together
    expr <- expression
    # expr <- unique(expression, by = "unique") # Keep only unique combis
    # expr <- expr[p_val < 0.1]
    
    ### Melt table
    ids <- colnames(expr)[6:27]
    expr <- melt(expr, id.vars = c("unique", "cluster_id", "marker_id", "p_val"), 
                 measure.vars = ids, 
                 value.name = "median_expr", 
                 variable.name = "sample_id")
    
    # Add condition
    expr[, condition := ifelse(grepl("control", sample_id), "control",
                               ifelse(grepl("ms", sample_id), "ms", 
                                      ifelse(grepl("ad", sample_id), "ad", NA)))]
    
    # Keep only significant (before correction)
    expr_for_plot <- expr[p_val < 0.1]
    expression <- unique(expr, by = c("unique", "sample_id")) # Keep only unique combis

    
    p_vals_table <- rbindlist(list(
        ms[, comparison := "ms_con"],
        ad[, comparison := "ad_con"],
        adms[, comparison := "ad_ms"]))
    
    return(list(expr_for_plot = expr_for_plot, expression = expression, p_vals_table = p_vals_table))
}

### Plot
plot_expr <- function(expr, p_vals_table, settings){
    colors = c("#043741", "#189cb3", "#e79d24")
    linecol = "gray30"
    black = c("black", "black", "black")
    
    # Plot
    markers <- unique(expr$marker_id)
    for (m in markers) {
        ex <- expr[marker_id == m] # keep only one marker
        ex$cluster_id <- reorder(ex$cluster_id, -ex$median_expr) # re-order by abundance
        ex_m <- ex[, 
                   .(median_expr = mean(median_expr, na.rm = T)), 
                   by = c("condition", "marker_id", "cluster_id")]
        ex$condition <- factor(ex$condition, levels = c("control", "ms", "ad"))
        ex_m$condition <- factor(ex_m$condition, levels = c("control", "ms", "ad"))
        w <- length(unique(ex$cluster_id))
        print(m)
        
        ex$condition <- factor(ex$condition, levels = c("control", "ad", "ms"))
        p <- ggplot(ex, aes(x = cluster_id, y = median_expr, fill = condition)) +
            geom_boxplot(width = 0.75, color = linecol, lwd = 0.25, 
                         outlier.shape = NA, alpha = 0.5) +
            geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                                       jitter.height = 0,
                                                       dodge.width = 0.75), 
                       size = 1, shape = 21, stroke = 0.1,
                       aes(color = condition, fill = condition)
                       # fill = "white"
            ) +
            xlab(NULL) +
            ylab(paste0(m, "expression")) + # bars show mean, dots are individual median expression
            expand_limits(y=0) +
            scale_fill_manual(values=colors, labels = c('Control','Dementia','MS'), name = NULL) +
            scale_color_manual(values=colors, labels = c('Control','Dementia','MS'), name = NULL) +
            theme_classic() + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)
                  , text=element_text(size=7)
                  # , legend.position = "none"
            )
       
        max_height <- 0
        for (cluster in ex[, unique(cluster_id)]) {
            
            height = ex[cluster_id == cluster, max(median_expr, na.rm = T)]
            for (comp in c('ad_con', 'ad_ms', 'ms_con')) {
                cond_a <- unlist(str_split(comp, "_"))[1]
                cond_b <- unlist(str_split(comp, "_"))[2]
                if (cond_a == "con") {cond_a <- "control"}
                if (cond_b == "con") {cond_b <- "control"}
                
                cond_pos <- data.table(cond = c("control", "ad", "ms"), pos = c(-0.33, 0, +0.33))
                pos_a <- cond_pos[cond == cond_a, pos]
                pos_b <- cond_pos[cond == cond_b, pos]
                
                p_label <- data.table(ref_p = c(0.1, 0.05, 0.01), label = c("*", "**", "***"))
                my_p_value <- p_vals_table[comparison == comp & cluster_id == cluster & marker_id == m, p_adj]
                label <- p_label[my_p_value <= ref_p][.N][, label]
                if (length(label) == 0) {
                    print(paste("Skipping", m, comp, cluster))
                    next
                } else {
                    height <- height * 1.1
                }
                
                p <- p + geom_segment(
                    x=as.numeric(factor(cluster, levels = levels(ex$cluster_id))) + pos_a,
                    y=height * 1.1,
                    xend=as.numeric(factor(cluster, levels = levels(ex$cluster_id))) + pos_b,
                    yend=height * 1.1
                ) + geom_text(
                    x = (as.numeric(factor(cluster, levels = levels(ex$cluster_id))) + pos_a + 
                             as.numeric(factor(cluster, levels = levels(ex$cluster_id))) + pos_b) / 2,
                    y = height * 1.1 + 0.02,
                    label = label
                )
            }
            max_height <- max(max_height, height)
        }
        p <- p + scale_y_continuous(expand = c(0,0), limits = c(0, max_height*1.5))
        
        print(p)
    }
}

