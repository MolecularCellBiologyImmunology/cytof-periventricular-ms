#' Exploring and visualising the data for diagnostics 

plot_marker_expression_distributions <- function(sce, settings, do_log = FALSE) {
    # Plot with per-sample marker expression distributions, colored by condition
    
    plot_name <- paste0("distributions", ifelse(do_log, "_xlog2", ""))
    
    set_dev(plot_path(settings, plot_name), width = 8)
    p <- plotExprs(sce, color_by = "condition")
    p$facet$params$ncol <- 6
    
    if (do_log) {
        p <- p + 
            scale_x_continuous(trans = 'log2', labels = scales::scientific) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    print(p)
    unset_dev()
    
    return(p)
}


plot_n_cells_per_sample <- function(sce, settings) {
    # Check number of cells per sample
    n_cells(sce) # or, equivalently, `metadata(sce)$experiment_info$n_cells`
    counts <- as.data.table(ei(sce))
    counts[, sample_id := as.character(sample_id)]
    counts[, condition := factor(condition, levels = c("control", "ad", "ms"))]
    counts[, sample_id := factor(sample_id, levels = sample_id[order(condition, sample_id)])]
    ### Plot
    colors = c("#043741", "#189cb3", "#e79d24")
    
    set_dev(plot_path(settings, "n_cells"), height = 5, width = 5)
    p <- ggplot(counts, aes(x = sample_id, y = n_cells, fill = condition)) +
        geom_col() +
        # geom_text(aes(label = n_cells, vjust = -0.5)) +
        # geom_hline(yintercept = 111065) +
        scale_fill_manual(values=colors, labels = c('NC','NNC','MS'), name = NULL) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
    unset_dev()
    
    return(p)
}


plot_mds <- function(sce, settings) {
    # MDS plot (similar to PCA)
    set_dev(plot_path(settings, "MDS"))
    p <- CATALYST::pbMDS(sce, color_by = "condition") # or plotMDS(sce)
    print(p)
    unset_dev()
    
    return(p)
}


plot_heatmap <- function(sce, settings) {
    # Heatmap
    set_dev(plot_path(settings, "heatmap"))
    p <- plotExprHeatmap(sce, bin_anno = FALSE, row_anno = FALSE)
    print(p)
    unset_dev()
    
    return(p)
}



plot_redundancy_score <- function(sce, settings) {
    ### MARKER RANKING BASED ON THE NON-REDUNDANCY SCORE ###
    set_dev(plot_path(SETTINGS, "redundancy"))
    p <- plotNRS(sce, features = type_markers(sce), color_by = "condition")
    print(p)
    unset_dev()
    
    return(p)
}
