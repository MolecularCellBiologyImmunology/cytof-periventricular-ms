#' Uses `cluster()` function from CATALYST to perform `FlowSOM` clustering 
#' & `ConsensusClusterPlus` metaclustering.

run_clustering <- function(sce, settings, skip_checkpoint = FALSE, features = type_markers(sce)) {
    ### Use CATALYST function cluster() to build SOM with FlowSOM and do
    ### metaclustering with ConsensusClusterPlus()
    # Check if RDS already there, open if so
    checkpoint_path <- object_path(settings, "sce_with_clustering.RDS")
    if (!skip_checkpoint && file.exists(checkpoint_path)) {
        cat(paste0("Loading result from checkpoint ", checkpoint_path))
        sce <- readRDS(checkpoint_path)
        return(sce)
    }
    
    # Create new folder and RDS if there's none from before
    dir.create(object_path(settings, NULL))
    sce <- CATALYST::cluster(sce, features = features,
                             xdim = 10, ydim = 10, maxK = 20, seed = 1234)
    
    saveRDS(sce, checkpoint_path)
    
    # max n of clusters as an overestimation
    return(sce)
}

plot_elbow_plot <- function(sce, settings) {
    # Delta plot or elbow plot to find the optimal number of clusters 
    set_dev(plot_path(settings, "elbow_plot"), height=10, width=5)
    p <- metadata(sce)$delta_area
    print(p)
    unset_dev()
    
    return(p)
}

plot_heatmap_cellpops <- function(sce, settings, clusters, features = "type", bin_anno = F, 
                                  name = "heatmap_cellpops") {
    ### Heatmap
    set_dev(plot_path(settings, name))
    # p <- plotExprHeatmap(sce, #hm2 = NULL, cluster_anno = TRUE, draw_freqs = TRUE,
    #                      k = clusters, m = NULL)
    p <- plotExprHeatmap(sce, features = features, fun = "median",
                    by = "cluster_id", k = clusters, 
                    bars = TRUE, perc = TRUE, bin_anno = bin_anno)
    
    print(p)
    unset_dev()
    
    return(p)
}

plot_distribution_markers <- function(sce, clusters, settings) {
    # Distributions of marker intensities
    set_dev(plot_path(settings, "distrib_cellpops"))
    p <- plotClusterExprs(sce, k = clusters, features = "type")
    print(p)
    unset_dev()
    
    return(p)
}


plot_specific_marker_heatmap <- function(sce, markers, clusters, settings) {
    # Visualize median marker expression in the 100 SOM codes
    p <- plotClusterHeatmap(sce,
                            hm2 = markers, k = "som100", m = clusters,
                            cluster_anno = FALSE, draw_freqs = TRUE)
    
    return(p)
}

run_tsne_and_umap <- function(sce, settings, skip_checkpoint = FALSE) {
    
    checkpoint_path <- object_path(settings, "sce_with_clustering_with_tsne_and_umap.RDS")
    if (!skip_checkpoint && file.exists(checkpoint_path)) {
        cat(paste0("Loading result from checkpoint ", checkpoint_path))
        sce <- readRDS(checkpoint_path)
        return(sce)
    }
    
    # Run t-SNE/UMAP on at most 500/1000 cells per sample
    sce <- runDR(sce, dr = "TSNE", cells = 500, features = "type")
    sce <- runDR(sce, dr = "UMAP", cells = 1e3, features = "type")
    
    saveRDS(sce, checkpoint_path)
    
    return(sce)
}

plot_tsne_and_umap <- function(sce, clusters, settings) {
    # Compare tSNE to uMAP
    p1 <- plotDR(sce, "TSNE", color_by = clusters) +
        theme(legend.position = "none")
    p2 <- plotDR(sce, "UMAP", color_by = clusters)
    lgd <- get_legend(p2 + guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))))
    p2 <- p2 + theme(legend.position = "none")
    
    set_dev(plot_path(settings, "tsne_umap"))
    p <- plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))
    print(p)
    unset_dev()
    
    return(p)
}

plot_umap_per_sample <- function(sce, clusters, settings) {
    # Facet per sample
    set_dev(plot_path(settings, "umap_sample"))
    p <- plotDR(sce, "UMAP", color_by = clusters) + 
        facet_wrap("sample_id") +
        guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
    print(p)
    unset_dev()
    
    return(p)
}

plot_umap_per_condition <- function(sce, clusters, settings) {
    # Facet per condition
    set_dev(plot_path(settings, "umap_condition"))
    p <- plotDR(sce, "UMAP", color_by = clusters) + 
        facet_wrap("condition") +
        guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
    print(p)
    unset_dev()
    
    return(p)
}

plot_umap_markers <- function(sce, features = type_markers(sce), settings) {
    set_dev(plot_path(settings, "umap_per_marker"))
    p <- plotDR(sce, "UMAP", color_by = features, ncol = 7, a_pal = rev(brewer.pal(11, "RdYlBu"))) # UMAP
    print(p)
    unset_dev()
    
    return(p)
}

apply_manual_merging <- function(sce, clusters, settings, name = "CP_CyTOF-cluster-merging-1.csv") {
    merge1 <- fread(object_path(settings, name))[, 1:2]
    
    # convert to factor with merged clusters in desired order
    cluster_order <- unique(merge1$new_cluster)
    merge1$new_cluster <- factor(merge1$new_cluster,
                                 levels = cluster_order)
    merge1$new_cluster <- factor(merge1$new_cluster)
    
    # apply manual merging
    sce <- mergeClusters(sce, k = clusters,
                         table = merge1, id = "merge1")
    
    return(sce)
}

plot_umap_merged <- function(sce, settings, name = "umap_merged_clusters") {
    set_dev(plot_path(settings, name = name), height=5, width=5)
    p <- plotDR(sce, "UMAP", color_by = "merge1") # UMAP
    print(p)
    unset_dev()
    
    return(p)
}

plot_heatmap_cellpops_2layers <- function(sce, clusters, settings) {
    ### Heatmap cell pops
    # clustering to use for computing cluster medians is specified with k
    # for visualization, we  specify a second layer of cluster annotations with m
    set_dev(plot_path(settings, "heatmap_2layers"))
    p <- plotClusterHeatmap(sce, k = clusters, m = "merge1")
    print(p)
    unset_dev() 
    
    return(p)
}

plot_heatmap_cellpops_2layers_merged <- function(sce, clusters, 
                                                 features = NULL, settings,
                                                 name = "heatmap_2layers_merged") {
    # with merged clusters and freqs
    set_dev(plot_path(settings, name), height = 5)
    p <- plotExprHeatmap(sce, features = features, 
                    by = "cluster_id", k = "merge1", 
                    bars = TRUE, perc = TRUE)
    print(p)
    unset_dev()
    
    return(p)
}

plot_umap_manual_vs_auto <- function(sce, settings, meta = "meta20") {
    p1 <- plotDR(sce, "UMAP", color_by = "merge1")
    p2 <- plotDR(sce, "UMAP", color_by = meta)
    
    set_dev(plot_path(settings, "manual_vs_auto"))
    p <- plot_grid(p1, p2, ncol = 1, align = "v", labels = c("A", "B"))
    print(p)
    unset_dev()
    
    return(p)
}

