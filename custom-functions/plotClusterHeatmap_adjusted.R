library(RColorBrewer)
plotClusterHeatmap_adjusted <- function (x, hm2 = NULL, k = "meta20", m = NULL, fun = c("median", 
                                                         "mean"), cluster_anno = TRUE, split_by = NULL, scale = TRUE, 
          draw_dend = TRUE, draw_freqs = FALSE, palette = rev(brewer.pal(11, 
                                                                         "RdYlBu"))) 
{
    CATALYST:::.check_sce(x)
    fun <- match.arg(fun)
    k <- CATALYST:::.check_validity_of_k(x, k)
    CATALYST:::.check_cd_factor(x, split_by)
    u <- c("abundances", "state_markers", rownames(x))
    if (!is.null(hm2)) 
        stopifnot(hm2 %in% u)
    nk <- nlevels(x$cluster_id <- cluster_ids(x, k))
    ms_by_k <- t(CATALYST:::.agg(x, "cluster_id", fun))
    d <- dist(ms_by_k[, type_markers(x)])
    row_clustering <- hclust(d, method = "average")
    if (cluster_anno) {
        anno <- levels(x$cluster_id)
        if (nk > 30) {
            cols <- colorRampPalette(CATALYST:::.cluster_cols)(nk)
        }
        else {
            cols <- CATALYST:::.cluster_cols[seq_len(nk)]
        }
        cols <- setNames(cols, anno)
        cluster_anno <- CATALYST:::.row_anno(anno, cols, "cluster_id", 
                                  row_clustering, draw_dend)
    }
    if (!is.null(m)) {
        CATALYST:::.check_validity_of_k(x, m)
        idx <- match(seq_len(nk), cluster_codes(x)[, k])
        anno <- factor(cluster_codes(x)[, m][idx])
        if (nlevels(anno) > 30) {
            cols <- colorRampPalette(CATALYST:::.cluster_cols)(nlevels(anno))
        }
        else {
            cols <- CATALYST:::.cluster_cols[seq_len(nlevels(anno))]
        }
        cols <- setNames(cols, levels(anno))
        merging_anno <- CATALYST:::.row_anno(anno, cols, "merging_id", 
                                  row_clustering, draw_dend)
    }
    many <- !is.null(split_by)
    cs <- seq_len(ncol(x))
    if (many) 
        groups <- split(cs, x[[split_by]])
    else groups <- list(cs)
    if (scale) 
        assay(x, "exprs") <- CATALYST:::.scale_exprs(assay(x, "exprs"), 
                                          1)
    hm_cols <- colorRampPalette(palette)(100)
    hms <- lapply(seq_along(groups), function(i) {
        idx <- groups[[i]]
        cs_by_k <- split(idx, x$cluster_id[idx])
        if (!many) {
            if (scale) {
                hm1_es <- t(CATALYST:::.agg(x, "cluster_id", fun))
            }
            else {
                hm1_es <- ms_by_k
            }
        }
        else {
            hm1_es <- t(CATALYST:::.agg(x[, idx], "cluster_id", fun))
        }
        hm1 <- Heatmap(matrix = hm1_es[, pop_markers], col = hm_cols, 
                       name = "expression", column_names_gp = gpar(fontsize = 8), 
                       rect_gp = gpar(col = "white"), na_col = "lightgrey", 
                       cluster_rows = row_clustering, cluster_columns = FALSE, 
                       show_row_dend = draw_dend, column_title = names(groups)[i][many])
        freq_bars <- freq_anno <- NULL
        if (draw_freqs) {
            fq <- round(tabulate(x$cluster_id[idx])/length(idx) * 
                            100, 2)
            freq_bars <- rowAnnotation(`Frequency [%]` = row_anno_barplot(x = fq, 
                                                                          axis = TRUE, border = FALSE, bar_with = 0.8, 
                                                                          gp = gpar(fill = "grey50", col = "white")), 
                                       width = unit(2, "cm"))
            labs <- paste0(levels(x$cluster_id), " (", fq, "%)")
            freq_anno <- rowAnnotation(text = row_anno_text(labs), 
                                       width = max_text_width(labs))
        }
        p <- hm1 + freq_bars + freq_anno
        if (is(cluster_anno, "Heatmap")) 
            p <- cluster_anno + p
        if (exists("merging_anno")) 
            p <- merging_anno + p
        if (!is.null(hm2)) {
            if (hm2 == "abundances") {
                cs <- table(x$cluster_id[idx], x$sample_id[idx])
                fq <- as.matrix(unclass(prop.table(cs, 2)))
                fq <- fq[, !is.na(colSums(fq)), drop = FALSE]
                p <- p + Heatmap(matrix = fq, name = "frequency", 
                                 na_col = "lightgrey", rect_gp = gpar(col = "white"), 
                                 show_row_names = FALSE, column_names_gp = gpar(fontsize = 8), 
                                 cluster_rows = row_clustering, cluster_columns = FALSE)
            }
            else if (hm2 == "state_markers") {
                p <- p + Heatmap(col = hm_cols, na_col = "lightgrey", 
                                 matrix = hm1_es[, state_markers(x)], rect_gp = gpar(col = "white"), 
                                 show_heatmap_legend = FALSE, cluster_rows = row_clustering, 
                                 cluster_columns = FALSE, column_names_gp = gpar(fontsize = 8))
            }
            else {
                for (ch in hm2) {
                    ms <- CATALYST:::.agg(x[ch, idx], c("cluster_id", "sample_id"), 
                               fun)
                    ms <- do.call("rbind", ms)
                    rownames(ms) <- levels(x$cluster_id)
                    p <- p + Heatmap(matrix = ms, col = hm_cols, 
                                     na_col = "lightgrey", rect_gp = gpar(col = "white"), 
                                     show_heatmap_legend = FALSE, show_row_names = FALSE, 
                                     cluster_rows = row_clustering, cluster_columns = FALSE, 
                                     column_title = ch, column_names_gp = gpar(fontsize = 8))
                }
            }
        }
        return(p)
    })
    for (i in seq_along(hms)) draw(hms[[i]])
    invisible(hms)
}