plotDiffHeatmap_adjusted <- function (x, y, top_n = 20, all = FALSE, order = TRUE, th = 0.1, 
                                      type = "DA", hm1 = TRUE, normalize = TRUE, 
                                      row_anno = TRUE, col_anno = TRUE) 
{
    CATALYST:::.check_sce(x)
    es <- assay(x, "exprs")
    stopifnot(is.numeric(top_n), length(top_n) == 1, is.logical(order), 
              length(order) == 1, is.numeric(th), length(th) == 1, 
              is.logical(hm1), length(hm1) == 1, is.logical(normalize), 
              length(normalize) == 1, is.logical(row_anno), length(row_anno) == 
                  1)
    stopifnot(!is.null(k <- metadata(y$res)$clustering_name))
    k <- CATALYST:::.check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
    factors <- dplyr::select(as.data.frame(colData(x)), -c("sample_id", 
                                                    "cluster_id"))
    y <- rowData(y$res)
    if (order) 
        y <- y[order(y$p_adj), , drop = FALSE]
    if (all | top_n > nrow(y)) 
        top_n <- nrow(y)
    top <- as.data.frame(y[seq_len(top_n), ])
    top <- mutate_if(top, is.factor, as.character)
    if (hm1) {
        ms_by_k <- t(CATALYST:::.agg(x[type_markers(x)], "cluster_id"))[top$cluster_id, 
                                                             ]
        qs <- quantile(ms_by_k, probs = c(0.01, 0.5, 0.99), 
                       na.rm = TRUE)
        hm_cols <- circlize::colorRamp2(qs, c("royalblue3", "white", "tomato2"))
        hm1 <- CATALYST:::.diff_hm(ms_by_k, hm_cols, "expression", cluster_rows = !order, 
                        xlab = "type_markers", row_title = "cluster_id"[!is.null(hm1)], 
                        row_names_side = "left")
    }
    else {
        hm1 <- NULL
    }
    if (col_anno) {
        # HERE YOU DEFINE THE ANNOTATIONS
        m <- match(levels(x$sample_id), x$sample_id)
        # m <- !is.na(m)
        df <- data.frame(factors[m, ], row.names = NULL)
        df <- df[, c("condition", "patient_id", "age",
                     "gender", "pmd_hours", "batch_id")] # subset df to get rid of age, gender and pmd
        col_anno <- CATALYST:::.anno_factors(df, "column")
    }
    else {
        col_anno <- NULL
    }
    if (type == "DA") {
        cnts <- table(x$cluster_id, x$sample_id)
        frqs <- prop.table(cnts, 2)
        frqs <- frqs[top$cluster_id, ]
        frqs <- as.matrix(unclass(frqs))
        # Remove unwanted samples:
        # ms <- c("ms01", "ms02", "ms03", "ms04", "ms05", "ms06",
        #         "ms07", "ms08", "ms09", "ms10", "ms11")
        # ad <- c("ad01", "ad02", "ad03", "ad04", "ad05", "ad06", "ad07")
        # frqs <- as.data.table(frqs)
        # frqs <- frqs[, colnames(frqs)[!(colnames(frqs) %in% ms)], with = F]
        # frqs <- as.matrix(frqs)
        
        if (normalize) {
            frqs <- CATALYST:::.z_normalize(asin(sqrt(frqs)))
            at <- seq(-2.5, 2.5, 0.5)
            labels <- at
            labels[-seq(2, length(at), 2)] <- ""
        }
        else {
            min <- floor(min(frqs)/0.1) * 0.1
            max <- ceiling(max(frqs)/0.1) * 0.1
            at <- seq(min, max, 0.1)
            labels <- at
        }
        # HERE YOU PLOT THE ANNOTATIONS
        hm2 <- CATALYST:::.diff_hm(matrix = frqs, cluster_rows = !order, 
                        col = c("skyblue", "cornflowerblue", "royalblue", 
                                "black", "orange3", "orange", "gold"), 
                        name = paste0("normalized\n"[normalize], "frequency"), 
                        show_row_names = is.null(hm1), 
                        row_names_side = "left", 
                        heatmap_legend_param = list(at = at, labels = labels), 
                        top_annotation = col_anno, 
                        xlab = "sample_id", row_title = "cluster_id")
    }
    else {
        cs_by_ks <- CATALYST:::.split_cells(x, c("cluster_id", "sample_id"))
        ms_by_ks <- t(mapply(function(k, g) vapply(cs_by_ks[[k]], 
                                                   function(cs) median(es[g, cs, drop = FALSE]), numeric(1)), 
                             k = top$cluster_id, g = top$marker_id))
        if (!is.null(hm1)) {
            rownames(ms_by_ks) <- top$marker_id
        }
        else {
            rownames(ms_by_ks) <- sprintf("%s(%s)", top$marker_id, 
                                          top$cluster_id)
        }
        if (normalize) 
            ms_by_ks <- CATALYST:::.z_normalize(ms_by_ks)
        hm2 <- CATALYST:::.diff_hm(matrix = ms_by_ks, cluster_rows = !order, 
                        name = paste0("normalized\n"[normalize], "expression"), 
                        col = c("skyblue", "cornflowerblue", "royalblue", 
                                "black", "orange3", "orange", "gold"), xlab = "sample_id", 
                        top_annotation = col_anno, row_names_side = c("right", 
                                                                      "left")[as.numeric(is.null(hm1)) + 1])
    }
    if (row_anno) {
        s <- top$p_adj <= th
        s[is.na(s)] <- FALSE
        s <- as.matrix(c("no", "yes")[as.numeric(s) + 1])
        rownames(s) <- format(top$p_adj, digits = 2)
        row_anno <- Heatmap(matrix = s, name = "significant", 
                            col = c(no = "lightgrey", yes = "limegreen"), width = unit(5, 
                                                                                       "mm"), rect_gp = gpar(col = "white"), show_row_names = TRUE, 
                            row_names_side = "right")
    }
    else {
        row_anno <- NULL
    }
    main <- switch(type, DA = "top DA clusters", DS = "top DS cluster-marker combinations")
    suppressWarnings(draw(hm1 + hm2 + row_anno, column_title = main, 
                          auto_adjust = FALSE, column_title_gp = gpar(fontface = "bold", 
                                                                      fontsize = 12)))
}

