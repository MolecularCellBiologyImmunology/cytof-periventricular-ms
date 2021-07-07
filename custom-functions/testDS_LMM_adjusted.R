testDS_LMM_adjusted <- function (d_counts, d_medians, formula, contrast, weights = TRUE, 
          markers_to_test = NULL, min_cells = 100, min_samples = NULL, 
          contrast_extra = F, contrast_extra_con = NULL) 
{
    # if (is.null(min_samples)) {
    #     min_samples <- ncol(d_counts)/3
    # }
    if (!is.null(markers_to_test)) {
        markers_to_test <- markers_to_test
    }
    else {
        markers_to_test <- metadata(d_medians)$id_state_markers
    }
    counts <- assays(d_counts)[["counts"]]
    cluster_id <- rowData(d_counts)$cluster_id
    ## Threshold of minimum number of cells per cluster:
    # tf <- counts >= min_cells
    # ix_keep <- apply(tf, 1, function(r) sum(r) >= min_samples)
    ix_keep <- apply(counts, 1, function(r) sum(r) > min_cells) # threshold per cluster, not per sample
    
    counts <- counts[ix_keep, , drop = FALSE]
    cluster_id <- cluster_id[ix_keep]
    if (is.logical(weights)) {
        if (weights) {
            weights <- colSums(counts)
        }
        else {
            weights <- NULL
        }
    }
    else if (is.numeric(weights)) {
        stopifnot(length(weights) == ncol(d_counts))
    }
    ## REMOVE this because we specify contrast below 
    # if (ncol(contrast) == 1 & nrow(contrast) > 1) { 
    #     contrast <- t(contrast)
    # }
    
    ### In order to see which cell pop and marker we're analysing:
    state_names <- names(assays(d_medians))[markers_to_test]
    state_names <- state_names[state_names != "CD39"]
    foo <- as.list(assays(d_medians)[state_names])
    bar <- foo
    meds <- rbindlist(lapply(names(bar), function(a) {
        aa <- as.data.table(bar[[a]], keep.rownames = T)
        aa[, marker := a]
        setnames(aa, 'rn', 'pops')
        aa
        }))
    meds <- meds[pops %in% cluster_id] # keep only pops with enough cells (>100 per cluster)
    meds[, unique := paste(pops, marker, sep = "_")]
    meds <- meds[unique %in% state_dt$unique] # Keep only selected combination marker-pop
    # Remove combinations marker-pop with zero expression 
    limit <- ncol(meds[, 2:ncol(meds)])*2/3 # set limit at third samples (at max 0 expression in 2/3 of samples)
    meds <- meds[rowSums(meds==0, na.rm=TRUE) < limit, ]
    
    combi <- unique(meds$unique)
    meds <- as.matrix(meds, rownames = meds$unique)
    meds_all <- do.call("rbind", as.list(assays(d_medians)[state_names]))
    
    p_vals <- c()
    for (i in combi) {
        p_vals[i] <- tryCatch({
            y <- meds[i,]
            print(y[length(y)]) # print cluster name and marker
            to_remove <- c(1, length(y), length(y)-1)
            y <- y[-to_remove]
            data_i <- cbind(cbind.data.frame(y, weights), formula$data)
            data_i$y <- as.numeric(data_i$y)
            if ("age" %in% colnames(data_i)) {
                data_i$age <- as.numeric(data_i$age) # if added to the model
            }
            if ("pmd_hours" %in% colnames(data_i)) {
                data_i$pmd_hours <- as.numeric(data_i$pmd_hours) #  if added to the model
            }
            if (formula$random_terms) {
                fit <- lmer(formula$formula, data = data_i, 
                            weights = weights)
            }
            else {
                fit <- lm(formula$formula, data = data_i, weights = weights)
            }
            print(summary(fit)) # print model
            ### Create contrast (from https://publicifsv.sund.ku.dk/~jufo/courses/rm2018/nlmePackage.pdf):
            name.coef <- names(coef(fit))#[["sample_id"]])
            n.coef <- length(name.coef)
            contr <- matrix(0, nrow = 1, ncol = n.coef,
                            dimnames = list("contrast", name.coef))
            contr[1, contrast] <- 1  # choose among: conditionms, conditionad, age, genderm, pmd_hours
            if (contrast_extra) {
                contr[1, contrast_extra_con] <- -1 
            }# C[1, "conditionms"] <- -1 # optional, if you want to compare to other than intercept (control)
            test <- glht(fit, contr)
            print(summary(test))
            print(contrast_extra)
            summary(test)$test$pvalues
        }, error = function(e) NA)
    }
    
    ### Reduce multiple testing adjustment
    #Adjust per pop + Keep only p_vals of state markers of interest per cell pop
    pvals_dt <- as.data.table(p_vals, keep.rownames = T) 
    setnames(pvals_dt, "rn", "unique")
    # Keep only selected markers:
    pvals_dt <- merge(pvals_dt, state_dt)[, c("unique", "p_vals", "pop", "marker")]
    padj_dt <- c()
    for (popu in unique(pvals_dt$pop)) {
        pop_vals <- pvals_dt[pop == popu]
        pop_vals[, p_adj := p.adjust(p_vals, method = "fdr")]
        padj_dt <- rbind(padj_dt, pop_vals)
    }
    setnames(padj_dt, c("pop", "p_vals", "marker"), c("cluster_id", "p_val", "marker_id"))
    # padj_dt <- padj_dt[, c("p_val", "cluster_id", "p_adj", "marker_id")]
    padj_dt <- as.data.frame(padj_dt)
    res <- padj_dt
    
    ### This is not necessary anymore, since we just return the list of pvals and padjs:
    # p_adj <- p.adjust(p_vals[!is.na(p_vals)], method = "fdr")
    # stopifnot(length(p_vals) == length(p_adj))
    # out <- data.frame(p_val = p_vals, p_adj = p_adj, stringsAsFactors = FALSE)
    # row_data <- as.data.frame(matrix(as.numeric(NA), nrow = length(p_adj), ncol = ncol(out)))
    # colnames(row_data) <- colnames(out)
    # cluster_id_nm <- as.numeric(cluster_id)
    # s <- seq(0, nlevels(cluster_id) * (length(state_names) - 
    #                                        1), by = nlevels(cluster_id))
    # r1 <- rep(cluster_id_nm, length(state_names))
    # r2 <- rep(s, each = length(cluster_id_nm))
    # stopifnot(length(s) == length(state_names))
    # stopifnot(length(r1) == length(r2))
    # rows <- r1 + r2
    # row_data[rows, ] <- out
    # clus <- factor(rep(levels(cluster_id), length(state_names)), 
    #                levels = levels(cluster_id))
    # stat <- factor(rep(state_names, each = length(levels(cluster_id))), 
    #                levels = state_names)
    # stopifnot(length(clus) == nrow(row_data), length(stat) == 
    #               nrow(row_data))
    # row_data <- cbind(data.frame(cluster_id = clus, marker_id = stat, 
    #                              stringsAsFactors = FALSE), row_data)
    # col_data <- colData(d_medians)
    # res <- SummarizedExperiment(meds_all, rowData = row_data, 
    #                             colData = col_data)
    return(list(res = res, meds = meds))
}
