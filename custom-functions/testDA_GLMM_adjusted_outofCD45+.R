testDA_GLMM_adjusted_outofCD45 <- function (d_counts, formula, contrast, min_cells = 100, 
                                  min_samples = NULL, normalize = FALSE, norm_factors = "TMM",
                                  contrast_extra = F, contrast_extra_con = NULL, settings) 
{
    # if (is.null(min_samples)) {
    #     min_samples <- ncol(d_counts)/3
    # }
    counts <- assays(d_counts)[["counts"]]
    counts <- counts[!rownames(counts) %in% c("Unknown", "CD45-", "Undefined"), ] # to reduce multiple testing correction
    cluster_id <- rowData(d_counts)$cluster_id
    cluster_id <- cluster_id[!cluster_id %in% c("Unknown", "CD45-", "Undefined")] # to reduce multiple testing correction
    ## Threshold of minimum number of cells per cluster:
    # tf <- counts >= min_cells
    # ix_keep <- apply(tf, 1, function(r) sum(r) >= min_samples)
    ix_keep <- apply(counts, 1, function(r) sum(r) > min_cells) # threshold per cluster, not per sample
    
    counts <- counts[ix_keep, , drop = FALSE]
    cluster_id <- cluster_id[ix_keep]
    if (normalize & norm_factors == "TMM") {
        norm_factors <- calcNormFactors(counts, method = "TMM")
    }
    if (normalize) {
        n_cells_smp <- colSums(counts)/norm_factors
    }
    else { # Import table with absolute counts of CD45+ immune cells:
        absolute_tab <- fread(paste0("~/Documents/Data analysis/CyTOF/CP_CyTOF_data/",
                                        settings$tissue, "/subset/original/", 
                                        settings$tissue, "_absolute_table.csv"))
        absolute_immune <- absolute_tab[pop != "CD45-", sum(n), by = "sample_id"]
        setnames(absolute_immune, "V1", "n_cells_smp")
    }
    # if (ncol(contrast) == 1 & nrow(contrast) > 1) {
    #     contrast <- t(contrast)
    # }
    p_vals <- rep(NA, length(cluster_id))
    for (i in seq_along(cluster_id)) {
        p_vals[i] <- tryCatch({ # Create table with proportions out of total CD45+ cells:
            i_counts <- as.data.table(counts[i, ], keep.rownames = T)
            setnames(i_counts, c("V1", "V2"), c("sample_id", "n"))
            data_i <- merge(i_counts, absolute_immune, by = "sample_id")
            data_i <- merge(data_i, formula$data, by = "sample_id")
            data_i[, y := n/n_cells_smp]
            data_i <- data_i[, c("y", "n_cells_smp", "condition", "sample_id")]
            
            if ("age" %in% colnames(data_i)) {
                data_i$age <- as.numeric(data_i$age) # if added to the model
            }
            if ("pmd_hours" %in% colnames(data_i)) {
                data_i$pmd_hours <- as.numeric(data_i$pmd_hours) #  if added to the model
            }
            fit <- tryCatch(
                glmer(formula$formula, data = data_i, family = "binomial",
                      weights = n_cells_smp)
                , warning = function(w) NULL
            ) 
            if (is.null(fit)) {
                fit <- tryCatch(
                    glmer(formula$formula, 
                          data = data_i,
                          family = "binomial", 
                          weights = n_cells_smp, 
                          control=glmerControl(optimizer="bobyqa", # to reach convergence (https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html; https://biologyforfun.wordpress.com/2018/04/09/help-i-have-convergence-warnings/)
                                               optCtrl=list(maxfun=2e6)))
                    , warning = function(w) NULL
                )
            }
            if (is.null(fit)) {
                fit <- glmer(formula$formula, 
                          data = data_i,
                          family = "binomial", 
                          weights = n_cells_smp)
                warn <- (paste0("*** Warning_", cluster_id[i]))
                print(warn)
            } 
            print(cluster_id[i]) # print cluster name
            print(summary(fit)) # print model
            ### Create contrast (from https://publicifsv.sund.ku.dk/~jufo/courses/rm2018/nlmePackage.pdf):
            name.coef <- names(coef(fit)[["sample_id"]])
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
            cat("\n")
            summary(test)$test$pvalues
        }, error = function(e) NA)
    }
    p_adj <- p.adjust(p_vals, method = "fdr")
    stopifnot(length(p_vals) == length(p_adj))
    out <- data.frame(p_val = p_vals, p_adj = p_adj, stringsAsFactors = FALSE)
    row_data <- as.data.frame(matrix(as.numeric(NA), nrow = nlevels(cluster_id), 
                                     ncol = ncol(out)))
    colnames(row_data) <- colnames(out)
    cluster_id_nm <- as.numeric(cluster_id)
    row_data[cluster_id_nm, ] <- out
    row_data <- cbind(cluster_id = rowData(d_counts)$cluster_id, 
                      row_data)
    res <- d_counts
    rowData(res) <- row_data
    if (normalize) {
        metadata(res)$norm_factors <- norm_factors
    }
    res
}
