#' Modified version of diffcyt that calls adjusted versions of the 
#' diffcyt-DA-GLMM and diffcyt-DS-LMM methods

source("./Custom_functions/testDA_GLMM_adjusted.R") 
source("./Custom_functions/testDA_GLMM_adjusted_outofCD45+.R") # For proportions out of CD45+ cells
source("./Custom_functions/testDS_LMM_adjusted.R") 

diffcyt_adjusted <- function (d_input, experiment_info = NULL, marker_info = NULL, 
          design = NULL, formula = NULL, contrast, 
          contrast_extra, contrast_extra_con,
          analysis_type = c("DA", "DS"), different_prop = FALSE, 
          method_DA = c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM"), settings, 
          method_DS = c("diffcyt-DS-limma", "diffcyt-DS-LMM"), 
          markers_to_test = NULL, clustering_to_use = NULL, 
          cols_to_include = NULL, subsampling = FALSE, n_sub = NULL, 
          seed_sub = NULL, transform = TRUE, cofactor = 5, cols_clustering = NULL, 
          xdim = 10, ydim = 10, meta_clustering = FALSE, meta_k = 40, 
          seed_clustering = NULL, min_cells = 100, min_samples = NULL, # adjusted min_cells to total 100
          normalize = FALSE, norm_factors = "TMM", trend_method = "none", 
          block_id = NULL, trend = TRUE, weights = TRUE, plot = TRUE, 
          path = ".", verbose = TRUE) 
{
    analysis_type <- match.arg(analysis_type)
    method_DA <- match.arg(method_DA)
    method_DS <- match.arg(method_DS)
    if (!is(d_input, "SingleCellExperiment")) {
        if (is.null(experiment_info) | is.null(marker_info)) {
            stop("'experiment_info' and 'marker_info' must be provided (unless using a SingleCellExperiment ", 
                 "object from CATALYST as input)")
        }
        if (verbose) 
            message("preparing data...")
        d_se <- prepareData(d_input, experiment_info, marker_info, 
                            cols_to_include, subsampling, n_sub, seed_sub)
        if (transform) {
            if (verbose) 
                message("transforming data...")
            d_se <- transformData(d_se, cofactor)
        }
        if (verbose) 
            message("generating clusters...")
        d_se <- generateClusters(d_se, cols_clustering, xdim, 
                                 ydim, meta_clustering, meta_k, seed_clustering)
    }
    else if (is(d_input, "SingleCellExperiment")) {
        if (verbose) 
            message("using SingleCellExperiment object from CATALYST as input")
        if (is.null(clustering_to_use)) {
            stopifnot("cluster_id" %in% colnames(colData(d_input)))
            if (verbose) 
                message("using cluster IDs stored in column named 'cluster_id' in 'colData' of ", 
                        "SingleCellExperiment object from CATALYST")
            clustering_name <- colnames(metadata(d_input)$cluster_codes)[1]
        }
        else if (!is.null(clustering_to_use)) {
            stopifnot(as.character(clustering_to_use) %in% colnames(metadata(d_input)$cluster_codes))
            stopifnot("cluster_id" %in% colnames(colData(d_input)))
            if (verbose) 
                message("using cluster IDs from clustering stored in column '", 
                        clustering_to_use, "' of 'cluster_codes' data frame in 'metadata' of SingleCellExperiment object from CATALYST")
            code_id <- colData(d_input)$cluster_id
            cluster_id <- metadata(d_input)$cluster_codes[, 
                                                          clustering_to_use][code_id]
            stopifnot(length(cluster_id) == nrow(colData(d_input)), 
                      length(code_id) == nrow(colData(d_input)))
            colData(d_input)$code_id <- code_id
            colData(d_input)$cluster_id <- cluster_id
            clustering_name <- clustering_to_use
        }
        stopifnot("sample_id" %in% colnames(colData(d_input)))
        stopifnot("experiment_info" %in% names(metadata(d_input)))
        stopifnot("cluster_id" %in% colnames(colData(d_input)))
        stopifnot("cluster_codes" %in% names(metadata(d_input)))
        cs_by_s <- split(seq_len(ncol(d_input)), colData(d_input)$sample_id)
        cs <- unlist(cs_by_s[as.character(metadata(d_input)$experiment_info$sample_id)])
        es <- t(assays(d_input)[["exprs"]])[cs, , drop = FALSE]
        d_se <- SummarizedExperiment(assays = list(exprs = es), 
                                     rowData = colData(d_input)[cs, ], colData = rowData(d_input), 
                                     metadata = metadata(d_input))
    }
    if (verbose) 
        message("calculating features...")
    d_counts <- calcCounts(d_se)
    d_medians <- calcMedians(d_se)
    d_medians_by_cluster_marker <- calcMediansByClusterMarker(d_se)
    d_medians_by_sample_marker <- calcMediansBySampleMarker(d_se)
    if (analysis_type == "DA" && method_DA == "diffcyt-DA-edgeR") {
        if (verbose) 
            message("calculating DA tests using method 'diffcyt-DA-edgeR'...")
        res <- testDA_edgeR(d_counts, design, contrast, trend_method, 
                            min_cells, min_samples, normalize, norm_factors)
    }
    if (analysis_type == "DA" && method_DA == "diffcyt-DA-voom") {
        if (verbose) 
            message("calculating DA tests using method 'diffcyt-DA-voom'...")
        res <- testDA_voom(d_counts, design, contrast, block_id, 
                           min_cells, min_samples, normalize, norm_factors, 
                           plot, path)
    }
    if (analysis_type == "DA" && method_DA == "diffcyt-DA-GLMM") {
        if (verbose) 
            message("calculating DA tests using method 'diffcyt-DA-GLMM'...")
        if (different_prop == T) {
            res <- testDA_GLMM_adjusted_outofCD45(d_counts, formula, contrast, min_cells, 
                                                  min_samples, normalize, norm_factors,
                                                  contrast_extra, contrast_extra_con, settings)
        } else {
            res <- testDA_GLMM_adjusted(d_counts, formula, contrast, min_cells, 
                                        min_samples, normalize, norm_factors,
                                        contrast_extra, contrast_extra_con)
        }
        
        meds <- NULL
    }
    if (analysis_type == "DS" && method_DS == "diffcyt-DS-limma") {
        if (verbose) 
            message("calculating DS tests using method 'diffcyt-DS-limma'...")
        res <- testDS_limma(d_counts, d_medians, design, contrast, 
                            block_id, trend, weights, markers_to_test, min_cells, 
                            min_samples, plot, path)
    }
    if (analysis_type == "DS" && method_DS == "diffcyt-DS-LMM") {
        if (verbose) 
            message("calculating DS tests using method 'diffcyt-DS-LMM'...")
        res_and_meds <- testDS_LMM_adjusted(d_counts, d_medians, formula, contrast, 
                          weights, markers_to_test, min_cells, min_samples, 
                          contrast_extra, contrast_extra_con)
        res <- res_and_meds$res
        meds <- res_and_meds$meds
    }
    if (!is(d_input, "SingleCellExperiment")) {
        return(list(res = res, d_se = d_se, d_counts = d_counts, 
                    d_medians = d_medians, d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
                    d_medians_by_sample_marker = d_medians_by_sample_marker))
    }
    else if (is(d_input, "SingleCellExperiment")) {
        if (analysis_type == "DA") {
            metadata(res) <- as.list(c(metadata(res), clustering_name = clustering_name))
        }
        return(list(res = res, d_counts = d_counts, d_medians = d_medians, 
                    meds = meds, 
                    d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
                    d_medians_by_sample_marker = d_medians_by_sample_marker))
    }
}
