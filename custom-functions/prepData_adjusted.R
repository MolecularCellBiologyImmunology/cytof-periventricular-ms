prepData_adjusted <- function (x, panel, md, features = NULL, cofactor = 5, panel_cols = list(channel = "fcs_colname", 
    antigen = "antigen", class = "marker_class"), md_cols = list(file = "file_name", 
    id = "sample_id", factors = c("condition", "patient_id", "age", "gender", "pmd_hours", "batch_id"))) 
{
    if (!is(panel, "data.frame")) 
        panel <- data.frame(panel, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is(md, "data.frame")) 
        md <- data.frame(md, check.names = FALSE, stringsAsFactors = FALSE)
    stopifnot(is.list(panel_cols), is.list(md_cols), c("channel", 
                                                       "antigen") %in% names(panel_cols), c("file", "id", "factors") %in% 
                  names(md_cols))
    if (!is.null(cofactor)) 
        stopifnot(is.numeric(cofactor), length(cofactor) == 
                      1, cofactor > 0)
    if (is(x, "flowSet")) {
        fs <- x
    }
    else if (is.character(x)) {
        stopifnot(dir.exists(x))
        fcs <- list.files(x, ".fcs$", full.names = TRUE, ignore.case = TRUE)
        if (length(fcs) < 2) 
            stop("The specified directory contains", " none or only a single FCS file.")
        stopifnot(all(vapply(fcs, isFCSfile, logical(1))))
        fs <- read.flowSet(fcs, transformation = FALSE, truncate_max_range = FALSE)
    }
    else {
        stop("Invalid argument 'x'; should be either a flowSet", 
             " or a character string specifying the path to", 
             " a directory containing a set of FCS files.")
    }
    stopifnot(panel[[panel_cols$channel]] %in% colnames(fs))
    if (is.null(features)) {
        features <- as.character(panel[[panel_cols$channel]])
    }
    else {
        chs <- colnames(fs)
        check1 <- is.logical(features) && length(features) == 
            length(chs)
        check2 <- is.integer(features) && all(features %in% 
                                                  seq_along(chs))
        check3 <- all(features %in% chs)
        if (!any(check1, check2, check3)) 
            stop("Invalid argument 'features'. Should be either", 
                 " a logial vector,\n  a numeric vector of indices, or", 
                 " a character vector of column names.")
    }
    ids <- c(keyword(fs, "FILENAME"))
    if (is.null(unlist(ids))) 
        ids <- c(fsApply(fs, identifier))
    stopifnot(all(ids %in% md[[md_cols$file]]))
    # idx <- match(ids, md[[md_cols$file]])
    # fs <- fs[idx]
    if (!is.null(cofactor)) 
        fs <- fsApply(fs, function(ff) {
            exprs(ff) <- asinh(exprs(ff)/cofactor)
            return(ff)
        })
    k <- c(md_cols$id, md_cols$factors)
    md <- data.frame(md)[, k] %>% mutate_all(factor) %>% dplyr::rename(sample_id = md_cols$id)
    o <- order(md[[md_cols$factors[1]]])
    md$sample_id <- factor(md$sample_id, levels = md$sample_id[o])
    antigens <- panel[[panel_cols$antigen]]
    antigens <- gsub("-", "_", antigens)
    antigens <- gsub(":", ".", antigens)
    fs <- fs[, features]
    chs0 <- colnames(fs)
    m1 <- match(panel[[panel_cols$channel]], chs0, nomatch = 0)
    m2 <- match(chs0, panel[[panel_cols$channel]], nomatch = 0)
    flowCore::colnames(fs)[m1] <- antigens[m2]
    chs <- colnames(fs)
    es <- matrix(fsApply(fs, exprs), byrow = TRUE, nrow = length(chs), 
                 dimnames = list(chs, NULL))
    md$n_cells <- as.numeric(fsApply(fs, nrow))
    valid_mcs <- c("type", "state", "none")
    if (is.null(panel_cols$class)) {
        mcs <- factor("none", levels = valid_mcs)
    }
    else {
        mcs <- factor(panel[[panel_cols$class]], levels = valid_mcs)
        mcs <- mcs[match(chs0, panel[[panel_cols$channel]])]
        if (any(is.na(mcs))) 
            stop("Invalid marker classes detected.", " Valid classes are 'type', 'state', and 'none'.")
    }
    rd <- DataFrame(row.names = chs, channel_name = chs0, marker_name = chs, 
                    marker_class = mcs)
    k <- setdiff(names(md), "n_cells")
    cd <- DataFrame(lapply(md[k], function(u) {
        v <- as.character(rep(u, md$n_cells))
        factor(v, levels = levels(u))
    }), row.names = NULL)
    SingleCellExperiment(assays = list(exprs = es), rowData = rd, 
                         colData = cd, metadata = list(experiment_info = md, 
                                                       cofactor = cofactor))
}
