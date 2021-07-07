### PLOTS ###
stack_and_boxplots <- function(sce, settings) {
    ### Stacked plot cell type composition per sample
    set_dev(plot_path(settings, "de_abundance_stack"))
    p1 <- plotAbundances(sce, k = "merge1", by = "sample_id")
    print(p1)
    unset_dev()
    
    ### Boxplots
    set_dev(plot_path(settings, "de_abundance_boxplot"))
    p2 <- plotAbundances(sce, k = "merge1", by = "cluster_id", shape = "patient_id")
    print(p2)
    unset_dev()
    
    return(p1)
    return(p2)
}

### PROP TABLE ###
export_prop_table <- function(sce, merge = "merge1", name = "prop_table", settings) {
    ### Export proportion data (out of CD45+ immune cells) ###
    # Exclude CD45- cells from sce:
    equi <- as.data.table(metadata(sce)$cluster_codes)
    cd45_som <- equi[merge1 == "CD45-", som100]
    if (length(cd45_som) != 0) {
        sce_immune <- sce[, sce$cluster_id != cd45_som]
    } else {
        sce_immune <- sce
    }
    
    # Exclude Undefined to calculate proportions:
    exclusion <- cluster_ids(sce_immune, merge) != "Undefined"
    final_pops <- cluster_ids(sce_immune, merge)[exclusion] 
    final_sample_ids <- sample_ids(sce_immune)[exclusion]
    # Create proportion table from sce_immune:
    prop_table <- prop.table(table(final_pops, final_sample_ids), 2)
    prop_table <- as.data.table(prop_table)
    setnames(prop_table, c("final_pops", "final_sample_ids", "N"), c("pop", "sample_id", "prop"))
    # Clarify that Undefined are not counted for calculating proportion by changing values to NA:
    if ("Undefined" %in% prop_table$pop) {
        prop_table[pop == "Undefined", prop := NA]    
    }
    print(table_path(settings, name = name))
    exp_tab(table = prop_table, tab_path = table_path(settings, name = name))
}


### ABSOLUTE TABLE ###
### Create absolute tables and export them:
export_absolute_table <- function(sce, merge = "merge1", name = "absolute_table", settings) {
    absolute_table <- as.data.table(table(cluster_ids(sce, merge), sample_ids(sce)))
    setnames(absolute_table, c("V1", "V2", "N"), c("pop", "sample_id", "n"))
    exp_tab(table = absolute_table, tab_path = table_path(settings, name = name))
}

export_absolute_props <- function(sce, merge = "merge1", name = "abs_props", settings,
                                  import_tab = "absolute_table") {
    abs_zoomed <- fread(paste0("./", settings$tissue,
                               "/subset/original/", settings$zoom_pop, "/", settings$tissue,
                               "_", import_tab, ".csv" )) 
    abs_immune <- fread(paste0("./", settings$tissue,
                              "/subset/original/", settings$tissue, "_", import_tab, ".csv" ))
    abs_immune <- abs_immune[, sum(n), by = "sample_id"]
    setnames(abs_immune, "V1", "n_immune")
    abs_props <- merge(abs_zoomed, abs_immune, by = "sample_id")
    abs_props[, prop := n/n_immune]
    abs_props <- abs_props[, c("sample_id", "pop", "prop")]
    exp_tab(table = abs_props, tab_path = table_path(settings, name = name))
}


### DIFFERENTIAL ABUNDANCE ###
### CREATE FORMULAS
create_da_formulas <- function(sce, settings) {
    ### Compare proportions with GLMM
    ei <- metadata(sce)$experiment_info
    
    ### Create formulas (two depending if we include PMD or not)
    # General
    (da_formula1 <- createFormula(ei,
                                  cols_fixed = c("condition"
                                                 # , "age"
                                                 # , "gender"
                                                 # , "pmd_hours"
                                                 # , "batch_id"
                                  ),
                                  cols_random = "sample_id")) 
    
    # For septum, PMD is not different in MS vs con and AD vs con
    if (settings$tissue == "septum") {
        (da_formula2 <- createFormula(ei,
                                      cols_fixed = c("condition"
                                                     # ,"age"
                                                     # , "gender"
                                                     # # , "pmd_hours" # not significant vs control
                                                     # , "batch_id"
                                      ),
                                      cols_random = "sample_id")) 
    } else {
        (da_formula2 <- createFormula(ei,
                                      cols_fixed = c("condition"
                                                     # ,"age"
                                                     # , "gender"
                                                     # , "pmd_hours"
                                                     # , "batch_id"
                                      ),
                                      cols_random = "sample_id")) 
    }
    return(list(da_formula1 = da_formula1, da_formula2 = da_formula2))
}


### DIFFERENTIAL ABUNDANCE
differential_abundance <- function(sce, settings, contrast = "conditionms", formula = "da_formula1", 
                      contrast_extra = F, contrast_extra_con = NULL, clustering_to_use = "merge1",
                      group1 = "control", group2 = "ms", different_prop = F,
                      name = "da_res_ms") {
    # Start writing to an output file
    sink(table_path(settings, name = name, format = ".txt"))
    da_res <- diffcyt_adjusted(sce,
                               formula = formula,
                               contrast = contrast, # choose among: conditionms, conditionad, age, genderm, pmd_hours
                               contrast_extra = contrast_extra, # optional, if you want to compare to other than intercept (control)
                               contrast_extra_con = contrast_extra_con, # optional, if you want to compare to other than intercept (control)
                               analysis_type = "DA", method_DA = "diffcyt-DA-GLMM", settings = settings, 
                               different_prop = different_prop, 
                               clustering_to_use = clustering_to_use, verbose = TRUE)
    sink()
    
    rowData(da_res$res)
    table(rowData(da_res$res)$p_adj < FDR_cutoff)
    topTable <- as.data.frame(topTable(da_res,show_props = TRUE,
                                       format_vals = TRUE, digits = 4))
    exp_tab(table = topTable, tab_path = table_path(settings, name = paste0("topTable_", name)))

    # Plot heatmap
    sce_comp <- sce[,sce$condition %in% c(group1, group2)]
    set_dev(plot_path(settings, paste0("diff_heatmap_", name)))
    hm <- plotDiffHeatmap(sce_comp, rowData(da_res$res), all = TRUE, fdr = FDR_cutoff)
    print(hm)
    unset_dev()
    
}

### NON-PARAMETRIC TESTS ###
### GET FREQUENCIES
get_freqs <- function(sce, merge = "merge1") {
    ## Create table with cell frequencies of each cluster:
    fq <- prop.table(table(cluster_ids(sce, "merge1"), sample_ids(sce)), 
                     2) * 100
    df <- melt(fq, value.name = "freq", varnames = c("cluster_id", 
                                                     "sample_id"))
    m <- match(df$sample_id, ei(sce)$sample_id)
    cols <- setdiff(names(ei(sce)), names(df))
    df <- data.frame(df, ei(sce)[m, cols])
    df <- as.data.table(df)
    return(df)
}

### KRUSKALL-WALLIS
## Perform Kruskall-Wallis for each cluster:
kruskal_per_cluster <- function(df, settings, name = "kruskall") {
    cluster_id <- unique(df$cluster_id)
    kruskall <- as.data.table(cluster_id)
    kruskall[, pvalue := -1]
    p <- c()
    for (cluster in cluster_id) {
        res.kruskal <- kruskal.test(freq ~ condition, data = df[cluster_id==cluster])
        print(cluster)
        print(res.kruskal)
        # p <- append(p, res.kruskal$p.value)
        p <- res.kruskal$p.value
        kruskall[cluster_id == cluster, pvalue := p]
    }
    
    kruskall[, padjust := p.adjust(kruskall$pvalue, method = "fdr")]
    exp_tab(table = kruskall, tab_path = table_path(settings, name = name))
}

### POSTHOC
posthoc_pairwise <- function(df) {
    ## Post hoc pairwise control, MS, AD
    cluster_id <- unique(df$cluster_id)
    p.pairwise <- c()
    for (cluster in cluster_id) {
        res.pairwise <- pairwise.wilcox.test(df[cluster_id==cluster]$freq,
                                             df[cluster_id==cluster]$condition,
                                             p.adjust.method = "BH")
        print(cluster)
        print(res.pairwise)
        p.pairwise <- append(p.pairwise, res.pairwise$p.value)
    }
}


### DIFFERENTIAL EXPRESSION ###
### Differential analysis of marker expression stratified by cell population ###
plot_median_expression <- function(sce, settings) {
    # Plot median expression of all markers in each cluster for each sample, by condition
    set_dev(plot_path(settings, "marker_expression"))
    p <- plotPbExprs(sce, k = "merge1",
                     facet = "cluster_id", shape_by = "patient_id")
    p$facet$params$ncol <- 2
    print(p)
    unset_dev()
}

### CREATE FORMULAS
create_ds_formulas <- function(sce) {
    ei <- metadata(sce)$experiment_info
    ### Compare expression with LMM linear model
    # Formula
    ds_formula1 <- createFormula(ei, cols_fixed = c("condition"
                                                    # ,"age"
                                                    # , "gender"
                                                    # , "pmd_hours"
                                                    # , "batch_id"
    )
    # ,cols_random = "" # no random effect
    )
    return(list(ds_formula1 = ds_formula1))
}

### Select markers of interest
select_markers <- function(markers_of_interest = "CP_CyTOF_state_markers_per_pop.csv") {
    ### To only study DE in markers of interest for each cell pop:
    # Import table with state markers of interest per cell pop
    markers_test <- fread(paste0("./", markers_of_interest),
                          na.strings = c("", NULL, "NA"))
    state <- colnames(markers_test)[2:17] # Create list with all state marker names
    state_dt <- melt(markers_test, id.vars = c("pop"), measure.vars = state, 
                     variable.name = "marker", value.name = "include")
    state_dt <- state_dt[include == "x"] # Keep only those of interest
    state_dt[, unique := paste(pop, marker, sep = "_")] # Create single unique identifier per combi
    return(state_dt)
}

### DIFFERENTIAL EXPRESSION
differential_expression <- function(sce, settings, contrast = "conditionms", formula = "ds_formula1", 
                                    contrast_extra = F, contrast_extra_con = NULL, 
                                    clustering_to_use = "merge1",
                                    group1 = "control", group2 = "ms",
                                    name = "ds_res_ms") {
    sink(table_path(settings, name = name, format = ".txt"))
    ds_res <- diffcyt_adjusted(sce,
                                  formula = formula, 
                                  contrast = contrast,
                                  contrast_extra = contrast_extra,
                                  contrast_extra_con = contrast_extra_con, 
                                  analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                                  clustering_to_use = clustering_to_use, verbose = FALSE)
    sink()
    
    ds_pvals <- as.data.table(ds_res$res)
    ds_medians <- ds_res$meds
    ds_results <- merge(ds_pvals, ds_medians, by = "unique")
    ds_results <- ds_results[, !c("pops", "marker")]
    table(ds_results$p_adj < FDR_cutoff)
    exp_tab(table = ds_results, tab_path = table_path(settings, name = paste0(name, "_results")))
    
    ###
    # Heatmap not working because I removed 
    # metadata(res) <- as.list(c(metadata(res), clustering_name = clustering_name))
    # from diffcyt_adjusted()
    # Heatmap to report the differential signals
    # plotDiffHeatmap(sce, ds_res_ms, top_n = 50, order = TRUE,
    #                 th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)
    # set_dev(plot_path(SETTINGS, "ds_heatmap_ms"))
    # plotDiffHeatmap(sce_ms, ds_res_ms, top_n = 50, order = TRUE,
    #                 th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)
    # unset_dev()
}


plot_heatmap_expresssion <- function(sce, fdr = 0.1, settings,
                                     name1 = "ds_res_ms_results.csv", name2 = "ds_res_advsms_results.csv") {
    # Load latest version of heatmap.3 function:
    source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
    
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
    } else {
        zoom_folder <- NULL
    }
    tissue <- settings$tissue
    
    populations <- as.list(unique(metadata(sce)$cluster_codes["merge1"]))$merge1
    populations <- populations[!populations %in% c("CD45-", "Undefined")]
    
    for (population in populations) {
        ### Import table with median expressions
        # MS
        dt_ms <- as.data.table(
            fread(paste0(settings$experiment_dir, "/", tissue, "/", 
                         subset_folder, normalised_folder, zoom_folder, tissue, "_", name1)))
        dt2_ms <- dt_ms[cluster_id == population]
        dt2_ms <-  unique(dt2_ms, by = "marker_id") # remove duplicated rows
        dt3 <- dt2_ms[, !(c("p_val", "p_adj"))] # remove the p-value for heatmap
        
        ### Define order of heatmap:
        sample_order <- c("control01", "control02", "control03", "control04", "control05", 
                          "control06", "control07", "control08", "control09", "control10",  
                          "control11", "control12",    
                          "ad01", "ad02", "ad03", "ad04", "ad05", "ad06","ad07", "ad08",
                          "ms01","ms02", "ms03", "ms04", "ms05", "ms06", 
                          "ms07", "ms08", "ms09", "ms10", "ms11", "ms12", "ms13" )
        sample_order <- sample_order[sample_order %in% colnames(dt3)] # keep only existing samples 
        setcolorder(dt3, sample_order) # reorder columns for heatmap
        
        # AD (load for p-values only)
        dt_ad <- as.data.table(
            fread(paste0(settings$experiment_dir, "/", tissue, "/", 
                         subset_folder, normalised_folder, zoom_folder, tissue, "_", name2)))
        dt2_ad <- dt_ad[cluster_id == population]
        dt2_ad <-  unique(dt2_ad, by = "marker_id") # remove duplicated rows
        
        ### Get adjusted p-values
        padj_ms <- dt2_ms[, p_adj] # MS adjusted p-value for annotations
        padj_ad <- dt2_ad[, p_adj] # AD adjusted p-value for annotations
        padj <- cbind(padj_ad, padj_ms)
        padj_colors <- ifelse(padj > fdr, "gray", 
                              ifelse(padj < 0.05 & padj > 0.01, "olivedrab3", "olivedrab2")) # assign colors
        padj_colors <- t(padj_colors) # transpose for the heatmap
        
        ### Create condition groups for annotating the heatmap
        condition_colors <- ifelse(grepl("control", sample_order), "#043741",
                                   ifelse(grepl("ms", sample_order), "#e79d24", 
                                          ifelse(grepl("ad", sample_order), "#189cb3", NA)))
        condition_colors <- as.matrix(condition_colors)
        # 
        # ### Create marker groups
        # marker_order <- c("CD49d", "CD54", "CD31/PECAM_1", 
        #                   "Tbet", "CD45RA", "CD45RO",
        #                   "PDL1", 
        #                   "CD25", "CD69")
        # foo <- data.table(marker_order)
        # foo[, marker_type := ifelse(marker_order %in% c("CD49d", "CD54", "CD31/PECAM_1"), "adhesion", 
        #                             ifelse(marker_order %in% c("Tbet", "CD45RA", "CD45RO"), "memory/maturation", 
        #                                    ifelse(marker_order %in% c("PDL1"), "inhibition", "residency/activation")))]
        
        ### Create matrix
        mat <- as.matrix(dt3[, !c("unique", "cluster_id", "marker_id")])
        rownames(mat) <- dt3$marker_id
        
        palette = rev(brewer.pal(11, "RdYlBu"))
        # palette = colorRampPalette(c("deepskyblue4","white","darkorange2"))(200)
        
        ### Create heatmap
        
        main_title <- paste(population, tissue)
        set_dev(plot_path(settings, paste0("expr_heatmap_", population)))
        heatmap.3(mat, #hclustfun=myclust, distfun=mydist, na.rm = TRUE,
                  scale="row",
                  dendrogram="none", 
                  # margins=c(9,9),
                  Rowv=TRUE, Colv=F, #ColSideColors=clab,
                  RowSideColors=padj_colors, RowSideColorsSize=2,
                  ColSideColors=condition_colors, ColSideColorsSize=1,
                  symbreaks=T,
                  key=TRUE, keysize = 0.75, symkey=T,
                  # density.info="histogram",
                  main=main_title, 
                  cexRow=2, cexCol=2,
                  col=palette
                  )
        legend("topright",legend=c("padj > 0.1","0.05 < padj < 0.1", "padj < 0.05"),
               fill=c("gray","olivedrab3", "olivedrab2"), border=FALSE, bty="n", 
               y.intersp = 1, cex=1.25)
        unset_dev()
        
              }
}



