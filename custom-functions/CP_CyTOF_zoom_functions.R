#' Contains functions for zooming into different cell subpopulations.

sce_for_zoom <- function(settings){
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
    # sce <- readRDS(paste0(settings$experiment_dir, settings$tissue, "/", subset_folder, normalised_folder, "sce.RDS"))
    
    sce <- readRDS(object_path(settings, "sce_with_clustering.RDS"))
    
    clusters <- "meta20"
    
    sce <- apply_manual_merging(sce, clusters, settings)
    
    return(sce)
}


### SELECT CLUSTERS
select_clusters <- function(settings){
    if (settings$zoom_pop == 'Tcells') {
        sub_panel <- "CP_CyTOF_panel_Tcells_state_PD1type.csv"
        pop_selec <- c("CD4+ T cells", "CD8+ T cells", "gd T cells")
    }
    if (settings$zoom_pop == 'CD4') {
        sub_panel <- "CP_CyTOF_panel_Tcells_state_PD1type.csv"
        pop_selec <- c("CD4+ T cells")
    }
    if (settings$zoom_pop == 'CD8') {
        sub_panel <- "CP_CyTOF_panel_Tcells_state_PD1type.csv"
        pop_selec <- c("CD8+ T cells")
    }
    if (settings$zoom_pop == 'Myeloid') {
        sub_panel <- "CP_CyTOF_panel_myeloid.csv"
        pop_selec <- c("Myeloid", "Monocytes", "Microglia")
    }
    if (settings$zoom_pop == 'BtoPlasma') {
        sub_panel <- "CP_CyTOF_panel_Bplasma.csv"
        pop_selec <- c("B cells", "B to plasma cells", "B cell lineage")
    }
    if (settings$zoom_pop == 'NKB') {
        sub_panel <- "CP_CyTOF_panel_NKB.csv"
        pop_selec <- c("B cell lineage", "B cells", "Plasma cells", "NK cells")
    }
    if (settings$zoom_pop == 'NK') {
        sub_panel <- "CP_CyTOF_panel_NK.csv"
        pop_selec <- c("NK cells")
    }
    return(list(sub_panel = sub_panel, pop_selec = pop_selec))
}


### CREATE SCE_POP ###
create_sce_pop <- function(sce) {
    ### OPTION A: "manually" selecting the cells with data.table
    # # If you want to check how many clusterings the SCE object contains:
    # show_clusters <- colnames(metadata(sce)$cluster_codes)
    # 
    # # Create dt with som and meta equivalences:
    # equi <- metadata(sce)$cluster_codes
    # equi <- as.data.table(equi)
    # someta <- equi[merge1 %in% pop_selec]
    # 
    # # Create list with som100 ids for the selected cell pops:
    # pop_som <- someta[[1]]
    # 
    # # Subset sce based on som100 clusters:
    # sce_pop <- sce[, sce$cluster_id %in% pop_som]
    
    ### OPTION B: use elegant filterSCE()
    sce_pop <- filterSCE(sce, k = "merge1", cluster_id %in% pop_selec)
    
    return(sce_pop = sce_pop)
}



