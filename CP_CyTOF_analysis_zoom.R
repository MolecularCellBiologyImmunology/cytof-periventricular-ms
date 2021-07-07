### Load packages
library(HDCytoData)
library(readxl)
library(CATALYST)
library(ggplot2)
library(cowplot)
library(diffcyt)
library(SummarizedExperiment)
library(data.table)
library(stringr)
library(tidyverse)
library(SingleCellExperiment)
library(lme4)
library(multcomp)
library(ComplexHeatmap) # for DA heatmaps
library(dplyr)
library(RColorBrewer)
library(devtools)

setwd("~/Documents/Data analysis/CyTOF/CP_CyTOF_data/")

############################
### CHOOSE YOUR ANALYSIS ###
# Settings for parent clustering (choose tissue, no zoom)
PRESETTINGS <- list(
    experiment_dir = normalizePath("."),
    SAVE_PLOTS = T, # Export plots
    SAVE_TABLES = T, # Export tables
    norm = F, # Normalisation
    subsetting = T, # Subsetting samples by cell number
    zoom = F, # Zoom into cell pops?
    zoom_pop = NULL, # Zoom into cell pops?
    tissue = 'cp' # Choose tissue
)

# Settings for zoom clustering
SETTINGS <- list(
    experiment_dir = normalizePath("."),
    SAVE_PLOTS = T, # Export plots
    SAVE_TABLES = T, # Export tables
    norm = F, # Normalisation
    subsetting = T, # Subsetting samples by cell number
    zoom = T, # Zoom into cell pops? 
    zoom_pop = "NKB", # Zoom into cell pops? (Tcells/CD4/CD8/Myeloid/BtoPlasma/NKB/NK)
    tissue = 'cp' # Choose tissue
)

source("./Custom_functions/CP_CyTOF_functions.R") # Load created functions
source("./Custom_functions/CP_CyTOF_zoom_functions.R")
source("./Custom_functions/clustering.R")

### LOAD SCE ###
sce <- sce_for_zoom(PRESETTINGS)

### SELECT CLUSTERS
selec_clus <- select_clusters(SETTINGS)
sub_panel <- selec_clus$sub_panel
pop_selec <- selec_clus$pop_selec

### CREATE SCE_POP ###
sce_pop <- create_sce_pop(sce)

####################
###  CLUSTERING  ###
####################
source("./Custom_functions/clustering.R")

set.seed(1234) # set random seed for reproducibility

### Adjust type markers for the specific subpopulation
pop_markers <- fread(paste0("./", 
                            sub_panel))
pop_markers <- pop_markers[marker_class == "type", antigen]

sce_pop <- run_clustering(sce = sce_pop, settings = SETTINGS, features = pop_markers, skip_checkpoint = T)

plot_elbow_plot(sce_pop, SETTINGS)

# ADJUST k:
clusters <- "meta8"

plot_heatmap_cellpops(sce_pop, SETTINGS, clusters, features = pop_markers, name = "heatmap_cellpops")

plot_heatmap_cellpops(sce_pop, SETTINGS, clusters, features = pop_markers, bin_anno = T, 
                      name = "heatmap_cellpops_anno")

plot_distribution_markers(sce_pop, clusters, SETTINGS)

plot_specific_marker_heatmap(sce_pop, "CD45RA", clusters, SETTINGS)

sce_pop <- run_tsne_and_umap(sce_pop, SETTINGS, skip_checkpoint = T)

plot_tsne_and_umap(sce_pop, clusters, SETTINGS)

plot_umap_per_sample(sce_pop, clusters, SETTINGS)

plot_umap_per_condition(sce_pop, clusters, SETTINGS)

plot_umap_markers(sce_pop, type_markers(sce_pop), SETTINGS)

# Compare tSNE and PCA (good practice)
plotCodes(sce_pop, k = clusters)

##################################### 
### CLUSTER MERGING AND ANOTATION ###
#####################################

### Manual cluster merging and annotation based on heatmaps ###

sce_checkpoint <- object_path(SETTINGS, "sce_with_clustering_with_tsne_and_umap.RDS")
# saveRDS(sce, sce_checkpoint)

if (!exists("sce_pop")) {
    sce_pop <- readRDS(sce_checkpoint)
}

clusters <- "meta8"

sce_pop <- apply_manual_merging(sce_pop, clusters, SETTINGS, name = "CP_CyTOF-cluster-merging-1.csv")

plot_umap_merged(sce_pop, SETTINGS)

# If you want to check how many clusterings the SCE object contains:
colnames(metadata(sce_pop)$cluster_codes)

plot_heatmap_cellpops_2layers(sce_pop, clusters, SETTINGS)

plot_heatmap_cellpops_2layers_merged(sce_pop, clusters, features = pop_markers, SETTINGS)

plot_umap_manual_vs_auto(sce_pop, SETTINGS, meta = clusters)


### ### ### ### ### ### ### ### 
### DIFFERENTIAL ANALYSIS ###
### ### ### ### ### ### ### ### 
FDR_cutoff <- 0.1

# Import adjusted functions:
source("./Custom_functions/differential_analysis.R")
source("./Custom_functions/diffcyt_adjusted.R") 
source("./Custom_functions/plotDiffHeatmap_adjusted.R") 

### ### ### ### ### ### ### ### ### ### ### ### 
### Differential cell population abundance ###
stack_and_boxplots(sce_pop, SETTINGS)

export_prop_table(sce = sce_pop, merge = "merge1", name = "prop_table", SETTINGS)
export_absolute_table(sce = sce_pop, merge = "merge1", name = "absolute_table", SETTINGS)
export_absolute_props(sce = sce_pop, merge = "merge1", name = "total_prop_table", SETTINGS)

da_formulas <- create_da_formulas(sce_pop, SETTINGS) # We don't include age, gender, etc as covariates
da_formula1 <- da_formulas$da_formula1
da_formula2 <- da_formulas$da_formula2

if (SETTINGS$zoom_pop == "Tcells") {
    different = F
    abs = F
}
if (SETTINGS$zoom_pop == "NKB") {
    different = T
    abs = T
}

### MS vs control:
differential_abundance(sce_pop, SETTINGS, contrast = "conditionms", formula = da_formula1, 
                       clustering_to_use = "merge1", group1 = "control", group2 = "ms",
                       different_prop = different, # if T, proportions are calculated out of total CD45+ cells
                       name = "da_res_ms")

### AD vs control:
differential_abundance(sce_pop, SETTINGS, contrast = "conditionad", formula = da_formula2, 
                       clustering_to_use = "merge1", group1 = "control", group2 = "ad",
                       different_prop = different, # if T, proportions are calculated out of total CD45+ cells
                       name = "da_res_ad")

### AD vs MS
differential_abundance(sce_pop, SETTINGS, contrast = "conditionad", formula = da_formula1, 
                       clustering_to_use = "merge1", group1 = "ms", group2 = "ad",
                       contrast_extra = T, contrast_extra_con = "conditionms", 
                       different_prop = different, # if T, proportions are calculated out of total CD45+ cells
                       name = "da_res_advsms")

### Non-parametric tests
df <- get_freqs(sce_pop, merge = "merge1") # Get frequencies

kruskal_per_cluster(df, SETTINGS, name = "kruskall") ## Do Kruskall-Wallis per cluster and export

posthoc_pairwise(df) # Post hoc pairwise control, MS, AD; not exporting

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Differential analysis of marker expression stratified by cell population ###
plot_median_expression(sce_pop, SETTINGS)

ds_formulas <- create_ds_formulas(sce_pop) # We don't include age, gender, etc as covariates
ds_formula1 <- ds_formulas$ds_formula1

state_dt <- select_markers(
    markers_of_interest = paste0("CP_CyTOF_state_markers_", SETTINGS$zoom_pop, ".csv"))


### MS vs control
differential_expression(sce_pop, SETTINGS, contrast = "conditionms", formula = ds_formula1, 
                        clustering_to_use = "merge1", group1 = "control", group2 = "ms",
                        name = "ds_res_ms")

### AD vs control
differential_expression(sce_pop, SETTINGS, contrast = "conditionad", formula = ds_formula1, 
                        clustering_to_use = "merge1", group1 = "control", group2 = "ad",
                        name = "ds_res_ad")

### AD vs MS
differential_expression(sce_pop, SETTINGS, contrast = "conditionad", formula = ds_formula1, 
                        clustering_to_use = "merge1", group1 = "ms", group2 = "ad",
                        contrast_extra = T, contrast_extra_con = "conditionms",
                        name = "ds_res_advsms")

plot_heatmap_expresssion(sce_pop, fdr = FDR_cutoff, SETTINGS)

###########################
### PLOTTING DA RESULTS ###
###########################
source("./Custom_functions/plot_percentages.R")

### Create frequency file
freq_file <- create_freq_file(abs_prop = abs, SETTINGS)
freqs <- freq_file$freqs
med_sum <- freq_file$med_sum
p_vals_table <- freq_file$p_vals_table

### Organise frequency file by abundance, per tissue:
org_freqs <- organise_freqs(freqs, SETTINGS)
abundant <- org_freqs$abundant
mid_ab <- org_freqs$mid_ab
low_ab <- org_freqs$low_ab
non_ab <- org_freqs$non_ab

### Plot
w <- length(unique(abundant$pop))
set_dev(plot_path(SETTINGS, "perc_abundant"), height=3, width = 1+w*0.5)
plot_perc(abundant, p_vals_table, settings = SETTINGS, abs_prop = abs)
unset_dev()

w <- length(unique(non_ab$pop))
set_dev(plot_path(SETTINGS, "perc_non_abundant"), height=3, width = 1+w*0.5)
plot_perc(non_ab, p_vals_table, settings = SETTINGS, abs_prop = abs)
unset_dev()

if (!is.null(mid_ab)) {
    w <- length(unique(mid_ab$pop))
    set_dev(plot_path(SETTINGS, "perc_mid_abundant"), height=3, width = 1+w*0.5)
    plot_perc(mid_ab, p_vals_table, settings = SETTINGS, abs_prop = abs)
    unset_dev()
}

# Plot percentages for each population individually/separately:
plot_perc_sep(abundant, p_vals_table, settings = SETTINGS, accuracy = 0.1)
plot_perc_sep(non_ab, p_vals_table, settings = SETTINGS, accuracy = 0.01)

###########################
### PLOTTING DE RESULTS ###
###########################
source("./Custom_functions/plot_exprs.R")

### Create expression file (p-value < 0.1)
expression_file <- create_expr_file(SETTINGS)
expr_for_plot <- expression_file$expr_for_plot # for plotting below
expression <- expression_file$expression # for obtaining median values if you're curious
p_vals_table <- expression_file$p_vals_table

### Plot
set_dev(plot_path(SETTINGS, "plots_expr"), height=3, width = 3)
plot_expr(expr_for_plot, p_vals_table, SETTINGS) # adj_p_value; ref_p = c(0.1, 0.05, 0.01), label = c("*", "**", "***")
unset_dev()

