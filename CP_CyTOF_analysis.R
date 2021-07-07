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
library(gplots)
library(RColorBrewer)
library(devtools)

setwd("~/Documents/Data analysis/CyTOF/CP_CyTOF_data/")

############################
### CHOOSE YOUR ANALYSIS ###
############################

# Make list with all options
SETTINGS <- list(
    experiment_dir = normalizePath("."),
    SAVE_PLOTS = T, # Export plots
    SAVE_TABLES = T, # Export tables
    norm = F, # Normalisation
    subsetting = T, # Subsetting samples by cell number
    zoom = F, # Zoom into cell pops?
    zoom_pop = NULL, # Zoom into cell pops?
    tissue = 'septum' # Choose tissue
)

source("./Custom_functions/CP_CyTOF_functions.R") # Load created functions

###################
### DATA IMPORT ###
###################

source("./Custom_functions/data_import.R")

experiment <- load_experiment(SETTINGS)

sce <- construct_sce(experiment)

########################
### DIAGNOSTIC PLOTS ###
########################

source("./Custom_functions/diagnostic_plots.R")

plot_marker_expression_distributions(sce, SETTINGS)

plot_marker_expression_distributions(sce, SETTINGS, do_log = TRUE)

plot_n_cells_per_sample(sce, SETTINGS)

plot_mds(sce, SETTINGS)

plot_heatmap(sce, SETTINGS)

plot_redundancy_score(sce, SETTINGS)

####################
###  CLUSTERING  ###
####################

source("./Custom_functions/clustering.R")

set.seed(1234) # set random seed for reproducibility

sce <- run_clustering(sce, SETTINGS, skip_checkpoint = T)

plot_elbow_plot(sce, SETTINGS)

# ADJUST k:
clusters <- "meta20"

plot_heatmap_cellpops(sce, SETTINGS, clusters, name = "heatmap_cellpops")

plot_heatmap_cellpops(sce, SETTINGS, clusters, bin_anno = T, name = "heatmap_cellpops_anno")

plot_distribution_markers(sce, clusters, SETTINGS)

plot_specific_marker_heatmap(sce, "CD45RA", clusters, SETTINGS)

sce <- run_tsne_and_umap(sce, SETTINGS, skip_checkpoint = T)

plot_tsne_and_umap(sce, clusters, SETTINGS)

plot_umap_per_sample(sce, clusters, SETTINGS)

plot_umap_per_condition(sce, clusters, SETTINGS)

plot_umap_markers(sce, type_markers(sce), SETTINGS)

##################################### 
### CLUSTER MERGING AND ANOTATION ###
#####################################

### Manual cluster merging and annotation based on heatmaps ###

sce_checkpoint <- object_path(SETTINGS, "sce_with_clustering_with_tsne_and_umap.RDS")
# saveRDS(sce, sce_checkpoint)

if (!exists("sce")) {
    sce <- readRDS(sce_checkpoint)
}

clusters <- "meta20"

sce <- apply_manual_merging(sce, clusters, SETTINGS)

plot_umap_merged(sce, SETTINGS)
if (SETTINGS$tissue == "septum") {
    sce_immune <- filterSCE(sce, cluster_id != "CD45-", k = "merge1")
    plot_umap_merged(sce_immune, SETTINGS, name = "umap_merged_clusters_immune")
}


# If you want to check how many clusterings the SCE object contains:
colnames(metadata(sce)$cluster_codes)

plot_heatmap_cellpops_2layers(sce, clusters, SETTINGS)

plot_heatmap_cellpops_2layers_merged(sce, clusters, features = "type", SETTINGS)
if (SETTINGS$tissue == "septum") {
    plot_heatmap_cellpops_2layers_merged(sce_immune, clusters, 
                                         features = "type", SETTINGS,
                                         name = "heatmap_2layers_merged_immune")
}

plot_umap_manual_vs_auto(sce, SETTINGS)

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
stack_and_boxplots(sce, SETTINGS)

# Proportion table (out of CD45+ immune cells):
export_prop_table(sce = sce, merge = "merge1", name = "prop_table", SETTINGS)
export_absolute_table(sce = sce, merge = "merge1", name = "absolute_table", SETTINGS)

da_formulas <- create_da_formulas(sce, SETTINGS) # We don't include age, gender, etc as covariates
da_formula1 <- da_formulas$da_formula1
da_formula2 <- da_formulas$da_formula2

### MS vs control:
differential_abundance(sce, SETTINGS, contrast = "conditionms", formula = da_formula1, 
                       clustering_to_use = "merge1", group1 = "control", group2 = "ms",
                       name = "da_res_ms")

### AD vs control:
differential_abundance(sce, SETTINGS, contrast = "conditionad", formula = da_formula2, 
                       clustering_to_use = "merge1", group1 = "control", group2 = "ad",
                       name = "da_res_ad")

### AD vs MS
differential_abundance(sce, SETTINGS, contrast = "conditionad", formula = da_formula1, 
                       clustering_to_use = "merge1", group1 = "ms", group2 = "ad",
                       contrast_extra = T, contrast_extra_con = "conditionms", 
                       name = "da_res_advsms")

### Non-parametric tests
df <- get_freqs(sce, merge = "merge1") # Get frequencies

kruskal_per_cluster(df, SETTINGS, name = "kruskall") ## Do Kruskall-Wallis per cluster and export

posthoc_pairwise(df) # Post hoc pairwise control, MS, AD; not exporting

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Differential analysis of marker expression stratified by cell population ###
plot_median_expression(sce, SETTINGS)

ds_formulas <- create_ds_formulas(sce) # We don't include age, gender, etc as covariates
ds_formula1 <- ds_formulas$ds_formula1

state_dt <- select_markers(markers_of_interest = "CP_CyTOF_state_markers_per_pop.csv")

### MS vs control
differential_expression(sce, SETTINGS, contrast = "conditionms", formula = ds_formula1, 
                       clustering_to_use = "merge1", group1 = "control", group2 = "ms",
                       name = "ds_res_ms")

### AD vs control
differential_expression(sce, SETTINGS, contrast = "conditionad", formula = ds_formula1, 
                        clustering_to_use = "merge1", group1 = "control", group2 = "ad",
                        name = "ds_res_ad")

### AD vs MS
differential_expression(sce, SETTINGS, contrast = "conditionad", formula = ds_formula1, 
                        clustering_to_use = "merge1", group1 = "ms", group2 = "ad",
                        contrast_extra = T, contrast_extra_con = "conditionms",
                        name = "ds_res_advsms")

plot_heatmap_expresssion(sce, fdr = FDR_cutoff, SETTINGS)


###########################
### PLOTTING DA RESULTS ###
###########################
source("./Custom_functions/plot_percentages.R")

### Create frequency file
freq_file <- create_freq_file(abs_prop = T, SETTINGS)
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
plot_perc(abundant, p_vals_table, settings = SETTINGS)
unset_dev()

w <- length(unique(low_ab$pop))
set_dev(plot_path(SETTINGS, "perc_low_abundant"), height= 3, width = 1+w*0.5)
plot_perc(low_ab, p_vals_table, settings = SETTINGS)
unset_dev()

if (!is.null(mid_ab)) {
    w <- length(unique(mid_ab$pop))
    set_dev(plot_path(SETTINGS, "perc_mid_abundant"), height=3, width = 1+w*0.5)
    plot_perc(mid_ab, p_vals_table, settings = SETTINGS)
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
set_dev(plot_path(SETTINGS, "plots_expr"), height=5, width = 5)
plot_expr(expr_for_plot, p_vals_table, SETTINGS) # adj_p_value; ref_p = c(0.1, 0.05, 0.01), label = c("*", "**", "***")
unset_dev()


