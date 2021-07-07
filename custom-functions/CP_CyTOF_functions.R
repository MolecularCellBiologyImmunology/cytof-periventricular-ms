### ### ### ### ### ### 
### CREATE FUNCTIONS ###
### ### ### ### ### ### 


####################
### EXPORT PLOTS ###

object_path <- function(settings, filename) {
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
    file_path <- paste0(
        settings$experiment_dir, "/", settings$tissue, "/", subset_folder, normalised_folder, zoom_folder, filename
    )
    
    return(file_path)
}

# Set directory 
plot_path <- function(settings, name) {
    return(object_path(settings, paste0(settings$tissue, "_", name, ".pdf")))
}

# Open
set_dev <- function(file_path, settings=SETTINGS,height=10, width=10, useDingbats = F) {
    if (settings$SAVE_PLOTS) {
        dir.create(dirname(file_path), showWarnings = F, recursive = T)
        pdf(file_path, height = height, width = width, useDingbats = useDingbats)
    }
}
# Close: 
unset_dev <- function(settings=SETTINGS) {
    if (settings$SAVE_PLOTS) {
        dev.off()
    } 
}


#####################
### EXPORT TABLES ###
table_path <- function(settings, dir_ = "./",
                       name, format = ".csv"
                      ) {
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
    tab_path <- paste0(
        dir_, settings$tissue, "/", subset_folder, normalised_folder, zoom_folder, settings$tissue, "_", name, format
    )
    return(tab_path)
}

# Export table
exp_tab <- function(table, tab_path, settings = SETTINGS) {
    if (settings$SAVE_TABLES) {
        dir.create(dirname(tab_path), showWarnings = F, recursive = T)
        fwrite(x = table, file = tab_path)
    }
}

# # Sink table
# # Start writing to an output file
# sink_open <- function(sink_path, settings = SETTINGS){
#     sink(sink_path)
# }


###############################
### RPHENOGRAPH ###

### Heatmap ###
# Calculate mean expression per marker per cluster:
pheno_hm <- function(dt) {
    hm <- melt(dt, id.vars = "phenograph_cluster")
    hm <- as.data.table(hm)
    hm <- hm[,
             mean(value),
             by = c("phenograph_cluster", "variable")]
    # Build matrix for heatmap:
    hm <- dcast(hm, phenograph_cluster ~ variable, value.var = "V1")
    # Plot heatmap:
    pal <- colorRampPalette(rev(brewer.pal(n=8, name="RdYlBu")))
    heatmap.2(as.matrix(hm[,2:22]),
              # Rowv = TRUE,
              Colv = FALSE,
              scale = "col",
              col = pal,
              # dendrogram = "row",
              trace = "none",
              key.title = NA,
              symm = F,
              keysize = 1.3
              # margins = c(5, 5.5),
    )
}

### t-SNE ###
### tSNE (option 1)
pheno_tsne_1 <- function(dt, perplexity = 30, max_iter = 3000) {
    # Set seed for reproducibility:
    set.seed(42)
    # Calculate tSNE, adjusting hyperparameters (https://distill.pub/2016/misread-tsne/):
    tsne_out <- Rtsne(dt[,1:22], 
                      perplexity = perplexity, # range 5-50
                      max_iter = max_iter # iterations
    )
    # Plot the result
    res <- as.data.table(tsne_out$Y)
    res <- res[, phenograph_cluster := dt$phenograph_cluster]
    res <- res[sample(1:nrow(dt), 1000), ] # sample random subset
    ggplot(res, aes(x = V1, y = V2, color = phenograph_cluster)) +
        geom_point() +
        theme_classic()
}



### ### ###
### tSNE (option 2)
# tsne with wrapper function that uses the Rtsne package and ggplot2
# from "Quick and easy t-SNE analysis in R"
# https://www.r-bloggers.com/quick-and-easy-t-sne-analysis-in-r/
pheno_tsne_2 <- function(dt, n_subsample = 1000, k = 30, perplex = 30, marker = 'CD3') {
    dt2 <- assay(sce, "exprs")
    dt2 <- dt2[, sample(1:ncol(dt2), n_subsample)] # sample random subset
    tdt2 <- t(as.matrix(dt2))
    Rphenograph_out2 <- Rphenograph(tdt2, k = k)
    tdt2$phenograph_cluster <- factor(membership(Rphenograph_out2[[2]]))
    
    # Visualise clusters
    tsne(dt2,
         labels=as.factor(tdt2$phenograph_cluster),
         perplex = perplex,
         dotsize = 1,
         axistextsize = 12,
         legendtextsize = 12
    )
}

# Visualise one marker
pheno_tsne_2_marker <- function(dt, marker = 'CD3') {
    tsne(dt,
         labels = scale(dt2[row.names(dt)== marker]),
         controlscale = TRUE,
         dotsize = 1,
         axistextsize = 12,
         legendtextsize = 12,
         scale=2)
}






