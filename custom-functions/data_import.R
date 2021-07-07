#' Loading the panel, the single-cell data and exploring normalisation options

# Load adjusted prepData function to allow for more metadata (age, gender, etc)
source("./Custom_functions/prepData_adjusted.R") 

### Load panel (it has the metal names swaped so we need to change them):
load_panel <- function(settings) {  
    panel <- fread(file.path(settings$experiment_dir, "CP_CyTOF_panel.csv"))
    panel[, name1 := sub('^([[:digit:]]+)([[:alpha:]]+)$', '\\1', fcs_colname)]
    panel[, name2 := sub('^([[:digit:]]+)([[:alpha:]]+)$', '\\2', fcs_colname)]
    
    pnames <- fread(file.path(settings$experiment_dir, "CP_CyTOF_panel_names.csv"))
    pnames[, name1 := sub('^([[:alpha:]]+)([[:digit:]]+)([[:alpha:]]+)$', '\\1', fcs_colname)]
    pnames[, name2 := sub('^([[:alpha:]]+)([[:digit:]]+)([[:alpha:]]+)$', '\\2', fcs_colname)]
    pnames[, name3 := sub('^([[:alpha:]]+)([[:digit:]]+)([[:alpha:]]+)$', '\\3', fcs_colname)]
    pnames <- pnames[fcs_colname != name1]
    
    panel <- merge(panel, pnames, by.x = "name1", by.y = "name2")
    panel <- panel[, c("fcs_colname.y", "antigen", "marker_class")]
    
    setnames(panel, "fcs_colname.y", "fcs_colname")
    return(panel)
}


#################
### LOAD DATA ###

load_data <- function(settings) {
    if (settings$tissue == 'cp') {
        if (settings$norm) {
            normalised_folder <- "normalised"
        } else {
            normalised_folder <- "original"
        }
        files <- list.files(file.path(settings$experiment_dir, "Pregated_fcs/Cp",
                                      normalised_folder), 
                            full.names = T)
        files <- files[grepl(".fcs$", files)] # keep only fcs files
        files <- files[!grepl("ref", files)] # exclude ref
        files <- files[!grepl("2019-110", files)] # exclude control6 (no cells)
        if (settings$subsetting == T) {
            fs <- list()
            ncells <- c()
            for (file in files) {
                cat("Reading File", file, "\n", sep = " ")
                temp <- read.FCS(filename = file)
                ncells <- c(ncells, dim(temp@exprs)[1])
            }
            n_subset <- quantile(ncells, 0.75)
            for (file in files) {   
                temp <- read.FCS(filename = file)
                if(nrow(temp) > n_subset) { 
                    lines <- sample(1:nrow(temp), n_subset)
                    temp <- temp[lines]
                } 
                print(dim(temp))
                fs[[file]] <- temp
            }
            fs <- as(fs, "flowSet")
        } else {
            fs <- read.flowSet(files = files,
                               transformation = F, truncate_max_range = F)
        }
        # fs <- flowCore::rbind2(fs1, fs2)
        
        ## Metadata:
        md <- fread(file.path(settings$experiment_dir, "Clinical_info/CP_CyTOF_metadata_cp.csv"))
        md <- md[patient_id != "2019-110"] # exclude control6 (no cells)
        for (i in 1:nrow(md)) {
            md[i, file_name := files[grepl(md[i]$patient_id, files)]]
        }
        batch <- fread(file.path(settings$experiment_dir, "CP_CyTOF_batches_cp.csv"))
        setnames(batch, "barcode", "batch_id")
        info <- fread(file.path(settings$experiment_dir, "Clinical_info/CP_CyTOF_clinical_info_selection_cp.csv"))
        info <- info[, c("donor", "age", "gender", "pmd_hours")]
        # Combine metadata of fixed effects
        md <- merge(md, batch, by = "patient_id")
        md <- merge(md, info, by.x = "patient_id", by.y = "donor")
        return(list(fs = fs, md = md))
    }
    if (settings$tissue == 'septum') {
        if (settings$norm) {
            normalised <- "normalised/"
        } else {
            normalised <- "original/"
        }
        files <- list.files(paste0(settings$experiment_dir, "/Pregated_fcs/Spt/",
                                   normalised), 
                            full.names = T)
        files <- files[grepl(".fcs$", files)] # keep only fcs files
        files <- files[!grepl("ref", files)] # exclude ref
        if (settings$subsetting == T) {
            fs <- list()
            ncells <- c()
            for (file in files) {
                cat("Reading File", file, "\n", sep = " ")
                temp <- read.FCS(filename = file)
                ncells <- c(ncells, dim(temp@exprs)[1])
            }
            n_subset <- quantile(ncells, 0.75)
            for (file in files) {   
                temp <- read.FCS(filename = file)
                if(nrow(temp) > n_subset) { 
                    lines <- sample(1:nrow(temp), n_subset)
                    temp <- temp[lines]
                } 
                print(dim(temp))
                fs[[file]] <- temp
            }
            fs <- as(fs, "flowSet")
        } else {
            fs <- read.flowSet(files = files,
                               transformation = F, truncate_max_range = F)
        }
        # fs <- flowCore::rbind2(fs1, fs2)
        
        md <- fread(file.path(settings$experiment_dir, "Clinical_info/CP_CyTOF_metadata_septum.csv"))
        # Keep only md from those samples that have septum
        md[, file_name := ""] 
        for (i in 1:nrow(md)) {
            md[i, file_name := {
                my_value <- files[grepl(md[i]$patient_id, files)]
                if (is_empty(my_value)) {
                    'no_septum'
                } else {
                    my_value
                }
            }
            ]
        }
        md <- md[file_name != "no_septum"]
        batch <- fread(file.path(settings$experiment_dir, "CP_CyTOF_batches_spt.csv"))
        setnames(batch, "barcode", "batch_id")
        info <- fread(file.path(settings$experiment_dir, "Clinical_info/CP_CyTOF_clinical_info_selection_septum.csv"))
        info <- info[, c("donor", "age", "gender", "pmd_hours")]
        # Combine metadata of fixed effects
        md <- merge(md, batch, by = "patient_id")
        md <- merge(md, info, by.x = "patient_id", by.y = "donor")
        return(list(fs = fs, md = md))    }
    if (settings$tissue == 'blood') {
        if (settings$norm) {
            normalised <- "normalised/"
        } else {
            normalised <- "original/"
        }
        files <- list.files(paste0(settings$experiment_dir, "/Pregated_fcs/Bld/",
                                   normalised), 
                            full.names = T)
        files <- files[grepl(".fcs$", files)] # keep only fcs files
        files <- files[!grepl("ref", files)] # exclude ref
        if (settings$subsetting == T) {
            fs <- list()
            ncells <- c()
            for (file in files) {
                cat("Reading File", file, "\n", sep = " ")
                temp <- read.FCS(filename = file)
                ncells <- c(ncells, dim(temp@exprs)[1])
            }
            n_subset <- quantile(ncells, 0.75)
            for (file in files) {   
                temp <- read.FCS(filename = file)
                if(nrow(temp) > n_subset) { 
                    lines <- sample(1:nrow(temp), n_subset)
                    temp <- temp[lines]
                } 
                print(dim(temp))
                fs[[file]] <- temp
            }
            fs <- as(fs, "flowSet")
        } else {
            fs <- read.flowSet(files = files,
                               transformation = F, truncate_max_range = F)
        }
        # fs <- flowCore::rbind2(fs1, fs2)
        
        md <- fread(file.path(settings$experiment_dir, "Clinical_info/CP_CyTOF_metadata_blood.csv"))
        # Keep only md from those samples that have septum
        md[, file_name := ""] 
        for (i in 1:nrow(md)) {
            md[i, file_name := {
                my_value <- files[grepl(md[i]$patient_id, files)]
                if (is_empty(my_value)) {
                    'no_bld'
                } else {
                    my_value
                }
            }
            ]
        }
        md <- md[file_name != "no_bld"]
        # for (i in 1:nrow(md)) {
        #     md[i, file_name := files[grepl(md[i]$patient_id, files)]]
        # }
        batch <- fread(file.path(settings$experiment_dir, "CP_CyTOF_batches_bld.csv"))
        setnames(batch, "barcode", "batch_id")
        info <- fread(file.path(settings$experiment_dir, "Clinical_info/CP_CyTOF_clinical_info_selection_blood.csv"))
        info <- info[, c("donor", "age", "gender", "pmd_hours")]
        # Combine metadata of fixed effects
        md <- merge(md, batch, by = "patient_id")
        md <- merge(md, info, by.x = "patient_id", by.y = "donor")
        return(list(fs = fs, md = md))   
    }
}


#####################
### NORMALISATION ###

normalise <- function(settings, exclude = c('CD39'), method="95p", transformation=FALSE) {
    source("./Custom_functions/BatchAdjust.R")
    
    # Create list of channels to adjust:
    panel_dt <- as.data.table(panel)
    exclude <- exclude # list antigens/channels to exclude
    panel_dt <- panel_dt[!(antigen %in% exclude),] # Exclude unreliable markers
    fwrite((panel_dt[, .(fcs_colname)]), 
           paste0(settings$experiment_dir, "/channels_to_adjust.txt"), 
           col.names = F)
    
    # Choose tissue folder name
    tissue <- settings$tissue
    folder_name <- if (tissue == 'cp') {
        'Cp'
    } else if (tissue == 'septum') {
        'Spt'
    } else if (tissue == 'blood') {
        'Bld'
    }
    
    BatchAdjust(
        basedir= paste0(settings$experiment_dir, "/Pregated_fcs/", 
                        folder_name, "/"),
        outdir= paste0(settings$experiment_dir, "/Pregated_fcs/", 
                       folder_name, "/normalised"),
        channelsFile = paste0(settings$experiment_dir, "/channels_to_adjust.txt"),
        batchKeyword=paste0("_", folder_name),
        anchorKeyword = "ref",
        method=method,
        transformation=transformation,
        addExt=NULL,
        plotDiagnostics=TRUE)
}


load_experiment <- function(settings) {
    #' Loads panel, fcs files and metadata.
    
    ### Load panel (marker_class varies: use latest version; or use cellpop-specific panel)
    panel <- load_panel(settings)
    
    ### NORMALISATION ###
    if (settings$norm == T) {
        normalise(exclude = c('CD39'), method="95p", transformation=FALSE)
    }
    
    ### Load fcs data into a flowset and metadata:
    # read.flowSet(), by default, may transform the marker intensities and
    # remove cells with extreme positive values. This behavior can be controlled
    # with arguments transformation and truncate_max_range, respectively.
    # Apply function to each tissue:
    data <- load_data(settings)
    fs <- data$fs
    md <- data$md
    
    # Check that all panel columns are in the flowSet object:
    all(panel$fcs_colname %in% colnames(fs))
    
    panel <- as.data.frame(panel)
    
    return(list(
        panel = panel,
        flowset = fs,
        metadata = md
    ))
}


construct_sce <- function(experiment) {
    md <- experiment$metadata
    fs <- experiment$flowset
    panel <- experiment$panel
    
    ### Specify levels for conditions & sample IDs to assure desired ordering:
    md$condition <- factor(md$condition, levels = c("control", "ms", "ad"))
    md$sample_id <- factor(md$sample_id,
                           levels = md$sample_id[order(md$condition, md$sample_id)])
    
    ### Construct SingleCellExperiment (default arcsinh transformation of
    # marker expressions with a cofactor of 5
    # Do the following manually because of error in function
    # when "merging" fs and md, otherwise sample names are wrong:
    md <- md[match(c(keyword(fs, "FILENAME")), md$file_name)]
    
    ### Construct sce:
    sce <- prepData_adjusted(fs, panel, md)#, features = panel$fcs_colname)
    
    return(sce)
}
