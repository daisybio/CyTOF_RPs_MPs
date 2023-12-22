# RP MP Clustering

##################### Normalization of RPs/MPs Data #####################

######## ----------------- Install packages ----------------- ########

#required_packages <- c("data.table", "BiocManager", "ggplot2", "stringr")
#for(package in required_packages){
#  if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
#  library(package, character.only = TRUE, quietly = TRUE)
#}
#BiocManager::install("CATALYST", force = TRUE)
#library(CATALYST)

source("functions.R")

######## ----------------- Prepare Data ----------------- ########

path_to_data <- "data/"

panel <- fread(file.path(path_to_data, "panel.csv"))

#built meta data file
files <- list.files(file.path(path_to_data, "CCS_baseline"))
md <- data.table(file_name = files)
md$sample_id <- sapply(strsplit(md$file_name,"[.]"), "[", 1)
md$sample_id <- str_replace(md$sample_id, "_platlet_specific", "")
md$sample_id <- str_replace(md$sample_id, "RPS_", "RPs") # inconsistent naming of files
md$type <- sapply(strsplit(md$sample_id,"_"), "[", 2)
md$patient_id <- sapply(strsplit(md$sample_id,"_"), "[", 1) 

sce <- prepData(file.path(path_to_data, "CCS_baseline"), 
                panel, 
                md, 
                transform = TRUE,
                cofactor = 5,
                md_cols = list(file = "file_name", id = "sample_id", factors = c("type", "patient_id")))

saveRDS(sce, file.path(path_to_data, "sce_original_RPs_MPs_rest.rds"))


######## ----------------- Should we really normalize MPs by RPs?  ----------------- ########

sce <- normalize_patient_wise_sce(sce, column = "type", ain = "exprs", aout = "exprs")

saveRDS(sce, file.path(path_to_data, "sce_CD42b_RPs_MPs_rest.rds"))

medians <- get_patient_FCs(sce, column = "type", ain = "exprs")
medians[marker == "CD42b",]

######## ----------------- Add Cell Annotation  ----------------- ########

annotation <- colData(sce)$type

######## ----------------- Perform Clustering  ----------------- ########

clusterSCE <-
  function (x,
            assayType = "exprs",
            features = "type",
            xdim = 10,
            ydim = 10,
            maxK = 20,
            verbose = TRUE,
            seed = 1){
    # set seed
    if (!is.null(seed)) {
      set.seed(seed)
    }
    # check input
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(
      is.logical(verbose),
      length(verbose) == 1,
      vapply(list(xdim,
                  ydim, maxK, seed), function(arg)
                    is.numeric(arg) &&
               length(arg) == 1, logical(1))
    )
    # extract markers to cluster
    features <- CATALYST:::.get_features(x, features)
    if (is.null(CATALYST::marker_classes(x))) {
      SummarizedExperiment::rowData(x)$marker_class <-
        factor(c("state", "type")[as.numeric(rownames(x) %in%
                                               features) + 1], levels = c("type", "state", "none"))
    }
    SummarizedExperiment::rowData(x)$used_for_clustering <- rownames(x) %in% features
    if (verbose)
      message("o running FlowSOM clustering...")
    fsom <-
      FlowSOM::ReadInput(flowCore::flowFrame(t(SummarizedExperiment::assay(x, assayType))))
    som <-
      FlowSOM::BuildSOM(
        fsom,
        colsToUse = features,
        silent = FALSE,
        xdim = xdim,
        ydim = ydim
      )
    som <- FlowSOM::BuildMST(som)
    if (verbose)
      message("o running ConsensusClusterPlus metaclustering...")
    pdf(NULL)
    mc <-
      suppressWarnings(
        ConsensusClusterPlus::ConsensusClusterPlus(
          t(som$map$codes),
          maxK = maxK,
          reps = 100,
          distance = "euclidean",
          seed = seed,
          plot = NULL
        )
      )
    dev.off()
    # PAC
    Kvec <- 2:maxK
    x1 <- 0.1
    x2 <- 0.9 # threshold defining the intermediate sub-interval
    PAC <- rep(NA, length(Kvec))
    names(PAC) <- Kvec
    ecdf_list <- list()
    for(i in Kvec){
      consensus_matrix <- mc[[i]]$consensusMatrix
      consensus_values <- consensus_matrix[lower.tri(consensus_matrix)]
      ecdf_data <- ecdf(consensus_values)
      PAC[i-1] <- ecdf_data(x2) - ecdf_data(x1)
      ecdf_list[[as.character(i)]] <- data.frame(ConsensusIndex = consensus_values, CDF = ecdf_data(consensus_values), Cluster = rep(i, length(consensus_values)), PAC = rep(PAC[i-1], length(consensus_values)))
      
    }
    # the optimal K
    optK <- Kvec[which.min(PAC)]
    message(paste0("The optimal K (according to PAC): ", optK))
    
    # plot
    # Combine all the data frames into one
    ecdf_data_combined <- do.call(rbind, ecdf_list)
    
    # Convert the 'Clusters' column to a factor for plotting
    ecdf_data_combined$Cluster <- factor(ecdf_data_combined$Cluster, levels = sort(unique(ecdf_data_combined$Cluster)))
    ecdf_data_combined$Legend <- interaction(ecdf_data_combined$Cluster, round(ecdf_data_combined$PAC, 3), sep = " - PAC: ", lex.order = TRUE)
    
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual",]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    # Plot the ECDFs for each k
    ggplot(ecdf_data_combined, aes(x = ConsensusIndex, y = CDF, color = Legend)) +
      geom_line() +
      labs(x = "Consensus Index Value", y = "CDF",
           title = "Consensus matrix CDFs") +
      theme_minimal() + scale_color_manual(name = "Cluster", values = col_vector) +
      geom_vline(xintercept = x1, color = "black", linetype = "dotdash") + 
      geom_vline(xintercept = x2, color = "black", linetype = "dotdash") 
    
    ggplot(data.frame(Cluster = factor(names(PAC), levels = names(PAC)), PAC = PAC), aes( x = Cluster, y = PAC, group= 1)) + geom_point() + geom_path()

    k <- xdim * ydim
    mcs <- seq_len(maxK)[-1]
    codes <-
      data.frame(seq_len(k), purrr::map(mc[-1], "consensusClass"))
    codes <-
      dplyr::mutate_all(codes, function(u)
        factor(u, levels = sort(unique(u))))
    colnames(codes) <- c(sprintf("som%s", k), sprintf("meta%s",
                                                      mcs))
    x$cluster_id <- factor(som$map$mapping[, 1])
    S4Vectors::metadata(x)$cluster_codes <- codes
    S4Vectors::metadata(x)$SOM_codes <- som$map$codes
    S4Vectors::metadata(x)$SOM_medianValues <- som$map$medianValues
    S4Vectors::metadata(x)$SOM_MST <- som$MST
    S4Vectors::metadata(x)$delta_area <- CATALYST:::.plot_delta_area(mc)
    x <-
      CATALYST::mergeClusters(
        x,
        k = sprintf("meta%s", maxK),
        id = "all",
        table = data.frame(old_cluster = seq_len(maxK), new_cluster = "all")
      )
    return(x)
  }


sce <- clusterSCE(sce, features = "state", seed = 1234, maxK = 40)

CATALYST::delta_area(sce)

saveRDS(sce, file.path(path_to_data, "sce_clustered_RPs_MPs_rest.rds"))

