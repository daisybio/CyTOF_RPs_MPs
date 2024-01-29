######## ----------------- Subset SummarizedExperiment ----------------- ########

subset_sce <- function(sce, coldata_column, coldata_value){
  sceEI <- sceEI <- CATALYST::ei(sce)
  # subset sce
  subset_sce <- sce[,sce[[coldata_column]] == coldata_value]
  # droplevels coldata
  coldata <- as.data.frame(SummarizedExperiment::colData(subset_sce))
  coldata[[coldata_column]] <- droplevels(coldata[, coldata_column]) 
  coldata[["sample_id"]] <- droplevels(coldata[, "sample_id"]) 
  SummarizedExperiment::colData(subset_sce) <- S4Vectors::DataFrame(as.data.frame(coldata))
  # subset experimental info
  sceEI_subset <- sceEI[sceEI$activation == "baseline",]
  metadata(subset_sce)$experiment_info <- sceEI_subset
  return(subset_sce)
}


######## ----------------- Functions for CD42b normalization ----------------- ########

get_patient_FCs <- function(sce, column = "subgroup2", ain = "exprs"){
  exprsData <- assays(sce)[[ain]]
  exprsData <- data.frame(t(exprsData), colData(sce), check.names = FALSE)
  gg_df <- reshape2::melt(exprsData, value.name = "expression", variable.name = "marker", id.vars = names(colData(sce)))
  gg_df <- as.data.table(gg_df)
  medians <- merge(gg_df[get(column) == "RP", median(expression), by=c("marker", "patient_id")], gg_df[get(column) == "MP", median(expression), by=c("marker", "patient_id")], by=c("marker", "patient_id"))
  colnames(medians) <- c("marker", "patient_id", "RP_median", "MP_median")
  medians[, foldchange := RP_median/MP_median]
  return(medians)
}

normalize_patient_wise_sce <- function(sce, marker = "CD42b", column = "subgroup2", ain = "exprs", aout = "norm_exprs"){
  exprs <- assays(sce)[[ain]]
  col_dt <- as.data.frame(colData(sce))
  patient_ids <- unique(col_dt$patient_id)
  
  norm_exprs <- exprs
  for (id in patient_ids){
    cols_MP <-  col_dt$patient_id == id & col_dt[[column]] == "MP"
    cols_RP <-  col_dt$patient_id == id & col_dt[[column]] == "RP"
    median_RP <- median(exprs[row.names(exprs) == marker, cols_RP], na.rm=TRUE)
    median_MP <- median(exprs[row.names(exprs) == marker, cols_MP], na.rm=TRUE)
    fold_change <- median_RP/median_MP
    norm_exprs[,cols_MP] <- exprs[,cols_MP] * fold_change
  }
  assay(sce, aout, TRUE) <- norm_exprs
  return(sce)
}


######## ----------------- Violin Plots ----------------- ########

plotViolinMod <- function (x, k = "meta20", features = "state", assay = "exprs", 
                           fun = c("mean", "sum", "median"), facet_by = "antigen", 
                           color_by = "condition", shape_by = NULL, size_by = FALSE, 
                           geom = "violins", jitter = TRUE, ncol = NULL) 
{
  
  fun <- match.arg(fun)
  geom <- match.arg(geom)
  facet_by <- match.arg(facet_by)
  stopifnot(is.logical(jitter), length(jitter) == 1)
  if (!is.null(ncol)) 
    stopifnot(is.numeric(ncol), length(ncol) == 1, ncol %% 1 == 0)
  
  if (facet_by == "cluster_id") {
    CATALYST:::.check_sce(x, TRUE)
    k <- CATALYST:::.check_k(x, k)
  } else {
    CATALYST:::.check_sce(x)
  }
  
  CATALYST:::.check_assay(x, assay)
  CATALYST:::.check_cd_factor(x, color_by)
  CATALYST:::.check_cd_factor(x, shape_by)
  shapes <- CATALYST:::.get_shapes(x, shape_by)
  if (is.null(shapes)) 
    shape_by <- NULL
  
  x <- x[CATALYST:::.get_features(x, features), ]
  
  if (facet_by == "cluster_id") {
    x$cluster_id <- cluster_ids(x, k)
    by <- c("cluster_id", "sample_id")
  } else {
    by <- "sample_id"
  }
  
  ms <- CATALYST:::.agg(x, by, fun, assay)
  
  df <- reshape2::melt(ms, varnames = c("antigen", by[length(by)]))
  df <- reshape2::melt(ms, varnames = c("antigen", by[length(by)]))
  df[[by[length(by)]]] <- factor(df[[by[length(by)]]], levels = unique(colData(x)[[by[length(by)]]]))
  
  if (length(by) == 2) 
    names(df)[ncol(df)] <- "cluster_id"
  
  x_var <- ifelse(facet_by == "antigen", color_by, "antigen")
  
  
  if (!is.null(df$cluster_id)) 
    df$cluster_id <- factor(df$cluster_id, levels = x$cluster_id)
  
  i <- match(df$sample_id, x$sample_id)
  j <- setdiff(names(colData(x)), c(names(df), "cluster_id"))
  df <- cbind(df, colData(x)[i, j, drop = FALSE])
  ncs <- table(as.list(colData(x)[by]))
  ncs <- rep(c(t(ncs)), each = nrow(x))
  if (size_by) {
    size_by <- "n_cells"
    df$n_cells <- ncs
  } else {
    size_by <- NULL
  }
  
  df <- df[ncs > 0, , drop = FALSE]
  
  ggplot(df, aes_string(x_var, "value", col = color_by, fill = color_by)) + 
    facet_wrap(facet_by, ncol = ncol, scales = "free_y") + 
    geom_violin(alpha = 0.3) +
    geom_point(alpha = 0.8, position = (if (jitter) {
      position_jitterdodge(jitter.width = 0.2, jitter.height = 0)
    } else "identity"), aes_string(fill = color_by, size = size_by, shape = shape_by)) +
    scale_shape_manual(values = shapes) + 
    scale_size_continuous(range = c(0.5, 3)) + 
    guides(size = guide_legend(order = 4), 
           shape = guide_legend(order = 3, override.aes = list(size = 3)), 
           col = guide_legend(order = 2, override.aes = list(alpha = 1, size = 3))) + 
    ylab(paste(fun, ifelse(assay == "exprs", "expression", assay))) + 
    theme_bw() + 
    theme(legend.key.height = unit(0.8, "lines"), 
          axis.text = element_text(color = "black"), 
          strip.text = element_text(face = "bold"), 
          strip.background = element_rect(fill = NA, color = NA), 
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_line(color = "grey", size = 0.2)) + 
    if (length(unique(c(x_var, color_by, facet_by))) == 1) {
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    } else {
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    }
}


######## ----------------- DE Results Table ----------------- ########

print_DE_top_table <- function(sce, de_res, method){
  out <- de_res[[method]]
  topTableOut <- data.frame(out)
  eff_r <- de_res$effect_size
  if(!is.null(eff_r)){
    eff_r[, marker_id := sapply(strsplit(eff_r$group2,'::'), "[", 1)]
    topTableOut <- merge(topTableOut, eff_r[, c("cluster_id", "marker_id", "overall_group","effsize", "magnitude")], by = c("cluster_id", "marker_id"), all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE)
    colnames(topTableOut) <- c(colnames(data.frame(out)), "overall_group","cohens_d", "magnitude")
    topTableOut$cohens_d <- formatC(topTableOut$cohens_d)
  }
  
  topTableOut$p_val <- formatC(topTableOut$p_val)
  topTableOut$p_adj <- formatC(topTableOut$p_adj)
  
  DT::datatable(topTableOut, rownames = FALSE, 
                options = list(pageLength = 10, searching = FALSE, 
                               columnDefs = list(list( targets = c(1,2), 
                                                       render = JS("function(data, type, row, meta) {","return data === null ? 'NA' : data;","}")))))
}

######## ----------------- Clustering Function ----------------- ########

clusterSCENew <-
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


