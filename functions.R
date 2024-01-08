######## ----------------- Subset SummarizedExperiment ----------------- ########

subset_sce <- function(sce, coldata_column, coldata_value){
  sceEI <- sceEI <- CATALYST::ei(sce)
  # subset sce
  subset_sce <- sce[,sce[[coldata_column]] == coldata_value]
  # droplevels coldata
  coldata <- as.data.frame(SummarizedExperiment::colData(subset_sce))
  coldata <- sapply(coldata, droplevels)
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
