######## ----------------- Functions ----------------- ########

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

# transform SingleCellExperiment
transformData <- function (sce,
                           cf = 5,
                           ain = "counts",
                           aout = "exprs") {
  y <- assay(sce, ain)
  chs <- CATALYST::channels(sce)
  stopifnot(is.numeric(cf), cf > 0)
  if (length(cf) == 1) {
    int_metadata(sce)$cofactor <- cf
    cf <- rep(cf, nrow(sce))
  }
  else {
    stopifnot(!is.null(names(cf)), chs %in% names(cf))
    cf <- cf[match(chs, names(cf))]
    int_metadata(sce)$cofactor <- cf
  }
  fun <- asinh
  op <- "/"
  y <- fun(sweep(y, 1, cf, op))
  assay(sce, aout, FALSE) <- y
  sce
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
