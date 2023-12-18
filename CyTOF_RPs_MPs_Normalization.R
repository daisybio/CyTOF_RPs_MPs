##################### Normalization of RPs/MPs Data #####################

######## ----------------- Install packages ----------------- ########

required_packages <- c("data.table", "BiocManager", "ggplot2")
for(package in required_packages){
  if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
  library(package, character.only = TRUE, quietly = TRUE)
}
BiocManager::install("CATALYST", force = TRUE)
library(CATALYST)

BiocManager::install("flowCore", force = TRUE)
library(flowCore)

######## ----------------- Read files ----------------- ########
path_to_data <- "data/"

panel <- fread(file.path(path_to_data, "panel.csv"))
md <- fread(file.path(path_to_data, "metafile RP MP all.csv"))
sce <- prepData(file.path(path_to_data, "RPs_MPs_files"), 
                panel, 
                md, 
                md_cols = list(file = "file_name", id = "sample_id", factors = c("activation", "patient_id", "subgroup", "subgroup2")))

saveRDS(sce, file.path(path_to_data, "sce_original_RPs_MPs.rds"))


######## ----------------- Functions ----------------- ########

get_patient_FCs <- function(sce, ain = "exprs"){
  exprsData <- assays(sce)[[ain]]
  exprsData <- data.frame(t(exprsData), colData(sce), check.names = FALSE)
  gg_df <- reshape2::melt(exprsData, value.name = "expression", variable.name = "marker", id.vars = names(colData(sce)))
  gg_df <- as.data.table(gg_df)
  medians <- merge(gg_df[subgroup2 == "RP", median(expression), by=c("marker", "patient_id")], gg_df[subgroup2 == "MP", median(expression), by=c("marker", "patient_id")], by=c("marker", "patient_id"))
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

normalize_patient_wise_sce <- function(sce, marker = "CD42b", ain = "exprs", aout = "norm_exprs"){
  exprs <- assays(sce)[[ain]]
  col_dt <- as.data.frame(colData(sce))
  patient_ids <- unique(col_dt$patient_id)
  
  norm_exprs <- exprs
  for (id in patient_ids){
    cols_MP <-  col_dt$patient_id == id & col_dt$subgroup2 == "MP"
    cols_RP <-  col_dt$patient_id == id & col_dt$subgroup2 == "RP"
    median_RP <- median(exprs[row.names(exprs) == marker, cols_RP], na.rm=TRUE)
    median_MP <- median(exprs[row.names(exprs) == marker, cols_MP], na.rm=TRUE)
    fold_change <- median_RP/median_MP
    norm_exprs[,cols_MP] <- exprs[,cols_MP] * fold_change
  }
  assay(sce, aout, TRUE) <- norm_exprs
  return(sce)
}

######## ----------------- Fold Change Plot ----------------- ########

medians <- get_patient_FCs(sce)
ggplot(medians[!marker %in% c("CD40", "CD154", "CD3", "CD107a"),], aes(x = reorder(marker, foldchange), y = foldchange))+
  geom_boxplot()+
  geom_jitter(aes(shape=patient_id, color = patient_id))+
  scale_shape_manual(values=1:nlevels(medians$patient_id)) + labs(x = "Marker", y = "Fold Change", color = "Patient ID", shape = "Patient ID") 

######## ----------------- Normalize MP Expression Values (based on CD42b) ----------------- ########

sce <- normalize_patient_wise_sce(sce, ain = "exprs", aout = "exprs")

medians_norm_new <- get_patient_FCs(sce, ain = "exprs")
medians_norm_new[marker == "CD42b",]

saveRDS(sce, file.path(path_to_data, "sce_CD42b_RPs_MPs.rds"))

#assays(sce)$exprs[1:10,1:10]

# try normalizing on counts + arcsinh

#sce <- normalize_patient_wise_sce(sce, ain = "counts", aout = "norm_counts")
#assays(sce)$norm_counts[1:10,1:10]

#medians_norm_new <- get_patient_FCs(sce, ain = "norm_counts")
#medians_norm_new[marker == "CD42b",]

#sce <- transformData(sce, ain = "norm_counts", aout = "norm_counts_exprs")
#assays(sce)$norm_counts_exprs[1:10,1:10]

######## ----------------- Check normalization values ----------------- ########

#sce_norm_old <- readRDS(file.path(path_to_data, "sce_CD42b_RPs_MPs_old.rds"))
#old_exprs <- assays(sce_norm_old)$exprs
#old_counts <- assays(sce_norm_old)$counts

#old_exprs[1:10,1:10]

#old_medians <- get_patient_FCs(sce_norm_old)
#old_medians[marker == "CD42b",]


######## ----------------- Write FCS Files ----------------- ########
norm_out_dir <- file.path(path_to_data, "RPs_MPs_files_normalized")
dir.create(norm_out_dir)

for(id in unique(sce$sample_id)){
  print(id)
  # extract the normalized expr
  norm_exprs <- assay(sce[, sce$sample_id == id], "exprs")
  # first order the rows
  norm_exprs <- norm_exprs[panel$antigen,]
  row.names(norm_exprs) <- panel$fcs_colname
  data_df <- t(norm_exprs)
  ff <- flowFrame(data_df)
  write.FCS(ff, paste0(norm_out_dir,"/",id, ".fcs"))
}
