##################### Normalization of RPs/MPs Data #####################
path_to_data <- "Cyanus_RPsMPs/"

######## ----------------- Install packages ----------------- ########

required_packages <- c("data.table", "BiocManager", "ggplot2")
for(package in required_packages){
  if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
  library(package, character.only = TRUE, quietly = TRUE)
}
BiocManager::install("CATALYST")
library(CATALYST)

BiocManager::install("flowCore")
library(flowCore)

######## ----------------- Read files ----------------- ########

panel <- fread(paste0(path_to_data, "panel.csv"))
md <- fread(paste0(path_to_data, "metafile RP MP all.csv"))
sce <- prepData(paste0(path_to_data, "files"), 
                panel, 
                md, 
                md_cols = list(file = "file_name", id = "sample_id", factors = c("activation", "patient_id", "subgroup", "subgroup2")))


######## ----------------- Plot ----------------- ########

exprsData <- assays(sce)$exprs
exprsData <- data.frame(t(exprsData), colData(sce), check.names = FALSE)
gg_df <- reshape2::melt(exprsData, value.name = "expression", variable.name = "marker", id.vars = names(colData(sce)))
gg_df <- as.data.table(gg_df)
medians <- merge(gg_df[subgroup2 == "RP", median(expression), by=c("marker", "patient_id")], gg_df[subgroup2 == "MP", median(expression), by=c("marker", "patient_id")], by=c("marker", "patient_id"))
colnames(medians) <- c("marker", "patient_id", "RP_median", "MP_median")
medians[, foldchange := RP_median/MP_median]


ggplot(medians[!marker %in% c("CD40", "CD154", "CD3", "CD107a"),], aes(x = reorder(marker, foldchange), y = foldchange))+
  geom_boxplot()+
  geom_jitter(aes(shape=patient_id, color = patient_id))+
  scale_shape_manual(values=1:nlevels(medians$patient_id)) + labs(x = "Marker", y = "Fold Change", color = "Patient ID", shape = "Patient ID") 

######## ----------------- Normalize MP Expression Values (based on CD42b) ----------------- ########

exprs <- assays(sce)$exprs
col_dt <- as.data.frame(colData(sce))
patient_ids <- unique(col_dt$patient_id)

norm_exprs <- exprs
for (id in patient_ids){
  cols_MP <-  col_dt$patient_id == id & col_dt$subgroup2 == "MP"
  cols_RP <-  col_dt$patient_id == id & col_dt$subgroup2 == "RP"
  median_RP <- median(exprs[row.names(exprs) == "CD42b",][cols_RP], na.rm=TRUE)
  median_MP <- median(exprs[row.names(exprs) == "CD42b",][cols_MP], na.rm=TRUE)
  fold_change <- median_RP/median_MP
  norm_exprs[,cols_MP] <- exprs[,cols_MP] * fold_change
}

assays(sce)$exprs <- norm_exprs
saveRDS(sce, paste0(path_to_data, "sce_norm_CD42b.rds"))


######## ----------------- Write FCS Files ----------------- ########
norm_out_dir <- file.path(path_to_data, "files_normalized")
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


sce_norm_old <- readRDS("/Users/lisiarend/Desktop/sce_norm_CD42b.rds")
old_exprs <- assays(sce_norm_old)$exprs

sce_norm_new <- readRDS("/Users/lisiarend/Desktop/RPs_MPs_Normalization/Cyanus_RPsMPs/sce_norm_CD42b.rds")
new_exprs <- assays(sce)$exprs

