##################### Normalization of RPs/MPs Data #####################
library(data.table)
library(flowCore)
path_to_data <- "data"

######## ----------------- Source Functions ----------------- ########

source("functions.R")

######## ----------------- Read files ----------------- ########

panel <- fread(file.path(path_to_data, "panel.csv"))
md <- fread(file.path(path_to_data, "metafile RP MP all.csv"))
sce <- prepData(file.path(path_to_data, "RPs_MPs_files"), 
                panel, 
                md, 
                md_cols = list(file = "file_name", id = "sample_id", factors = c("activation", "patient_id", "subgroup", "subgroup2")))

saveRDS(sce, file.path(path_to_data, "sce_original_RPs_MPs.rds"))


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
  file_name <- md[md$sample_id == id,]$file_name
  write.FCS(ff, paste0(norm_out_dir,"/",file_name))
}

######## ----------------- Check Normalized FCS Files ----------------- ########

sce <- prepData(file.path(path_to_data, "RPs_MPs_files_normalized"), 
                panel, 
                md, 
                md_cols = list(file = "file_name", id = "sample_id", factors = c("activation", "patient_id", "subgroup", "subgroup2")))

medians <- get_patient_FCs(sce, ain = "exprs")
medians[marker == "CD42b",]

