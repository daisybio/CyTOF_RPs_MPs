############# ------------- Prepare Data ------------- #############

# Load packages
suppressMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(CATALYST)
  source("functions.R")
})

path_to_data <- "/nfs/data/Bongiovanni-KrdIsar-platelets/Cyanus_RPsMPs/data/"

out_path <- "/nfs/data/Bongiovanni-KrdIsar-platelets/Cyanus_RPsMPs/data/sce_objects/"


######## Original Data ########

# add DNA channels in panel
panel <- fread(file.path(path_to_data, "panel.csv"))
panel <- rbind(panel, data.table(fcs_colname = c("Ir191Di", "Ir193Di"), antigen = c('DNA1', 'DNA2'), marker_class = c('none', 'none')))

# built meta data file
files <- list.files(file.path(path_to_data, "CCS_baseline"))
md <- data.table(file_name = files)
md$sample_id <- sapply(strsplit(md$file_name,"[.]"), "[", 1)
md$sample_id <- str_replace(md$sample_id, "_platlet_specific", "")
md$sample_id <- str_replace(md$sample_id, "RPS_", "RPs") # inconsistent naming of files
md$type <- sapply(strsplit(md$sample_id,"_"), "[", 2)
md$patient_id <- sapply(strsplit(md$sample_id,"_"), "[", 1) 

# TODO: once we have the stimulated files -> add column for activation

sce <- prepData(file.path(path_to_data, "CCS_baseline"), 
                panel, 
                md, 
                transform = TRUE,
                cofactor = 5,
                md_cols = list(file = "file_name", id = "sample_id", factors = c("type", "patient_id")))
saveRDS(sce, file.path(out_path, paste0("sce_original_RPs_MPs_rest.rds")))

######## Data Normalized by CD42b ########

normalization_marker <- "CD42b"
sce_CD42b <- normalize_patient_wise_sce(sce, column = "type", ain = "exprs", aout = "exprs", marker = normalization_marker)
saveRDS(sce_CD42b, file.path(out_path, paste0("sce_", normalization_marker, "_RPs_MPs_rest.rds")))

######## Data Normalized by DNA2 ########

normalization_marker <- "DNA2"
sce_DNA2 <- normalize_patient_wise_sce(sce, column = "type", ain = "exprs", aout = "exprs", marker = normalization_marker)
saveRDS(sce_DNA2, file.path(out_path, paste0("sce_", normalization_marker, "_RPs_MPs_rest.rds")))

