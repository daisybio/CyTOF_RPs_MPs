
##################### Clustering state markers of RPs/MPs/Rest Data #####################

normalize <- TRUE

######## ----------------- Source functions and install packages ----------------- ########

source("functions.R")

library(data.table)
library(SummarizedExperiment)
library(RColorBrewer)
library(ggplot2)

######## ----------------- Prepare Data ----------------- ########
path_to_data <- "data/"


if(file.exists(file.path(path_to_data, "sce_original_RPs_MPs_rest.rds"))){
  sce <- readRDS(file.path(path_to_data, "sce_original_RPs_MPs_rest.rds"))
} else {
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
}

######## ----------------- Should we really normalize MPs by RPs?  ----------------- ########

if(normalize){
  sce <- normalize_patient_wise_sce(sce, column = "type", ain = "exprs", aout = "exprs")
}


######## ----------------- Perform Clustering  ----------------- ########

sce <- clusterSCENew(sce, features = "state", seed = 1234, maxK = 40)

CATALYST::delta_area(sce)

if(normalize){
  file_name <- paste0("sce_CD42b_clustered_RPs_MPs_rest.rds")
} else {
  file_name <- paste0("sce_original_clustered_RPs_MPs_rest.rds")
}

saveRDS(sce, file.path(path_to_data, file_name))


######## ----------------- Compare Clustering  ----------------- ########

gated_annotation <- colData(sce)$type

cluster_annotation <- 

rand_index <- rand.index(cluster_annotation, gated_annotation)
adj_rand_index <- adj.rand.index(cluster_annotation, gated_annotation)

