
##################### Clustering state markers of RPs/MPs/Rest Data #####################

######## ----------------- Source functions ----------------- ########
library(data.table)
library(CATALYST)
library(stringr)
library(RColorBrewer)
source("functions.R")

normalize <- TRUE
normalization_marker <- 'DNA2'
add_DNA <- TRUE

######## ----------------- Prepare Data ----------------- ########
path_to_data <- "/nfs/data/Bongiovanni-KrdIsar-platelets/Cyanus_RPsMPs/data/"

if(add_DNA){
  outfile <- paste0("sce_", normalization_marker, "_RPs_MPs_rest.rds")
}else{
  outfile <- paste0("sce_", normalization_marker, "_RPs_MPs_rest_DNA.rds")
}

if(file.exists(file.path(path_to_data, outfile))){
  sce <- readRDS(file.path(path_to_data, outfile))
} else {
  panel <- fread(file.path(path_to_data, "panel.csv"))
  if(add_DNA){
    panel <- rbind(panel, data.table(fcs_colname = c("Ir191Di", "Ir193Di"), antigen = c('DNA1', 'DNA2'), marker_class = c('none', 'none')))
  }
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
  
  #saveRDS(sce, file.path(path_to_data, outfile))
}

######## ----------------- Should we really normalize MPs by RPs?  ----------------- ########
if(normalize){
  sce <- normalize_patient_wise_sce(sce, column = "type", ain = "exprs", aout = "exprs", marker = normalization_marker)
  saveRDS(sce, file.path(path_to_data, paste0("sce_", normalization_marker, "_RPs_MPs_rest", ifelse(add_DNA, "_DNA", ""), ".rds")))
}

######## ----------------- Add Cell Annotation  ----------------- ########

annotation <- colData(sce)$type

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

