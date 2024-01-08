##################### Differential Expression Analysis of RPs/MPs Data #####################

######## ----------------- Read data ----------------- ########

sce_file <- "data/sce_original_RPs_MPs.rds"
type <- "original"

#sce_file <- "data/sce_CD42b_RPs_MPs.rds"
#type <- "CD42b"


sce <- readRDS(sce_file)

######## ----------------- Source Functions ----------------- ########

sapply(list.files("functions/", full.names = TRUE), source)

######## ----------------- Clustering ----------------- ########

if(is.null(sce$cluster_id)){
  sce <- clusterSCE(sce, features = "type", seed = 1234, maxK = 40)
  
  saveRDS(sce, sce_file)
}

######## ----------------- DE Analysis ----------------- ########

condition <- "subgroup2"
random_effect <- "patient_id"
sceEI <- CATALYST::ei(sce)

# baseline RPs vs. MPs.
message(".... baseline RPs vs. MPs ....")

sce_b <- sce[,sce$activation == "baseline"]
sceEI_b <- sceEI[sceEI$activation == "baseline",]
metadata(sce_b)$experiment_info <- sceEI_b

de_res_b <- runDS(sce_b,
                clustering_to_use = "all",
                contrast_vars = condition,
                markers_to_test = c("state", "type"),
                ds_methods = c("CyEMD", "t_test"),
                design_matrix_vars = c(random_effect, condition),
                fixed_effects = condition,
                random_effects = random_effect,
                cyEMD_nperm = 500)

# lets compute and save the effect size:
de_res_b$effect_size <- effectSize(sce_b, condition, random_effect, k = "all")
saveRDS(de_res_b, paste0("data/de_results_baseline_", type, ".rds"))


# stimulated RPs vs. MPs
message(".... stimulated RPs vs. MPs ....")

sce_a <- sce[,sce$activation == "stimulated"]
sceEI_a <- sceEI[sceEI$activation == "stimulated",]
metadata(sce_a)$experiment_info <- sceEI_a

de_res_a <- runDS(sce_a,
                  clustering_to_use = "all",
                  contrast_vars = condition,
                  markers_to_test = c("state", "type"),
                  ds_methods = c("CyEMD", "t_test"),
                  design_matrix_vars = c(random_effect, condition),
                  fixed_effects = condition,
                  random_effects = random_effect,
                  cyEMD_nperm = 500)

de_res_a$effect_size <- effectSize(sce_a, condition, random_effect, k = "all")
saveRDS(de_res_a, paste0("data/de_results_stimulated_", type, ".rds"))

