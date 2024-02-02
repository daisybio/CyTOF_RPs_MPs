##################### Differential Expression Analysis of RPs/MPs Data #####################

# run the script like this for example:
# Rscript run_DE.R /nfs/data/Bongiovanni-KrdIsar-platelets/Cyanus_RPsMPs/data/sce_original_RPs_MPs.rds original

# first argument -> path of the sce object
# second argument -> if it is the "original" sce object or the "CD42b" object

#source("renv/activate.R")

######## ----------------- Source Functions ----------------- ########

sapply(list.files("functions/", full.names = TRUE), source)


######## ----------------- Save Arguments ----------------- ########
args <- commandArgs(TRUE)
sceFile <- args[1]
type <- args[2]

outputPath <- strsplit(sceFile, "/")[[1]]
outputPath <- outputPath[1:(length(outputPath)-1)]
outputPath <- paste(outputPath, collapse = "/")


######## ----------------- Read data ----------------- ########

if(file.exists(scePath)){
  sce <- readRDS(sceFile)
} else {
  stop("There is no file or directory with the given name!")
}  


######## ----------------- Clustering ----------------- ########

if(is.null(sce$cluster_id)){
  sce <- clusterSCE(sce, features = "type", seed = 1234, maxK = 40)
}

######## ----------------- DE Analysis ----------------- ########

condition <- "subgroup2"
random_effect <- "patient_id"
sceEI <- CATALYST::ei(sce)

# baseline RPs vs. MPs.
message(".... baseline RPs vs. MPs ....")

sce_b <- sce[,sce$activation == "baseline"]
sceEI_b <- sceEI[sceEI$activation == "baseline",]
S4Vectors::metadata(sce_b)$experiment_info <- sceEI_b

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
saveRDS(de_res_b, file.path(outputPath, paste("de_results_baseline_", type, ".rds")))


# stimulated RPs vs. MPs
message(".... stimulated RPs vs. MPs ....")

sce_a <- sce[,sce$activation == "stimulated"]
sceEI_a <- sceEI[sceEI$activation == "stimulated",]
S4Vectors::metadata(sce_a)$experiment_info <- sceEI_a

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
saveRDS(de_res_b, file.path(outputPath, paste("de_results_stimulated_", type, ".rds")))

