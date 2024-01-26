suppressMessages({
  library(data.table)
  library(ggplot2)
  library(SingleCellExperiment)
  library(CATALYST)
  library(ggrepel)
  library(ggpubr)
  library(RColorBrewer)
  library(ggplotify)
  library(pheatmap)
  library(patchwork)
  source("functions/prep_functions.R")
})

theme_set(theme_bw(base_size = 16))

sce <- readRDS("/nfs/data/Bongiovanni-KrdIsar-platelets/Cyanus_RPsMPs/data/sce_original_RPs_MPs_rest_DNA.rds")
sce_norm <- readRDS("/nfs/data/Bongiovanni-KrdIsar-platelets/Cyanus_RPsMPs/data/sce_CD42b_RPs_MPs_rest_DNA.rds")

paired_boxes <- function(sce) {
  df <- as.data.table(t(assays(sce)$exprs))
  df[, group := colData(sce)$type]
  df[, patient_id := colData(sce)$patient_id]
  
  df <- melt(df, id.vars = c('group', 'patient_id'), variable.name = 'marker', value.name = 'expression')
  
  df_medians <- df[!group == 'rest', median(expression), by=c("patient_id", "group", "marker")]
  colnames(df_medians) <- c('patient_id', 'group', 'marker', 'expression')
  df_medians[, patient_id := as.factor(patient_id)]
  
  ggpaired(df_medians, x = "group", y = "expression",
           color = "group", line.color = "gray", line.size = 0.4, id="patient_id",
           palette = "jco")+
    #stat_compare_means(paired = TRUE, method = "t.test")+
    facet_wrap(~marker, scales = 'free')
}

paired_boxes(sce)
ggsave('plots/paired_boxes.png', height=8, width=12)
paired_boxes(sce_norm)
ggsave('plots/paired_boxes_norm.png', height=8, width=12)


dna_boxes <- function(sce) {
  df <- as.data.table(t(assays(sce)$exprs))
  df[, patient_id := colData(sce)$patient_id]
  
  df_melt <- melt(df, id.vars = c('patient_id', 'DNA1', 'DNA2'), variable.name = 'marker', value.name = 'expression')
  DNA_1 <- df_melt[, cor(x=DNA1, y=expression), by=c('patient_id', 'marker')]
  DNA_1[, median_cor := median(V1), by = marker]
  colnames(DNA_1) <- c('patient_id', 'marker', 'correlation', 'median_cor')
  
  DNA_2 <- df_melt[, cor(x=DNA2, y=expression), by=c('patient_id', 'marker')]
  DNA_2[, median_cor := median(V1), by = marker]
  DNA_2 <- DNA_2[order(median_cor)]
  DNA_2[, marker := factor(marker, levels = unique(DNA_2$marker))]
  DNA_1[, marker := factor(marker, levels = unique(DNA_2$marker))]
  colnames(DNA_2) <- c('patient_id', 'marker', 'correlation', 'median_cor')
  
  dna_df <- merge(DNA_1[, -c('median_cor')], DNA_2[, -c('median_cor')], by = c('patient_id', 'marker'))
  colnames(dna_df) <- c('patient_id', 'marker', 'DNA_1', 'DNA_2')
  dna_df <- melt(dna_df, id.vars = c('patient_id', 'marker'), variable.name = 'DNA_marker', value.name = 'correlation')
  
  ggplot(dna_df, aes(x = marker, y = correlation, fill = DNA_marker))+
    geom_boxplot()
}

dna_boxes(sce)
ggsave('plots/dna_boxes.png', height = 5, width=17)
dna_boxes(sce_norm)
ggsave('plots/dna_boxes_norm.png', height = 5, width=17)
