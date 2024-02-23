############# ------------- Overview Figure CyTOF ------------- #############

# Load packages
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
  source("functions/de_functions.R")
  source("functions.R")
})

path_to_data <- "/nfs/data/Bongiovanni-KrdIsar-platelets/Cyanus_RPsMPs/data/sce_objects"

analysis_state <- "baseline" # or "stimulated"


# For reproducibility
set.seed(1234)

# Fix colors
colors <- c("#FF1F5B", "#009ADE", "#C4C4C4")
names(colors) <- c("RP", "MP", "rest")


# Read data
mapping <- c(paste0("sce_", analysis_state, "_original_RPs_MPs_rest.rds"), 
             paste0("sce_", analysis_state, "_CD42b_RPs_MPs_rest.rds"), 
             paste0("sce_", analysis_state, "_DNA2_RPs_MPs_rest.rds"))
names(mapping) <- c('Original', 'CD42b', 'DNA2')

sce <- readRDS(file.path(path_to_data, mapping['Original']))
sce_CD42b <- readRDS(file.path(path_to_data, mapping['CD42b'])) 
sce_DNA2 <- readRDS(file.path(path_to_data, mapping['DNA2']))


######## Dimensionality Reduction ########

# tSNE colored by type (RPs, MPs, rest)
sce <- runDR(sce,
             dr = c("TSNE"),
             cells = 1000,
             features = "type",
             assay = "exprs")

tsne_RPs_MPs_rest_plot <- plotDR(sce,
                    dr = "TSNE",
                    color_by = "type") + 
  scale_color_manual(name = "", values = colors) + 
  theme(legend.position = "bottom") +
  guides(color = guide_legend(direction = "horizontal", title = "", override.aes = list(size = 5)))

ggsave(paste0("plots/TSNE_", analysis_state, "_RPs_MPs_rest.png"), width = 4, height = 4, dpi = 300)

# UMAP colored by type (RPs, MPs, rest)
sce <- runDR(sce,
             dr = c("UMAP"),
             cells = 1000,
             features = "type",
             assay = "exprs")

umap_RPs_MPs_rest_plot <- plotDR(sce,
                                 dr = "UMAP",
                                 color_by = "type") + 
  scale_color_manual(name = "", values = colors) + 
  theme(legend.position = "bottom") +
  guides(color = guide_legend(direction = "horizontal", title = "", override.aes = list(size = 5)))

ggsave(paste0("plots/UMAP_", analysis_state, "_RPs_MPs_rest.png"), width = 4, height = 4, dpi = 300)


# tSNE colored by expression, separated by RP and MP

sce_RPs_MPs <- sce[,sce$type != "rest"]

markers1 <- c("CD62P", "CD63", "GPVI", "PAR1", "CD40") 
markers2 <- c("CD42a", "CD42b", "PEAR", "CD31", "PAC1")

tsne_expression_plot_1 <- plotDR(sce_RPs_MPs,
                               dr = c("TSNE"),
                               color_by = markers1,
                               facet_by = "type",
                               ncol = 4) + 
  theme(legend.position = "bottom", legend.title = element_text(hjust = 1)) +
  guides(color = guide_colorbar(title = "Scaled Expression", direction = "horizontal", title.position = "left", title.vjust = 0.8))

tsne_expression_plot_2 <- plotDR(sce_RPs_MPs,
                                 dr = "TSNE",
                                 color_by = markers2,
                                 facet_by = "type",
                                 ncol = 4) + 
  theme(legend.position = "bottom", legend.title = element_text(hjust = 1)) +
  guides(color = guide_colorbar(title = "Scaled Expression", direction = "horizontal", title.position = "left", title.vjust = 0.8))

tsne_expression_plots <- ggarrange(tsne_expression_plot_1, tsne_expression_plot_2, ncol = 1, common.legend = TRUE, labels = NULL, legend = "bottom")

ggsave(paste0("plots/TSNE_", analysis_state, "_marker_expression.png"), width = 12, height = 8, dpi = 300)

# UMAP colored by expression, separated by RP and MP
umap_expression_plot_1 <- plotDR(sce_RPs_MPs,
                                 dr = c("UMAP"),
                                 color_by = markers1,
                                 facet_by = "type",
                                 ncol = 4) + 
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(title = "Scaled Expression", direction = "horizontal", title.position = "left", title.vjust = 0.8))

umap_expression_plot_2 <- plotDR(sce_RPs_MPs,
                                 dr = "UMAP",
                                 color_by = markers2,
                                 facet_by = "type",
                                 ncol = 4) + 
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(title = "Scaled Expression", direction = "horizontal", title.position = "left", title.vjust = 0.8))

umap_expression_plots <- ggarrange(umap_expression_plot_1, umap_expression_plot_2, ncol = 1, common.legend = TRUE, labels = NULL, legend = "bottom")

ggsave(paste0("plots/UMAP_", analysis_state, "_marker_expression.png"), width = 12, height = 8, dpi = 300)


# Paired (patient-wise) analysis of marker expressions

# On original data
df_medians_original <- paired_boxes(sce, 'Original', paste0("tables/median_table_with_paired_results_", analysis_state, "_original.csv"))

df <- df_medians_original[!df_medians_original$marker %in% c("CD45", "DNA1", "DNA2", "CD47"),]
violins_original <- ggplot(df[signif != ""], aes(x = group, y = Expression, color = group, fill = group))+
  geom_violin(alpha = 0.3)+
  geom_point()+
  geom_line(aes(group = patient_id), color = '#C4C4C4')+
  scale_color_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  scale_fill_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  facet_wrap(~marker_title, scales = 'free', ncol = 8)+
  theme_minimal()+
  theme(legend.position = 'none', axis.title.x=element_blank(), axis.text.x = element_blank(), strip.text.x = element_text(face = "bold"))
ggsave(paste0("plots/paired_boxes_", analysis_state, "_original.png"), width = 12, height = 4, dpi = 300)



# Normalized by size
df_medians_CD42b <- paired_boxes(sce_CD42b, 'CD42b', paste0("tables/median_table_with_paired_results_", analysis_state, "_CD42b.csv"))

df <- df_medians_CD42b[!df_medians_CD42b$marker %in% c("CD45", "DNA1", "DNA2"),]
violins_CD42b <- ggplot(df[signif != ""], aes(x = group, y = Expression, color = group, fill = group))+
  geom_violin(alpha = 0.3)+
  geom_point()+
  geom_line(aes(group = patient_id), color = '#C4C4C4')+
  scale_color_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  scale_fill_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  facet_wrap(~marker_title, scales = 'free', ncol = 4)+
  theme_minimal()+
  theme(legend.position = 'none', axis.title.x=element_blank(), axis.text.x = element_blank(), strip.text.x = element_text(face = "bold"))
ggsave(paste0("plots/paired_boxes_", analysis_state, "_CD42b.png"), width = 6, height = 5, dpi = 300)


# Normalized by RNA
df_medians_DNA2 <- paired_boxes(sce_DNA2, 'DNA2', paste0("tables/median_table_with_paired_results_", analysis_state, "_DNA2.csv"))

df <- df_medians_DNA2[!df_medians_DNA2$marker %in% c("CD45"),]
violins_DNA2 <- ggplot(df[signif != ""], aes(x = group, y = Expression, color = group, fill = group))+
  geom_violin(alpha = 0.3)+
  geom_point()+
  geom_line(aes(group = patient_id), color = '#C4C4C4')+
  scale_color_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  scale_fill_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  facet_wrap(~marker_title, scales = 'free', ncol = 4)+
  theme_minimal()+
  theme(legend.position = 'none', axis.title.x=element_blank(), axis.text.x = element_blank(), strip.text.x = element_text(face = "bold"))
ggsave(paste0("plots/paired_boxes_", analysis_state, "_DNA2.png"), width = 6, height = 5, dpi = 300)


# Overall Figure
#plot_A_B <- ggarrange(tsne_RPs_MPs_rest_plot, tsne_expression_plots, ncol = 2, widths = c(0.5, 1))
#plot_C <- ggarrange(violins_original, legend = NULL)
#plot_D_E <- ggarrange(violins_CD42b, violins_DNA2, legend = NULL, labels = NULL)
#ggarrange(NULL, plot_A_B, NULL, plot_C, NULL, plot_D_E, nrow = 6, heights = c(0.05, 0.6, 0.05, 0.5, 0.05, 0.7), labels = NULL)
