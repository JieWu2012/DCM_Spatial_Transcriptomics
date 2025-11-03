################################################

R scritps for analyzing single cell RNA-seq data

                  By Jie Wu

##############################################




#load r package
#some of the package needs to be installed in the terminal

library(Seurat)
library(ggplot2)
library(purrr)
library(dplyr)
library(data.table)
library(R.utils)
library(viridis)
library(fields)
library(dplyr)
library(EnhancedVolcano)


# set the working directory

setwd("/oak/stanford/groups/larsms/Users/jwu/SpatialTranscriptomics/data1/")

# read RDS
all <- readRDS("CustomVariables_results_all_default/all_default.RDS")


# sub-sampling cells to make UMAP

set.seed(1234)
sel = sample(ncol(all), 5e4)
seu = all[,sel]

DimPlot(seu, reduction = "umap", label = FALSE, repel=TRUE,pt.size = 0.1,raster=FALSE, group.by=c("celltype_new"),
        cols=colors, label.size = 5) + coord_fixed(ratio = 1)
ggsave(file="CustomVariables_results_all_default/all_umap_Celltype_anno_sampling5e4_new.pdf", dpi = "print", height = 7, width = 11)

DimPlot(seu, reduction = "umap", label = FALSE, repel=TRUE,pt.size = 0.1,raster=FALSE, split.by = c("DiseaseStatus_edited"),group.by=c("celltype_new"),
    cols=colors, label.size = 5) + coord_fixed(ratio = 1) + theme(legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 5), nrow = 2, byrow = TRUE))
ggsave(file="CustomVariables_results_all_default/all_umap_Celltype_anno_byDisease_sampling5e4_new.pdf", dpi = "print", height = 7, width = 15)








color1 <- viridis(10, option = "D")
FeaturePlot(
  seu,
  reduction = "umap", # pick tsne or umap
  features = c("F830016B08Rik"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90",
  raster=FALSE,
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)
ggsave(file=paste0("CustomVariables_results_all_default/seu_umap_", "FRik", ".pdf"), dpi = "print", height = 6, width = 14)


color1 <- viridis(10, option = "D")
FeaturePlot(
  seu,
  reduction = "umap", # pick tsne or umap
  features = c("Ttn-N2A"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  min.cutoff = "q1",
  max.cutoff = "q99",
  raster=FALSE, 
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)
ggsave(file=paste0("CustomVariables_results_all_default/seu_umap_", "Ttn_N2A", ".pdf"), dpi = "print", height = 6, width = 14)


FeaturePlot(
  seu,
  reduction = "umap", # pick tsne or umap
  features = c("Ttn-N2B"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90",
  raster=FALSE, 
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)
ggsave(file=paste0("CustomVariables_results_all_default/seu_umap_", "Ttn_N2B", ".pdf"), dpi = "print", height = 6, width = 14)

FeaturePlot(
  seu,
  reduction = "umap", # pick tsne or umap
  features = c("Nppa"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90",
  raster=FALSE, 
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)
ggsave(file=paste0("CustomVariables_results_all_default/seu_umap_", "Nppa", ".pdf"), dpi = "print", height = 6, width = 14)


FeaturePlot(
  seu,
  reduction = "umap", # pick tsne or umap
  features = c("Myh7"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90",
  raster=FALSE, 
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)
ggsave(file=paste0("CustomVariables_results_all_default/seu_umap_", "Myh7", ".pdf"), dpi = "print", height = 6, width = 14)


FeaturePlot(
  seu,
  reduction = "umap", # pick tsne or umap
  features = c("Nppb"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90",
  raster=FALSE, 
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)
ggsave(file=paste0("CustomVariables_results_all_default/seu_umap_", "Nppb", ".pdf"), dpi = "print", height = 6, width = 14)


