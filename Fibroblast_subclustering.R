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



# extract fibroblasts
all_Fibro <- subset(all, subset = celltype=="Fibroblasts")


# same scRNA-seq workflow for only fibroblasts
panel_gene <- read.table("CustomVariables_results/panel_gene_new.txt", sep="\t")
panel_gene <- panel_gene[grepl("ENSMUSG", panel_gene$V1),]
all_Fibro <- FindVariableFeatures(all_Fibro,nfeatures = 200)
gene <- VariableFeatures(all_Fibro)
valgene <- gene[gene %in% panel_gene$V2]
VariableFeatures(all_Fibro)<-valgene
all_Fibro <- RunPCA(all_Fibro,approx = FALSE)
ElbowPlot(all_Fibro,ndims = 30)
all_Fibro <- FindNeighbors(all_Fibro, dims = 1:20)
all_Fibro <- FindClusters(all_Fibro , resolution = 0.1)
all_Fibro <- RunUMAP(all_Fibro, dims = 1:20, n.neighbors = 100, min.dist = 0.05, spread=2)


DimPlot(all_Fibro, reduction = "umap", label = FALSE, repel=TRUE,pt.size = 0.1, cols=colors, raster=FALSE)+coord_fixed(ratio = 1)
ggsave(file="CustomVariables_results_all_default/Fibroblast_umap_subclusters.pdf", dpi = "print", height = 7, width = 11)

DimPlot(all_Fibro, reduction = "umap",raster=FALSE, split.by = c("DiseaseStatus_edited"),cols = colors)+coord_fixed(ratio = 1)
ggsave(file="CustomVariables_results_all_default/Fibroblast_umap_subclusters_byDisease.pdf", dpi = "print", height = 7, width = 15)

all_markers_only_positive <-FindAllMarkers(all_Fibro,only.pos = TRUE, features = valgene)
gen_marker_table <- function(x){
    all_markers_only_positive[all_markers_only_positive$cluster == x, ] %>%
    head(n=5)
}
all_top5_markers_positiv <- map_dfr(unique(all_markers_only_positive$cluster), gen_marker_table)
write.table(all_top5_markers_positiv, "CustomVariables_results_all_default/all_fibroblasts_top5_markers_positiv.txt", sep = "\t",quote = FALSE)

Fibro_subcluster <- c(
  "0" = "Fibroblasts", 
  "1" = "Fibroblasts",
  "2" = "Myofibroblasts", 
  "3" = "Macrophages-like Fibroblasts"         
)


all_Fibro@meta.data$fibroblast_subcluster = Fibro_subcluster[as.character(all_Fibro$seurat_clusters)]

color_fibro <- c(
  "#aa6e28", # Brown
  "#D73027",  # Bold Red
  "#008080" # Teal
)

DimPlot(all_Fibro, reduction = "umap",raster=FALSE, split.by = c("DiseaseStatus_edited"),cols = color_fibro, group.by = c("fibroblast_subcluster"))+coord_fixed(ratio = 1)
ggsave(file="CustomVariables_results_all_default/Fibroblast_umap_subclusters_byDisease_anno.pdf", dpi = "print", height = 7, width = 15)


celltypes <- all_Fibro$fibroblast_subcluster
names(celltypes) <- colnames(all_Fibro)
all$celltype_new <- all$celltype
all$celltype_new[names(celltypes)] <- celltypes

df <- as.data.frame(all@meta.data[all@meta.data$orig.ident=="BE_rep3",])

ggplot(df, aes(x=x_centroid, y=y_centroid))+geom_point(size=0.1, color="gray") + 
  geom_point(data = df[df$celltype_new=="Stressed VCMs",], aes(x=x_centroid, y=y_centroid, color=celltype_new), size=0.2) + 
  scale_color_manual(values = "#911eb4") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "top") +  
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  coord_fixed(ratio = 1) + 
  scale_y_reverse()+
  ggtitle("BE_rep3")+
  ggplot(df, aes(x=x_centroid, y=y_centroid)) +
  geom_point(size=0.1, color="gray") + 
  geom_point(data = df[df$celltype_new=="Myofibroblasts",], aes(x=x_centroid, y=y_centroid, color=celltype_new), size=0.2) + 
  scale_color_manual(values = "#008080") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  theme(legend.position = "top") +  
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  coord_fixed(ratio = 1) + 
  scale_y_reverse()+
  ggplot(df, aes(x=x_centroid, y=y_centroid)) +
  geom_point(size=0.1, color="gray") + 
  geom_point(data = df[df$celltype_new=="Fibroblasts",], aes(x=x_centroid, y=y_centroid, color=celltype_new), size=0.2) + 
  scale_color_manual(values = "#aa6e28") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  theme(legend.position = "top") +  
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  coord_fixed(ratio = 1) + 
  scale_y_reverse()+
  ggplot(df, aes(x=x_centroid, y=y_centroid)) +
  geom_point(size=0.1, color="gray") + 
  geom_point(data = df[df$celltype_new=="Macrophages-like Fibroblasts",], aes(x=x_centroid, y=y_centroid, color=celltype_new), size=0.2) + 
  scale_color_manual(values = "#D73027") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  theme(legend.position = "top") +  
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  coord_fixed(ratio = 1) + 
  scale_y_reverse()

ggsave(file="CustomVariables_results_all_default/Fibroblast_stressedCM_BE_rep3.pdf", dpi = "print", height = 12, width = 10)

saveRDS(object = all, file="CustomVariables_results_all_default/all_default_fibrosubcluster.RDS")





# fibroblasts subclustering with high resolution

all <- readRDS("CustomVariables_results_all_default/all_default.RDS")
panel_gene <- read.table("CustomVariables_results/panel_gene_new.txt", sep="\t")
all_Fibro <- subset(all, subset = celltype=="Fibroblasts")
panel_gene <- panel_gene[grepl("ENSMUSG", panel_gene$V1),]
all_Fibro <- FindVariableFeatures(all_Fibro,nfeatures = 200)
gene <- VariableFeatures(all_Fibro)
valgene <- gene[gene %in% panel_gene$V2]
VariableFeatures(all_Fibro)<-valgene
all_Fibro <- RunPCA(all_Fibro,approx = FALSE)
ElbowPlot(all_Fibro,ndims = 30)
all_Fibro <- FindNeighbors(all_Fibro, dims = 1:20)
all_Fibro <- FindClusters(all_Fibro , resolution = 0.5)

DimPlot(all_Fibro, reduction = "umap", label = FALSE, repel=TRUE,pt.size = 0.1, raster=FALSE)+coord_fixed(ratio = 1)
ggsave(file="CustomVariables_results_all_default/Fibroblast_umap_subclusters_resol0.5.pdf", dpi = "print", height = 7, width = 11)


colors <- c(
"#1F77B4", # blue
"#AED6F1", # light blue
"#FF7F0E", # orange
"#FFBB78", # light orange
"#2CA02C", # green
"#98DF8A", # light green
"#D62728", # red
"#FF9896", # pink
"#9467BD", # purple
"#8C564B", # brown
"#C49C94", # light brown
"#E377C2", # magenta
"#F7B6D2", # light magenta
"#7F7F7F", # gray
"#C7C7C7", # light gray
"#BCBD22", # olive
"#DBDB8D", # light olive
"#17BECF"  # teal
)


DimPlot(all_Fibro, reduction = "umap", label = FALSE, repel=TRUE,pt.size = 0.1, raster=FALSE, cols = colors)+coord_fixed(ratio = 1)
ggsave(file="CustomVariables_results_all_default/Fibroblast_umap_subclusters_resol0.5.pdf", dpi = "print", height = 7, width = 11)

DimPlot(all_Fibro, reduction = "umap", label = FALSE, repel=TRUE,pt.size = 0.1, raster=FALSE, cols = colors, split.by = "DiseaseStatus_edited")+coord_fixed(ratio = 1)
ggsave(file="CustomVariables_results_all_default/Fibroblast_umap_subclusters_resol0.5_splitbyDisease.pdf", dpi = "print", height = 7, width = 15)


all_markers_only_positive <-FindAllMarkers(all_Fibro,only.pos = TRUE, features = valgene)

gen_marker_table <- function(x){
    all_markers_only_positive[all_markers_only_positive$cluster == x, ] %>%
    head(n=10)
}

all_top10_markers_positiv <- map_dfr(unique(all_markers_only_positive$cluster), gen_marker_table)
write.table(all_top10_markers_positiv, "CustomVariables_results_all_default/all_fibroblasts_deepclustering_top10_markers_positiv.txt", sep = "\t",quote = FALSE)

##cell annotation explanation: 

##________
#| Subcluster | Top Markers                                                                    | Likely Identity / Notes                                                                                                                                                   |
#  | ---------- | ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
#  | 0          | Pln, Atp2a2, Myh6, Myl2, mt-Nd5, mt-Atp6, Tpm1, Fhl2, Pdgfra, Hopx             | Likely **cardiomyocyte contamination** or **myofibroblast with contractile features**, since Pln, Atp2a2, Myh6, Myl2 are cardiac muscle genes.                            |
#  | 1          | Pdgfra1, Egfr1, Fgfr1, Gas1, Csf1, Col1a1, Tgfbr2, Col6a2, Osmr, Il6st         | Classic **fibroblasts / mesenchymal cells** (Pdgfra, Col1a1, Col6a2, Tgfbr2). Possibly **proliferative or ECM-producing fibroblasts**.                                    |
#  | 2          | Pln1, Hopx1, Vegfb1, Fgf1, Ppargc1a1, Fhl21, Ttn, Mybpc3, S100a11, Actn2       | Could be **fibroblast subset with some cardiomyocyte-like features**; maybe **activated fibroblasts or pericyte-like**.                                                   |
#  | 3          | Nppa, Ankrd1, Ttn1, Flnc, Actn21, Dsp1, Mybpc31, Obscn1, Des1, S100a12         | Markers like Nppa, Ankrd1 indicate **stress-responsive fibroblasts or early cardiomyocyte contamination**. Could represent **fibroblasts responding to injury / stress**. |
#  | 4          | Prg4, Dkk3, Npr11, Postn, Tbx201, Tgfbr11, Runx1, Sox91, Itga10, Sdc41         | **Synovial-like fibroblasts / ECM-modulating fibroblasts**. Prg4 is a lubricin gene; Dkk3 is Wnt-modulator.                                                               |
#  | 5          | Postn2, Meox13, Prg41, Nppa1, Tgfb12, Fn12, Ankrd11, Col1a12, Itga52, Tgfbr12  | **Myofibroblast / activated fibroblasts** (Postn, Tgfb1, Fn1, Col1a1).                                                                                                    |
 # | 6          | Egfl7, Cdh5, Pecam11, Adgrg1, Pdgfd1, Rgs51, Epas11, Itga12, Npr12, Sema6a1    | **Endothelial contamination** (Cdh5, Pecam1, Egfl7) and possibly pericytes.                                                                                               |
#  | 7          | Postn2, Sdc42, Chl1, Itga13, Actb4, Lmna3, Il6st2, Fos3, Mfn13, Erbb23         | Likely **activated fibroblasts / signaling fibroblasts**, some stress-response genes (Fos).                                                                               |
#  | 8          | F13a1, Cd163, Mrc1, Adgre1, Ms4a71, Cd681, Ptprc1, Cd742, Cd861, Actb5         | **Macrophages / immune cells** (Cd163, Mrc1, Adgre1).                                                                                                                     |
#  | 9          | Chad, Prg42, Dkk31, Itga102, Fn13, Runx13, Cacna1g4, Col1a13, Itga83, Col1a24  | Likely **perivascular fibroblasts / ECM-modulating fibroblasts**.                                                                                                         |
#  | 10         | Ccl2, Il61, Junb3, Itga53, Sdc45, Csf13, Actb6, Meox14, Sod24, Osmr3           | **Inflammatory fibroblasts** (Ccl2, Il6) and stress-response cells.                                                                                                       |
#  | 11         | mt-Atp63, mt-Nd53, Myl24, Ppargc1a2, Pln4, Tpm13, Myh63, Fhl23, Hopx2, Atp2a23 | **Cardiomyocyte contamination / mitochondrial-rich fibroblasts**.                                                                                                         |
#  | 12         | mt-Atp64, mt-Nd54, Myl25, Ppargc1a3, Ttn3, Fgf13, Obscn4, Pln5, Vegfb4, Fhl24  | Similar to **subcluster 11**: cardiomyocyte-like or mitochondrial-rich fibroblasts.                                                                                       |
#  | 13         | Wt11, Gja14, Npr14, Myl73, Sdc46, Csf14, Gas13, Wnt5a4, Il6st5, Calr6          | Likely **epicardial fibroblasts / developmental fibroblasts** (Wt1, Gja1, Wnt5a).                                                                                         |

##_________
cluster_annotations <- c(
  "0"  = "Cardiomyocyte-like fibroblast",
  "1"  = "Fibroblast",
  "2"  = "Activated fibroblast / pericyte-like",
  "3"  = "Stress-responsive fibroblast",
  "4"  = "Synovial-like / ECM-modulating fibroblast",
  "5"  = "Myofibroblast",
  "6"  = "Endothelial",
  "7"  = "Activated fibroblast",
  "8"  = "Macrophage",
  "9"  = "Perivascular fibroblast",
  "10" = "Inflammatory fibroblast",
  "11" = "Cardiomyocyte-like / mitochondrial-rich fibroblast",
  "12" = "Cardiomyocyte-like / mitochondrial-rich fibroblast",
  "13" = "Epicardial fibroblast"
)


cluster_annotations <- c(
  "0"  = "Cardiomyocyte-like fibroblast",
  "1"  = "Fibroblast",
  "2"  = "Cardiomyocyte-like fibroblast",
  "3"  = "Stress-responsive fibroblast",
  "4"  = "Synovial-like / ECM-modulating fibroblast",
  "5"  = "Myofibroblast",
  "6"  = "Endothelial/pericyte contamination",
  "7"  = "Activated / signaling fibroblast",
  "8"  = "Macrophage-like fibroblast",
  "9"  = "Perivascular fibroblast",
  "10" = "Inflammatory fibroblast",
  "11" = "Cardiomyocyte-like / mitochondrial-rich fibroblast",
  "12" = "Cardiomyocyte-like / mitochondrial-rich fibroblast",
  "13" = "Epicardial fibroblast"
)

all_Fibro@meta.data$fibroblast_subcluster = cluster_annotations[as.character(all_Fibro$seurat_clusters)]

DimPlot(all_Fibro, reduction = "umap", label = FALSE, repel=TRUE,pt.size = 0.1, raster=FALSE, cols = colors, split.by = "DiseaseStatus_edited", group.by = c("fibroblast_subcluster"))+coord_fixed(ratio = 1) + theme(legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 5), nrow = 2, byrow = TRUE))
ggsave(file="CustomVariables_results_all_default/Fibroblast_umap_subclusters_resol0.5_splitbyDisease_celltype.pdf", dpi = "print", height =9, width = 15)

saveRDS(object = all_Fibro, file="CustomVariables_results_all_default/all_default_deep_fibrosubcluster.RDS")