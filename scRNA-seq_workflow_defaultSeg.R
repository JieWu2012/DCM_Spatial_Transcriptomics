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




#set the working directory

setwd("/oak/stanford/groups/larsms/Users/jwu/SpatialTranscriptomics/data1/")


#find the samples with directory name with "rep"

ids <- list.files(path = "./", pattern = "rep")
ids <- ids[!grepl("txt", ids)]

ids

#read the cell feature matrix from output of resegmentation by default, exp=0,3,5
#read all the samples to one matrix

all.data <- sapply(ids, function(i){
  #d10x <- Read10X(file.path(".", i, "resegment_expdis0_formal/outs/cell_feature_matrix")) 
  d10x <- Read10X(file.path(".", i, "cell_feature_matrix"))
  d10x <- d10x$`Gene Expression`
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})


#read the information of cells for the x_centroid and y_centroid information. 

all.cells<- lapply(ids, function(i){
  
  #cells <- fread(file.path(".", i, "resegment_expdis0_formal/outs/cells.csv.gz"))
  cells <- fread(file.path(".", i, "cells.csv.gz"))
  cells$sample <- i
  cells$cell_id <- gsub("1$", i, cells$cell_id)
  cells
})

all.data_cb <- do.call("cbind", all.data)
all.cells_cb <- do.call("rbind", all.cells)
rownames(all.cells_cb) <- all.cells_cb$cell_id



#Create Seurat object
all <- CreateSeuratObject(counts = all.data_cb, project = "all", min.cells = 3, min.features = 3)

#Add cell information to the metadata of seurat object

all <- AddMetaData(all, metadata = all.cells_cb)

#extract the cell_id without the 
all$orig.ident <- sapply(strsplit(colnames(all), "-"), `[`, 2)


#QC
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^mt-")
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)


VlnPlot(all, features = c("nFeature_RNA"),pt.size = 0, group.by = "orig.ident")+NoLegend() +xlab("") + ylab("Number of genes") + scale_fill_manual(values = rep("steelblue", length(unique(all$orig.ident))))
ggsave(file="CustomVariables_results_all_default/nFeature.pdf", dpi = "print", height = 5, width = 7)

VlnPlot(all, features = c("nCount_RNA"),pt.size = 0, group.by = "orig.ident")+NoLegend() +xlab("") + ylab("Number of transcripts") + scale_fill_manual(values = rep("#E7298A", length(unique(all$orig.ident))))
ggsave(file="CustomVariables_results_all_default/nCounts.pdf", dpi = "print", height = 3, width = 7)
VlnPlot(all, features = c("nCount_RNA"),pt.size = 0, group.by = "orig.ident")+NoLegend() +xlab("") + ylab("Number of transcripts") + scale_fill_manual(values = rep("#E7298A", length(unique(all$orig.ident)))) + coord_cartesian(ylim = c(0,1000))
ggsave(file="CustomVariables_results_all_default/nCounts_lt1000.pdf", dpi = "print", height = 4, width = 7)


#rename mutant samples; order the samples to mutant, wt, BE; Rename samples with clear names
all@meta.data$ident_edited <- apply(all@meta.data, 1, function(x){a=x[1];a=gsub("PBS", "Mutant", a);return(a)})
samples <- c("Mutant_rep1", "Mutant_rep2", "Mutant_rep3", "Mutant_rep4", "WT_rep1", "WT_rep2", "WT_rep3", "WT_rep4", "WT_rep5", "BE_rep1", "BE_rep2", "BE_rep3", "BE_rep4")
all$ident_edited <- factor(all$ident_edited, levels = samples)
all$DiseaseStatus_edited <- apply(all@meta.data, 1, function(x){a=x[1];a=gsub("_rep[0-9]","", a);if(a=="BE"){return("Base-edited")}else if(a=="PBS"){return("Mutant")}else if(a=="WT"){return("WT")}})
all@meta.data$DiseaseStatus_edited=factor(all@meta.data$DiseaseStatus_edited, levels = c("Mutant", "WT", "Base-edited"))


#QC
VlnPlot(all, features = c("nFeature_RNA"),pt.size = 0, group.by = "ident_edited") + aes(fill=all@meta.data$DiseaseStatus_edited)+NoLegend() +xlab("Samples") + ylab("Number of genes") + scale_fill_manual(values = c("#f4a582", "#74add1", "#a6d96a")) + geom_boxplot(width = 0.3, outlier.shape = NA, color = "black")
ggsave(file="CustomVariables_results_all_default/nFeature_byDisease.pdf", dpi = "print", height = 5, width = 7)

VlnPlot(all, features = c("nCount_RNA"),pt.size = 0, group.by = "ident_edited") + aes(fill=all@meta.data$DiseaseStatus_edited)+NoLegend() +xlab("Samples") + ylab("Number of transcripts") + scale_fill_manual(values = c("#f4a582", "#74add1", "#a6d96a")) + geom_boxplot(width = 0.3, outlier.shape = NA, color = "black")
ggsave(file="CustomVariables_results_all_default/nCount_byDisease.pdf", dpi = "print", height = 5, width = 7)

VlnPlot(all, features = c("nCount_RNA"),pt.size = 0, group.by = "ident_edited") + aes(fill=all@meta.data$DiseaseStatus_edited)+NoLegend() +xlab("Samples") + ylab("Number of transcripts") + scale_fill_manual(values = c("#f4a582", "#74add1", "#a6d96a")) + geom_boxplot(width = 0.2, outlier.shape = NA, color = "black")  + coord_cartesian(ylim = c(0,1000))
ggsave(file="CustomVariables_results_all_default/nCount_byDisease_lt500.pdf", dpi = "print", height = 5, width = 7)

#filter genes
subset_all <- subset(all, subset = nCount_RNA >=50)
#all <- subset(all, subset = nCount_RNA >=50)

#SCTransform is equal to NormalizeData, FindVariableFeatures, ScaleData
all <- SCTransform(all)

#read panel genes
panel_gene <- read.table("CustomVariables_results/panel_gene_new.txt", sep="\t")

#filter splicing genes
panel_gene <- panel_gene[grepl("ENSMUSG", panel_gene$V1),]


#Use top 200 variable features
all <- FindVariableFeatures(all,nfeatures = 200)
gene <- VariableFeatures(all)

#only select variables genes with "ENSMUSG"
valgene <- gene[gene %in% panel_gene$V2]


#reset the variable genes of seruat object
VariableFeatures(all)<-valgene

# Run PCA
all <- RunPCA(all,approx = FALSE)

# Generate elbow plto
ElbowPlot(all,ndims = 30)
ggsave(file="CustomVariables_results_all_default/elbow_plot_all.pdf", dpi = "print", height = 5, width = 8)

# use dim=1:25 according to elbow plot
all <- FindNeighbors(all, dims = 1:25)

# set resolution to multiple values, 0.1, 0.3, 0.5
all <- FindClusters(all , resolution = 0.3)

# can change n.neighbors, min.dist, spread to make a clean UMAP
all <- RunUMAP(all, dims = 1:25, n.neighbors = 100, min.dist = 0.05, spread=2)
#all <- RunTSNE(all, dims = 1:25, check_duplicates = FALSE)

colors_cluster <- c(
  "#e6194b", # Red
  "#3cb44b", # Green
  "#ffe119", # Yellow
  "#0082c8", # Blue
  "#f58231", # Orange
  "#911eb4", # Purple
  "#46f0f0", # Cyan
  "#f032e6", # Magenta
  "#d2f53c", # Lime
  "#fabebe", # Pink
  "#008080", # Teal
  "#e6beff", # Lavender
  "#aa6e28", # Brown
  "#800000", # Maroon
  "#aaffc3", # Mint
  "#808000", # Olive
  "#ffd8b1", # Apricot
  "#000080", # Navy
  "#808080"  # Gray
)



# Generate UMAP for all cells. 
DimPlot(all, reduction = "umap", label = FALSE, repel=TRUE,pt.size = 0.1, cols=colors_cluster, raster=FALSE)+coord_fixed(ratio = 1) 
ggsave(file="CustomVariables_results_all_default/all_umap_nneibor100_mdist0.05.pdf", dpi = "print", height = 6, width = 9)


# Generate UMAP for all cells splitted by seurat clusters. 
DimPlot(all, reduction = "umap",raster=FALSE, split.by = c("seurat_clusters"),ncol = 4)+coord_fixed(ratio = 1) 
ggsave(file="CustomVariables_results_all_default/all_umap_nneibor100_mdist0.05_splitbycluster.pdf", dpi = "print", height = 10, width = 14)


# Generate UMAP for all cells splitted by Disease status. 
DimPlot(all, reduction = "umap",raster=FALSE, split.by = c("DiseaseStatus_edited"),cols = colors)+coord_fixed(ratio = 1) 
ggsave(file="CustomVariables_results_all_default/all_umap_nneibor100_mdist0.05_splitbyDiseaseStatus.pdf", dpi = "print", height = 6, width = 16)




## Find only positive markers for each clusters

all_markers_only_positive <-FindAllMarkers(all,only.pos = TRUE, features = valgene)

# function to select top 5 positive markers
gen_marker_table <- function(x){
  all_markers_only_positive[all_markers_only_positive$cluster == x, ] %>%
    head(n=5)
}

all_top5_markers_positiv <- map_dfr(unique(all_markers_only_positive$cluster), gen_marker_table)
write.table(all_top5_markers_positiv, "CustomVariables_results_all_default/all_top5_markers_positiv.txt", sep = "\t",quote = FALSE)

celltype.annotations <- c(
  "0" = "Endothelial",                 
  "1" = "Fibroblasts",                 
  "2" = "Ventricular Cardiomyocytes", 
  "3" = "Cardiomyocytes (Myh6+)",     
  "4" = "Stressed Cardiomyocytes",    
  "5" = "Macrophages",                
  "6" = "Atrial Cardiomyocytes",      
  "7" = "Fibroblast-like SMCs",       
  "8" = "Endocardial",                
  "9" = "Smooth Muscle Cells",        
  "10" = "Immature Cardiomyocytes",   
  "11" = "Cardiomyocytes (Myh7+)",    
  "12" = "Cardiomyocyte Subtype",     
  "13" = "T cells",                   
  "14" = "Fibroblast-like Macrophages",
  "15" = "Adipocytes"                 
)

celltype.annotations <- c(
  "0" = "Endothelial",                 
  "1" = "Fibroblasts",                 
  "2" = "Stressed Cardiomyocytes", 
  "3" = "Ventricular Cardiomyocytes",     
  "4" = "Ventricular Cardiomyocytes",    
  "5" = "Macrophages",                
  "6" = "Atrial Cardiomyocytes",      
  "7" = "Pericytes",       
  "8" = "Endocardial",                
  "9" = "Smooth Muscle Cells",        
  "10" = "Cardiomyocytes (F830016B08Rik+)",   
  "11" = "Stressed Cardiomyocytes",    
  "12" = "Cardiomyocytes (Myh7+)",     
  "13" = "T cells",                   
  "14" = "Fibroblasts",
  "15" = "Adipocytes"                 
)


colors <-  c(
  "#4477AA",  # Vivid Blue
  "#CC6677",  # Bold Red
  "#66CCAA",   # Fresh Green
  "#984EA3",  # Rich Purple
  "#FC8D59",  # Deep Coral
  "#A6761D",  # Warm Brown
  "#d2f53c", # Lime
  "#87CEEB", # Blue
  "#f032e6",  # Bright Scarlet
  "#E6AB02",  
  "#fabebe",  # Strong Green
  "#008080",  # Strong Teal
  "#999999"  # Mid Gray
)


colors <- c(
  "#800000", # Maroon
  "#e6194b", # Red
  "#ffe119", # Yellow
  "#3cb44b", # Green
  "#87CEEB", # Blue
  "#f58231", # Orange
  "#911eb4", # Purple
  "#f032e6", # Magenta
  "#d2f53c", # Lime
  "#008080", # Teal
  "#fabebe", # Pink
  "#aa6e28", # Brown
  "#377EB8"  # Clear Blue

)

color <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#6A3D9A")


# assign cell types for each clusters
all@meta.data$celltype = celltype.annotations[as.character(all$seurat_clusters)]


# make UMAPs
DimPlot(all, reduction = "umap", label = FALSE, repel=TRUE,pt.size = 0.1,raster=FALSE, group.by=c("celltype"),
        cols=colors, label.size = 5) + coord_fixed(ratio = 1) 
ggsave(file="CustomVariables_results_all_default/all_umap_Celltype_anno.pdf", dpi = "print", height = 7, width = 11)


DimPlot(all, reduction = "umap", label = FALSE, repel=TRUE,pt.size = 0.1,raster=FALSE, split.by = c("DiseaseStatus_edited"),group.by=c("celltype"),
        cols=colors, label.size = 5) + coord_fixed(ratio = 1) + theme(legend.position = "bottom") 
ggsave(file="CustomVariables_results_all_default/all_umap_Celltype_anno_byDisease.pdf", dpi = "print", height = 7, width = 15)

#Save RDS data. 
saveRDS(object = all, file="CustomVariables_results_all_default/all_default.RDS")


### downstream analysis

# change names of the cell type

all@meta.data$celltype_new <- celltype_new[as.character(all$celltype)]

orders_new <- c("VCMs", "Stressed VCMs", "VCMs (F830016B08Rik+)", "VCMs (Myh7+)", "ACMs", "Endothelial", "Endocardial", "SMCs", "Fibroblasts", "Macrophages", "Pericytes", "T cells", "Adipocytes")
#orders_new <- c("Normal vCMs", "Stress-responsive vCMs", "Ifgga4+ vCMs", "Myh7+ vCMs", "aCMs", "Endothelial cells", "Endocardial cells", "SMCs", "Fibroblasts", "Macrophages", "Pericytes", "T cells", "Adipocytes")
all$celltype_new <- factor(all$celltype_new, levels = orders_new)



# proportion of cell types in each sample. 

ggData = data.frame(prop.table(table(all$ident_edited, all$celltype_new), margin = 1))
colnames(ggData) = c("Samples", "Celltype", "value")
ggplot(ggData, aes(Samples, value, fill = Celltype)) + geom_col() + xlab("Samples") + ylab("Proportion of Cells (%)") + theme_bw() + theme(panel.grid = element_blank(), axis.text.x  = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = color_rmACM) + xlab("")


#use table(all$ident_edited, all$celltype_new) and check which column the ACMs is.
# remove ACMs. 
table(all$ident_edited, all$celltype_new)
color_rmACM <- colors[-5]
ggData = data.frame(prop.table(table(all$ident_edited, all$celltype_new)[,-5], margin = 1))
colnames(ggData) = c("Samples", "Celltype", "value")
ggplot(ggData, aes(Samples, value, fill = Celltype)) + geom_col() + xlab("Samples") + ylab("Proportion of Cells (%)") + theme_bw() + theme(panel.grid = element_blank(), axis.text.x  = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = color_rmACM) + xlab("")
ggsave(file="CustomVariables_results_all_default/all_prop_cell_samples_exclACMs_new.pdf", dpi = "print", height = 4, width = 8)





# generate violin plot for marker genes in each cell type


# find marker genes by "celltype_new"
all_markers_only_positive_celltype <-FindAllMarkers(all,only.pos = TRUE, group.by = "celltype_new")

gen_marker_table <- function(x){
  all_markers_only_positive_celltype[all_markers_only_positive_celltype$cluster == x, ] %>%
    head(n=5)
}

all_top5_markers_positive_celltype <- map_dfr(unique(all_markers_only_positive_celltype$cluster), gen_marker_table)


#make violin plot for marker genes for each cell type. 
genes_for_markerVln <- c("Rbm20", "Nppa", "F830016B08Rik", "Myh7", "Myl4", "Pecam1", "Vwf", "Acta2", "Dcn", "F13a1", "Rgs5", "Ccl5", "Adipoq")
Idents(all) <- "celltype_new"
VlnPlot(all, features = genes_for_markerVln, stack = TRUE, fill.by = "ident", cols = colors, flip = TRUE, pt.size = 0)
ggsave(file="CustomVariables_results_all_default/all_marker_celltype_vlnplot.pdf", dpi = "print", height = 7, width = 12)




## generate feature plot for genes

#feature plot for all cells 
#feature plot for AAV related genes

p <- FeaturePlot(
  all,
  reduction = "umap", # pick tsne or umap
  features = c("Rad50"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  raster=FALSE,
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)

LabelClusters(p, id = "ident",  color = "red")
ggsave(file=paste0("CustomVariables_results_all_default/all_umap_", "Rad50", ".pdf"), dpi = "print", height = 6, width = 14)



p <- FeaturePlot(
  all,
  reduction = "umap", # pick tsne or umap
  features = c("Mre11a"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  raster=FALSE,
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)

LabelClusters(p, id = "ident",  color = "red")
ggsave(file=paste0("CustomVariables_results_all_default/all_umap_", "Mre11a", ".pdf"), dpi = "print", height = 6, width = 14)


p <- FeaturePlot(
  all,
  reduction = "umap", # pick tsne or umap
  features = c("Zfp638"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  raster=FALSE,
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)

LabelClusters(p, id = "ident",  color = "red")
ggsave(file=paste0("CustomVariables_results_all_default/all_umap_", "Zfp638", ".pdf"), dpi = "print", height = 6, width = 14)

p <- FeaturePlot(
  all,
  reduction = "umap", # pick tsne or umap
  features = c("Gpr108"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  raster=FALSE,
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)

LabelClusters(p, id = "ident",  color = "red")
ggsave(file=paste0("CustomVariables_results_all_default/all_umap_", "Gpr108", ".pdf"), dpi = "print", height = 6, width = 14)


p <- FeaturePlot(
  all,
  reduction = "umap", # pick tsne or umap
  features = c("Itga7"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  raster=FALSE,
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)

LabelClusters(p, id = "ident",  color = "red")
ggsave(file=paste0("CustomVariables_results_all_default/all_umap_", "Itga7", ".pdf"), dpi = "print", height = 6, width = 14)


p <- FeaturePlot(
  all,
  reduction = "umap", # pick tsne or umap
  features = c("Sdc4"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 0.2,
  label = FALSE,
  repel = TRUE,
  raster=FALSE,
  split.by = "DiseaseStatus_edited"
)  + theme(legend.position = c(0.95,0.2)) & coord_fixed(ratio = 1) & scale_colour_gradientn(colours =color1)

LabelClusters(p, id = "ident",  color = "red")
ggsave(file=paste0("CustomVariables_results_all_default/all_umap_", "Sdc4", ".pdf"), dpi = "print", height = 6, width = 14)




## DE analysis between mutant and wt for specific cell typ. 


Idents(all) <- "celltype"
Stressed_CM_comparison <-
  FindMarkers(
    all,
    subset.ident = "Stressed Cardiomyocytes", 
    group.by = "DiseaseStatus_edited", #
    ident.1 = "Mutant", 
    ident.2 = "WT", 

    only.pos = FALSE,
    test.use = "MAST"
  )



# generate volcano plot for specific cell type. 

EnhancedVolcano(Stressed_CM_comparison,
                lab = rownames(Stressed_CM_comparison),
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'DEG in Mutant vs WT',
                subtitle = "Stressed Cardiomyocytes",
                caption = bquote(~Log[2]~ "fold change cutoff, 0.25; p-value cutoff, 10e-10"),
                pCutoff = 10e-10, # cut off line for p-value
                FCcutoff = 0.25, # cut off line for fold-change
                pointSize = 2.0,
                labSize = 5.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1,
                legendPosition = 'none',
                drawConnectors = T,
                widthConnectors = 0.3)




#spatial of cells with Cardiomyocytes (F830016B08Rik+)

df <- as.data.frame(all@meta.data[all@meta.data$orig.ident=="BE_rep2",])
ggplot(df, aes(x=x_centroid, y=y_centroid))+geom_point(size=0.1, color="gray") + 
  geom_point(data = df[df$celltype_new=="VCMs (F830016B08Rik+)",], aes(x=x_centroid, y=y_centroid, color=celltype_new), size=0.1) + 
  scale_color_manual(values = "#3cb44b", name= NULL) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "top") +  
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  coord_fixed(ratio = 1) + 
  scale_y_reverse()+
  ggplot(df, aes(x=x_centroid, y=y_centroid)) +
  geom_point(size=0.1, color="gray") + 
  geom_point(data = df[df$celltype_new=="Macrophages",], aes(x=x_centroid, y=y_centroid, color=celltype_new), size=0.1) + 
  scale_color_manual(values = "#f032e6", name= NULL) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  theme(legend.position = "top") +  
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  coord_fixed(ratio = 1) + 
  scale_y_reverse() +
  ggplot(df, aes(x=x_centroid, y=y_centroid)) +
  geom_point(size=0.1, color="gray") + 
  geom_point(data = df[df$celltype_new=="T cells",], aes(x=x_centroid, y=y_centroid, color=celltype_new), size=0.1) + 
  scale_color_manual(values = "#377EB8" , name= NULL) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  theme(legend.position = "top") +  
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  coord_fixed(ratio = 1) + 
  scale_y_reverse()

ggsave(file="CustomVariables_results_all_default/Spatial_BE_rep2_F830016B08Rik+_macrophages_new.pdf", dpi = "print", height = 5, width = 12)









# test the cell type proportion between mutant and wt. 

#extract mutant and wt
seu <- subset(all, subset= DiseaseStatus_edited!="Base-edited")
seu$group <- seu$DiseaseStatus_edited
seu$group <- factor(seu$group, levels = c("Mutant", "WT"))

#use "propeller" to test the significant differential proportion of cell types between mutant and wt
prop_result <- propeller(seu)
prop_df <- as.data.frame(prop_result)
prop_df$cluster <- factor(rownames(prop_df), levels = rownames(prop_df)[order(prop_df$Tstatistic)])
prop_df <- prop_df %>%
  mutate(sig_label = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01  ~ "**",
    FDR < 0.05  ~ "*",
    TRUE            ~ ""
  ))

#generate bar plot
ggplot(prop_df, aes(x = cluster, y = Tstatistic, fill = Tstatistic > 0)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(
    aes(label = sig_label, 
        vjust = ifelse(Tstatistic < 0, 1, 0.1)),  # below bar if negative, above if positive
    size = 5
  ) +
  scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue")) +
  theme_minimal(base_size = 12) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  labs(
    x = "",
    y = "T-statistic",
    title = "Differential Cell Type Proportions: Mutant vs WT "
  ) + 
  geom_hline(yintercept = 0, color = "black")



ggsave(file="CustomVariables_results_all_default/Proportion_celltype_Mutant_vs_WT.pdf", dpi = "print", height = 5, width = 6)

