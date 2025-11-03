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

#make sure ids is only the directory with "rep", otherwise need to remove other files
ids <- list.files(path = "./", pattern = "rep")


cell_types <- unique(all@meta.data$celltype)


# k is 13 samples
for(k in 1:13){
print(ids[k])
df <- as.data.frame(all@meta.data[all@meta.data$orig.ident==ids[k],])


cross_moran <- data.frame()

# i and j are 13 cell types
for(i in 1:13){
  for(j in 1:13){
    if(i != j ){
    print(i)
    print(j)
    
    meta <- df %>%
      select(celltype, x_centroid,  y_centroid) %>%  mutate(
        is_A = as.numeric(celltype == cell_types[i]),
        is_B = as.numeric(celltype == cell_types[j])
      )
    coords <- meta[, c("x_centroid", "y_centroid")]
    knn <- knearneigh(coords, k = 10)
    nb <- knn2nb(knn)
    lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
    A <- meta$is_A
    B <- meta$is_B
    numerator <- sum(unlist(lapply(1:length(A), function(i) {
      sum(lw$weights[[i]] * (A[i] - mean(A)) * (B[lw$neighbours[[i]]] - mean(B)))
    })))
    denominator <- sqrt(sum((A - mean(A))^2) * sum((B - mean(B))^2))
    cross_moran_I <- length(A) / sum(unlist(lw$weights)) * numerator / denominator
    a = data.frame(cell1 = cell_types[i], cell2 = cell_types[j], cross_moran_I)
    cross_moran <- rbind(cross_moran, a)
    }
  }
}

#write.table(cross_moran, "./CustomVariables_results_all_default/cross_moran_I.txt", sep = "\t", quote = F)
write.table(cross_moran, paste0("./CustomVariables_results_all_default/cross_moran_I_", ids[k], ".txt"), sep = "\t", quote = F)


gggplot(cross_moran, aes(x = cell1, y = cell2, fill = cross_moran_I)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(ids[k])
ggsave(file=paste0("CustomVariables_results_all_default/Cross_Moran_I_all_", ids[k], ".pdf"), dpi = "print", height = 8, width = 10)



colnames(cross_moran) <- c("cell1", "cell2", "Cross_moran_I")
test <- cross_moran[cross_moran$cell1!=cross_moran$cell2,]
ggplot(test[test$cell2=="Cardiomyocytes (F830016B08Rik+)" | test$cell2=="Stressed Cardiomyocytes"  | test$cell2=="Cardiomyocytes (Myh7+)",], aes(x = cell1, y = cell2, fill = Cross_moran_I)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  xlab("") + ylab("") + ggtitle(ids[k])

#ggsave(file="CustomVariables_results_all_default/Cross_Moran_I.pdf", dpi = "print", height = 5, width = 8)
ggsave(file=paste0("CustomVariables_results_all_default/Cross_Moran_I_", ids[k], ".pdf"), dpi = "print", height = 5, width = 8)


}


#cross_moran_I for Fibroblast subclusters to other cell types. 

ids <- list.files(path = "./", pattern = "rep")


cell_types_new <- unique(all@meta.data$celltype_new)


# k is 13 samples
for(k in 1:13){
  print(ids[k])
  df <- as.data.frame(all@meta.data[all@meta.data$orig.ident==ids[k],])
  
  cross_moran <- data.frame()

# i and j are 13 cell types
  for(i in 1:13){
    for(j in c(5,12,14)){
      if(i != j ){
        print(i)
        print(j)
        
        meta <- df %>%
          select(celltype_new, x_centroid,  y_centroid) %>%  mutate(
            is_A = as.numeric(celltype_new == cell_types_new[i]),
            is_B = as.numeric(celltype_new == cell_types_new[j])
          )
        coords <- meta[, c("x_centroid", "y_centroid")]
        knn <- knearneigh(coords, k = 10)
        nb <- knn2nb(knn)
        lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
        A <- meta$is_A
        B <- meta$is_B
        numerator <- sum(unlist(lapply(1:length(A), function(i) {
          sum(lw$weights[[i]] * (A[i] - mean(A)) * (B[lw$neighbours[[i]]] - mean(B)))
        })))
        denominator <- sqrt(sum((A - mean(A))^2) * sum((B - mean(B))^2))
        cross_moran_I <- length(A) / sum(unlist(lw$weights)) * numerator / denominator
        a = data.frame(cell1 = cell_types_new[i], cell2 = cell_types_new[j], cross_moran_I)
        cross_moran <- rbind(cross_moran, a)
      }
    }
  }
  
  #write.table(cross_moran, "./CustomVariables_results_all_default/cross_moran_I.txt", sep = "\t", quote = F)
  write.table(cross_moran, paste0("./CustomVariables_results_all_default/cross_moran_I_Fibroblast2others", ids[k], ".txt"), sep = "\t", quote = F)
  

  
  
  
  colnames(cross_moran) <- c("cell1", "cell2", "Cross_moran_I")
  test <- cross_moran[cross_moran$cell1!=cross_moran$cell2,]
  ggplot(test, aes(x = cell1, y = cell2, fill = Cross_moran_I)) +
    geom_tile(color = "black") +
    scale_fill_viridis_c() +
    coord_fixed() +
    theme_minimal() + 
    theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    xlab("") + ylab("") + ggtitle(ids[k])
  
  #ggsave(file="CustomVariables_results_all_default/Cross_Moran_I.pdf", dpi = "print", height = 5, width = 8)
  ggsave(file=paste0("CustomVariables_results_all_default/Cross_Moran_I_Fibroblast2others", ids[k], ".pdf"), dpi = "print", height = 5, width = 8)
  
  
}






#draw cross Moran I index 
ids <- list.files(path = "./CustomVariables_results_all_default/", pattern = "^cross_moran_I.*.txt$")
ids_filtered <- ids[!grepl("Fibroblast", ids)]
ids_filtered <- ids_filtered[grepl("rep", ids_filtered)]
ids_filtered
cmi_all <- data.frame()
for(i in 1:13){
  cmi <- read.table(file.path("./CustomVariables_results_all_default",ids_filtered[i]), sep="\t", row.names = 1, header = TRUE)
  cmi$cell1_new <- celltype_new[cmi$cell1]
  cmi$cell2_new <- celltype_new[cmi$cell2]
  sample_name <- gsub("cross_moran_I_", "", ids_filtered[i])
  sample_name <- gsub(".txt", "", sample_name)
  sample_name <- gsub("PBS", "Mutant", sample_name)
  cmi_all <- rbind(cmi_all, data.frame(cmi, sample=sample_name))
}

samples <- c("Mutant_rep1", "Mutant_rep2", "Mutant_rep3", "Mutant_rep4", "WT_rep1", "WT_rep2", "WT_rep3", "WT_rep4", "WT_rep5", "BE_rep1", "BE_rep2", "BE_rep3", "BE_rep4")
cmi_all$cell1_new <- factor(cmi_all$cell1_new, levels = rev(orders_new))
cmi_all$sample <- factor(cmi_all$sample, levels = samples)
ggplot(cmi_all[cmi$cell2_new=="VCMs (F830016B08Rik+)", ], aes(y = cell1_new, x = sample, fill = cross_moran_I)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("") + ylab("") + ggtitle("cell types close to VCMs (F830016B08Rik+)")

ggsave(file="CustomVariables_results_all_default/cross_moran_I_celltoF830016B08Rik.pdf", dpi = "print", height = 5, width = 8)


ggplot(cmi_all[cmi$cell2_new=="VCMs (Myh7+)", ], aes(y = cell1_new, x = sample, fill = cross_moran_I)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("") + ylab("") + ggtitle("cell types close to VCMs (Myh7+)")
ggsave(file="CustomVariables_results_all_default/cross_moran_I_celltoMyh7.pdf", dpi = "print", height = 5, width = 8)

ggplot(cmi_all[cmi$cell2_new=="Stressed VCMs", ], aes(y = cell1_new, x = sample, fill = cross_moran_I)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("") + ylab("") + ggtitle("cell types close to Stressed VCMs")

ggsave(file="CustomVariables_results_all_default/cross_moran_I_celltoStressedVCM.pdf", dpi = "print", height = 5, width = 8)




##for all samples
## cell types to Ifgga4+ vCMs

ids <- list.files(path = "./CustomVariables_results_all_default/", pattern = "^cross_moran_I.*.txt$")
ids
ids_filtered <- ids[!grepl("Fibroblast", ids)]
ids_filtered <- ids_filtered[grepl("rep", ids_filtered)]
ids_filtered
cmi_all <- data.frame()

for(i in 1:13){
  cmi <- read.table(file.path("./CustomVariables_results_all_default",ids_filtered[i]), sep="\t", row.names = 1, header = TRUE)
  cmi$cell1_new <- celltype_new[cmi$cell1]
  cmi$cell2_new <- celltype_new[cmi$cell2]
  sample_name <- gsub("cross_moran_I_", "", ids_filtered[i])
  sample_name <- gsub(".txt", "", sample_name)
  sample_name <- gsub("PBS", "Mutant", sample_name)
  cmi_all <- rbind(cmi_all, data.frame(cmi, sample=sample_name))
}

samples <- c("Mutant_rep1", "Mutant_rep2", "Mutant_rep3", "Mutant_rep4", "WT_rep1", "WT_rep2", "WT_rep3", "WT_rep4", "WT_rep5", "BE_rep1", "BE_rep2", "BE_rep3", "BE_rep4")
cmi_all$cell1_new <- factor(cmi_all$cell1_new, levels = rev(orders_new))
cmi_all$sample <- factor(cmi_all$sample, levels = samples)

## cell types to Ifgga4+ vCMs
ggplot(cmi_all[cmi$cell2_new=="Ifgga4+ vCMs", ], aes(y = cell1_new, x = sample, fill = cross_moran_I)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text=element_text(size=14)) + 
  xlab("") + ylab("") + ggtitle("cell types close to Ifgga4+ vCMs")
ggsave(file="CustomVariables_results_all_default/cross_moran_I_celltoIfgga4vCMs.pdf", dpi = "print", height = 6, width = 8)

## cell types to Myh7+ vCMs

ggplot(cmi_all[cmi$cell2_new=="Myh7+ vCMs", ], aes(y = cell1_new, x = sample, fill = cross_moran_I)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text=element_text(size=14)) + 
  xlab("") + ylab("") + ggtitle("cell types close to Myh7+ vCMs")
ggsave(file="CustomVariables_results_all_default/cross_moran_I_cellto_myh7vCMs.pdf", dpi = "print", height = 6, width = 8)


## cell types to Stress-responsive vCMs

ggplot(cmi_all[cmi$cell2_new=="Stress-responsive vCMs", ], aes(y = cell1_new, x = sample, fill = cross_moran_I)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text=element_text(size=14)) + 
  xlab("") + ylab("") + ggtitle("cell types close to Stress-responsive vCMs")
ggsave(file="CustomVariables_results_all_default/cross_moran_I_cellto_Stress-responsive_vCMs.pdf", dpi = "print", height = 6, width = 8)


