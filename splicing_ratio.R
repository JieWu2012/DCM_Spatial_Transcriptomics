#splicing ratio
library(Seurat)
library(ggplot2)
all <- readRDS("CustomVariables_results_all_default/all_default.RDS")


celltype_selected <- unique(all$celltype)[grepl("Cardiomyocytes", unique(all$celltype))]

celltype_selected <- celltype_selected[!grepl("Atrial", celltype_selected)]

all <- subset(all, subset = celltype %in% celltype_selected)

#pair Ttn-N2A (Mutant) vs Ttn-N2B (WT)
gene_exp_Ttn_N2A <- FetchData(all, vars = "Ttn-N2A", assay = "SCT", slot = "data")
gene_exp_Ttn_N2B <- FetchData(all, vars = "Ttn-N2B", assay = "SCT", slot = "data")
all[["Expr_Ttn_N2A"]] <- gene_exp_Ttn_N2A
all[["Expr_Ttn_N2B"]] <- gene_exp_Ttn_N2B
df <- all@meta.data
#df$M_W_ratio <- log2((df$`Expr_Ttn_N2A`+0.01)/(df$`Expr_Ttn_N2B`+0.01))
ggplot(df) + geom_histogram(aes(x=M_W_ratio), bins = 100, fill = "steelblue", color = "black") + facet_wrap(~ident_edited)

df$M_W_ratio <- (df$`Expr_Ttn_N2A`+0.01)/(df$`Expr_Ttn_N2B`+0.01 + df$`Expr_Ttn_N2A`+0.01)
df_new <- as.data.frame( df %>% group_by(M_W_ratio, DiseaseStatus_edited) %>% summarise(count = n()))
ggplot(df_new, aes(x=M_W_ratio,y=log2(count), color=DiseaseStatus_edited))+geom_point()+xlab("PSI")
ggsave("CustomVariables_results_all_default/Ttn_N2A_N2B_splicing_percentage_distribution_all.pdf", width = 13, height = 5)


df_new$bin <- cut(df_new$M_W_ratio, breaks = seq(0, 1, by = 0.2), include.lowest = F)
ggplot(df_new, aes(x=M_W_ratio,y=count, color=DiseaseStatus_edited))+geom_point()+facet_wrap(~bin, scales = "free")+theme(panel.grid.minor = element_blank())+theme_bw()+xlab("PSI")
ggsave("CustomVariables_results_all_default/Ttn_N2A_N2B_splicing_percentage_distribution.pdf", width = 10, height = 5)


ggplot(df_new, aes(x=M_W_ratio,y=log2(count), color=DiseaseStatus_edited))+geom_point()+facet_wrap(~bin, scales = "free")+theme(panel.grid.minor = element_blank())+theme_bw()+xlab("PSI")
ggsave("CustomVariables_results_all_default/Ttn_N2A_N2B_splicing_percentage_distribution_log2.pdf", width = 10, height = 5)

# generate raster heatmap plot of PSI score
df_full <- df_new %>% complete(M_W_ratio, DiseaseStatus_edited, fill = list(count = 0))
head(df_full)
df_full <- df_full %>% mutate(log2_count = log2(count+1))

ggplot(df_full, aes(x=factor(M_W_ratio),y=DiseaseStatus_edited, fill=log2_count))+geom_tile()+theme_bw()+
  scale_fill_viridis_c(option = "magma") +
  labs(x = "PSI", y = "Disease status", fill = "log2(count)")+scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + theme(axis.text.x = element_blank())

ggsave("CustomVariables_results_all_default/Ttn_N2A_N2B_splicing_percentage_distribution_log2_raster.pdf", width = 10, height = 2.5)

ggplot(df_full, aes(x=factor(M_W_ratio),y=DiseaseStatus_edited, fill=count))+geom_tile()+theme_bw()+
  scale_fill_viridis_c(option = "magma") +
  labs(x = "PSI", y = "Disease status", fill = "count")+scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + theme(axis.text.x = element_blank())

ggsave("CustomVariables_results_all_default/Ttn_N2A_N2B_splicing_percentage_distribution_raster.pdf", width = 10, height = 2.5)



#pair Ttn-alt-exon (Mutant) vs Ttn-flanking-exons (WT)
gene_exp_Ttn_alt <- FetchData(all, vars = "Ttn-alt-exon", assay = "SCT", slot = "data")
gene_exp_Ttn_flanking <- FetchData(all, vars = "Ttn-flanking-exons", assay = "SCT", slot = "data")
all[["Expr_Ttn_alt"]] <- gene_exp_Ttn_alt
all[["Expr_Ttn_flanking"]] <- gene_exp_Ttn_flanking
df <- all@meta.data

df$M_W_ratio <- (df$`Expr_Ttn_alt`+0.01)/(df$`Expr_Ttn_flanking`+0.01 + df$`Expr_Ttn_alt`+0.01)
df_new <- as.data.frame( df %>% group_by(M_W_ratio, DiseaseStatus_edited) %>% summarise(count = n()))
ggplot(df_new, aes(x=M_W_ratio,y=log2(count), color=DiseaseStatus_edited))+geom_point()+xlab("PSI")
ggsave("CustomVariables_results_all_default/Ttn_alt_flanking_splicing_percentage_distribution_all.pdf", width = 13, height = 5)


df_new$bin <- cut(df_new$M_W_ratio, breaks = seq(0, 1, by = 0.1), include.lowest = F)
ggplot(df_new, aes(x=M_W_ratio,y=count, color=DiseaseStatus_edited))+geom_point()+facet_wrap(~bin, scales = "free")+theme(panel.grid.minor = element_blank())+theme_bw()+xlab("PSI")
ggsave("CustomVariables_results_all_default/Ttn_alt_flanking_splicing_percentage_distribution.pdf", width = 10, height = 5)


ggplot(df_new, aes(x=M_W_ratio,y=log2(count), color=DiseaseStatus_edited))+geom_point()+facet_wrap(~bin, scales = "free")+theme(panel.grid.minor = element_blank())+theme_bw()+xlab("PSI")
ggsave("CustomVariables_results_all_default/Ttn_alt_flanking_splicing_percentage_distribution_log2.pdf", width = 10, height = 5)

# generate raster heatmap plot of PSI score
df_full <- df_new %>% complete(M_W_ratio, DiseaseStatus_edited, fill = list(count = 0))
head(df_full)
df_full <- df_full %>% mutate(log2_count = log2(count+1))

ggplot(df_full, aes(x=factor(M_W_ratio),y=DiseaseStatus_edited, fill=log2_count))+geom_tile()+theme_bw()+
scale_fill_viridis_c(option = "magma") +
labs(x = "PSI", y = "Disease status", fill = "log2(count)")+scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + theme(axis.text.x = element_blank())

ggsave("CustomVariables_results_all_default/Ttn_alt_flanking_splicing_percentage_distribution_log2_raster.pdf", width = 10, height = 2.5)

ggplot(df_full, aes(x=factor(M_W_ratio),y=DiseaseStatus_edited, fill=count))+geom_tile()+theme_bw()+
  scale_fill_viridis_c(option = "magma") +
  labs(x = "PSI", y = "Disease status", fill = "count")+scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + theme(axis.text.x = element_blank())

ggsave("CustomVariables_results_all_default/Ttn_alt_flanking_splicing_percentage_distribution_raster.pdf", width = 10, height = 2.5)



# set threshold and check the amount of mutant and wt cells in each sample
df$M_W_ratio_Ttn_alt_fla <- (df$`Expr_Ttn_alt`+0.01)/(df$`Expr_Ttn_flanking`+0.01 + df$`Expr_Ttn_alt`+0.01)
df$M_W_ratio_Ttn_N2A_N2B <- (df$`Expr_Ttn_N2A`+0.01)/(df$`Expr_Ttn_N2B`+0.01 + df$`Expr_Ttn_N2A`+0.01)
#df$cell_state <- apply(df, 1, function(x){if(x[32]>0.9925 | x[33]>0.65){return("mutant_cell")}else{return("wt_cell")}})
df$cell_state <- apply(df, 1, function(x){if(x[32]>0.9925 | x[33]>0.65){return("mutant_cell")}else if(x[32]<0.51 | x[33]<0.2){return("wt_cell")}else{return("unknown")}})
table(df$cell_state)
table(df$cell_state, df$DiseaseStatus_edited)
table(df$cell_state, df$celltype)

df_samples <- as.data.frame( df %>% group_by(cell_state, DiseaseStatus_edited) %>% summarise(count = n()))
ggplot(df_samples, aes(x=DiseaseStatus_edited, y=count, fill=cell_state))+geom_bar(stat = "identity",position = "dodge")
ggsave("CustomVariables_results_all_default/mutant_wt_unknown_cell_Disease_status.pdf", width = 6, height = 4)

df_samples <- as.data.frame( df %>% group_by(cell_state, celltype) %>% summarise(count = n()))
ggplot(df_samples, aes(x=celltype, y=count, fill=cell_state))+geom_bar(stat = "identity",position = "dodge")+theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(t = 5, r = 5, b = 5, l = 50))
ggsave("CustomVariables_results_all_default/mutant_wt_unknown_cell_celltype.pdf", width = 6, height = 4)


df$cell_state <- apply(df, 1, function(x){if(x[32]>0.9925 | x[33]>0.65){return("mutant_cell")}else{return("wt_cell")}})
table(df$cell_state)
table(df$cell_state, df$DiseaseStatus_edited)
table(df$cell_state, df$celltype)

df_samples <- as.data.frame( df %>% group_by(cell_state, DiseaseStatus_edited) %>% summarise(count = n()))
ggplot(df_samples, aes(x=DiseaseStatus_edited, y=count, fill=cell_state))+geom_bar(stat = "identity",position = "dodge")
ggsave("CustomVariables_results_all_default/mutant_wt_cell_Disease_status.pdf", width = 6, height = 4)

df_samples <- as.data.frame( df %>% group_by(cell_state, celltype) %>% summarise(count = n()))
ggplot(df_samples, aes(x=celltype, y=count, fill=cell_state))+geom_bar(stat = "identity",position = "dodge")+theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(t = 5, r = 5, b = 5, l = 50))
ggsave("CustomVariables_results_all_default/mutant_wt_cell_celltype.pdf", width = 6, height = 4)

#new cut off
df$cell_state <- apply(df, 1, function(x){if(as.numeric(x[32])>0.9950){return("mutant_cell")}else{return("wt_cell")}})
df_samples <- as.data.frame( df %>% group_by(cell_state, DiseaseStatus_edited) %>% summarise(count = n()))
ggplot(df_samples, aes(x=DiseaseStatus_edited, y=count, fill=cell_state))+geom_bar(stat = "identity",position = "dodge")
ggsave("CustomVariables_results_all_default/mutant_wt_cell_Disease_status_Ttn_alt_fla.pdf", width = 6, height = 4)

df$cell_state <- apply(df, 1, function(x){if(as.numeric(x[33])>0.8){return("mutant_cell")}else{return("wt_cell")}})
df_samples <- as.data.frame( df %>% group_by(cell_state, DiseaseStatus_edited) %>% summarise(count = n()))
ggplot(df_samples, aes(x=DiseaseStatus_edited, y=count, fill=cell_state))+geom_bar(stat = "identity",position = "dodge")
ggsave("CustomVariables_results_all_default/mutant_wt_cell_Disease_status_Ttn_N2A_N2B.pdf", width = 6, height = 4)

df$cell_state <- apply(df, 1, function(x){if(as.numeric(x[32])>0.9950 | as.numeric(x[33])>0.8){return("mutant_cell")}else{return("wt_cell")}})
df_samples <- as.data.frame( df %>% group_by(cell_state, DiseaseStatus_edited) %>% summarise(count = n()))
ggplot(df_samples, aes(x=DiseaseStatus_edited, y=count, fill=cell_state))+geom_bar(stat = "identity",position = "dodge")
ggsave("CustomVariables_results_all_default/mutant_wt_cell_Disease_status_bothpairs.pdf", width = 6, height = 4)

df$cell_state <- apply(df, 1, function(x){if(as.numeric(x[32])>0.9950 & as.numeric(x[33])>0.8){return("mutant_cell")}else{return("wt_cell")}})
df_samples <- as.data.frame( df %>% group_by(cell_state, DiseaseStatus_edited) %>% summarise(count = n()))
ggplot(df_samples, aes(x=DiseaseStatus_edited, y=count, fill=cell_state))+geom_bar(stat = "identity",position = "dodge")
ggsave("CustomVariables_results_all_default/mutant_wt_cell_Disease_status_bothpairs_bothpairs_and.pdf", width = 6, height = 4)
df_samples <- as.data.frame( df %>% group_by(cell_state, celltype) %>% summarise(count = n()))

# check the amount of mutant and wt cells in each cell type. 
ggplot(df_samples, aes(x=celltype, y=count, fill=cell_state))+geom_bar(stat = "identity",position = "dodge")+theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(t = 5, r = 5, b = 5, l = 50))
ggsave("CustomVariables_results_all_default/mutant_wt_cell_celltype_bothpairs_and.pdf", width = 6, height = 4)



#for all pairs:
Mu <- c("Ttn-alt-exon", "Ttn-N2A", "Ldb3-MUT-alt-exon", "Camk2d-isoformB", "Pdlim5-flanking-exons", "Immt-alt-exon-junction")
WT <- c("Ttn-flanking-exons", "Ttn-N2B", "Ldb3-WT-alt-exon", "Camk2d-isoformA", "Pdlim5-double-junction", "Immt-flanking-exons")

df <- all@meta.data


pdf("CustomVariables_results_all_default/splicing_pari_ratio_histogram.pdf", width = 8, height = 5)
for(i in 1:6){
  df <- all@meta.data
  Mutant_expr <- FetchData(all, vars = Mu[i], assay = "SCT", slot = "data")[[1]]
  WT_expr <- FetchData(all, vars = WT[i], assay = "SCT", slot = "data")[[1]]
  df$Mutant_expr <- Mutant_expr
  df$WT_expr <- WT_expr
  df$M_W_ratio <- log2((df$Mutant_expr+0.01)/(df$WT_expr+0.01))
  p <- ggplot(df) + 
    geom_histogram(aes(x=M_W_ratio), bins = 100, fill = "steelblue", color = "black") + 
    facet_wrap(~ident_edited) + 
    theme(panel.grid.minor = element_blank()) + 
    theme_bw() + 
    xlab("log2(ratio)") +
    ggtitle(paste(Mu[i], "_vs_", WT[i], sep = ""))
  print(p)
}

dev.off()



pdf("CustomVariables_results_all_default/splicing_pair_percentage_histogram.pdf", width = 8, height = 5)
for(i in 1:6){
  df <- all@meta.data
  Mutant_expr <- FetchData(all, vars = Mu[i], assay = "SCT", slot = "data")[[1]]
  WT_expr <- FetchData(all, vars = WT[i], assay = "SCT", slot = "data")[[1]]
  df$Mutant_expr <- Mutant_expr
  df$WT_expr <- WT_expr
  df$M_W_ratio <- log2((df$Mutant_expr+0.01)/(df$WT_expr+df$Mutant_expr))
  p <- ggplot(df) + 
    geom_histogram(aes(x=M_W_ratio), bins = 100, fill = "steelblue", color = "black") + 
    facet_wrap(~ident_edited) + 
    theme(panel.grid.minor = element_blank()) + 
    theme_bw() + 
    xlab("log2(ratio)") +
    ggtitle(paste(Mu[i], "_vs_", WT[i], sep = ""))
  print(p)
}

dev.off()



for(i in 1:6){
  df <- all@meta.data
  Mutant_expr <- FetchData(all, vars = Mu[i], assay = "SCT", slot = "data")[[1]]
  WT_expr <- FetchData(all, vars = WT[i], assay = "SCT", slot = "data")[[1]]
  df$Mutant_expr <- Mutant_expr
  df$WT_expr <- WT_expr
  df$M_W_ratio <- log2((df$Mutant_expr+0.01)/(df$WT_expr+0.01))
  samples <- unique(df$ident_edited)
  samples=samples[!grepl("WT_rep3", samples)]
  rot_degrees <- data.frame(ident_edited=samples, degree = c(90, -30, 90, 90, 30, -60, -20, 90, 180, 0, 180, 0))
  df$ident_edited <- factor(df$ident_edited, levels=c("WT_rep1", "WT_rep2", "WT_rep4", "WT_rep5", "Mutant_rep1", "Mutant_rep2", "Mutant_rep3", "Mutant_rep4", "BE_rep1", "BE_rep2", "BE_rep3", "BE_rep4"))
  df_merge <- merge(df, rot_degrees, by=c("ident_edited"))
  new_coord <- t(apply(df_merge, 1, function(x){theta <- as.numeric(x[30]) * pi / 180; x_rot <- as.numeric(x[6]) * cos(theta) + as.numeric(x[7]) * sin(theta); y_rot <- -as.numeric(x[6]) * sin(theta) + as.numeric(x[7]) * cos(theta);c(x_rot, y_rot)}))
  colnames(new_coord) <- c("x_rot", "y_rot")
  df_merge <- cbind(df_merge, new_coord)
  
  png(filename = paste0("CustomVariables_results_all_default/FeaturePlot/FeaturePlot_splicingRatio_", Mu[i], "_vs_", WT[i], ".png"), width = 1400, height = 1000)
  
  p <- ggplot(df_merge[df_merge$ident_edited!="WT_rep3",], aes(x=x_rot, y=y_rot, color=M_W_ratio)) + 
    geom_point(size=0.1) + 
    facet_wrap(~ident_edited, scales = "free") + 
    theme_minimal() + 
    theme(panel.grid = element_blank(), axis.title = element_blank(),axis.ticks = element_blank(), axis.text = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 25,colour = "white", face = "bold"), plot.background  = element_rect(fill = "black"), legend.text = element_text(color = "white"), legend.title = element_text(color="white"))  + 
    scale_color_gradient(name = paste(Mu[i], "_vs_", WT[i], sep = ""), low = "white", high = "darkred")
  
  
  print(p)
  
  dev.off()
}

cell_type <- unique(all$celltype)

for(i in 1:13){
  df <- all@meta.data
  
  samples <- unique(df$ident_edited)
  samples=samples[!grepl("WT_rep3", samples)]
  rot_degrees <- data.frame(ident_edited=samples, degree = c(90, -30, 90, 90, 30, -60, -20, 90, 180, 0, 180, 0))
  df$ident_edited <- factor(df$ident_edited, levels=c("WT_rep1", "WT_rep2", "WT_rep4", "WT_rep5", "Mutant_rep1", "Mutant_rep2", "Mutant_rep3", "Mutant_rep4", "BE_rep1", "BE_rep2", "BE_rep3", "BE_rep4"))
  df_merge <- merge(df, rot_degrees, by=c("ident_edited"))
  new_coord <- t(apply(df_merge, 1, function(x){theta <- as.numeric(x[27]) * pi / 180; x_rot <- as.numeric(x[6]) * cos(theta) + as.numeric(x[7]) * sin(theta); y_rot <- -as.numeric(x[6]) * sin(theta) + as.numeric(x[7]) * cos(theta);c(x_rot, y_rot)}))
  colnames(new_coord) <- c("x_rot", "y_rot")
  df_merge <- cbind(df_merge, new_coord)
  df_merge$cell_type_new <- apply(df_merge, 1, function(x){if(x[26]==cell_type[i]){return(x[26])}else{return("")}})
  pdf(file = paste0("CustomVariables_results_all_default/FeaturePlot/FeaturePlot_", cell_type[i], ".pdf"), width = 14, height = 10)
  df_merge[df_merge$ident_edited!="WT_rep3",]
  p <- ggplot(df_merge, aes(x=x_rot, y=y_rot)) + 
    geom_point(data = subset(df_merge, cell_type_new != cell_type[i]), size=0.1, color="white") + 
    geom_point(data = subset(df_merge, cell_type_new == cell_type[i]), size=0.1, color="darkred") + 
    facet_wrap(~ident_edited, scales = "free") + 
    theme_minimal() + 
    theme(panel.grid = element_blank(), plot.title = element_text(color = "darkred", size = 16, face = "bold"), axis.title = element_blank(),axis.ticks = element_blank(), axis.text = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 15,colour = "white", face = "bold"), plot.background  = element_rect(fill = "black"), legend.text = element_text(color = "darkred")) + 
    ggtitle(cell_type[i])
  #scale_color_manual(name = , values = c( "white", "darkred"))
  
  
  print(p)
  
  dev.off()
}



#PSI distribution for all pairs

Mu <- c("Ttn-alt-exon", "Ttn-N2A", "Ldb3-MUT-alt-exon", "Camk2d-isoformB", "Pdlim5-flanking-exons", "Immt-alt-exon-junction")
WT <- c("Ttn-flanking-exons", "Ttn-N2B", "Ldb3-WT-alt-exon", "Camk2d-isoformA", "Pdlim5-double-junction", "Immt-flanking-exons")


pdf("CustomVariables_results_all_default/all_splicing_pairs_PSI_distribution_raster.pdf", width = 10, height = 2.5)
for(i in 1:6){
  df <- all@meta.data
  Mutant_expr <- FetchData(all, vars = Mu[i], assay = "SCT", slot = "data")[[1]]
  WT_expr <- FetchData(all, vars = WT[i], assay = "SCT", slot = "data")[[1]]
  df$Mutant_expr <- Mutant_expr
  df$WT_expr <- WT_expr
  df$PSI <- (df$Mutant_expr+0.01)/(df$WT_expr+0.01 + df$Mutant_expr+0.01)
  df_new <- as.data.frame( df %>% group_by(PSI, DiseaseStatus_edited) %>% summarise(count = n()))
  df_full <- df_new %>% complete(PSI, DiseaseStatus_edited, fill = list(count = 0))
  df_full <- df_full %>% mutate(log2_count = log2(count+1))

  p=ggplot(df_full, aes(x=factor(PSI),y=DiseaseStatus_edited, fill=log2_count))+geom_tile()+theme_bw()+
    scale_fill_viridis_c(option = "magma") +
    labs(x = "PSI", y = "Disease status", fill = "log2(count)")+scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + theme(axis.text.x = element_blank()) +
    ggtitle(paste(Mu[i], WT[i],sep=":"))
  print(p)
  }
dev.off()
