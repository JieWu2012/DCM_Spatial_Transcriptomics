library(Seurat)
library(ggplot2)
library(randomForest)

setwd("/oak/stanford/groups/larsms/Users/jwu/SpatialTranscriptomics/data1/")

all <- readRDS("CustomVariables_results_all_default/all_default.RDS")


celltype_selected <- unique(all$celltype)[grepl("Cardiomyocytes", unique(all$celltype))]

celltype_selected <- celltype_selected[!grepl("Atrial", celltype_selected)]

all <- subset(all, subset = celltype %in% celltype_selected)

expr_mat <- Seurat::GetAssayData(all, assay = "SCT", slot = "data")

expr_mat <- expr_mat[grepl("-", rownames(expr_mat)),]
expr_mat <- expr_mat[!grepl("^mt-", rownames(expr_mat)),]
meta <- all@meta.data

table(meta$DiseaseStatus_edited)

train_cells <- rownames(meta[meta$DiseaseStatus_edited %in% c("WT", "Mutant"),])
test_cells <- rownames(meta[meta$DiseaseStatus_edited %in% c("Base-edited"),])

train_data <- t(as.matrix(expr_mat[, train_cells])) 
train_label <- meta[train_cells, "DiseaseStatus_edited"]
train_label <- droplevels(as.factor(train_label))

test_data <- t(as.matrix(expr_mat[, test_cells])) 

set.seed(123)


rf_model <- randomForest(
  x = train_data,
  y = as.factor(train_label),
  ntree = 500,
  importance = TRUE
)

rf_pred <- predict(rf_model, newdata = test_data, type = "prob")

meta$wt_prob <- NA
meta$mutant_prob <- NA

meta[test_cells, "wt_prob"] <- rf_pred[, "WT"]
meta[test_cells, "mutant_prob"] <- rf_pred[, "Mutant"]
meta[test_cells, "predicted_class"] <- ifelse(rf_pred[, "WT"] > 0.5, "wt_like", "mutant_like")
meta$predicted_class <- ifelse(
  is.na(meta$predicted_class),
  meta$celltype,
  meta$predicted_class
)


all@meta.data <- meta


DimPlot(all, reduction = "umap",pt.size = 1, group.by = "predicted_class", split.by = "DiseaseStatus_edited")+coord_fixed(ratio = 1)+ ggtitle(NULL)  + guides(color = guide_legend(override.aes = list(size = 5), nrow = 2, byrow = TRUE))+theme(legend.position = "bottom")
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt.pdf", width = 9, height = 5)
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT.pdf", width = 9, height = 5)

saveRDS(rf_model, file = "CustomVariables_results_all_default/rf_model.rds")
saveRDS(rf_model, file = "CustomVariables_results_all_default/rf_model_noMT.rds")




#use sub samples of WT and Mutant as test data

set.seed(123)
train_idx <- sample(seq_len(nrow(train_data)), 0.8 * nrow(train_data))

train_x <- train_data[train_idx, ]
train_y <- train_label[train_idx]

test_x  <- train_data[-train_idx, ]
test_y  <- train_label[-train_idx]

rf_model <- randomForest(x = train_x, y = train_y, ntree = 500, importance = TRUE)

rf_pred_test <- predict(rf_model, test_x, type = "prob")
roc_obj <- roc(test_y, rf_pred_test[, "Mutant"], levels=c("WT", "Mutant"), direction="<")
par(pty = "s")
pdf("CustomVariables_results_all_default/rf_prediction_model_noMT_crossvalidation0.2_ROC.pdf", height=5, width = 5)
plot.roc(roc_obj, col = "red",print.auc = TRUE, print.auc.x = 0.9, print.auc.y = 0.8)
dev.off()

meta_test_x <- meta[meta$cell_id %in% rownames(test_x),]
head(meta_test_x)
meta_test_x$wt_prob <- NA
meta_test_x$mutant_prob <- NA

meta_test_x[, "wt_prob"] <- rf_pred_test[, "WT"]
meta_test_x[, "mutant_prob"] <- rf_pred_test[, "Mutant"]
opt_thresh <- coords(roc_obj, "best", best.method = "youden")$threshold
opt_thresh
meta_test_x[, "predicted_class"] <- ifelse(rf_pred_test[, "Mutant"] > opt_thresh, "mutant_like", "wt_like")
meta_test_x$predicted_class <- ifelse(
  is.na(meta_test_x$predicted_class),
  meta_test_x$celltype,
  meta_test_x$predicted_class
)
df_cell_predict <- as.data.frame(table(meta_test_x$celltype, meta_test_x$predicted_class))
ggplot(df_cell_predict, aes(x=Var1, y=Freq, fill=Var2))+ geom_bar(stat="identity", position="dodge")+theme_bw()+coord_flip()+theme(panel.grid = element_blank())+ylab("Counts")+xlab("")+scale_fill_manual(name="Prediction", values=c("#ff7f00", "#0868ac"))
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2_test_x_celltypes.pdf", width = 8, height = 4)


#box plot of prediction score of mutant cells in each cell types 
ggplot(meta_test_x, aes(x=celltype, y=mutant_prob)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  xlab("") +
  ylab("Prediction score for mutant cells")
ggsave("CustomVariables_results_all_default/rf_prediction_score_mutant_crossvalidation0.2.pdf", width = 8, height = 4)

write.table(meta_test_x, "CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2_meta_test_x.txt",sep = "\t", row.names = F,quote = F)



test_cells <- rownames(meta[meta$DiseaseStatus_edited %in% c("Base-edited"),])
test_data <- t(as.matrix(expr_mat[, test_cells])) 
rf_pred_test_be <- predict(rf_model, test_data, type = "prob")
meta$wt_prob <- NA
meta$mutant_prob <- NA

meta[test_cells, "wt_prob"] <- rf_pred_test_be[, "WT"]
meta[test_cells, "mutant_prob"] <- rf_pred_test_be[, "Mutant"]

#Youden's J statistic to select cutoff for probability
opt_thresh <- coords(roc_obj, "best", best.method = "youden")$threshold

meta[test_cells, "predicted_class"] <- ifelse(rf_pred_test_be[, "Mutant"] > opt_thresh, "mutant_like", "wt_like")
meta$predicted_class <- ifelse(
  is.na(meta$predicted_class),
  meta$celltype,
  meta$predicted_class
)

write.table(meta, "CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2_meta.txt",sep = "\t", row.names = F,quote = F)
all@meta.data <- meta


DimPlot(all, reduction = "umap",pt.size = 1, group.by = "predicted_class", split.by = "DiseaseStatus_edited")+coord_fixed(ratio = 1)+ ggtitle(NULL)  + guides(color = guide_legend(override.aes = list(size = 5), nrow = 2, byrow = TRUE))+theme(legend.position = "bottom")
#ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt.pdf", width = 9, height = 5)
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2.pdf", width = 9, height = 5)

saveRDS(rf_model, file = "CustomVariables_results_all_default/rf_model_noMT_crossvalidation0.2.rds")

par(pty = "m")
rf_importance = importance(rf_model)
rf_importance_df <- data.frame(rf_importance)
rf_importance_df$gene <- rownames(rf_importance_df)
top30_features$gene <- factor(top30_features$gene,levels = rev(top30_features$gene))
ggplot(top30_features, aes(x=gene, y=MeanDecreaseGini))+geom_bar(stat = "identity", fill="darkgreen")+theme_bw()+coord_flip()+theme(panel.grid = element_blank())+ylab("Feature importance")+xlab("")
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2_featureimportance.pdf", width = 4, height = 5)

# BE samples cell types vs prediction
meta_BE <- meta[meta$DiseaseStatus_edited=="Base-edited",]
df_cell_predict <- as.data.frame(table(meta_BE$celltype, meta_BE$predicted_class))
ggplot(df_cell_predict, aes(x=Var1, y=Freq, fill=Var2))+ geom_bar(stat="identity", position="dodge")+theme_bw()+coord_flip()+theme(panel.grid = element_blank())+ylab("Counts")+xlab("")+scale_fill_manual(name="Prediction", values=c("#ff7f00", "#0868ac"))
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2_BE_celltypes.pdf", width = 8, height = 4)

# change prediction cut off
meta$wt_prob <- NA
meta$mutant_prob <- NA

meta[test_cells, "wt_prob"] <- rf_pred_test_be[, "WT"]
meta[test_cells, "mutant_prob"] <- rf_pred_test_be[, "Mutant"]
meta[test_cells, "predicted_class"] <- ifelse(rf_pred_test_be[, "Mutant"] > 0.8, "mutant_like", "wt_like")
meta$predicted_class <- ifelse(
  is.na(meta$predicted_class),
  meta$celltype,
  meta$predicted_class
)

meta_BE <- meta[meta$DiseaseStatus_edited=="Base-edited",]
df_cell_predict <- as.data.frame(table(meta_BE$celltype, meta_BE$predicted_class))
ggplot(df_cell_predict, aes(x=Var1, y=Freq, fill=Var2))+ geom_bar(stat="identity", position="dodge")+theme_bw()+coord_flip()+theme(panel.grid = element_blank())+ylab("Counts")+xlab("")+scale_fill_manual(name="Prediction", values=c("#ff7f00", "#0868ac"))
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2_BE_celltypes_cutoff0.8.pdf", width = 8, height = 4)


meta$wt_prob <- NA
meta$mutant_prob <- NA

meta[test_cells, "wt_prob"] <- rf_pred_test_be[, "WT"]
meta[test_cells, "mutant_prob"] <- rf_pred_test_be[, "Mutant"]
meta[test_cells, "predicted_class"] <- ifelse(rf_pred_test_be[, "Mutant"] > 0.9, "mutant_like", "wt_like")
meta$predicted_class <- ifelse(
  is.na(meta$predicted_class),
  meta$celltype,
  meta$predicted_class
)

meta_BE <- meta[meta$DiseaseStatus_edited=="Base-edited",]
df_cell_predict <- as.data.frame(table(meta_BE$celltype, meta_BE$predicted_class))
ggplot(df_cell_predict, aes(x=Var1, y=Freq, fill=Var2))+ geom_bar(stat="identity", position="dodge")+theme_bw()+coord_flip()+theme(panel.grid = element_blank())+ylab("Counts")+xlab("")+scale_fill_manual(name="Prediction", values=c("#ff7f00", "#0868ac"))
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2_BE_celltypes_cutoff0.9.pdf", width = 8, height = 4)


meta$wt_prob <- NA
meta$mutant_prob <- NA

meta[test_cells, "wt_prob"] <- rf_pred_test_be[, "WT"]
meta[test_cells, "mutant_prob"] <- rf_pred_test_be[, "Mutant"]
meta[test_cells, "predicted_class"] <- ifelse(rf_pred_test_be[, "Mutant"] > 0.7, "mutant_like", "wt_like")
meta$predicted_class <- ifelse(
  is.na(meta$predicted_class),
  meta$celltype,
  meta$predicted_class
)

meta_BE <- meta[meta$DiseaseStatus_edited=="Base-edited",]
df_cell_predict <- as.data.frame(table(meta_BE$celltype, meta_BE$predicted_class))
ggplot(df_cell_predict, aes(x=Var1, y=Freq, fill=Var2))+ geom_bar(stat="identity", position="dodge")+theme_bw()+coord_flip()+theme(panel.grid = element_blank())+ylab("Counts")+xlab("")+scale_fill_manual(name="Prediction", values=c("#ff7f00", "#0868ac"))
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2_BE_celltypes_cutoff0.7.pdf", width = 8, height = 4)


meta$wt_prob <- NA
meta$mutant_prob <- NA

meta[test_cells, "wt_prob"] <- rf_pred_test_be[, "WT"]
meta[test_cells, "mutant_prob"] <- rf_pred_test_be[, "Mutant"]
meta[test_cells, "predicted_class"] <- ifelse(rf_pred_test_be[, "Mutant"] > 0.8, "mutant_like", ifelse(rf_pred_test_be[, "Mutant"] < 0.5, "wt_like", "unknown"))
meta$predicted_class <- ifelse(
  is.na(meta$predicted_class),
  meta$celltype,
  meta$predicted_class
)

meta_BE <- meta[meta$DiseaseStatus_edited=="Base-edited",]
df_cell_predict <- as.data.frame(table(meta_BE$celltype, meta_BE$predicted_class))
ggplot(df_cell_predict, aes(x=Var1, y=Freq, fill=Var2))+ geom_bar(stat="identity", position="dodge")+theme_bw()+coord_flip()+theme(panel.grid = element_blank())+ylab("Counts")+xlab("")+scale_fill_manual(name="Prediction", values=c("#ff7f00", "gray" , "#0868ac"))
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2_BE_celltypes_cutoff0.8mutatn_0.5wt.pdf", width = 8, height = 4)


hist(meta$mutant_prob, breaks = 100)

#use cutoff 0.3 and 0.9

meta$wt_prob <- NA
meta$mutant_prob <- NA

meta[test_cells, "wt_prob"] <- rf_pred_test_be[, "WT"]
meta[test_cells, "mutant_prob"] <- rf_pred_test_be[, "Mutant"]
meta[test_cells, "predicted_class"] <- ifelse(rf_pred_test_be[, "Mutant"] > 0.9, "mutant_like", ifelse(rf_pred_test_be[, "Mutant"] < 0.3, "wt_like", "unknown"))
meta$predicted_class <- ifelse(
  is.na(meta$predicted_class),
  meta$celltype,
  meta$predicted_class
)

meta_BE <- meta[meta$DiseaseStatus_edited=="Base-edited",]
df_cell_predict <- as.data.frame(table(meta_BE$celltype, meta_BE$predicted_class))
ggplot(df_cell_predict, aes(x=Var1, y=Freq, fill=Var2))+ geom_bar(stat="identity", position="dodge")+theme_bw()+coord_flip()+theme(panel.grid = element_blank())+ylab("Counts")+xlab("")+scale_fill_manual(name="Prediction", values=c("#ff7f00", "gray" , "#0868ac"))
ggsave("CustomVariables_results_all_default/rf_prediction_mutant_wt_noMT_crossvalidation0.2_BE_celltypes_cutoff0.9mutatn_0.3wt.pdf", width = 8, height = 4)




FeaturePlot(
  subset_all,
  reduction = "umap", # pick tsne or umap
  features = c("mutant_prob"),
  keep.scale = "feature",
  order = T,
  label.size = 3,
  pt.size = 1.5,
  label = FALSE,
  repel = TRUE,
  #min.cutoff = "q10",
  #max.cutoff = "q90",
  raster=T
)  + theme(legend.position = c(0.95,0.2))  + scale_color_viridis_c(option = "D") & coord_fixed(ratio = 1)

ggsave(file="CustomVariables_results_all_default/rf_prediction_mutant_prob_featureplot.pdf", dpi = "print", height = 6, width = 14)


p1 <- ggplot(meta[meta$ident_edited=="BE_rep1",], aes(x=x_centroid, y=y_centroid))+geom_point(size=0.1, color="gray") + 
  geom_point(data = meta[meta$celltype=="Stressed Cardiomyocytes" & meta$ident_edited=="BE_rep1",], aes(x=x_centroid, y=y_centroid, color=celltype), size=0.1) + 
  scale_color_manual(values = "darkgreen", name= NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = "top") +  
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  coord_fixed(ratio = 1) + 
  scale_y_reverse()+
  facet_wrap(~ident_edited) +
  xlab("")+
  ylab("")


p2 <- ggplot(meta[meta$ident_edited=="BE_rep1",], aes(x=x_centroid, y=y_centroid))+geom_point(size=0.1, color="gray") + 
  geom_point(data = meta[meta$predicted_class=="mutant_like" & meta$ident_edited=="BE_rep1",], aes(x=x_centroid, y=y_centroid, color=predicted_class), size=0.1) + 
  scale_color_manual(values = "darkgreen", name= NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = "top") +  
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  coord_fixed(ratio = 1) + 
  scale_y_reverse()+
  facet_wrap(~ident_edited) +
  xlab("")+
  ylab("")

pdf("CustomVariables_results_all_default/BE_rep1_StressedVCMs_vs_mutant_like.pdf", height = 6, width = 10)
grid.arrange(p1, p2, ncol = 2)
dev.off()






