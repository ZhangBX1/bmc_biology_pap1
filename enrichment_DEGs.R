library(ggplot2)
library(dplyr)
setwd("D://paper//single_cell//aa//marker")
library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)
library(DESeq2)
df=read.csv("all_clusters_markers.csv")

clusters <- unique(df$cluster)
kegg_results <- list()
go_results <- list()


for (cluster in clusters) {
  cluster_genes <- df %>% 
    dplyr::filter(cluster == !!cluster & avg_log2FC > 0.5 & pct.1 > 0.25  & p_val_adj < 0.05) %>% 
    dplyr::pull(gene)
  print(cluster)
  if (length(cluster_genes) > 0) {
    # GO富集分析
    go_res <- enrichGO(gene = cluster_genes,
                       OrgDb = org.At.tair.db,
                       keyType = "SYMBOL",
                       ont = "ALL",
                       pvalueCutoff = 0.05)
    if (!is.null(go_res) && nrow(go_res) > 0) {
      go_results[[as.character(cluster)]] <- go_res
    }
  }
}




kegg_df <- do.call(rbind,lapply(names(kegg_results),function(cluster) {
  kegg_res <- kegg_results[[cluster]]
  kegg_res_df <- as.data.frame(kegg_res)
  kegg_res_df$Cluster <- cluster
  return(kegg_res_df)
}))


kegg_df$Description <- gsub(" - Arabidopsis thaliana \\(thale cress\\)$","",kegg_df$Description)

# 绘制修改后的图
ggplot(kegg_df,aes(x = Cluster,y = Description,size = Count,color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#F6631C",  high = "#6D65A3") +
  theme_bw() +
  labs(title = "KEGG Enrichment",x = "Cluster",y = "KEGG Pathway",
       color = "Adjusted p-value",size = "Gene Count") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


# 处理GO结果
top_go_terms <- data.frame()
for (cluster in names(go_results)) {
  go_res <- go_results[[cluster]]
  go_df <- as.data.frame(go_res) %>%
    dplyr::arrange(p.adjust) %>%
    head(30)
  go_df$cluster <- cluster  # 使用小写的cluster
  top_go_terms <- rbind(top_go_terms,go_df)
}

# 确保数据框包含所需的列
print(colnames(top_go_terms))  # 检查列名

top_go_terms <- top_go_terms %>%
  dplyr::select(cluster,Description,Count,p.adjust) %>%  # 使用小写的cluster
  dplyr::mutate(Description = gsub(" - Arabidopsis thaliana \\(thale cress\\)","",Description))

# 绘制GO气泡图
ggplot(top_go_terms,aes(x = cluster,y = Description,size = Count,color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#C5272D",high = "#037F77") +
  theme_bw() +
  labs(title = "Top 5 GO Terms per Cluster",
       x = "Cluster",
       y = "GO Term",
       color = "Adjusted p-value",
       size = "Gene Count") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


kegg_df <- do.call(rbind,lapply(names(kegg_results),function(cluster) {
  kegg_res <- kegg_results[[cluster]]
  kegg_res_df <- as.data.frame(kegg_res)
  kegg_res_df$Cluster <- cluster
  return(kegg_res_df)
}))

# 导出KEGG富集分析结果
write.csv(kegg_df,file = "KEGG_all_clusters.csv",row.names = FALSE)

# 合并所有GO富集分析结果到一个数据框
go_df <- do.call(rbind,lapply(names(go_results),function(cluster) {
  go_res <- go_results[[cluster]]
  go_res_df <- as.data.frame(go_res)
  go_res_df$Cluster <- cluster
  return(go_res_df)
}))

# 导出GO富集分析结果
write.csv(go_df,file = "avg_log2FC_0.5_GO_all_clusters.csv",row.names = FALSE)



for (cluster in clusters) {
  cluster_genes <- df %>% 
    dplyr::filter(cluster == !!cluster & avg_log2FC >1 & pct.1 > 0.3 & p_val_adj > 0 & p_val_adj < 0.05) %>% 
    dplyr::pull(gene)
  print(cluster)
  if (length(cluster_genes) > 0) {
    # KEGG富集分析
    kegg_res <- enrichKEGG(gene = cluster_genes,
                           organism = 'ath',
                           keyType = "kegg",
                           pvalueCutoff = 0.05)
    if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
      kegg_results[[as.character(cluster)]] <- kegg_res
    }
    
    # GO富集分析
    go_res <- enrichGO(gene = cluster_genes,
                       OrgDb = org.At.tair.db,
                       keyType = "TAIR",
                       ont = "ALL",
                       pvalueCutoff = 0.05)
    if (!is.null(go_res) && nrow(go_res) > 0) {
      go_results[[as.character(cluster)]] <- go_res
    }
  }
}


kegg_df <- do.call(rbind,lapply(names(kegg_results),function(cluster) {
  kegg_res <- kegg_results[[cluster]]
  kegg_res_df <- as.data.frame(kegg_res)
  kegg_res_df$Cluster <- cluster
  return(kegg_res_df)
}))


kegg_df$Description <- gsub(" - Arabidopsis thaliana \\(thale cress\\)$","",kegg_df$Description)

# 绘制修改后的图
ggplot(kegg_df,aes(x = Cluster,y = Description,size = Count,color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#F6631C",  high = "#6D65A3") +
  theme_bw() +
  labs(title = "KEGG Enrichment",x = "Cluster",y = "KEGG Pathway",
       color = "Adjusted p-value",size = "Gene Count") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


# 处理GO结果
top_go_terms <- data.frame()
for (cluster in names(go_results)) {
  go_res <- go_results[[cluster]]
  go_df <- as.data.frame(go_res) %>%
    dplyr::arrange(p.adjust) %>%
    head(30)
  go_df$cluster <- cluster  # 使用小写的cluster
  top_go_terms <- rbind(top_go_terms,go_df)
}

# 确保数据框包含所需的列
print(colnames(top_go_terms))  # 检查列名

top_go_terms <- top_go_terms %>%
  dplyr::select(cluster,Description,Count,p.adjust) %>%  # 使用小写的cluster
  dplyr::mutate(Description = gsub(" - Arabidopsis thaliana \\(thale cress\\)","",Description))

# 绘制GO气泡图
ggplot(top_go_terms,aes(x = cluster,y = Description,size = Count,color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#C5272D",high = "#037F77") +
  theme_bw() +
  labs(title = "Top 5 GO Terms per Cluster",
       x = "Cluster",
       y = "GO Term",
       color = "Adjusted p-value",
       size = "Gene Count") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


kegg_df <- do.call(rbind,lapply(names(kegg_results),function(cluster) {
  kegg_res <- kegg_results[[cluster]]
  kegg_res_df <- as.data.frame(kegg_res)
  kegg_res_df$Cluster <- cluster
  return(kegg_res_df)
}))

# 导出KEGG富集分析结果
write.csv(kegg_df,file = "KEGG_all_clusters_up.csv",row.names = FALSE)

# 合并所有GO富集分析结果到一个数据框
go_df <- do.call(rbind,lapply(names(go_results),function(cluster) {
  go_res <- go_results[[cluster]]
  go_res_df <- as.data.frame(go_res)
  go_res_df$Cluster <- cluster
  return(go_res_df)
}))

# 导出GO富集分析结果
write.csv(go_df,file = "GO_all_clusters.csv",row.names = FALSE)


go_terms_of_interest <- unique(c("4-coumarate-CoA ligase activity",
                                 "4-coumarate-CoA ligase activity",
                                 "4-coumarate-CoA ligase activity",
                                 "CoA-ligase activity",
                                 "CoA-ligase activity",
                                 "daphnetin 3-O-glucosyltransferase activity",
                                 "daphnetin 3-O-glucosyltransferase activity",
                                 "enoyl-CoA hydratase activity",
                                 "flavonoid binding",
                                 "flavonol 3-O-glucosyltransferase activity",
                                 "flavonol 3-O-glucosyltransferase activity",
                                 "glucosyltransferase activity",
                                 "glucosyltransferase activity",
                                 "glucosyltransferase activity",
                                 "glucosyltransferase activity",
                                 "lignin biosynthetic process",
                                 "lignin biosynthetic process",
                                 "lignin metabolic process",
                                 "lignin metabolic process",
                                 "myricetin 3-O-glucosyltransferase activity",
                                 "myricetin 3-O-glucosyltransferase activity",
                                 "phenylpropanoid biosynthetic process",
                                 "phenylpropanoid biosynthetic process",
                                 "phenylpropanoid metabolic process",
                                 "photosynthesis",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 3-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercetin 7-O-glucosyltransferase activity",
                                 "quercitrin binding",
                                 "UDP-glucosyltransferase activity",
                                 "UDP-glucosyltransferase activity",
                                 "UDP-glucosyltransferase activity",
                                 "UDP-glucosyltransferase activity",
                                 "UDP-glucosyltransferase activity",
                                 "UDP-glucosyltransferase activity",
                                 "UDP-glucosyltransferase activity",
                                 "UDP-glucosyltransferase activity",
                                 "UDP-glycosyltransferase activity",
                                 "UDP-glycosyltransferase activity",
                                 "UDP-glycosyltransferase activity",
                                 "UDP-glycosyltransferase activity",
                                 "UDP-glycosyltransferase activity"
))
                                 
go_terms_of_interest <- unique(c(
  "photosynthesis",
  "photosynthesis, light reaction",
  "photosynthesis, light harvesting in photosystem I",
  "photosynthesis, light harvesting",
  "photorespiration",
  "carbon fixation",
  "photosynthesis, light harvesting in photosystem II",
  "fructose 1,6-bisphosphate metabolic process",
  "photosystem",
  "photosystem I",
  "photosystem I reaction center",
  "photosystem II",
  "photosynthesis, light harvesting in photosystem I",
  "photosynthesis, light harvesting",
  "photosynthesis, light reaction",
  "photosynthesis",
  "photosynthesis, light harvesting in photosystem II",
  "photosystem I",
  "photosystem",
  "photosystem I reaction center",
  "photosynthesis, light harvesting in photosystem I",
  "photosystem",
  "photosynthesis",
  "photosynthesis, light reaction",
  "photosynthesis, light harvesting in photosystem I",
  "carbon fixation",
  "photosynthesis, light harvesting",
  "photorespiration",
  "photosystem",
  "photosystem I",
  "photosystem I reaction center",
  "photosystem II",
  "photosynthesis, light harvesting in photosystem I",
  "photosynthesis, light harvesting",
  "photosynthesis",
  "photosynthesis, light harvesting in photosystem I",
  "photosynthesis, light reaction",
  "photosynthesis, light harvesting",
  "carbon fixation",
  "phloem loading",
  "vascular transport",
  "phloem transport",
  "photosystem",
  "photosystem I",
  "photosystem I reaction center",
  "photosystem II",
  "photosynthesis",
  "photosynthesis, light reaction",
  "photosynthesis, light harvesting in photosystem I",
  "photosynthesis, light harvesting",
  "carbon fixation",
  "photorespiration",
  "photosynthesis, light harvesting in photosystem II",
  "photosystem",
  "photosystem I",
  "photosystem I reaction center",
  "photosystem II"
))


# 创建结果数据框
plot_data <- data.frame()

# 创建cluster到细胞类型的映射
cluster_to_celltype <- c("0"="Phloem parenchyma",
                         "1"="Mesophyll",
                         "2"="Phloem parenchyma",
                         "3"="Mesophyll",
                         "4"="Mesophyll",
                         "5"="Phloem parenchyma",
                         "6"="Mesophyll",
                         "7"="Mesophyll",
                         "8"="Epidermis",
                         "9"="Vasculature",
                         "10"="Mesophyll",
                         "11"="Xylem",
                         "12"="Phloem",
                         "13"="Vasculature",
                         "14"="Mesophyll",
                         "15"="Epidermis",
                         "16"="Companion",
                         "17"="Guard")
# 设定细胞类型的显示顺序
cell_type_order <- c('Mesophyll','Vasculature','Epidermis','Companion cells','Guard cells','Others')

plot_data <- data.frame()

# 只关注特定的cluster
clusters_of_interest <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
clusters_of_interest <- c(1,3,4,6,7,10,14)
# 遍历每个感兴趣的cluster（修改后的循环）
for(cluster in as.character(clusters_of_interest)) {  # 确保cluster是字符型
  if(cluster %in% names(go_results)) {
    go_res <- go_results[[cluster]]
    
    for(term in go_terms_of_interest) {
      term_row <- go_res@result[go_res@result$Description == term, ]
      
      if(nrow(term_row) > 0) {
        plot_data <- rbind(plot_data, data.frame(
          Cluster = cluster,
          CellType = cluster_to_celltype[cluster],
          Term = term,
          p.adjust = term_row$p.adjust,
          GeneRatio = as.numeric(sub("/\\d+", "", term_row$GeneRatio)) / 
            as.numeric(sub("^\\d+/", "", term_row$GeneRatio))
        ))
      } else {
        plot_data <- rbind(plot_data, data.frame(
          Cluster = cluster,
          CellType = cluster_to_celltype[cluster],
          Term = term,
          p.adjust = 1,
          GeneRatio = 0
        ))
      }
    }
  } else {
    # 如果该cluster没有GO结果，添加空行
    for(term in go_terms_of_interest) {
      plot_data <- rbind(plot_data, data.frame(
        Cluster = cluster,
        CellType = cluster_to_celltype[cluster],
        Term = term,
        p.adjust = 1,
        GeneRatio = 0
      ))
    }
  }
}

# 按细胞类型对cluster重新排序
plot_data$CellType <- factor(plot_data$CellType, levels = cell_type_order)
plot_data <- plot_data[order(plot_data$CellType), ]
plot_data$Cluster <- factor(plot_data$Cluster, levels = unique(plot_data$Cluster))
plot_data$Term <- factor(plot_data$Term, levels = rev(go_terms_of_interest))
# 创建气泡图
p <- ggplot(plot_data, aes(x = Cluster, y = Term)) +
  geom_point(aes(size = GeneRatio, color = -log10(p.adjust)), alpha = 0.7) +
  scale_color_gradient(low = "darkblue", high = "red", 
                       breaks = seq(0, max(-log10(plot_data$p.adjust)), by = 1)) +
  coord_cartesian(ylim = c(0.5, length(unique(plot_data$Term)) + 1)) +
  ggtitle("Enriched pathways of genes (avg_log2FC>0.5) related to photosynthesis") +
  scale_size_continuous(range = c(0, 8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, "cm")
  ) +
  labs(
    x = "Cluster",
    y = "GO Term",
    color = "-log10(p.adjust)",
    size = "Gene Ratio"
  )

# 找出每个细胞类型的最后一个cluster的位置
cluster_positions <- seq_along(levels(plot_data$Cluster))
names(cluster_positions) <- levels(plot_data$Cluster)
cell_type_data <- unique(plot_data[, c("Cluster", "CellType")])
cell_type_changes <- cluster_positions[as.character(cell_type_data$Cluster[which(diff(as.integer(factor(cell_type_data$CellType))) != 0)])]

# 添加垂直分隔线
for(pos in cell_type_changes) {
  p <- p + geom_vline(xintercept = pos + 0.5, 
                      linetype = "dashed", 
                      color = "grey70", 
                      alpha = 0.5)
}

# 计算每个细胞类型标签的位置
cell_type_labels <- tapply(cluster_positions, 
                           cell_type_data$CellType, 
                           function(x) mean(range(x)))

# 添加细胞类型标签
p <- p + annotate("text",
                  x = unname(cell_type_labels),
                  y = length(unique(plot_data$Term)) + 1.2,
                  label = names(cell_type_labels),
                  size = 3)

# 添加cluster标签
p <- p + annotate("text",
                  x = cluster_positions,
                  y = length(unique(plot_data$Term)) + 0.7,
                  label = paste("C", levels(plot_data$Cluster), sep=""),
                  size = 2.5)

print(p)
