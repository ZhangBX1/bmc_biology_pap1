# 加载必要的包
library(readxl)
library(limma)
library(splines)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci) 
library(ggrepel) 
library(patchwork) 
library(viridis)  
library(pheatmap)  

# ============== 1 ==============

data <- read_excel("meta.xlsx")
data <- as.data.frame(data)

metabolite_names <- data[, ncol(data)]

name_counts <- table(metabolite_names)
duplicate_names <- names(name_counts[name_counts > 1])

for (dup_name in duplicate_names) {
  indices <- which(metabolite_names == dup_name)
  for (i in 1:length(indices)) {
    metabolite_names[indices[i]] <- paste0(dup_name, "-", i)
  }
}

rownames(data) <- metabolite_names
data <- data[, -ncol(data)]

edesign <- read_excel("group.xlsx", col_names = TRUE)
edesign <- as.data.frame(edesign)

row_names <- edesign[[1]] 
edesign <- edesign[, -1] 
rownames(edesign) <- row_names 

if(!all(rownames(edesign) == colnames(data))) {
  warning("edesign的行名与data的列名不匹配，请检查数据")
  print("edesign行名：")
  print(rownames(edesign))
  print("data列名：")
  print(colnames(data))
  stop("样本名不匹配，分析终止")
}

if(any(data <= 0, na.rm = TRUE)) {
  # 添加一个小的常数以避免log2(0)
  min_non_zero <- min(data[data > 0], na.rm = TRUE)
  offset <- min_non_zero / 10
  cat("数据包含0或负值，添加偏移量", offset, "后进行log2转换\n")
  data_log2 <- log2(data + offset)
} else {
  data_log2 <- log2(data)
}
data_log2=data
# ============== 第2步：时间序列差异分析 ==============
# 从edesign提取时间点和分组信息
times <- edesign$Time
groups <- factor(ifelse(edesign$Control == 1, "Control", "Mutant"))

# 创建基于自然样条的设计矩阵
# 使用3个自由度的自然样条（可以根据时间点数量调整）
design <- model.matrix(~ ns(times, df=3) * groups)
colnames(design)

# 使用limma拟合模型
fit <- lmFit(data_log2, design)
fit <- eBayes(fit)

# 提取差异表达结果
# 获取组间差异随时间变化的结果
interaction_cols <- grep("ns.*groupsMutant", colnames(design))
results <- topTable(fit, coef=interaction_cols, n=Inf)

# 筛选显著差异的代谢物
sig_genes <- rownames(results[results$adj.P.Val < 0.01, ])
cat(paste("找到", length(sig_genes), "个显著差异的代谢物\n"))
sig_genes <- rownames(data_log2)

# ============== 第3步：聚类分析 ==============
# 提取显著代谢物的表达数据
sig_data <- data_log2[sig_genes, ]

# 计算相关性距离
gene_dist <- as.dist(1 - cor(t(sig_data)))

# 层次聚类
gene_hclust <- hclust(gene_dist, method = "complete")

# 确定聚类数
k <- 6


# 切割树形图
clusters <- cutree(gene_hclust, k=k)
names(clusters) <- sig_genes

cat(paste("将显著差异代谢物分为", k, "个聚类\n"))

# ============== 第4步：准备可视化数据 ==============
# 准备绘图数据
plot_data <- data.frame()
for(i in 1:length(sig_genes)) {
  gene_data <- data.frame(
    Gene = sig_genes[i],
    Time = times,
    Group = groups,
    Expression = as.numeric(sig_data[i, ]),
    Cluster = paste("Cluster", clusters[sig_genes[i]])
  )
  plot_data <- rbind(plot_data, gene_data)
}

# 计算每个时间点和组别的平均值和标准误
metabolite_means <- plot_data %>%
  group_by(Gene, Time, Group, Cluster) %>%
  summarize(
    Mean = mean(Expression),
    .groups = "drop"
  )

# 获取每个聚类在每个时间点、每个组的平均值和95%置信区间
cluster_stats <- plot_data %>%
  group_by(Time, Group, Cluster) %>%
  summarize(
    Mean = mean(Expression),
    SE = sd(Expression)/sqrt(n()),
    CI_lower = Mean - 1.96*SE,
    CI_upper = Mean + 1.96*SE,
    n_genes = n_distinct(Gene),
    .groups = "drop"
  )

# ============== 第5步：为每个聚类中的代谢物分配颜色 ==============
# 获取每个聚类中的唯一代谢物
cluster_genes <- split(sig_genes, clusters[sig_genes])

# 为每个聚类中的代谢物分配颜色
gene_colors <- list()
for(cluster_name in names(cluster_genes)) {
  genes_in_cluster <- cluster_genes[[cluster_name]]
  n_genes <- length(genes_in_cluster)
  
  # 使用viridis为每个聚类创建不同的颜色方案
  viridis_option <- switch(as.character(cluster_name),
                           "1" = "magma",
                           "2" = "plasma",
                           "3" = "viridis",
                           "4" = "cividis",
                           "5" = "rocket",
                           "6" = "mako",
                           "magma")  # 默认
  
  # 创建颜色向量
  colors <- viridis::viridis_pal(option = viridis_option)(n_genes)
  gene_colors[[cluster_name]] <- setNames(colors, genes_in_cluster)
}

# 将颜色信息添加到数据中
metabolite_means$Color <- sapply(1:nrow(metabolite_means), function(i) {
  gene <- metabolite_means$Gene[i]
  cluster <- gsub("Cluster ", "", metabolite_means$Cluster[i])
  return(gene_colors[[cluster]][gene])
})

# ============== 第6步：创建增强版聚类图 ==============
p2_enhanced <- ggplot() +
  # 1. 添加每个代谢物的线条 (半透明)
  geom_line(data = metabolite_means, 
            aes(x = Time, y = Mean, group = interaction(Gene, Group), 
                color = Gene, linetype = Group),
            size = 0.8, alpha = 0.4) +
  
  # 2. 添加聚类平均线 (粗线)
  geom_line(data = cluster_stats, 
            aes(x = Time, y = Mean, group = Group, linetype = Group),
            size = 2.5, color = "black", alpha = 0.7) +
  
  # 3. 添加置信区间
  geom_ribbon(data = cluster_stats,
              aes(x = Time, ymin = CI_lower, ymax = CI_upper, 
                  fill = Cluster, group = Group),
              alpha = 0.2) +
  
  # 4. 添加平均值点
  geom_point(data = cluster_stats,
             aes(x = Time, y = Mean, fill = Cluster, shape = Group),
             size = 5, color = "black", stroke = 1.2) +
  
  # 5. 添加代谢物数量标签
  geom_text(data = cluster_stats %>% 
              filter(Time == max(Time)) %>%
              group_by(Cluster) %>%
              slice(1),
            aes(x = max(Time) + 0.2, y = Mean, 
                label = paste0("n=", n_genes)),
            hjust = 0, size = 4, fontface = "bold") +
  
  # 6. 设置分面
  facet_wrap(~Cluster, scales = "free_y", ncol = 2) +
  
  # 7. 设置主题和标签
  labs(title = "Metabolite Expression Patterns by Cluster",
       subtitle = "Individual metabolites shown with cluster average (thick line) and 95% CI (shaded area)",
       x = "Time Point",
       y = "Expression Level (log2)",
       color = "Metabolite",
       linetype = "Group",
       shape = "Group",
       fill = "Cluster") +
  
  # 8. 设置颜色和形状
  scale_shape_manual(values = c("Control" = 21, "Mutant" = 24)) +
  scale_fill_manual(values = scales::hue_pal()(length(unique(cluster_stats$Cluster)))) +
  scale_linetype_manual(values = c("Control" = "solid", "Mutant" = "dashed")) +
  
  # 9. 设置x轴刻度
  scale_x_continuous(breaks = sort(unique(times))) +
  
  # 10. 设置主题
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "none",  # 隐藏图例，太多代谢物会使图例过大
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "black"),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
  )

# 保存增强版聚类图
ggsave("enhanced_cluster_trends.pdf", p2_enhanced, width = 12, height = 10, dpi = 300)

#============= 6改 =====================

for(group_name in c("Control", "Mutant")) {
  # 筛选该组别的数据
  group_metabolite_means <- metabolite_means %>% filter(Group == group_name)
  group_cluster_stats <- cluster_stats %>% filter(Group == group_name)
  
  p2_group <- ggplot() +
    # 1. 添加每个代谢物的线条
    geom_line(data = group_metabolite_means, 
              aes(x = Time, y = Mean, group = Gene, color = Gene),
              size = 0.8, alpha = 0.4) +
    
    # 2. 添加聚类平均线
    geom_line(data = group_cluster_stats, 
              aes(x = Time, y = Mean),
              size = 2.5, color = "black", alpha = 0.7) +
    
    # 3. 添加置信区间
    geom_ribbon(data = group_cluster_stats,
                aes(x = Time, ymin = CI_lower, ymax = CI_upper, fill = Cluster),
                alpha = 0.2) +
    
    # 4. 添加平均值点
    geom_point(data = group_cluster_stats,
               aes(x = Time, y = Mean, fill = Cluster),
               size = 5, color = "black", stroke = 1.2, shape = 21) +
    
    # 5. 添加代谢物数量标签
    geom_text(data = group_cluster_stats %>% 
                filter(Time == max(Time)),
              aes(x = max(Time) + 0.2, y = Mean, 
                  label = paste0("n=", n_genes)),
              hjust = 0, size = 4, fontface = "bold") +
    
    # 6. 设置分面
    facet_wrap(~Cluster, scales = "free_y", ncol = 2) +
    
    # 7. 设置主题和标签
    labs(title = paste(group_name, "Metabolite Expression Patterns by Cluster"),
         subtitle = "Individual metabolites shown with cluster average (thick line) and 95% CI (shaded area)",
         x = "Time Point",
         y = "Expression Level (log2)",
         color = "Metabolite",
         fill = "Cluster") +
    
    # 8. 设置颜色
    scale_fill_manual(values = scales::hue_pal()(length(unique(cluster_stats$Cluster)))) +
    
    # 9. 设置x轴刻度
    scale_x_continuous(breaks = sort(unique(times))) +
    
    # 10. 设置主题
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      strip.text = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "black"),
      plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
    )
  
  # 保存每个组别的聚类图
  ggsave(paste0(group_name, "_cluster_trends.pdf"), p2_group, width = 12, height = 10, dpi = 300)
}


# ============== 第7步：创建每个聚类的独立详细图表 ==============
cluster_plots <- list()

for(cluster_id in unique(gsub("Cluster ", "", cluster_stats$Cluster))) {
  # 获取该聚类的数据
  cluster_data <- metabolite_means %>% 
    filter(grepl(paste0("Cluster ", cluster_id), Cluster))
  
  cluster_avg <- cluster_stats %>% 
    filter(grepl(paste0("Cluster ", cluster_id), Cluster))
  
  # 获取该聚类中的代谢物数量
  n_metabolites <- length(unique(cluster_data$Gene))
  
  # 创建该聚类的图表
  p_cluster <- ggplot() +
    # 背景色带表示置信区间
    geom_ribbon(data = cluster_avg,
                aes(x = Time, ymin = CI_lower, ymax = CI_upper, 
                    fill = Group),
                alpha = 0.15) +
    
    # 个体代谢物线条
    geom_line(data = cluster_data, 
              aes(x = Time, y = Mean, group = interaction(Gene, Group), 
                  color = Gene),
              size = 0.9, alpha = 0.7) +
    
    # 平均趋势线
    geom_line(data = cluster_avg, 
              aes(x = Time, y = Mean, group = Group, color = Group),
              size = 2, alpha = 0.9) +
    
    # 平均值点
    geom_point(data = cluster_avg,
               aes(x = Time, y = Mean, fill = Group),
               size = 4, shape = 21, stroke = 1.2, color = "black") +
    
    # 在最后一个时间点添加代谢物标签
    geom_text_repel(data = cluster_data %>% 
                      filter(Time == max(Time) & Group == "Control"),
                    aes(x = Time, y = Mean, label = Gene, color = Gene),
                    size = 3, fontface = "bold", box.padding = 0.5,
                    direction = "y", hjust = 0, segment.alpha = 0.5) +
    
    # 设置标题和标签
    labs(title = paste0("Cluster ", cluster_id, " (", n_metabolites, " metabolites)"),
         subtitle = "Individual metabolites with group averages",
         x = "Time Point",
         y = "Expression Level (log2)",
         color = "Metabolite",
         fill = "Group") +
    
    # 设置颜色
    scale_color_manual(values = gene_colors[[cluster_id]]) +
    scale_fill_manual(values = c("Control" = "#3366CC", "Mutant" = "#CC3366")) +
    
    # 设置x轴刻度
    scale_x_continuous(breaks = sort(unique(times))) +
    
    # 设置主题
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  # 保存到列表
  cluster_plots[[cluster_id]] <- p_cluster
  
  # 保存单个聚类图
  ggsave(paste0("cluster_", cluster_id, "_detail.pdf"), p_cluster, 
         width = 8, height = 6, dpi = 300)
}
#==================== 7 改 - 只显示趋势线 ===========================
for(cluster_id in unique(gsub("Cluster ", "", cluster_stats$Cluster))) {
  # 获取该聚类的数据 - 使用精确匹配
  cluster_data <- metabolite_means %>% 
    filter(Cluster == paste0("Cluster ", cluster_id))
  
  cluster_avg <- cluster_stats %>% 
    filter(Cluster == paste0("Cluster ", cluster_id))
  
  # 确保cluster_avg没有重复行
  cluster_avg <- cluster_avg %>% distinct(Time, Group, .keep_all = TRUE)
  
  # 获取该聚类中的代谢物数量
  n_metabolites <- length(unique(cluster_data$Gene))
  
  # 创建只有趋势线的对比图
  p_cluster_trends <- ggplot() +
    # 背景色带表示置信区间
    geom_ribbon(data = cluster_avg,
                aes(x = Time, ymin = CI_lower, ymax = CI_upper, fill = Group),
                alpha = 0.2) +
    
    # 平均趋势线
    geom_line(data = cluster_avg, 
              aes(x = Time, y = Mean, group = Group, color = Group),
              size = 2, alpha = 1) +
    
    # 平均值点
    geom_point(data = cluster_avg,
               aes(x = Time, y = Mean, fill = Group),
               size = 4, shape = 21, stroke = 1.2, color = "black") +
    
    # 设置标题和标签
    labs(title = paste0("Cluster ", cluster_id, " - Trend Comparison (", n_metabolites, " metabolites)"),
         subtitle = "Control vs Mutant average trends",
         x = "Time Point",
         y = "Expression Level (log2)",
         color = "Group",
         fill = "Group") +
    
    # 设置颜色
    scale_color_manual(values = c("Control" = "#3366CC", "Mutant" = "#CC3366")) +
    scale_fill_manual(values = c("Control" = "#3366CC", "Mutant" = "#CC3366")) +
    
    # 设置x轴刻度
    scale_x_continuous(breaks = sort(unique(times))) +
    
    # 设置主题
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  # 保存只有趋势线的图
  ggsave(paste0("cluster_", cluster_id, "_trends_only.pdf"), p_cluster_trends, 
         width = 8, height = 6, dpi = 300)
  
  ggsave(paste0("cluster_", cluster_id, "_trends_only.tiff"), p_cluster_trends, 
         width = 8, height = 6, dpi = 300, bg = "white")
  
  # 您仍然可以保留原来的完整图作为参考
  p_cluster_compare <- ggplot() +
    # 个体代谢物线条
    geom_line(data = cluster_data, 
              aes(x = Time, y = Mean, group = interaction(Gene, Group), 
                  color = Gene, linetype = Group),
              size = 0.7, alpha = 0.5) +
    
    # 平均趋势线
    geom_line(data = cluster_avg, 
              aes(x = Time, y = Mean, group = Group, color = Group),
              size = 1.5, alpha = 1) +
    
    # 平均值点
    geom_point(data = cluster_avg,
               aes(x = Time, y = Mean, fill = Group),
               size = 3, shape = 21, stroke = 1, color = "black") +
    
    # 设置标题和标签
    labs(title = paste0("Cluster ", cluster_id, " - Detailed Comparison (", n_metabolites, " metabolites)"),
         subtitle = "Control vs Mutant with individual metabolites",
         x = "Time Point",
         y = "Expression Level (log2)",
         color = "Group",
         linetype = "Group",
         fill = "Group") +
    
    # 设置颜色和线型
    scale_color_manual(values = c("Control" = "#3366CC", "Mutant" = "#CC3366")) +
    scale_fill_manual(values = c("Control" = "#3366CC", "Mutant" = "#CC3366")) +
    scale_linetype_manual(values = c("Control" = "solid", "Mutant" = "dashed")) +
    
    # 设置x轴刻度
    scale_x_continuous(breaks = sort(unique(times))) +
    
    # 设置主题
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  # 保存详细对比图
  ggsave(paste0("cluster_", cluster_id, "_detailed_comparison.tiff"), p_cluster_compare, 
         width = 8, height = 12, dpi = 300, bg = "white")
}
# ============== 第8步：创建组合图表 ==============
if(length(cluster_plots) <= 6) {
  # 使用patchwork组合所有聚类图
  combined_plot <- wrap_plots(cluster_plots, ncol = 2)
  
  # 添加整体标题
  combined_plot <- combined_plot + 
    plot_annotation(
      title = "Detailed Metabolite Expression Patterns by Cluster",
      subtitle = "Individual metabolites shown with group averages and 95% CI",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  # 保存组合图
  ggsave("all_clusters_detailed.pdf", combined_plot, 
         width = 16, height = 12, dpi = 300)
}
# ============== 第9步：创建热图风格的时间序列可视化 ==============
# 准备热图数据
heatmap_data_long <- metabolite_means %>%
  pivot_wider(
    id_cols = c(Gene, Cluster),
    names_from = c(Group, Time),
    values_from = Mean
  ) %>%
  mutate(
    # 计算Control和Mutant之间的差异
    Diff_2 = Mutant_2 - Control_2,
    Diff_3 = Mutant_3 - Control_3,
    Diff_4 = Mutant_4 - Control_4,
    Diff_5 = Mutant_5 - Control_5
  ) %>%
  pivot_longer(
    cols = c(starts_with("Control_"), starts_with("Mutant_"), starts_with("Diff_")),
    names_to = c("Type", "Time"),
    names_pattern = "(.+)_(.+)",
    values_to = "Value"
  ) %>%
  mutate(
    Time = as.numeric(Time),
    Type = factor(Type, levels = c("Control", "Mutant", "Diff"))
  )

# 为每个聚类排序代谢物
ordered_genes <- list()
for(cluster_id in unique(gsub("Cluster ", "", cluster_stats$Cluster))) {
  # 获取该聚类的代谢物
  genes_in_cluster <- names(clusters[clusters == cluster_id])
  
  # 根据差异模式排序
  gene_patterns <- heatmap_data_long %>%
    filter(Gene %in% genes_in_cluster, Type == "Diff") %>%
    pivot_wider(
      id_cols = Gene,
      names_from = Time,
      values_from = Value
    )
  
  # 使用主成分分析排序
  if(length(genes_in_cluster) >= 3) {
    gene_matrix <- as.matrix(gene_patterns[, -1])
    rownames(gene_matrix) <- gene_patterns$Gene
    pca_result <- prcomp(gene_matrix)
    ordered_genes[[cluster_id]] <- rownames(gene_matrix)[order(pca_result$x[,1])]
  } else {
    ordered_genes[[cluster_id]] <- genes_in_cluster
  }
}

# 创建一个因子变量来控制代谢物排序
heatmap_data_long$Gene <- factor(
  heatmap_data_long$Gene,
  levels = unlist(ordered_genes)
)

# 创建热图风格的时间序列可视化
p_heatmap <- ggplot(heatmap_data_long, aes(x = Time, y = Gene, fill = Value)) +
  geom_tile(color = "white", size = 0.2) +
  facet_grid(Cluster ~ Type, scales = "free_y", space = "free_y") +
  scale_fill_viridis_c() +
  scale_x_continuous(breaks = sort(unique(times))) +
  labs(
    title = "Metabolite Expression Patterns Over Time",
    subtitle = "Control vs Mutant with Difference",
    x = "Time Point",
    y = "Metabolite",
    fill = "Expression\n(log2)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "gray95"),
    panel.spacing = unit(0.3, "lines")
  )

# 保存热图
ggsave("metabolite_heatmap_timeseries.pdf", p_heatmap, width = 14, height = 16, dpi = 300)

# ============== 第10步：创建每个聚类的单独热图 ==============
# 为每个聚类创建单独的热图
for(cluster_id in unique(gsub("Cluster ", "", cluster_stats$Cluster))) {
  # 获取该聚类的数据
  cluster_heatmap_data <- heatmap_data_long %>% 
    filter(grepl(paste0("Cluster ", cluster_id), Cluster))
  
  # 获取该聚类中的代谢物数量
  n_metabolites <- length(unique(cluster_heatmap_data$Gene))
  
  # 创建热图
  p_cluster_heatmap <- ggplot(cluster_heatmap_data, aes(x = Time, y = Gene, fill = Value)) +
    geom_tile(color = "white", size = 0.2) +
    facet_grid(. ~ Type, scales = "free_x", space = "free_x") +
    scale_fill_viridis_c(option = "plasma") +
    scale_x_continuous(breaks = sort(unique(times))) +
    labs(
      title = paste0("Cluster ", cluster_id, " - ", n_metabolites, " Metabolites"),
      subtitle = "Expression Patterns Over Time",
      x = "Time Point",
      y = "Metabolite",
      fill = "Expression\n(log2)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      strip.text = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = "gray95")
    )
  
  # 根据代谢物数量调整图表高度
  height_value <- max(6, min(n_metabolites * 0.3, 20))
  
  # 保存热图
  ggsave(paste0("cluster_", cluster_id, "_heatmap.pdf"), p_cluster_heatmap, 
         width = 10, height = height_value, dpi = 300)
}



# ============== 第10步：创建放射状时间序列图 ==============
# 函数：将笛卡尔坐标转换为极坐标
cart2pol <- function(x, y) {
  r <- sqrt(x^2 + y^2)
  theta <- atan2(y, x)
  return(list(r = r, theta = theta))
}

# 函数：将极坐标转换为笛卡尔坐标
pol2cart <- function(r, theta) {
  x <- r * cos(theta)
  y <- r * sin(theta)
  return(list(x = x, y = y))
}

# 为每个聚类创建放射状图
for(cluster_id in unique(gsub("Cluster ", "", cluster_stats$Cluster))) {
  # 获取该聚类的数据
  cluster_data <- metabolite_means %>% 
    filter(grepl(paste0("Cluster ", cluster_id), Cluster))
  
  # 获取该聚类中的代谢物
  metabolites <- unique(cluster_data$Gene)
  n_metabolites <- length(metabolites)
  
  # 创建角度
  angles <- seq(0, 2*pi, length.out = n_metabolites + 1)[-(n_metabolites + 1)]
  
  # 创建放射状数据
  radial_data <- data.frame()
  
  for(i in 1:length(metabolites)) {
    gene_data <- cluster_data %>% filter(Gene == metabolites[i])
    
    # 对每个时间点
    for(t in unique(gene_data$Time)) {
      # 获取Control和Mutant的值
      ctrl_val <- gene_data %>% filter(Time == t, Group == "Control") %>% pull(Mean)
      mut_val <- gene_data %>% filter(Time == t, Group == "Mutant") %>% pull(Mean)
      
      # 如果数据存在
      if(length(ctrl_val) > 0 && length(mut_val) > 0) {
        # 计算极坐标
        r_ctrl <- ctrl_val - min(cluster_data$Mean) + 1  # 确保r为正
        r_mut <- mut_val - min(cluster_data$Mean) + 1
        
        # 转换为笛卡尔坐标
        pos_ctrl <- pol2cart(r_ctrl, angles[i])
        pos_mut <- pol2cart(r_mut, angles[i])
        
        # 添加到数据框
        radial_data <- rbind(radial_data, data.frame(
          Gene = metabolites[i],
          Time = t,
          Group = "Control",
          x = pos_ctrl$x,
          y = pos_ctrl$y,
          value = ctrl_val,
          angle = angles[i],
          stringsAsFactors = FALSE
        ))
        
        radial_data <- rbind(radial_data, data.frame(
          Gene = metabolites[i],
          Time = t,
          Group = "Mutant",
          x = pos_mut$x,
          y = pos_mut$y,
          value = mut_val,
          angle = angles[i],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # 创建放射状图
  p_radial <- ggplot() +
    # 添加同心圆
    geom_path(
      data = data.frame(
        x = cos(seq(0, 2*pi, length.out = 100)) * 1:5,
        y = sin(seq(0, 2*pi, length.out = 100)) * 1:5,
        r = rep(1:5, each = 100)
      ),
      aes(x = x, y = y, group = r),
      color = "gray80", size = 0.3, linetype = "dashed"
    ) +
    
    # 添加放射线
    geom_segment(
      data = data.frame(
        angle = angles,
        Gene = metabolites
      ),
      aes(x = 0, y = 0, 
          xend = 6 * cos(angle), 
          yend = 6 * sin(angle)),
      color = "gray70", size = 0.3
    ) +
    
    # 添加代谢物标签
    geom_text(
      data = data.frame(
        angle = angles,
        Gene = metabolites,
        x = 6.5 * cos(angles),
        y = 6.5 * sin(angles)
      ),
      aes(x = x, y = y, label = Gene, angle = (angle * 180/pi + 90) %% 180 - 90),
      size = 3, hjust = 0.5, 
      color = sapply(metabolites, function(g) gene_colors[[cluster_id]][g])
    ) +
    
    # 添加时间点连线
    geom_polygon(
      data = radial_data %>% 
        group_by(Time, Group) %>% 
        arrange(angle) %>%
        ungroup(),
      aes(x = x, y = y, fill = factor(Time), group = interaction(Time, Group),
          alpha = Group),
      color = "black", size = 0.3
    ) +
    
    # 设置填充颜色和透明度
    scale_fill_viridis_d(option = "plasma") +
    scale_alpha_manual(values = c("Control" = 0.4, "Mutant" = 0.7)) +
    
    # 添加标题
    labs(
      title = paste0("Cluster ", cluster_id, " - Radial Time Series Plot"),
      subtitle = paste0(n_metabolites, " metabolites"),
      fill = "Time Point",
      alpha = "Group"
    ) +
    
    # 设置坐标系和主题
    coord_equal() +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom",
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  
  # 保存放射状图
  ggsave(paste0("cluster_", cluster_id, "_radial.pdf"), p_radial, 
         width = 10, height = 10, dpi = 300)
}

               

# 保存结果表格
results_table <- data.frame(
  Metabolite = sig_genes,
  Cluster = clusters[sig_genes],
  results[sig_genes, ]
)

write.csv(results_table, "significant_metabolites.csv", row.names = FALSE)

# 打印分析完成消息
cat("\n分析完成! 已生成以下文件:\n")
cat("1. enhanced_cluster_trends.pdf - 增强版聚类趋势图\n")
cat("2. all_clusters_detailed.pdf - 所有聚类的详细图\n")
cat("3. metabolite_heatmap_timeseries.pdf - 时间序列热图\n")
cat("4. metabolite_heatmap.pdf - 传统热图\n")
cat("5. significant_metabolites.csv - 显著差异代谢物列表\n")
cat("6. 各聚类的单独详细图和放射状图\n")
               

# ============== 添加部分：创建比较数据并导出为CSV ==============
# 准备比较数据
comparison_data <- metabolite_means %>%
  # 转换为宽格式，每个代谢物在每个时间点有Control和Mutant的值
  pivot_wider(
    id_cols = c(Gene, Cluster, Time),
    names_from = Group,
    values_from = Mean
  ) %>%
  # 计算差异和比率
  mutate(
    Difference = Mutant - Control,
    Log2FoldChange = log2((Mutant + 0.01) / (Control + 0.01)),  # 添加小偏移量避免除以零
    PercentChange = (Mutant - Control) / abs(Control) * 100
  ) %>%
  # 添加聚类信息
  mutate(
    ClusterID = gsub("Cluster ", "", Cluster)
  )

# 添加统计显著性信息（如果需要）
# 这里我们可以从之前的结果中提取调整后的p值
if(exists("results") && "adj.P.Val" %in% colnames(results)) {
  # 为每个代谢物添加调整后的p值
  comparison_data$AdjustedPValue <- NA
  for(gene in unique(comparison_data$Gene)) {
    if(gene %in% rownames(results)) {
      comparison_data$AdjustedPValue[comparison_data$Gene == gene] <- results[gene, "adj.P.Val"]
    }
  }
  
  # 添加显著性标记
  comparison_data <- comparison_data %>%
    mutate(
      Significance = case_when(
        is.na(AdjustedPValue) ~ "Unknown",
        AdjustedPValue < 0.001 ~ "***",
        AdjustedPValue < 0.01 ~ "**",
        AdjustedPValue < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

# 导出比较数据
write.csv(comparison_data, "metabolite_comparison_data.csv", row.names = FALSE)

# 创建聚类摘要数据
cluster_summary <- comparison_data %>%
  group_by(ClusterID, Time) %>%
  summarize(
    MeanControl = mean(Control, na.rm = TRUE),
    MeanMutant = mean(Mutant, na.rm = TRUE),
    MeanDifference = mean(Difference, na.rm = TRUE),
    MeanLog2FC = mean(Log2FoldChange, na.rm = TRUE),
    MeanPercentChange = mean(PercentChange, na.rm = TRUE),
    MedianControl = median(Control, na.rm = TRUE),
    MedianMutant = median(Mutant, na.rm = TRUE),
    MedianDifference = median(Difference, na.rm = TRUE),
    MedianLog2FC = median(Log2FoldChange, na.rm = TRUE),
    MedianPercentChange = median(PercentChange, na.rm = TRUE),
    SDControl = sd(Control, na.rm = TRUE),
    SDMutant = sd(Mutant, na.rm = TRUE),
    SDDifference = sd(Difference, na.rm = TRUE),
    SDLog2FC = sd(Log2FoldChange, na.rm = TRUE),
    SDPercentChange = sd(PercentChange, na.rm = TRUE),
    MetaboliteCount = n_distinct(Gene),
    .groups = "drop"
  )

# 导出聚类摘要
write.csv(cluster_summary, "cluster_comparison_summary.csv", row.names = FALSE)

# 加载必要的包
library(ggplot2)
library(dplyr)
library(readxl)

# 读取数据
data <- read_excel("enrichment.xlsx")

# 筛选p值<0.05的通路
sig_pathways <- data %>%
  filter(p < 0.05) %>%
  # 按cluster和p值排序
  arrange(cluster, p)
# 为每个cluster创建气泡图
plot_bubble_chart <- function(cluster_data, cluster_num) {
  if(nrow(cluster_data) == 0) {
    return(NULL)  # 如果该cluster没有显著通路，返回NULL
  }
  
  # 创建气泡图，使用现成的-log(p)和Impact
  p <- ggplot(cluster_data, aes(x = Impact, y = reorder(`Pathway Name`, `-log(p)`))) +
    geom_point(aes(size = `-log(p)`, color = `-log(p)`)) +
    scale_size_continuous(name = "-log10(p-value)", range = c(3, 10)) +
    scale_color_gradient(name = "-log10(p-value)", low = "blue", high = "red") +
    labs(
      title = paste("Cluster", cluster_num, "Enriched Pathways (p < 0.05)"),
      x = "Pathway Impact",
      y = ""
    ) +
    theme_bw() +
    theme(
      # 将Y轴文本（通路名称）设置为加粗
      axis.text.y = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # 保存图片，分辨率300dpi，使用TIFF格式
  ggsave(
    filename = paste0("Cluster_", cluster_num, "_Pathways.tiff"),
    plot = p,
    width = 10,
    height = 2 + nrow(cluster_data) * 0.3,  # 根据通路数量调整高度
    dpi = 300
  )
  
  return(p)
}

# 获取所有唯一的cluster值
clusters <- unique(sig_pathways$cluster)

# 为每个cluster创建并保存图表
plots <- list()
for(c in clusters) {
  cluster_data <- sig_pathways %>% filter(cluster == c)
  plots[[paste0("cluster_", c)]] <- plot_bubble_chart(cluster_data, c)
}
