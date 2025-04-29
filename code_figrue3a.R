# 加载必要的包
library(ggplot2)
library(dplyr)
library(reshape2)
library(scales)
library(stringr)

# 读取数据
df <- read.csv("F://paper//single_cell//aa//df//data.csv")

# 计算 -log10(p_val_adj)，用于圆圈大小
df$neg_log10_padj <- -log10(df$p_val_adj)
# 将 p_val_adj = 1 的值设置为一个很小的值，避免 -log10(1) = 0
df$neg_log10_padj[df$p_val_adj == 1] <- 0
# 限制最大值，避免极端值
df$neg_log10_padj <- pmin(df$neg_log10_padj, 10)

# 确保所有聚类都存在（0-17）
all_clusters <- 0:17
all_genes <- unique(df$gene)

# 创建完整的数据框，包括缺失的聚类
complete_df <- expand.grid(gene = all_genes, cluster = all_clusters)
df_full <- merge(complete_df, df, by = c("gene", "cluster"), all.x = TRUE)

# 填充缺失值
df_full$avg_log2FC[is.na(df_full$avg_log2FC)] <- 0
df_full$neg_log10_padj[is.na(df_full$neg_log10_padj)] <- 0

# 自定义基因排序
# 1. 提取符合特定模式的基因
at1g56650_genes <- all_genes[grep("AT1G56650", all_genes)]
pal_genes <- all_genes[grep("PAL", all_genes)]
c4h_genes <- all_genes[grep("C4H", all_genes)]
cl4_genes <- all_genes[grep("4CL", all_genes)]
chi_genes <- all_genes[grep("CHI", all_genes)]
adt_genes <- all_genes[grep("ADT", all_genes)]
ugt_genes <- all_genes[grep("UGT", all_genes)]

# 2. 找出其他基因
specified_genes <- c(at1g56650_genes, adt_genes, pal_genes, c4h_genes, cl4_genes, chi_genes, ugt_genes)
other_genes <- setdiff(all_genes, specified_genes)

# 3. 按指定顺序组合所有基因（反序）
ordered_genes <- c(other_genes, ugt_genes, chi_genes, cl4_genes, c4h_genes, pal_genes, adt_genes, at1g56650_genes)

# 将基因转换为因子并设置顺序
df_full$gene <- factor(df_full$gene, levels = ordered_genes)

# 创建热图
p <- ggplot(df_full, aes(x = factor(cluster), y = gene)) +
  # 添加方块，颜色代表 log2FC
  geom_tile(aes(fill = avg_log2FC), color = "black", linewidth = 0.5) +
  # 添加圆圈，大小代表 -log10(adj.p.value)，调整大小范围为更小的值
  geom_point(aes(size = neg_log10_padj), color = "black") +
  # 设置颜色方案
  scale_fill_gradient2(
    name = "log2FC",
    low = "#3b4cc0", 
    mid = "white", 
    high = "#a63603", 
    midpoint = 0,
    limits = c(-max(abs(df_full$avg_log2FC)), max(abs(df_full$avg_log2FC)))
  ) +
  # 设置圆圈大小 - 调小圆圈的范围
  scale_size_continuous(
    name = "-log10(adj.p.value)", 
    range = c(0, 3),  # 将最大值设为3，使圆圈适中
    breaks = c(0, 2, 4, 6, 8, 10),
    labels = c("0", "2", "4", "6", "8", "10+")
  ) +
  # 设置主题 - 加粗放大字体
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    plot.background = element_rect(fill = "white", color = NA),  # 添加白色背景
    panel.background = element_rect(fill = "white", color = NA)  # 添加白色面板背景
  ) +
  # 设置标签
  labs(
    x = "Cluster", 
    y = "Gene",
    title = "Gene Expression Differences (PAP1-D vs Col0) Across Clusters",
    subtitle = "Color: log2FC | Circle size: -log10(adj.p.value)"
  )
p

# 保存图片 - 增加分辨率
ggsave("gene_expression_heatmap_new_order.tiff", p, width = 10, height = 15, dpi = 400)

# 显示图形
print(p)
