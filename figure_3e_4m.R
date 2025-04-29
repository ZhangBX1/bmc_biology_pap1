setwd("F://paper//single_cell//aa//test//trac//mesophyll//pap1//point2")
setwd("E://bud//DEG")
# 加载必要的包
library(org.At.tair.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# 读取基因列表文件
genes <- read.table("cluster_1_genes.txt", header = FALSE, stringsAsFactors = FALSE)
genes <- genes$V1  # 转换为向量

genes <- read.table("1.txt", header = FALSE, stringsAsFactors = FALSE)
genes <- genes$V1  # 转换为向量
# 检查基因列表
head(genes)
length(genes)  # 查看基因数量

# 准备背景基因集（可选，如果不提供则使用org.At.tair.db中的所有基因）
# bg_genes <- keys(org.At.tair.db, keytype="TAIR")

# GO富集分析
go_enrich <- enrichGO(
  gene = genes,
  OrgDb = org.At.tair.db,
  keyType = "SYMBOL",  # 指定输入的ID类型为TAIR ID
  ont = "ALL",       # 分析所有GO类别: BP, CC, MF
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 查看GO富集结果
head(go_enrich,10)

# 可视化GO富集结果
# 条形图
barplot(go_enrich, showCategory=15, font.size=10)

# 点图
dotplot(go_enrich, showCategory=30)

# 基因-功能网络图
cnetplot(go_enrich, categorySize="pvalue", 
         foldChange=NULL, 
         showCategory=10)


cnetplot(go_enrich, 
         categorySize = "pvalue",
         showCategory = 30,
         node_label_size = 5,      # 增大标签字体
         cex_category = 1.2,       # 增大通路节点大小
         cex_gene = 0.8,           # 调整基因节点大小
         
         color_category = "#8A2BE2", # 紫色用于通路节点(BlueViolet)
         use_ggrepel = TRUE,
         color_gene = "#4682B4"      # 蓝色用于基因节点(SteelBlue)
)


cnetplot(go_enrich,
         categorySize = "pvalue",
         showCategory = 20,
         node_label_size = 5,
         cex_category = 1.2,
         cex_gene = 0.8,
         color_category = "#8A2BE2",
         use_ggrepel = TRUE,
         color_gene = "#4682B4",
         layout = "kk"           # 使用 Kamada-Kawai 布局
)




cnetplot(go_enrich,
         categorySize = "pvalue",
         showCategory = 25,
         node_label_size = 5,
         cex_category = 1.2,
         cex_gene = 0.8,
         color_category = "#8A2BE2",
         use_ggrepel = TRUE,
         color_gene = "#4682B4",
         layout = "kk"           # 使用 Kamada-Kawai 布局
)


edox2 <- pairwise_termsim(go_enrich)
p2 <- treeplot(edox2,hang=-1,cluster.params= list(method = "complete"),offset.params = list(bar_tree = 6))
p2=p2+theme(plot.title.position = "plot",plot.title = element_text(size = 20,hjust = 0.5))
p2
