library(Seurat)
library(WGCNA)
library(tidyverse)
library(ComplexHeatmap)

# 假设您的数据已经是Seurat对象，命名为seurat_obj
# 只选择PAP1-D的细胞
pap1d_cells <- subset(inte_ident, orig.ident == "PAP1-D")

# 准备表达矩阵
expr_matrix <- GetAssayData(pap1d_cells, layer = "data")
expr_matrix <- as.matrix(expr_matrix)

# 选择变异基因
var_genes <- VariableFeatures(pap1d_cells)
expr_matrix <- expr_matrix[var_genes,]

# 转置矩阵以符合WGCNA要求
input_mat <- t(expr_matrix)

# 检查是否有缺失值或0方差
gsg <- goodSamplesGenes(input_mat, verbose = 3)
if (!gsg$allOK) {
    input_mat <- input_mat[gsg$goodSamples, gsg$goodGenes]
}

# 选择软阈值(power)
powers <- c(1:30)
sft <- pickSoftThreshold(input_mat, powerVector = powers)

# 构建共表达网络
adjacency <- adjacency(input_mat, power = sft$powerEstimate)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM

# 识别模块
geneTree <- hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

# 模块颜色分配
moduleColors <- labels2colors(dynamicMods)

# 绘制树状图和模块颜色
plotDendroAndColors(geneTree, moduleColors, "Module colors",
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05)

# 计算模块特征值
MEs <- moduleEigengenes(input_mat, colors = moduleColors)$eigengenes

# 模块相关性热图
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                     marDendro = c(3,3,2,4),
                     marHeatmap = c(3,4,2,2))

# 获取细胞类型注释或cluster信息作为表型数据
traits <- pap1d_cells@meta.data[, c("seurat_clusters", "celltype")] # 可以同时使用多个注释列


# 如果cell_type是字符型,需要转换为数值型矩阵
traits_dummy <- model.matrix(~0 + celltype, data = pap1d_cells@meta.data)


# 检查和对齐样本
commonSamples <- intersect(rownames(MEs), rownames(traits_dummy))
MEs <- MEs[commonSamples, ]
traits_dummy <- traits_dummy[commonSamples, ]

# 计算相关性
moduleTraitCor <- cor(MEs, traits_dummy, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(input_mat))

# 绘制模块-细胞类型相关性热图
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

par(mar = c(6, 8.5, 3, 3))

# 方法1：使用tiff()函数
tiff(filename = "Module_CellType_Heatmap.tiff", 
     width = 14 * 300,    # 宽度(英寸) * DPI
     height = 10 * 300,   # 高度(英寸) * DPI
     units = "px",
     res = 300,
     compression = "lzw")

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traits_dummy),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-cell type relationships"))

dev.off()
