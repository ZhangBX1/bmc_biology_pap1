setwd("/home/bingxu/scrna_paper/aa/")

library(Seurat)
library(patchwork) 
library(dplyr)
library(future)
library(ggplot2)
library(DoubletFinder)

wt=Read10X("/home/bingxu/pap1_scrna/all/result/wt/output/filter_matrix/",gene.column=1)
mutant=Read10X("/home/bingxu/pap1_scrna/all/result/mutant/output/filter_matrix/",gene.column=1)
chl=c("accD", "ArthCp001", "ArthCp002", "ArthCp003", "ArthCp004", "ArthCp005", "ArthCp006", "ArthCp007", "ArthCp008", "ArthCp009", "ArthCp010", "ArthCp011", "ArthCp012", "ArthCp013", "ArthCp014", "ArthCp015", "ArthCp016", "ArthCp017", "ArthCp018", "ArthCp019", "ArthCp020", "ArthCp021", "ArthCp022", "ArthCp023", "ArthCp024", "ArthCp025", "ArthCp026", "ArthCp027", "ArthCp028", "ArthCp029", "ArthCp030", "ArthCp031", "ArthCp032", "ArthCp033", "ArthCp034", "ArthCp035", "ArthCp036", "ArthCp037", "ArthCp038", "ArthCp039", "ArthCp040", "ArthCp041", "ArthCp042", "ArthCp043", "ArthCp044", "ArthCp045", "ArthCp047", "ArthCp048", "ArthCp049", "ArthCp050", "ArthCp051", "ArthCp052", "ArthCp053", "ArthCp054", "ArthCp055", "ArthCp056", "ArthCp057", "ArthCp058", "ArthCp059", "ArthCp060", "ArthCp061", "ArthCp062", "ArthCp063", "ArthCp064", "ArthCp065", "ArthCp066", "ArthCp068", "ArthCp069", "ArthCp070", "ArthCp071", "ArthCp072", "ArthCp073", "ArthCp074", "ArthCp075", "ArthCp076", "ArthCp077", "ArthCp078", "ArthCp079", "ArthCp080", "ArthCp081", "ArthCp083", "ArthCp084", "ArthCp085", "ArthCp086", "Arthcp087", "ArthCp088", "ArthCr088", "ArthCr089", "ArthCr090", "ArthCr091", "ArthCt100", "ArthCt112", "atpA", "atpB", "atpE", "atpF", "atpH", "atpI", "ccsA", "cemA", "clpP", "matK", "ndhA", "ndhB", "ndhC", "ndhD", "ndhE", "ndhF", "ndhG", "ndhH", "ndhI", "ndhJ", "ndhK", "petA", "petB", "petD", "petG", "petL", "petN", "psaA", "psaB", "psaC", "psaI", "psaJ", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbN", "psbT", "psbZ", "rbcL", "rpl14", "rpl16", "rpl2", "rpl20", "rpl22", "rpl23", "rpl32", "rpl33", "rpl36", "rpoA", "rpoB", "rpoC1", "rpoC2", "rps11", "rps12", "rps14", "rps15", "rps16", "rps18", "rps19", "rps2", "rps3", "rps4", "rps7", "rps8", "rrn16S", "rrn23S", "rrn4.5S", "trnA", "trnC", "trnD", "trnE", "trnF", "trnfM", "trnG", "trnH", "trnI", "trnK", "trnL", "trnM", "trnN", "trnP", "trnQ", "trnR", "trnS", "trnT", "trnV", "trnW", "trnY", "ycf1", "ycf2", "ycf3", "ycf4"
)

wt <- wt[!rownames(wt) %in% chl, ]
mutant <- mutant[!rownames(mutant) %in% chl, ]

rownames(wt) <- gsub("_", "-", rownames(wt))  
rownames(mutant) <- gsub("_", "-", rownames(mutant))  

wt_object <- CreateSeuratObject(counts =wt, project = "Col-0", min.cells = 5, min.features = 1000)
mutant_object <- CreateSeuratObject(counts = mutant, project = "PAP1-D", min.cells =5, min.features = 1000)

wt_object <- subset(wt_object , features = setdiff(rownames(wt_object), chl))
mutant_object <- subset(mutant_object, features = setdiff(rownames(mutant_object), chl))
object.list <- list(wt_object, mutant_object)


data <- data.frame(orig.ident = inte$orig.ident, cluster = Idents(inte))

cell_counts <- data %>%
  group_by(orig.ident, cluster) %>%
  summarise(cell_count = n())

p <- ggplot(cell_counts, aes(x = cluster, y = cell_count, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = cell_count), vjust = -0.3, position = position_dodge(0.9), size = 3.5) +
  labs(title = "Cell Counts per Cluster for Each orig.ident", x = "Cluster", y = "Cell Count") +
  theme_minimal()

ggsave(filename = "cells_cluster_number.pdf", plot = p, width = 12, height = 8, dpi = 300)
ggsave(filename = "cells_cluster_number.tiff", plot = p, width = 12, height = 8, dpi = 300)

object.list <- lapply(X = object.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = object.list)

anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)

integrated_object <- IntegrateData(anchorset = anchors)

inte=integrated_object
DefaultAssay(inte) <- "integrated"

inte <- FindVariableFeatures(inte, selection.method = "vst", nfeatures = 1000)
top10 <- head(VariableFeatures(inte), 10)

plot1 <- VariableFeaturePlot(inte)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

inte <- ScaleData(inte, verbose = FALSE)
inte <- RunPCA(inte, npcs = 40)
inte <- FindNeighbors(inte, reduction = "pca", dims = 1:40)

inte$source <- inte$orig.ident

inte2 <- RunPCA(inte, npcs = 50)
inte2 <- FindNeighbors(inte2, reduction = "pca", dims = 1:50)
a=ElbowPlot(inte2, ndims = 50)
ggsave(filename = "PCs.pdf", plot = a, width = 8, height = 8, dpi = 300)
ggsave(filename = "PCs.tiff", plot = a, width = 8, height = 8, dpi = 300)

inte2 <- RunPCA(inte, npcs = 40)
inte2 <- FindNeighbors(inte2, reduction = "pca", dims = 1:40)
resolutions <- seq(0.1, 1.2, by = 0.1)
for(res in resolutions) {
    inte <- FindClusters(inte2, resolution = res, 
                        algorithm = 1, 
                        random.seed = 42,
                        verbose = FALSE)
}
library(clustree)
reso=clustree(inte@meta.data, prefix = "integrated_snn_res.")

ggsave(filename = "resolution.pdf", plot = reso, width = 10, height = 10, dpi = 300)
ggsave(filename = "resolution.tiff", plot = reso, width = 10, height = 10, dpi = 300)

inte <- RunUMAP(inte, reduction = "pca", dims = 1:40)
inte <- RunTSNE(inte, reduction = "pca", dims = 1:40)
inte <- FindClusters(inte, resolution = 0.7)
DimPlot(inte, 
        reduction = "umap", 
        split.by = "orig.ident", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("UMAP") +theme_minimal()
umap=DimPlot(inte, 
        reduction = "umap", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("UMAP")
ggsave(filename = "umap.pdf", plot = umap, width = 8, height = 8, dpi = 300)
ggsave(filename = "umap.tiff", plot = umap, width = 8, height = 8, dpi = 300)
umap_fen=DimPlot(inte, 
        reduction = "umap", 
        split.by = "orig.ident", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("UMAP")
ggsave(filename = "umap_fen.pdf", plot = umap_fen, width = 14, height = 8, dpi = 300)
ggsave(filename = "umap_fen.tiff", plot = umap_fen, width = 14, height = 8, dpi = 300)

DimPlot(inte, 
        reduction = "tsne", 
        split.by = "orig.ident", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("TSNE") +theme_minimal()
tsne=DimPlot(inte, 
        reduction = "tsne", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("TSNE")
ggsave(filename = "tsne.pdf", plot = tsne, width = 8, height = 8, dpi = 300)
ggsave(filename = "tsne.tiff", plot = tsne, width = 8, height = 8, dpi = 300)
tsne_fen=DimPlot(inte, 
        reduction = "tsne", 
        split.by = "orig.ident", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("TSNE")
ggsave(filename = "tsne_fen.pdf", plot = tsne_fen, width = 14, height = 8, dpi = 300)
ggsave(filename = "tsne_fen.tiff", plot = tsne_fen, width = 14, height = 8, dpi = 300)


library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)

df <- data.frame(
  ident = inte$orig.ident,
  nCount_RNA = inte$nCount_RNA,
  nFeature_RNA = inte$nFeature_RNA
)

df_long <- melt(df, id.vars = "ident")
df_long$type <- ifelse(df_long$variable == "nCount_RNA", "UMI counts per cell", "Genes per cell")
df_long$variable <- NULL

p1 <- ggplot(subset(df_long, type == "UMI counts per cell"), aes(x = ident, y = value, fill = ident)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_y_continuous(trans = 'log10') +
  theme_classic() +
  labs(title = "Distribution of UMI counts per cell",
       y = "UMI Counts (log10 scale)",
       x = "") +
  theme(legend.position = "none")+
theme(
    plot.title = element_text(face = "bold", size = 16) 
  )

p2 <- ggplot(subset(df_long, type == "Genes per cell"), aes(x = ident, y = value, fill = ident)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() +
  labs(title = "Distribution of Genes per cell",
       y = "Genes Counts",
       x = "") +
  theme(legend.position = "none")+
theme(
    plot.title = element_text(face = "bold", size = 16) 
  )

cell_count_df <- df %>%
  count(ident) %>%
  rename(cell_count = n)

p3 <- ggplot(cell_count_df, aes(x = ident, y = cell_count, fill = ident)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  theme_classic() +
  labs(title = "Cell Counts per Group",
       y = "Cell Counts",
       x = "") +
  theme(legend.position = "none")+
theme(
    plot.title = element_text(face = "bold", size = 16)
  )

ggsave(filename = "UMI.pdf", plot = p1, width = 8, height = 8, dpi = 300)
ggsave(filename = "UMI.tiff", plot = p1, width = 8, height = 8, dpi = 300)
ggsave(filename = "UMI.tiff", plot = p1, width = 8, height = 8, dpi = 300)
ggsave(filename = "genes.pdf", plot = p2, width = 8, height = 8, dpi = 300)
ggsave(filename = "genes.tiff", plot = p2, width = 8, height = 8, dpi = 300)
ggsave(filename = "genes.tiff", plot = p2, width = 8, height = 8, dpi = 300)
ggsave(filename = "cells.pdf", plot = p3, width = 8, height = 8, dpi = 300)
ggsave(filename = "cells.tiff", plot = p3, width = 8, height = 8, dpi = 300)
ggsave(filename = "cells.tiff", plot = p3, width = 8, height = 8, dpi = 300)


cell_metadata <- inte@meta.data

cell_counts <- cell_metadata %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarize(cell_count = n(), .groups = "drop")

p <- ggplot(cell_counts, aes(x = factor(seurat_clusters), y = cell_count, fill = orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = cell_count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 3) +
  labs(title = "Number of Cells per Cluster by orig.ident",
       x = "Cluster",
       y = "Number of Cells") +
  theme_minimal()

ggsave(filename = "cluster_cell_nu.pdf", plot = p, width = 14, height = 8, dpi = 300)
ggsave(filename = "cluster_cell_nu.tiff", plot = p, width = 14, height = 8, dpi = 300)

metadata <- inte@meta.data

average_umi_per_cell <- aggregate(metadata$nCount_RNA, by=list(metadata$orig.ident), FUN=mean)
colnames(average_umi_per_cell) <- c("orig.ident", "average_UMI_per_cell")

cell_counts <- table(metadata$orig.ident)
cell_counts <- as.data.frame(cell_counts)
colnames(cell_counts) <- c("orig.ident", "cell_count")

average_genes_per_cell <- aggregate(metadata$nFeature_RNA, by=list(metadata$orig.ident), FUN=mean)
colnames(average_genes_per_cell) <- c("orig.ident", "average_genes_per_cell")

print("每细胞平均UMI数:")
print(average_umi_per_cell)

print("\n每个orig.ident的细胞总数:")
print(cell_counts)

print("\n每细胞平均基因数:")
print(average_genes_per_cell)

all_clusters_markers <- FindAllMarkers(inte)
write.csv(all_clusters_markers, file = "markers.csv")

DefaultAssay(inte) <- "RNA"

inte=JoinLayers(inte)

library(dplyr)
all_clusters <- unique(Idents(inte))
diff_genes_list <- list()

clusters <- unique(Idents(inte))

for(cluster in clusters) {
    inte_subset <- subset(inte, idents = cluster)
    diff_genes <- FindMarkers(inte_subset, 
                            ident.1 = "PAP1-D", 
                            ident.2 = "Col-0",
                            group.by = "orig.ident",
                            min.pct = 0.25,
                            logfc.threshold = 0.25)

    diff_genes$cluster <- cluster
    diff_genes$gene <- rownames(diff_genes)

    up_genes <- diff_genes[diff_genes$p_val_adj < 0.05 & diff_genes$avg_log2FC > 0, ]

    diff_genes_list[[paste0("cluster_", cluster)]] <- up_genes
}

all_diff_genes <- do.call(rbind, diff_genes_list)

all_diff_genes <- all_diff_genes[, c("gene", "cluster", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]

write.csv(all_diff_genes, 
          file = "all_clusters_PAP1D_vs_Col0_up.csv", 
          row.names = FALSE)

cluster_stats <- table(all_diff_genes$cluster)
print("Each cluster up-regulated genes count:")
print(cluster_stats)

VlnPlot(inte, features = c("rna_AAP2"), split.by = "orig.ident", split.plot = TRUE,log=TRUE)

FeaturePlot(inte, features = c("rna_CHS"),split.by = "orig.ident", cols = c("grey", "red"))

VlnPlot(inte, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

celltype_detailed <- c("0"="Phloem parenchyma-1",
                       "1"="Mesophyll-1",
                       "2"="Phloem parenchyma-2",
                       "3"="Mesophyll-2",
                       "4"="Mesophyll-3",
                       "5"="Phloem parenchyma-3",
                       "6"="Mesophyll-4",
                       "7"="Mesophyll-5",
                       "8"="Epidermis-1",
                       "9"="Vasculature-1",
                       "10"="Mesophyll-6",
                       "11"="Xylem",
                       "12"="Phloem",
                       "13"="Vasculature-2",
                       "14"="Mesophyll-7",
                       "15"="Epidermis-2",
                       "16"="Companion cell",
                       "17"="Guard cell"
                       
)

celltype_levels <- c(
  "Mesophyll-1",
  "Mesophyll-2",
  "Mesophyll-3",
  "Mesophyll-4",
  "Mesophyll-5",
  "Mesophyll-6",
  "Mesophyll-7",
  "Mesophyll-8",
  "Mesophyll-9",
  
  "Phloem parenchyma-1",
  "Phloem parenchyma-2",
  "Phloem parenchyma-3",
  "Vasculature-1",
  "Vasculature-2",
  "Xylem",
  "Phloem",
  
  "Epidermis-1",
    "Epidermis-2" ,
  "Companion cell",
  "Guard cell"
)

inte_ident <- RenameIdents(inte, celltype_detailed)
Idents(inte_ident) <- factor(Idents(inte_ident), levels = celltype_levels)
inte_ident@meta.data$celltype <- Idents(inte_ident)



color_scheme <- c(
    "Mesophyll-1" = "#9CCC65",
    "Mesophyll-2" = "#8BC34A",
    "Mesophyll-3" = "#7CB342",
    "Mesophyll-4" = "#689F38",
    "Mesophyll-5" = "#558B2F",
    "Mesophyll-6" = "#33691E",
    "Mesophyll-7" = "#2E7D32",
    "Mesophyll-8" = "#1B5E20",
    "Mesophyll-9" = "#004D40",
    
    "Phloem parenchyma-1" = "#5C6BC0",
    "Phloem parenchyma-2" = "#3949AB",
    "Vasculature-1" = "#303F9F",
    "Vasculature-2" = "#1A237E",
    "Xylem"="#4387B5",
"Phloem"="#84CAC0",

    "Epidermis-1" = "#FFB74D",
    "Epidermis-2" = "#FF9800",
    "Companion cell" = "#F57C00",
    "Guard cell" = "#E65100",
    
    "Phloem parenchyma-3" = "#72ccff"
)

umap <- DimPlot(inte_ident, 
    reduction = "umap", 
    repel = TRUE,
    label = TRUE,
    label.size = 4,
    label.color = "black") +
    scale_color_manual(values = color_scheme) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))

umap_split <- DimPlot(inte_ident, 
    reduction = "umap", 
    split.by = "orig.ident",
    repel = TRUE,
    label = TRUE,
    label.size = 4,
    label.color = "black") +
    scale_color_manual(values = color_scheme) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(filename = "ident_fen.tiff", plot = umap_split , width = 14, height = 8, dpi = 300)
ggsave(filename = "ident.tiff", plot = umap , width = 10, height = 8, dpi = 300)

library(ggExtra)

inte$logUMI <- log10(inte$nCount_RNA)

umap_coords <- Embeddings(inte, "umap")

umap_df <- data.frame(UMAP_1 = umap_coords[,1], UMAP_2 = umap_coords[,2], logUMI = inte$logUMI, orig.ident = inte$orig.ident)

umap_df_col0 <- umap_df[umap_df$orig.ident == "Col-0", ]
umap_df_pap1d <- umap_df[umap_df$orig.ident == "PAP1-D", ]

min_umi_col0 <- umap_df_col0[which.min(umap_df_col0$logUMI), ]
min_umi_pap1d <- umap_df_pap1d[which.min(umap_df_pap1d$logUMI), ]

umap_plot_col0 <- ggplot(umap_df_col0, aes(x = UMAP_1, y = UMAP_2, color = logUMI)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "UMAP with log10(UMI) Density (Col-0)", x = "UMAP 1", y = "UMAP 2", color = "log10(UMI)") +
  geom_point(data = min_umi_col0, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 2) + 
  geom_text(data = min_umi_col0, aes(x = UMAP_1, y = UMAP_2, label = paste("Col-0 Min UMI:", round(logUMI, 2))), vjust = -1, color = "red", size = 7)+
 theme(
    plot.title = element_text(face = "bold", size = 16)
  )

umap_plot_pap1d <- ggplot(umap_df_pap1d, aes(x = UMAP_1, y = UMAP_2, color = logUMI)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "UMAP with log10(UMI) Density (PAP1-D)", x = "UMAP 1", y = "UMAP 2", color = "log10(UMI)") +
  geom_point(data = min_umi_pap1d, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 2) + 
  geom_text(data = min_umi_pap1d, aes(x = UMAP_1, y = UMAP_2, label = paste("PAP1-D Min UMI:", round(logUMI, 2))), vjust = -1, color = "red", size = 7)+
 theme(
    plot.title = element_text(face = "bold", size = 16) 
  )

ggsave(filename = "Col0_UMI.tiff", plot = umap_plot_col0 , width = 10, height = 8, dpi = 300)
ggsave(filename = "PAP1-D_UMI.tiff", plot = umap_plot_pap1d , width = 10, height = 8, dpi = 300)

p <- DotPlot(inte,
    features = genes_to_plot,
    cols = c("yellow", "purple"),
    dot.scale = 8,
    cluster.idents = TRUE) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")
    ) +
    labs(x = "Genes", y = "Clusters") +
    scale_y_discrete(limits = rev(cluster_order)) +
    coord_cartesian(ylim = c(0.5, length(unique(Idents(inte))) + 1))

cell_type_changes <- which(diff(as.integer(genes_df$cell_type)) != 0)
for(pos in cell_type_changes) {
    p <- p + geom_vline(xintercept = pos + 0.5, linetype = "dashed", color = "grey70", alpha = 0.5)
}

cell_type_labels <- tapply(1:nrow(genes_df), genes_df$cell_type, mean)
p <- p + annotate("text",
    x = unname(cell_type_labels),
    y = length(unique(Idents(inte))) + 1,
    label = names(cell_type_labels),
    size = 5, 
    fontface = "bold" 
)

p <- p + annotate("text",
    x = 1:nrow(genes_df),
    y = length(unique(Idents(inte))) + 0.5,
    label = paste("C", genes_df$cluster, sep=""),
    size = 5, 
    fontface = "bold" 
)

ggsave("markers_points.pdf", p, width = 24, height = 9, dpi = 300)
ggsave("markers_points.tiff", p, width = 24, height = 9, dpi = 300)


volcano_plots <- list()

clusters <- unique(inte@meta.data$seurat_clusters)

for(cluster in clusters) {
   
    cells_in_cluster <- WhichCells(inte, expression = seurat_clusters == cluster)
    cluster_obj <- subset(inte, cells = cells_in_cluster)
    
    markers <- FindMarkers(cluster_obj,
                          ident.1 = "PAP1-D",
                          ident.2 = "Col-0",
                          group.by = "orig.ident",
                          min.pct = 0.1,
                          logfc.threshold = 0)
    
    markers$gene <- rownames(markers)
    markers$significant <- ifelse(markers$p_val_adj < 0.05,
                                ifelse(markers$avg_log2FC > 0.25, "Up",
                                       ifelse(markers$avg_log2FC < -0.25, "Down", "NS")),
                                "NS")
    
    top_genes <- markers %>% 
        filter(p_val_adj < 0.00001 & abs(avg_log2FC) > 1) %>%
        top_n(10, wt = abs(avg_log2FC))
    
    p <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
        geom_point(size = 1, alpha = 0.7) +
        scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey")) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
        geom_text_repel(data = top_genes, 
                       aes(label = gene),
                       size = 3,
                       box.padding = 0.5,
                       max.overlaps = 20) +
        labs(title = paste("Cluster", cluster),
             x = "log2(Fold Change)",
             y = "-log10(adjusted p-value)") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "right")
    
    volcano_plots[[paste0("cluster_", cluster)]] <- p
}

n_plots <- length(clusters)
n_cols <- 3 
n_rows <- ceiling(n_plots/n_cols)

pdf("all_clusters_volcano_plots.pdf", width = 15, height = 5*n_rows)
plot_grid(plotlist = volcano_plots, ncol = n_cols)
dev.off()
