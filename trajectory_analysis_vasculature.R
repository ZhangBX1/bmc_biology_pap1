colsp <-c('#FED439FF','#709AE1FF','#8A9197FF','#D2AF81FF','#FD7446FF','#D5E4A2FF','#197EC0FF','#F05C3BFF','#46732EFF',
          '#71D0F5FF','#370335FF','#075149FF','#C80813FF','#91331FFF','#1A9993FF','#FD8CC1FF','#FF6700','#9370DB',
          '#F8D568','#00AD43','#89CFF0','#BA160C','#FF91AF','#A6A6A6','#006DB0','#C154C1','#D99A6C','#96C8A2','#FBEC5D')
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
  
  "mesophyll-1",
  "mesophyll-2",
  "Companion cell",
  "Guard cell"
)

inte_ident <- RenameIdents(inte, celltype_detailed)
Idents(inte_ident) <- factor(Idents(inte_ident), levels = celltype_levels)
inte_ident@meta.data$celltype <- Idents(inte_ident)

expr_matrix <- GetAssayData(inte_ident, assay = "RNA", layer = "counts")
pd <- new('AnnotatedDataFrame', data = inte_ident@meta.data)
fData <- data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

selected_cells <- which(grepl("Phloem parenchyma", pData(cds)$celltype) & pData(cds)$orig.ident == "PAP1-D"|
                          grepl("Vasculature", pData(cds)$celltype) & pData(cds)$orig.ident == "PAP1-D"|
                          grepl("Xylem", pData(cds)$celltype) & pData(cds)$orig.ident == "PAP1-D"|
                          grepl("Phloem", pData(cds)$celltype) & pData(cds)$orig.ident == "PAP1-D")
vasculature <- cds[, selected_cells]

vasculature <- estimateSizeFactors(vasculature)
vasculature <- estimateDispersions(vasculature)

disp_table <- dispersionTable(vasculature)
ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
vasculature <- setOrderingFilter(vasculature, ordering_genes)

plot_pc_variance_explained(vasculature, return_all = F) 

vasculature <- reduceDimension(vasculature, max_components = 2,num_dim = 25, method = 'DDRTree',cores=40)

vasculature <- orderCells(vasculature)
vasculature <- orderCells(vasculature, root_state = 8)

diff_test_res <- differentialGeneTest(vasculature,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime) * orig.ident")

# 7. 绘制轨迹图
p<-plot_cell_trajectory( vasculature, color_by = "State") +facet_wrap(~orig.ident)+
  ggtitle("State") + guides(color = guide_legend(override.aes = list(size = 10)))+  # 调整点的大小
  theme(legend.text = element_text(size = 25),  # 调整图例字体大小
        plot.title = element_text(size = 25)) 
ggsave(p,file = "state_trajectory.tiff", width = 14, height =10,limitsize = FALSE)

# 7. 绘制轨迹图
p=plot_cell_trajectory(vasculature, color_by = "orig.ident")+scale_color_manual(values=colsp)+facet_wrap(~orig.ident)+
  ggtitle("Density") + guides(color = guide_legend(override.aes = list(size = 10)))+
  geom_density_2d()+  # 调整点的大小
  theme(legend.text = element_text(size = 25),  # 调整图例字体大小
        plot.title = element_text(size = 25)) 
ggsave(p,file = "number_trajectory.tiff", width = 14, height =10,limitsize = FALSE)

# 可以根据需要添加其他可视化，例如：
plot_cell_trajectory(vasculature, color_by = "Pseudotime")


p1 <- plot_cell_trajectory(vasculature, color_by = "Pseudotime") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25),
    legend.position = "right"
  ) +
  guides(colour = guide_colorbar(
    barwidth = 2,           # 色条宽度
    barheight = 15,         # 色条高度
    ticks = TRUE,           # 显示刻度
    nbin = 50,              # 颜色分段数
    title = "Pseudotime",   # 标题
    title.position = "top",
    title.hjust = 0.5,
    label.theme = element_text(size = 12)  # 标签文字大小
  )) 
p1
ggsave(p1,file = "Pseudotime_trajectory.tiff",width = 8, height = 6,limitsize = FALSE)


p3 <-plot_cell_trajectory( vasculature, color_by = 'celltype')+facet_wrap(~orig.ident)+
  ggtitle("Celltype") + guides(color = guide_legend(override.aes = list(size = 10))  # 调整legend文字大小
  )+  # 调整点的大小
  theme(legend.text = element_text(size = 25),  # 调整图例字体大小
        plot.title = element_text(size = 25)) 
p3
ggsave(p3 ,file =("Celltype_trajectory.tiff"), width = 14, height =10,limitsize = FALSE)

p=plot_cell_trajectory(vasculature, color_by = "Pseudotime") +
  geom_density_2d()+ guides(color = guide_legend(override.aes = list(size = 10))  # 调整legend文字大小
  )+  # 调整点的大小
  theme(legend.text = element_text(size = 25),  # 调整图例字体大小
        plot.title = element_text(size = 25)) +
  ggtitle("Pseudotime") 
ggsave(p ,file =("number_trajectory.tiff"), width = 14, height =10,limitsize = FALSE)

BEAM_res <- BEAM(vasculature, 
                branch_point = 2, 
                cores = 4,
                progenitor_method = "sequential_split")

library(Cairo)

new_colors = colorRampPalette(c("#3A539B", "#FFFFFF", "#E67E22"))(62)

new_branch_colors = c("#6A6A6A", "#E67E22", "#3A539B")  # 灰色、橙色、蓝色

p=plot_genes_branched_heatmap(
  vasculature[row.names(subset(BEAM_res, qval < 1e-4)), ],
  branch_point = 2,
  num_clusters = 4,
  cores = 1,
  branch_labels = c("x", "p"),
  hmcols = new_colors,
  branch_colors = new_branch_colors,
  use_gene_short_name = T,
  show_rownames = F,
    return_heatmap = FALSE
)

a <- "CHI"

c <- paste(a, ".tiff", sep = "")
gly <- plot_genes_in_pseudotime(vasculature[a], color_by = "State")
gly <- gly + ggtitle(a) + theme(plot.title = element_text(size = 20, face = "bold"))
ggsave(c, gly, width = 10, height = 10, dpi = 300)
