colsp <-c('#FED439FF','#709AE1FF','#8A9197FF','#D2AF81FF','#FD7446FF','#D5E4A2FF','#197EC0FF','#F05C3BFF','#46732EFF',
          '#71D0F5FF','#370335FF','#075149FF','#C80813FF','#91331FFF','#1A9993FF','#FD8CC1FF','#FF6700','#9370DB',
          '#F8D568','#00AD43','#89CFF0','#BA160C','#FF91AF','#A6A6A6','#006DB0','#C154C1','#D99A6C','#96C8A2','#FBEC5D')


celltype_levels <- c(
  "Mesophyll-1",
  "Mesophyll-2",
  "Mesophyll-3",
  "Mesophyll-4",
  "Mesophyll-5",
  "Mesophyll-6",
  "Mesophyll-7",
  
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

selected_cells <- which(grepl("Mesophyll", pData(cds)$celltype))


mesophyll <- cds[, selected_cells]

mesophyll <- estimateSizeFactors(mesophyll)
mesophyll <- estimateDispersions(mesophyll)

disp_table <- dispersionTable(mesophyll)
ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mesophyll <- setOrderingFilter(mesophyll, ordering_genes)

plot_pc_variance_explained(mesophyll, return_all = F) 

mesophyll <- reduceDimension(mesophyll, 
                             max_components = 2, 
                             method = 'DDRTree',ndims=25,
                             cores = 40) 

mesophyll <- reduceDimension(mesophyll, max_components = 2, method = 'DDRTree',cores=40)

mesophyll <- orderCells(mesophyll)
mesophyll <- orderCells(mesophyll, root_state = 4)

plot_cell_trajectory(mesophyll, color_by = "orig.ident")+scale_color_manual(values=colsp)+facet_wrap(~orig.ident)+
  ggtitle("Density") + guides(color = guide_legend(override.aes = list(size = 10)))+
  geom_density_2d()+ 
  theme(legend.text = element_text(size = 25),
        plot.title = element_text(size = 25)) 
        
plot_cell_trajectory(mesophyll, color_by = "Pseudotime") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25),
    legend.position = "right"
  ) +
  guides(colour = guide_colorbar(
    barwidth = 2,         
    barheight = 15,      
    ticks = TRUE,          
    nbin = 50,              
    title = "Pseudotime",   
    title.position = "top",
    title.hjust = 0.5,
    label.theme = element_text(size = 12)  
  )) 


plot_cell_trajectory( mesophyll, color_by = "State") +facet_wrap(~orig.ident)+
  ggtitle("State") + guides(color = guide_legend(override.aes = list(size = 10)))+ 
  theme(legend.text = element_text(size = 25), 
        plot.title = element_text(size = 25)) 
plot_cell_trajectory(mesophyll, color_by = 'celltype')+facet_wrap(~orig.ident)+
  ggtitle("Celltype") + guides(color = guide_legend(override.aes = list(size = 5))  
  )+ 
  theme(legend.text = element_text(size = 15), 
        plot.title = element_text(size = 15)) 

        
library(RColorBrewer)
BEAM_res <- BEAM(mesophyll, 
                branch_point = 5, 
                cores = 4,
                progenitor_method = "sequential_split")

library(grid)
library(pheatmap)
new_colors = colorRampPalette(c("#3A539B", "#FFFFFF", "#E67E22"))(62)

new_branch_colors = c("#6A6A6A", "#E67E22", "blue")

heatmap_obj=plot_genes_branched_heatmap(
  mesophyll[row.names(subset(BEAM_res, qval < 1e-4)), ],
  branch_point = 5,
  num_clusters = 4,
  cores = 1,
  branch_labels = c("Spongy mesophyll", "Palisade mesophyll"),
  hmcols = new_colors,
  branch_colors = new_branch_colors,
  use_gene_short_name = T,
  show_rownames = F,
    return_heatmap = T
)

library(grid)
library(pheatmap)

tiff("branched_heatmap.tiff", width = 6, height = 6, units = "in", res = 300)
grid::grid.draw(heatmap_obj$ph_res)
dev.off()

cluster_info <- heatmap_obj$annotation_row$Cluster
head(cluster_info)

genes <- rownames(heatmap_obj$annotation_row)

cluster_genes <- list()
for(i in 1:4) { 
  genes_in_cluster <- genes[cluster_info == i]
  cluster_genes[[paste0("cluster_", i)]] <- genes_in_cluster
}
sapply(cluster_genes, length)

lapply(cluster_genes, head)

for(i in 1:4) {
  write.table(cluster_genes[[paste0("cluster_", i)]], 
              file = paste0("cluster_", i, "_genes.txt"), 
              quote = FALSE, 
              row.names = FALSE, 
              col.names = FALSE)
}

plot_genes_in_pseudotime(mesophyll[c("CYP710A1")], 
                         color_by = "State")
