library(Seurat)
library(patchwork) 
library(dplyr)
library(future)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(monocle)
library(ggplot2)
library(dplyr)
library(ggpubr)

load("single_cell.RData")

plot_gene_expression_comparison <- function(
  seurat_obj,               
  gene,                     
  group1 = "PAP1-D",      
  group2 = "Col-0",      
  group_col = "orig.ident",
  color1 = "purple",   
  color2 = "#4DBBD5",  
  save_plot = TRUE,     
  output_dir = ".",    
  width = 7,              
  height = 8,    
  dpi = 300,           
  layer = "data",   
  adj_p_threshold = 0.0001,  
  log2FC_threshold = 0.25  
) {
  
  library(ggplot2)
  library(dplyr)
  
  if(!gene %in% rownames(seurat_obj)) {
    warning(paste0("基因 '", gene, "' 在数据集中不存在"))
    return(NULL)
  }

  expr_data <- LayerData(seurat_obj, layer = layer)
  metadata <- seurat_obj[[]]
  
  cells_group1 <- rownames(metadata)[metadata[[group_col]] == group1]
  cells_group2 <- rownames(metadata)[metadata[[group_col]] == group2]
  
  mean_expr1 <- mean(expr_data[gene, cells_group1], na.rm = TRUE)
  mean_expr2 <- mean(expr_data[gene, cells_group2], na.rm = TRUE)
  
  pseudocount <- 1e-5
  log2FC <- log2((mean_expr1 + pseudocount) / (mean_expr2 + pseudocount))
  
  seurat_markers <- FindMarkers(seurat_obj, 
                               ident.1 = group1, 
                               ident.2 = group2, 
                               group.by = group_col,
                               features = gene,
                               test.use = "wilcox",
                               min.pct = 0)
  
  adj_p_val <- seurat_markers[1, "p_val_adj"]
  
  gene_expr <- data.frame(
    Expression = as.vector(expr_data[gene, ]),
    Group = as.character(metadata[[group_col]])
  )
  
  gene_expr <- gene_expr[gene_expr$Group %in% c(group1, group2), ]
  gene_expr$Group <- factor(gene_expr$Group, levels = c(group1, group2))
  
  ci_stats <- gene_expr %>%
    group_by(Group) %>%
    summarise(
      mean = mean(Expression, na.rm = TRUE),
      sd = sd(Expression, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      ci_lower = mean - qt(0.975, n-1) * se,
      ci_upper = mean + qt(0.975, n-1) * se
    )
  
  y_max <- max(gene_expr$Expression, na.rm = TRUE)
  y_min <- min(gene_expr$Expression, na.rm = TRUE)
  y_range <- max(0.1, y_max - y_min)
  
  is_significant <- adj_p_val < adj_p_threshold && abs(log2FC) > log2FC_threshold
  
  scatter_plot <- ggplot() +
    geom_jitter(data = gene_expr, 
               aes(x = Group, y = Expression, color = Group),
               position = position_jitter(width = 0.3, seed = 123),
               alpha = 0.5, size = 1.2) +
    
    geom_point(data = ci_stats,
              aes(x = Group, y = mean),
              size = 4, shape = 18, color = "black") +
    
    geom_errorbar(data = ci_stats,
                 aes(x = Group, ymin = ci_lower, ymax = ci_upper),
                 width = 0.2, size = 1, color = "black") +
    
    {if (is_significant)
      annotate("segment", x = 1, xend = 2, 
              y = y_max + y_range * 0.15, 
              yend = y_max + y_range * 0.15,
              size = 0.5) } +
    
    {if (is_significant && adj_p_val < 0.00000000000000000001)
      annotate("text", x = 1.5, y = y_max + y_range * 0.2, 
              label = "***", size = 6)
     else if (is_significant && adj_p_val < 0.000000000001)
      annotate("text", x = 1.5, y = y_max + y_range * 0.2, 
              label = "**", size = 6)
     else if (is_significant && adj_p_val < 0.0001)
      annotate("text", x = 1.5, y = y_max + y_range * 0.2, 
              label = "*", size = 6)
    } +
    
    scale_color_manual(values = setNames(c(color1, color2), c(group1, group2))) +
    
    ylim(NA, y_max + y_range * 0.3) +
    
    labs(title = paste0(gene, " Expression"),
        subtitle = paste0(group1, " vs ", group2, ", log2FC = ", round(log2FC, 2)),
        x = "",
        y = "Expression Level",
        color = "Genotype") +
    
    theme_classic() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      panel.grid.major.y = element_line(color = "gray90", linetype = "dashed")
    )
  
  print(scatter_plot)
  
  if(save_plot) {
    output_file <- file.path(output_dir, paste0(gene, "_expression_comparison.tiff"))
    ggsave(output_file, scatter_plot, width = width, height = height, dpi = dpi)
    message(paste0("图表已保存至: ", output_file))
  }
  
  return(list(
    plot = scatter_plot,
    stats = list(
      log2FC = log2FC,
      adj_p_val = adj_p_val,
      is_significant = is_significant
    )
  ))
}
