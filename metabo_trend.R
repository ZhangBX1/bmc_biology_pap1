# 加载必要的包
library(tidyverse)
library(readr)
library(ggplot2)

# 读取数据
data <- read_csv("pos_measure.csv")

data <- read_csv("neg_measure.csv")
# 删除第二列没有名字的行（即Accepted Description为空的行）
data <- data %>% filter(!is.na(`Accepted Description`) & `Accepted Description` != "")

# 将数据转换为长格式，便于分析和绘图
data_long <- data %>%
  pivot_longer(
    cols = -`Accepted Description`,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  # 从样本名称中提取信息
  mutate(
    genotype = case_when(
      substr(sample, 1, 2) == "WT" ~ "Wild Type",
      substr(sample, 1, 2) == "MU" ~ "Mutant",
      TRUE ~ "Unknown"
    ),
    # 提取W后面的数字作为周数
    week = as.numeric(str_extract(sample, "(?<=W)\\d")),
    # 提取最后一个数字作为重复数
    replicate = as.numeric(str_extract(sample, "\\d$"))
  )

# 检查提取是否正确
sample_check <- data_long %>% 
  select(sample, genotype, week, replicate) %>% 
  distinct()
print(head(sample_check, 10))



plot_metabolite_trends <- function(metabolite_name, save_path = NULL, data_source = data_long) {
  if (!is.character(metabolite_name)) {
    stop("metabolite_name must be a character string")
  }
  
  # 筛选目标代谢物
  metabolite_data <- data_source %>%
    filter(`Accepted Description` == metabolite_name)
  
  # 检查是否有匹配的代谢物
  if (nrow(metabolite_data) == 0) {
    warning(paste("No data found for metabolite:", metabolite_name))
    return(NULL)
  }
  
  # 计算显著性（t 检验）
  significance_data <- metabolite_data %>%
    group_by(week) %>%
    summarise(
      p_value = ifelse(
        length(unique(genotype)) > 1,  # 确保有多个基因型
        t.test(abundance ~ genotype)$p.value,
        NA
      ),
      .groups = "drop"
    ) %>%
    mutate(
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  # 合并显著性数据
  summary_data <- metabolite_data %>%
    group_by(`Accepted Description`, genotype, week) %>%
    summarise(
      mean_abundance = mean(abundance, na.rm = TRUE),
      se = sd(abundance, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    left_join(significance_data, by = "week")
  
  # 创建折线图
  p <- ggplot(summary_data, aes(x = week, y = mean_abundance, color = genotype, group = genotype)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = mean_abundance - se, ymax = mean_abundance + se), width = 0.3, linewidth = 1) +
    facet_wrap(~ `Accepted Description`, scales = "free_y") +
    labs(
      x = "Week",
      y = "Abundance",
      color = "Groups"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "bottom",
      legend.key.size = unit(1.2, "cm"),
      panel.grid.major = element_line(linewidth = 0.5),
      panel.grid.minor = element_line(linewidth = 0.25)
    ) +
    geom_text(
      data = summary_data %>% filter(genotype == "Mutant"),  # 仅对 Mutant 添加星标
      aes(x = week, y = mean_abundance + se + 0.05 * max(mean_abundance), label = significance),
      color = "black",
      size = 6,
      vjust = 0
    )
  
  # 如果提供了保存路径，则保存图片
  if (!is.null(save_path)) {
    file_path <- file.path(save_path, paste0(metabolite_name, ".tiff"))
    ggsave(file_path, plot = p, width = 6, height = 8, dpi = 300)
    message(paste("Plot saved to:", file_path))
  }
  
  return(p)
}


a="L-Tryptophan"

plot_metabolite_trends(a,"F://paper//wet_experiment//lcms//pos//")


a="Glucoiberin"

p <- plot_metabolite_trends(a)
print(p)

