# 1. 修改为相对路径（假设你把数据放到了脚本同级或下级的 data 文件夹中）
dir1 <- "./data/20260202/子勤"
dir2 <- "./data/20260210/子勤"
dir3 <- "./data/20260402/子勤"

# 定义统一保存合并结果的新文件夹（也放在当前目录下）
out_dir <- "./Merged_Results"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
# 2. 抓取两个文件夹下所有的 .txt 文件
all_txt_files <- c(
  list.files(path = dir1, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE),
  list.files(path = dir2, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE),
  list.files(path = dir3, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
)

# 排除掉原始的 Indel_histogram.txt
target_files <- all_txt_files[!grepl("Indel_histogram\\.txt$", all_txt_files)]

# 3. 智能提取样本组名并进行配对
file_names <- tools::file_path_sans_ext(basename(target_files))
group_names <- sub("-[123]$", "", file_names)

file_info <- data.frame(
  path = target_files,
  group = group_names,
  stringsAsFactors = FALSE
)

grouped_files <- split(file_info$path, file_info$group)

# 4. 遍历每一组进行合并与绘图
for (grp in names(grouped_files)) {
  paths <- grouped_files[[grp]]
  
  if (length(paths) != 3) {
    message(sprintf("⚠️ 跳过样本组 [%s] : 找到的文件数量不是2个 (当前为 %d 个)", grp, length(paths)))
    next
  }
  
  # 使用更稳定、专为 Tab 分割设计的 read.delim，并关闭 check.names 防止列名被 R 乱改
  data1 <- read.delim(paths[1], stringsAsFactors = FALSE, check.names = FALSE)
  data2 <- read.delim(paths[2], stringsAsFactors = FALSE, check.names = FALSE)
  data3 <- read.delim(paths[3], stringsAsFactors = FALSE, check.names = FALSE)
  
  # --- 新增的安检门：检查是否真的包含 indel_size 列 ---
  if (!"indel_size" %in% colnames(data1)) {
    message(sprintf("❌ 格式错误跳过: 文件 [%s] 中没有 indel_size 列！", basename(paths[1])))
    next
  }
  if (!"indel_size" %in% colnames(data2)) {
    message(sprintf("❌ 格式错误跳过: 文件 [%s] 中没有 indel_size 列！", basename(paths[2])))
    next
  }
  if (!"indel_size" %in% colnames(data3)) {
    message(sprintf("❌ 格式错误跳过: 文件 [%s] 中没有 indel_size 列！", basename(paths[3])))
    next
  }
  
  
  # ----------------------------------------------------
  
  # 5. 合并数据并计算 Mean 和 SEM
  combined_data <- bind_rows(data1, data2, data3) %>%
    filter(!is.na(indel_size), indel_size >= -15, indel_size <= 5)
  
  summary_data <- combined_data %>%
    group_by(indel_size) %>%
    summarise(
      mean_pct = mean(percentage, na.rm = TRUE),
      sd_pct = sd(percentage, na.rm = TRUE),
      n = n(),
      sem_pct = sd_pct / sqrt(n), 
      .groups = 'drop'
    )
  
  summary_data$sem_pct[is.na(summary_data$sem_pct)] <- 0
  
  # 6. 导出数据表
  txt_output_path <- file.path(out_dir, paste0(grp, "_merged.txt"))
  write.table(summary_data, txt_output_path, row.names = FALSE, sep = "\t", quote = FALSE)
  
  # 7. 绘图
  y_max_limit <- max(50, max(summary_data$mean_pct + summary_data$sem_pct) * 1.05)
  
  p <- ggplot(summary_data, aes(x = factor(indel_size), y = mean_pct)) +
    geom_bar(stat = "identity", fill = "#2ca0d3", color = NA, width = 0.8) +
    geom_errorbar(
      aes(
        ymin = pmax(0, mean_pct - sem_pct), 
        ymax = mean_pct + sem_pct
      ), 
      width = 0.3, 
      color = "black", 
      size = 0.6
    ) +
    scale_y_continuous(limits = c(0, y_max_limit), expand = c(0, 0)) + 
    scale_x_discrete(
      breaks = seq(-15, 5, by = 5),  # 只在这些位置显示刻度标签
      labels = seq(-15, 5, by = 5)   # 对应的标签
    ) + 
    labs(
      x = "Indel Size (bp)",
      y = "indel frequency (%)"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 0, size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.y = element_line(color = "grey90")
    )
  
  # 8. 保存 PDF
  pdf_output_path <- file.path(out_dir, paste0(grp, "_merged.pdf"))
  ggsave(pdf_output_path, plot = p, width = 4, height = 5)
  
  message(paste("✅ 成功合并:", grp))
}

print("--- 任务运行完毕 ---")

# 载入必须的包
library(ggplot2)
library(patchwork)
library(dplyr)
root_path <- out_dir
# 1. 基础配置
groups <- c("sg10", "sg6",  "sg123", "sg245","sg265")

# 存储所有绘图对象的列表
plot_list_top <- list()
plot_list_bottom <- list()

# 2. 核心处理逻辑函数 - 新增控制y轴标题显示的参数
get_merged_plot <- function(f_path, title_name, show_y_title = TRUE) {
  # 读取已经合并好的数据
  data <- read.table(f_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  data <- data %>% 
    mutate(indel_size = as.numeric(indel_size)) %>%
    filter(!is.na(indel_size), indel_size >= -15, indel_size <= 5)
  
  # 找出最高柱子的 indel_size（基于 mean_pct）
  max_indel_size <- data$indel_size[which.max(data$mean_pct)]
  
  # 创建一个新列，标记是否为最高柱子
  data$is_max <- ifelse(data$indel_size == max_indel_size, "max", "other")
  
  # 为了美观，如果有误差棒超过 50，动态调整 Y 轴上限
  y_max <- max(50, max(data$mean_pct + data$sem_pct, na.rm = TRUE) * 1.05)
  
  p <- ggplot(data, aes(x = factor(indel_size), y = mean_pct, fill = is_max)) +
    geom_bar(stat = "identity", color = NA, width = 0.8) +
    # 手动设置填充颜色：最高柱子红色，其他灰色
    scale_fill_manual(values = c("max" = "red", "other" = "#2ca0d3")) +
    # 增加误差棒
    geom_errorbar(
      aes(ymin = pmax(0, mean_pct - sem_pct), ymax = mean_pct + sem_pct), 
      width = 0.3, color = "black", size = 0.6
    ) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    scale_x_discrete(
      breaks = seq(-15, 5, by = 5),
      labels = seq(-15, 5, by = 5)
    ) + 
    labs(title = title_name, 
         x = NULL, 
         y = ifelse(show_y_title, "InDels frequency (%)", "")) + 
    theme_classic(base_size = 8) + 
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
      axis.title = element_text(size = 15, face = "bold"),              
      axis.text.y = element_text(size = 12, color = "black"),           
      axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, face = "bold"), 
      panel.grid.major.y = element_line(color = "grey95"),
      legend.position = "none"
    )
  
  # 如果不显示y轴标题，同时隐藏y轴标题的占位空间
  if (!show_y_title) {
    p <- p + theme(axis.title.y = element_blank())
  }
  
  return(p)
}

# 生成空白图的函数
get_empty_plot <- function(title_name) {
  ggplot() + 
    theme_void() + 
    labs(title = paste(title_name, "(Data Missing)")) +
    theme(plot.title = element_text(size = 12, hjust = 0.5, color = "red"))
}

# 3. 按照 5 列结构进行精准匹配
for (i in 1:length(groups)) {
  g <- groups[i]
  f_top <- file.path(root_path, paste0(g, "_merged.txt"))
  f_bottom <- file.path(root_path, paste0("out1-", g, "_merged.txt"))
  
  # 判断是否是最左边列（第1列 = 索引1，或第6列 = 索引1在第二行）
  # 对于第一行（top），第1列显示y轴标题；其他列不显示
  # 对于第二行（bottom），第1列也显示y轴标题
  show_y_top <- (i == 1)  # 第一行第一列
  show_y_bottom <- (i == 1)  # 第二行第一列
  
  if (file.exists(f_top)) {
    plot_list_top[[g]] <- get_merged_plot(f_top, g, show_y_title = show_y_top)
  } else {
    plot_list_top[[g]] <- get_empty_plot(g)
  }
  
  if (file.exists(f_bottom)) {
    plot_list_bottom[[g]] <- get_merged_plot(f_bottom, paste0("out1-", g), show_y_title = show_y_bottom)
  } else {
    plot_list_bottom[[g]] <- get_empty_plot(paste0("out1-", g))
  }
}

# 4. 拼图 (2行 5列)
combined_plots <- c(plot_list_top, plot_list_bottom)

final_plot <- wrap_plots(combined_plots, ncol = 5, nrow = 2) +
  plot_annotation(
    title = "Indel Distribution Comparison (Top: Original | Bottom: Out1)",
    subtitle = "Mean ± SEM (n=3 biological replicates). Red bar indicates the highest frequency.",
  )

# 5. 保存最终大图
output_pdf <- file.path(root_path, "Combined_Merged_Indel_Report_Highlighted.pdf")
ggsave(output_pdf, plot = final_plot, width = 12, height = 8)

message(paste("✅ 拼图完成！最高柱子已标红，文件存至:", output_pdf))
