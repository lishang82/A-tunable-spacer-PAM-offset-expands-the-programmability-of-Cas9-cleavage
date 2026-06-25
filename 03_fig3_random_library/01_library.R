#!/usr/bin/env Rscript

# 加载必要的包
library(tidyverse)
library(ggVennDiagram)

# ==============================================================================
# 0. 命令行参数接收与路径动态处理
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript process_replicates.R <sample1_path> <sample2_path> <sample3_path>")
}

sample1_path <- args[1]
sample2_path <- args[2]
sample3_path <- args[3]

# 检查文件是否存在
if (!all(file.exists(c(sample1_path, sample2_path, sample3_path)))) {
  stop("Error: One or more input files do not exist. Please check your paths.")
}

# 动态获取输出目录与样本文件夹名称
output_dir <- dirname(sample1_path)
folder_name <- basename(dirname(sample1_path))

# 读入数据
data2 <- read.table(sample1_path, header = FALSE, fill = TRUE, sep = "\t")
data3 <- read.table(sample2_path, header = FALSE, fill = TRUE, sep = "\t")
data4 <- read.table(sample3_path, header = FALSE, fill = TRUE, sep = "\t")

# 设置原始列名
orig_cols <- c("Seq", "Counts", "Cigar", "Start Pos", "Length", "Indel Times", 
               "Indel Len", "Cigar2", "IndelInN", "code1_1", "code1_2", "code1_3", 
               "randomfull", "code2_1", "code2_2", "code2_3")
colnames(data2) <- orig_cols
colnames(data3) <- orig_cols
colnames(data4) <- orig_cols

# ==============================================================================
# 【核心修改点】拼接标签 并 自动拆分 randomfull -> random + PAM
# ==============================================================================
# 假设 random 长度为 7 bp，PAM 在后面。如果是其他长度（例如random 8bp），请修改 sep = 7
process_mutations <- function(df) {
  df %>% 
    mutate(code2 = paste(code2_1, code2_2, code2_3, sep = "")) %>%
    # separate 逻辑：在第 7 个碱基后面切开，前面叫 random，后面叫 PAM
    separate(col = randomfull, into = c("random", "PAM"), sep = 7, remove = FALSE)
}

data2 <- process_mutations(data2)
data3 <- process_mutations(data3)
data4 <- process_mutations(data4)

# --- 1. 定义一个通用的清洗函数 ---
process_replicate <- function(df, rep_name) {
  
  # 预聚合 (防止原始数据里有重复行)
  df_agg <- df %>%
    group_by(code2, randomfull) %>%
    summarise(Counts = sum(Counts), .groups = 'drop')
  
  # 标记每个 code2 对应多少种 randomfull
  df_labeled <- df_agg %>%
    group_by(code2) %>%
    mutate(n_variants = n()) %>%
    ungroup()
  
  # --- 提取干净数据 (Unique: 1对1) ---
  df_clean <- df_labeled %>%
    filter(n_variants == 1) %>%
    select(code2, randomfull, Counts) %>%
    rename(!!paste0("Counts_", rep_name) := Counts)
  
  # --- 提取污染数据 (Multi: 1对多) ---
  df_multi_codes <- df_labeled %>%
    filter(n_variants > 1) %>%
    distinct(code2) %>%
    mutate(Source = rep_name)
  
  return(list("clean" = df_clean, "multi" = df_multi_codes))
}

# --- 2. 分别处理三个数据集 ---
res2 <- process_replicate(data2, "Rep1")
res3 <- process_replicate(data3, "Rep2")
res4 <- process_replicate(data4, "Rep3")

print(paste("Rep1 干净序列数:", nrow(res2$clean)))
print(paste("Rep2 干净序列数:", nrow(res3$clean)))
print(paste("Rep3 干净序列数:", nrow(res4$clean)))

# ==============================================================================
# 3. 干净数据的交集与合并
# ==============================================================================
merged_clean_data2 <- res2$clean %>%
  inner_join(res3$clean, by = c("code2", "randomfull")) %>%
  inner_join(res4$clean, by = c("code2", "randomfull"))

# 在最终输出的交集母表中，也把 random 和 PAM 拆分出来方便下游 Origin 使用
# 同样假设 random 长度为 7 bp
merged_clean_data2 <- merged_clean_data2 %>%
  separate(col = randomfull, into = c("random", "PAM"), sep = 7, remove = FALSE)

print("------------------------------------------------")
print(paste("【干净数据】三次重复的完美交集数量:", nrow(merged_clean_data2)))

# 动态保存到对应的样本输出目录下，带有样本前缀
output_csv <- file.path(output_dir, paste0(folder_name, "_Merged_Clean_Counts_3Reps.csv"))
write.csv(merged_clean_data2, output_csv, row.names = FALSE)
print(paste("结果已保存至:", output_csv))

# ==============================================================================
# 4. 污染数据的重叠情况 (查看)
# ==============================================================================
intersect_multi_codes <- intersect(res2$multi$code2, res3$multi$code2) %>%
  intersect(res4$multi$code2)

print("------------------------------------------------")
print(paste("【污染数据】三次重复中都表现为'一对多'的 code2 数量:", length(intersect_multi_codes)))

# ==============================================================================
# 5. 可视化检查 (韦恩图自动保存)
# ==============================================================================
venn_list <- list(
  Rep1 = paste0(res2$clean$code2, "_", res2$clean$randomfull),
  Rep2 = paste0(res3$clean$code2, "_", res3$clean$randomfull),
  Rep3 = paste0(res4$clean$code2, "_", res4$clean$randomfull)
)

p_venn <- ggVennDiagram(venn_list, label_alpha = 0) + 
  ggplot2::scale_fill_gradient(low="white", high = "red") +
  ggplot2::labs(title = paste0(folder_name, " - Overlap of Clean Sequences"))

# 自动保存韦恩图到本地，避免 Plots 窗口丢失
output_venn <- file.path(output_dir, paste0(folder_name, "_clean_venn.png"))
ggplot2::ggsave(output_venn, plot = p_venn, width = 7, height = 6, dpi = 300)
print(paste("韦恩图已保存至:", output_venn))
