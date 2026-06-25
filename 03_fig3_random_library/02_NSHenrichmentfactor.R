# ==========================================================================
# 载入依赖包
# ==========================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork) 
library(Biostrings)
library(stringr) 
library(scales)  

# ==========================================================================
# 1. 设置绝对路径与提取样本名
# ==========================================================================
# 读取库文件 (lib_file)
lib_file_path <- "D:/0 depth sequencing/20251020-5NRR(20250924_5NRR_uncut+cut)/results/repeat1/M12K_1/过滤code2对应多个random/Merged_Clean_Counts_3Reps.txt"
lib_data <- read.table(lib_file_path, header=TRUE, sep="\t")

# 读取三份重复数据
file1 <- "D:/0 depth sequencing/20260321-5NRR 胶回收 cut/results/cut/NSH/result.txt"
file2 <- "D:/0 depth sequencing/20260408-5NRR 胶回收 cut/results/cut/NSH/result.txt"
file3 <- "D:/0 depth sequencing/20260409-5NRR 胶回收 cut/results/cut/NSH/result.txt"

# 设置新的输出目录并自动提取样本名 (NSH)
output_dir <- "D:/0 depth sequencing/20260409-5NRR 胶回收 cut/merged/cut/NSH/20250527"
sample_name <- "NSH" 
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 设置列名并读取数据
col_names <- c("Seq", "Counts", "Cigar", "Start Pos", "Length", "Indel Times", 
               "Indel Len", "Cigar2", "IndelInN", "code2_1", "code2_2", "code2_3")

data1 <- read.table(file1, header=FALSE, sep="\t", stringsAsFactors=FALSE)
data2 <- read.table(file2, header=FALSE, sep="\t", stringsAsFactors=FALSE)
data3 <- read.table(file3, header=FALSE, sep="\t", stringsAsFactors=FALSE)

colnames(data1) <- col_names
colnames(data2) <- col_names
colnames(data3) <- col_names

# ==========================================================================
# 2. 定义函数：反向互补与序列汇总
# ==========================================================================
get_rc <- function(seqs) {
  seqs[is.na(seqs)] <- "" 
  as.character(reverseComplement(DNAStringSet(seqs)))
}

get_randomfull_counts <- function(df, lib) {
  df %>%
    filter(Counts > 3) %>%
    mutate(across(c(code2_1, code2_2, code2_3), get_rc)) %>%
    filter(nchar(code2_1) == 4, nchar(code2_2) == 4, nchar(code2_3) == 4) %>%
    mutate(code2 = paste0(code2_3, code2_2, code2_1)) %>%
    inner_join(lib, by = "code2") %>%
    group_by(randomfull) %>%
    summarise(Counts = sum(Counts, na.rm = TRUE), .groups = "drop")
}

# ==========================================================================
# 3. 核心计算与 Master 表格构建 (严格先合并，后计算频率)
# ==========================================================================
rf_counts1 <- get_randomfull_counts(data1, lib_data)
rf_counts2 <- get_randomfull_counts(data2, lib_data)
rf_counts3 <- get_randomfull_counts(data3, lib_data)

colnames(rf_counts1)[2] <- "NSH_Counts_Rep1"
colnames(rf_counts2)[2] <- "NSH_Counts_Rep2"
colnames(rf_counts3)[2] <- "NSH_Counts_Rep3"

# 3.1 整理 Library 数据：严格按 randomfull 分组，把对应的多行 count 累加 (sum)
lib_info <- lib_data %>%
  group_by(randomfull) %>%
  summarise(
    code2 = paste(unique(code2), collapse = ","),
    random = paste(unique(random), collapse = ","),
    PAM = paste(unique(PAM), collapse = ","),
    Lib_Counts_Rep1 = sum(Counts_Rep1, na.rm = TRUE),
    Lib_Counts_Rep2 = sum(Counts_Rep2, na.rm = TRUE),
    Lib_Counts_Rep3 = sum(Counts_Rep3, na.rm = TRUE),
    Lib_avg_counts = round(sum(avg_counts, na.rm = TRUE)),
    .groups = "drop"
  )

# 3.2 整理 NSH 实验数据：合并三次重复，并在行级保留 NA 为 0
merged_nsh <- rf_counts1 %>%
  full_join(rf_counts2, by = "randomfull") %>%
  full_join(rf_counts3, by = "randomfull")

merged_nsh[is.na(merged_nsh)] <- 0
merged_nsh <- merged_nsh %>%
  mutate(NSH_avg_counts = round((NSH_Counts_Rep1 + NSH_Counts_Rep2 + NSH_Counts_Rep3) / 3))

# 3.3 构建 Master 数据表：完全合并 Lib 和 NSH (锁定最终行数)
master_df <- lib_info %>%
  left_join(merged_nsh, by = "randomfull")

# 把 NSH 没测到、但 Lib 存在的靶点 Counts 补 0
master_df[is.na(master_df)] <- 0

# 3.4 全局频率计算：基于锁定后的 Master 表行数计算，保证分母 100% 正确
master_df <- master_df %>%
  mutate(
    # 文库总频率
    lib_Freq = Lib_avg_counts / sum(Lib_avg_counts, na.rm = TRUE),
    
    # 实验组各重复频率 (用于画散点图)
    Freq1 = NSH_Counts_Rep1 / sum(NSH_Counts_Rep1, na.rm = TRUE),
    Freq2 = NSH_Counts_Rep2 / sum(NSH_Counts_Rep2, na.rm = TRUE),
    Freq3 = NSH_Counts_Rep3 / sum(NSH_Counts_Rep3, na.rm = TRUE),
    
    # 实验组均值频率
    mean_Freq = NSH_avg_counts / sum(NSH_avg_counts, na.rm = TRUE),
    
    # 富集因子
    enrichment_factor = mean_Freq / lib_Freq
  )

# ==========================================================================
# 4. 数据拆分与导出 (3个文件底层数据绝对统一)
# ==========================================================================

# 4.1 导出中间文件 (按 Lib_avg_counts 降序)
intermediate_df <- master_df %>%
  select(-lib_Freq, -Freq1, -Freq2, -Freq3, -mean_Freq, -enrichment_factor) %>%
  arrange(desc(Lib_avg_counts))

intermediate_file <- file.path(output_dir, paste0(sample_name, "_Intermediate_Counts_Merged.txt"))
write.table(intermediate_df, intermediate_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste0("  [√] 中间合并文件已保存 (已按 Lib_avg_counts 降序排列): ", intermediate_file))

# 4.2 导出 Enrichment Factor 文件 (按 enrichment_factor 降序)
enrichment_df <- master_df %>%
  select(
    random, 
    PAM, 
    lib_counts = Lib_avg_counts, 
    non_offset_gRNA_counts = NSH_avg_counts, 
    lib_freq = lib_Freq, 
    non_offset_gRNA_freq = mean_Freq, 
    enrichment_factor
  ) %>%
  arrange(desc(enrichment_factor))

enrichment_file <- file.path(output_dir, paste0(sample_name, "_Enrichment_Factor.txt"))
write.table(enrichment_df, enrichment_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste0("  [√] Enrichment Factor 数据表已保存 (已按 enrichment_factor 降序排列): ", enrichment_file))

# ==========================================================================
# 5. 补充：计算并保存 Reads 保留率统计表 (Retention Rates)
# ==========================================================================
calc_stats <- function(raw_df, final_counts_col) {
  original <- sum(raw_df$Counts, na.rm = TRUE)
  denoised <- sum(raw_df$Counts[raw_df$Counts > 3], na.rm = TRUE)
  final_match <- sum(final_counts_col, na.rm = TRUE)
  return(c(original, denoised, final_match))
}

stats1 <- calc_stats(data1, rf_counts1$NSH_Counts_Rep1)
stats2 <- calc_stats(data2, rf_counts2$NSH_Counts_Rep2)
stats3 <- calc_stats(data3, rf_counts3$NSH_Counts_Rep3)

total_counts_df <- data.frame(
  Sample = c("Rep1_0321", "Rep2_0408", "Rep3_0409"),
  Original_Counts = c(stats1[1], stats2[1], stats3[1]),
  Denoised_Counts = c(stats1[2], stats2[2], stats3[2]),
  Final_Match_Counts = c(stats1[3], stats2[3], stats3[3])
) %>%
  mutate(
    Denoising_Retention = Denoised_Counts / Original_Counts,
    Total_Retention = Final_Match_Counts / Original_Counts,
    Denoising_Rate_Pct = paste0(round(Denoising_Retention * 100, 2), "%"),
    Total_Retention_Pct = paste0(round(Total_Retention * 100, 2), "%")
  )

retention_file <- file.path(output_dir, paste0(sample_name, "_Retention_Rates.txt"))
write.table(total_counts_df, retention_file, sep = "\t", quote = FALSE, row.names = FALSE)

# ==========================================================================
# 6. 绘图 1: 找回来的 三组重复性散点图 (p1 | p2 | p3)
# ==========================================================================
plot_scatter_rep <- function(df, x_col, y_col, x_label, y_label) {
  valid_df <- df %>% filter(.data[[x_col]] > 0, .data[[y_col]] > 0)
  cor_val <- cor(log10(valid_df[[x_col]]), log10(valid_df[[y_col]]), method = "pearson")
  cor_text <- sprintf("r = %.3f", cor_val)
  
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.3, size = 1, color = "dodgerblue4") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(x = paste(x_label, "Frequency"), y = paste(y_label, "Frequency"), title = paste(x_label, "vs", y_label)) +
    annotate("text", x = 10^-5, y = 10^-1, label = cor_text, size = 5, hjust = 0) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid.minor = element_blank())
}

p1 <- plot_scatter_rep(master_df, "Freq1", "Freq2", "Rep1 (0321)", "Rep2 (0408)")
p2 <- plot_scatter_rep(master_df, "Freq1", "Freq3", "Rep1 (0321)", "Rep3 (0409)")
p3 <- plot_scatter_rep(master_df, "Freq2", "Freq3", "Rep2 (0408)", "Rep3 (0409)")

rep_plot <- p1 | p2 | p3
ggsave(file.path(output_dir, paste0(sample_name, "_Replicates_Correlation.pdf")), plot = rep_plot, width = 15, height = 5)
message(paste0("  [√] 重复性散点图已保存"))

# ==========================================================================
# 7. 绘图数据准备 (直接基于 Master 构建 final_comparison_df)
# ==========================================================================
final_comparison_df <- master_df %>%
  filter(!str_detect(str_sub(PAM, 2, 3), "[TC]")) %>%
  mutate(
    NGG = ifelse(PAM %in% c("AGG", "TGG", "CGG", "GGG"), 1, 0),
    non_NGG = ifelse(PAM %in% c("AGG", "TGG", "CGG", "GGG"), 0, 1),
    Point_Type = ifelse(NGG == 1, "NGG points", "non-NGG points")
  ) %>%
  select(randomfull, random, PAM, lib_Freq, mean_Freq, enrichment_factor, NGG, non_NGG, Point_Type)

txt_filename <- file.path(output_dir, paste0(sample_name, "_processed_data.txt"))
write.table(final_comparison_df, txt_filename, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste0("  [√] 最终绘图数据表已保存"))

# ==========================================================================
# 8. 绘图 2: 实验频率 vs 文库频率 Scatter
# ==========================================================================
set.seed(123) 
common_limits <- c(10^-5.5, 10^-2.5)
p_scatter <- ggplot(final_comparison_df %>% sample_frac(1), aes(x = lib_Freq, y = mean_Freq, color = Point_Type)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("non-NGG points" = "#255C99", "NGG points" = "#F7B32B")) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  scale_x_log10(limits = common_limits, breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(limits = common_limits, breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  labs(title = sample_name, x = "Frequency in Library", y = paste("Frequency in", sample_name),
       subtitle = paste("NGG:", sum(final_comparison_df$NGG == 1), "| non-NGG:", sum(final_comparison_df$NGG == 0))) +  
  theme_minimal() +
  theme(legend.position = "none", panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), axis.title = element_text(size = 20), axis.text = element_text(size = 20))

ggsave(file.path(output_dir, paste0(sample_name, "_counts_frequency.pdf")), plot = p_scatter, width = 8, height = 8)

# ==========================================================================
# 9. 绘图 3: Enrichment Factor Histogram (真实最大最小值，50个bin)
# ==========================================================================
jiao_plot <- final_comparison_df %>%
  mutate(direction = ifelse(NGG == 1, "up", "down")) %>%
  filter(is.finite(enrichment_factor))

limit_min <- min(jiao_plot$enrichment_factor, na.rm = TRUE)
limit_max <- max(jiao_plot$enrichment_factor, na.rm = TRUE)

if (limit_max <= limit_min) limit_max <- limit_min + 1 

bin_count <- 50
breaks_seq <- seq(limit_min, limit_max, length.out = bin_count + 1)

hist_up <- hist(jiao_plot %>% filter(direction == "up") %>% pull(enrichment_factor), breaks = breaks_seq, plot = FALSE)
hist_down <- hist(jiao_plot %>% filter(direction == "down") %>% pull(enrichment_factor), breaks = breaks_seq, plot = FALSE)
plot_data <- data.frame(x_mid = hist_up$mids, up_count = hist_up$counts, down_count = -hist_down$counts, x_min = head(breaks_seq, -1), x_max = tail(breaks_seq, -1))

max_up_val <- max(plot_data$up_count, na.rm = TRUE); max_down_val <- max(abs(plot_data$down_count), na.rm = TRUE)
breaks_up_raw <- pretty(c(0, max_up_val), n = 4); breaks_down_raw <- pretty(c(0, max_down_val), n = 4)
max_up_tick <- max(breaks_up_raw); max_down_tick <- max(breaks_down_raw)
if(max_up_tick == 0) max_up_tick <- 10; if(max_down_tick == 0) max_down_tick <- 10

scale_ratio <- max_up_tick / max_down_tick
plot_data <- plot_data %>% mutate(scaled_down_count = down_count * scale_ratio)
breaks_down_mapped <- -breaks_down_raw * scale_ratio
all_breaks <- c(breaks_down_mapped[breaks_down_raw != 0], breaks_up_raw)
all_labels <- c(breaks_down_raw[breaks_down_raw != 0], breaks_up_raw)

p_hist <- ggplot() +
  geom_rect(data = plot_data, aes(xmin = x_min, xmax = x_max, ymin = 0, ymax = up_count), fill = "#F7B32B", alpha = 0.7, color = "black", linewidth = 0.2) +
  geom_rect(data = plot_data, aes(xmin = x_min, xmax = x_max, ymin = scaled_down_count, ymax = 0), fill = "#255C99", alpha = 0.7, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  scale_x_continuous(name = "Enrichment Factor", limits = c(limit_min, limit_max), breaks = pretty(c(limit_min, limit_max), n = 10), expand = c(0, 0)) +
  scale_y_continuous(name = "Frequency", breaks = all_breaks, labels = all_labels) +
  annotate("text", x = limit_max * 0.95, y = max_up_tick * 0.9, label = paste("NGG:", sum(hist_up$counts)), color = "#F7B32B", size = 6, fontface = "bold", hjust = 1) +
  annotate("text", x = limit_max * 0.95, y = -max_up_tick * 0.9, label = paste("non-NGG:", sum(hist_down$counts)), color = "#255C99", size = 6, fontface = "bold", hjust = 1) +
  labs(title = paste(sample_name, "gRNA")) + 
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = 14, face = "bold", color = "black"), axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(output_dir, paste0(sample_name, "_Enrichment_Factor_Histogram.pdf")), plot = p_hist, width = 6, height = 4)
message(paste0("  [√] 最终分析图表全部保存至: ", output_dir))
message(">>> Pipeline 跑完啦！\n")
