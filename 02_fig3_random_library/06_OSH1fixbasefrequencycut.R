library(tidyverse)
library(gridExtra)
library(ggseqlogo) 
library(patchwork)

# --- 1. 数据读取 ---
# 读取 cut 文件夹下的 OSH1 20nt 数据
intersection_data3 <- read.table("D:/0 depth sequencing/20260409-5NRR 胶回收 cut/merged/cut/OSH1/20260527/20nt_out1_M.txt", header=T, sep="\t")

# --- 2. 计算比率与对数转换 ---
# 使用你提供的归一化系数：0.308 和 0.418

# --- 3. 定义颜色方案 ---
my_color_scheme <- make_col_scheme(
  chars = c('A', 'T', 'C', 'G'),
  cols  = c('#109648', '#D62839', '#255C99', '#F7B32B') 
)


# --- 5. 数据过滤与预处理 ---
df_base <- intersection_data3 %>%
  # 仅排除后两位是 GG 的组合 (PAM2 形如 ?GG)
  filter(!str_detect(PAM2, ".GG")) %>% 
  # 生成 7 碱基全序列 (random + PAM) 用于校验
  mutate(full_seq = paste0(random, PAM)) %>%
  # 过滤掉非 7 位碱基的异常行
  filter(nchar(full_seq) == 7) %>%
  # 筛选 NSH_transformed >= 0.5
  filter(enrichment_factor_nonoffset >= 0.7) %>%
  # 按 OSH1_transformed 降序排列，取前 50 条 (此处使用1000)
  arrange(desc(enrichment_factor_out1)) %>%
  slice_head(n = 50) %>%
  # 新增：截取仅用于绘图的前 4 位碱基 (即 random 部分)
  mutate(plot_seq = substr(full_seq, 1, 5))

# --- 2. 提取 5 个位置的碱基因子 ---
df_positions <- df_base %>%
  mutate(
    P1 = factor(substr(plot_seq, 1, 1), levels = c("A", "C", "G", "T")),
    P2 = factor(substr(plot_seq, 2, 2), levels = c("A", "C", "G", "T")),
    P3 = factor(substr(plot_seq, 3, 3), levels = c("A", "C", "G", "T")),
    P4 = factor(substr(plot_seq, 4, 4), levels = c("A", "C", "G", "T")),
    P5 = factor(substr(plot_seq, 5, 5), levels = c("A", "C", "G", "T"))
  )

# --- 3. 核心单图函数：只画单行单列的无边框 Logo ---
plot_single_cell <- function(sequences) {
  # 如果该碱基下没有序列，返回一个完全空白的占位图
  if (length(sequences) == 0) {
    return(plot_spacer() + theme(plot.margin = margin(0,0,0,0)))
  }
  
  ggseqlogo(sequences, method = 'probability', col_scheme = my_color_scheme) +
    # 彻底移除所有多余组件，仅保留坐标轴线与 Major Ticks
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.line  = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(-0.15, "cm"), # 刻度线内朝
      legend.position = "none",
      # 消除所有外边距空间
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      panel.spacing = unit(0, "pt")
    ) +
    # 保证透明度
    scale_alpha_manual(values = c(0.5)) 
}

# --- 4. 循环生成 4x5 = 20 个子图 ---
plot_list <- list()
bases_order <- c("A", "C", "G", "T")
pos_cols <- c("P5", "P4", "P3", "P2", "P1")

# 按照 行(碱基) -> 列(位置) 的顺序生成
for (b in bases_order) {
  for (p in pos_cols) {
    # 筛选出当前位置(p)固定为当前碱基(b)的所有序列
    sub_seqs <- df_positions$plot_seq[df_positions[[p]] == b]
    # 生成单格图并存入列表
    plot_list[[length(plot_list) + 1]] <- plot_single_cell(sub_seqs)
  }
}

# --- 5. 使用 patchwork 强力拼成无缝网格 ---
final_plot <- wrap_plots(plot_list, ncol = 5) +
  plot_layout(widths = rep(1, 5), heights = rep(1, 4)) &
  theme(
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing.x = unit(0, "pt"),
    panel.spacing.y = unit(0, "pt")
  )

# 在图的上方加上你需要的列标题（Fixed Pos...）
# 如果需要加上标题，可以通过 patchwork 的 plot_annotation 完美对齐
final_plot <- final_plot + 
  plot_annotation(
    title = "Fixed Pos +1      Fixed Pos -2      Fixed Pos -3      Fixed Pos -4      Fixed Pos -5",
    theme = theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b=10)))
  )

print(final_plot)


# 打印拼接后的大图
ggsave("D:/0 depth sequencing/20260409-5NRR 胶回收 cut/merged/cut/OSH1/20260527/Sequence_Frequency_Logo50.pdf", final_plot, width = 20, height = 16)
