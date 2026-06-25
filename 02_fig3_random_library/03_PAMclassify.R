# ==========================================================================
# 载入依赖包
# ==========================================================================
library(dplyr)
library(stringr)

# ==========================================================================
# 1. 设置输入路径与输出目录
# ==========================================================================
# 输入文件绝对路径
input_file <- "D:/0 depth sequencing/20260409-5NRR 胶回收 cut/merged/cut/OSH1/20260527/20nt_in1_M.txt"

# 设置两个不同的输出目录
OSH1_out_dir <- "D:/0 depth sequencing/20260409-5NRR 胶回收 cut/merged/cut/OSH1/20260527"
# 如果输出目录不存在，则自动创建
if (!dir.exists(nsh_out_dir)) dir.create(nsh_out_dir, recursive = TRUE)
if (!dir.exists(osh1_out_dir)) dir.create(osh1_out_dir, recursive = TRUE)

# 读取源数据表
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ==========================================================================
# 2. 核心逻辑：基于 PAM 和 PAM2 打上全局统一标签
# ==========================================================================
categorized_df <- df %>%
  mutate(
    # 提取 PAM 和 PAM2 的第 2、3 个碱基
    PAM_suffix = str_sub(PAM, 2, 3),
    PAM2_suffix = str_sub(PAM2, 2, 3),
    
    # 判定 PAM2 是否为 NGG (即后两位是 GG)
    is_PAM2_NGG = ifelse(PAM2_suffix == "GG", TRUE, FALSE),
    
    # 判定 PAM 属于四类中的哪一类
    PAM_Category = case_when(
      PAM_suffix == "AA" ~ "NAA",
      PAM_suffix == "AG" ~ "NAG",
      PAM_suffix == "GA" ~ "NGA",
      PAM_suffix == "GG" ~ "NGG",
      TRUE ~ "Other" # 过滤掉不符合规律的杂项序列
    )
  ) %>%
  filter(PAM_Category != "Other")


# ==========================================================================
# 3. 任务 A: 导出 NSH 目录的 4 个文件
# ==========================================================================
message("\n=== 开始生成 NSH 目录下的 4 个文件 ===")
nsh_categories <- c("NAA", "NAG", "NGA", "NGG")
nsh_file_count <- 0

for (cat in nsh_categories) {
  final_nsh_df <- categorized_df %>%
    filter(PAM_Category == cat) %>%
    select(-PAM_suffix, -PAM2_suffix, -is_PAM2_NGG, -PAM_Category) # 剔除辅助列，保持源表干净
  
  if (nrow(final_nsh_df) > 0) {
    out_name <- file.path(nsh_out_dir, paste0("NSH_20nt_out1_M_PAM-", cat, ".txt"))
    write.table(final_nsh_df, out_name, sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("  [√] 已保存: %s (共 %d 行)", basename(out_name), nrow(final_nsh_df)))
    nsh_file_count <- nsh_file_count + 1
  }
}


# ==========================================================================
# 4. 任务 B: 导出 OSH1 目录的 6 个文件
# ==========================================================================
message("\n=== 开始生成 OSH1 目录下的 6 个文件 ===")
osh1_tasks <- list(
  # 集合 1：PAM2 为 NGG，按 PAM 提取 NGG 和 NGA
  list(pam2_condition = TRUE,  pam2_label = "PAM2-NGG", target_pams = c("NGG", "NGA")),
  # 集合 2：去除 PAM2 为 NGG 的情况，按 PAM 提取 NAA, NAG, NGA, NGG
  list(pam2_condition = FALSE, pam2_label = "PAM2-nonNGG", target_pams = c("NAA", "NAG", "NGA", "NGG"))
)
osh1_file_count <- 0

for (task in osh1_tasks) {
  # 根据 PAM2 逻辑筛选大池子
  sub_df <- categorized_df %>% filter(is_PAM2_NGG == task$pam2_condition)
  
  for (cat in task$target_pams) {
    final_osh1_df <- sub_df %>%
      filter(PAM_Category == cat) %>%
      select(-PAM_suffix, -PAM2_suffix, -is_PAM2_NGG, -PAM_Category) # 剔除辅助列
    
    if (nrow(final_osh1_df) > 0) {
      out_name <- file.path(osh1_out_dir, paste0("OSH1_20nt_out1_M_", task$pam2_label, "_PAM-", cat, ".txt"))
      write.table(final_osh1_df, out_name, sep = "\t", quote = FALSE, row.names = FALSE)
      message(sprintf("  [√] 已保存: %s (共 %d 行)", basename(out_name), nrow(final_osh1_df)))
      osh1_file_count <- osh1_file_count + 1
    }
  }
}

# ==========================================================================
# 5. 完成总结
# ==========================================================================
message("\n>>> 数据拆分分发完毕！")
message(sprintf("  -> 成功输出 NSH 目标文件: %d 个", nsh_file_count))
message(sprintf("  -> 成功输出 OSH1 目标文件: %d 个", osh1_file_count))
