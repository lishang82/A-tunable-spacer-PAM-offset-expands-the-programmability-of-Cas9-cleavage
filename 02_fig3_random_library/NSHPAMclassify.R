# ==========================================================================
# 载入依赖包
# ==========================================================================
library(dplyr)
library(stringr)

# ==========================================================================
# 1. 设置输入路径与输出目录
# ==========================================================================
# 输入文件绝对路径
input_file <- "D:/0 depth sequencing/20260409-5NRR 胶回收 cut/merged/cut/ISH1/20260527/20nt_in1_M.txt"

# 设置两个不同的输出目录
ish1_out_dir <- "D:/0 depth sequencing/20260409-5NRR 胶回收 cut/merged/cut/ISH1/20260527"

# 如果输出目录不存在，则自动创建
if (!dir.exists(ish1_out_dir)) dir.create(ish1_out_dir, recursive = TRUE)

# 读取源数据表
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ==========================================================================
# 2. 核心逻辑：基于 PAM 和 PAM2 打上全局统一标签
# ==========================================================================
categorized_df <- df %>%
  mutate(
    # 提取 PAM 和 PAM2 的第 2、3 个碱基
    PAM_suffix = str_sub(PAM, 2, 3),
    
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
ish1_categories <- c("NAA", "NAG", "NGA", "NGG")
ish1_file_count <- 0

for (cat in ish1_categories) {
  final_ish1_df <- categorized_df %>%
    filter(PAM_Category == cat) %>%
    select(-PAM_suffix, -PAM_Category) # 剔除辅助列，保持源表干净
  
  if (nrow(final_ish1_df) > 0) {
    out_name <- file.path(ish1_out_dir, paste0("NSH_20nt_in1_M_PAM-", cat, ".txt"))
    write.table(final_ish1_df, out_name, sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("  [√] 已保存: %s (共 %d 行)", basename(out_name), nrow(final_ish1_df)))
    ish1_file_count <- ish1_file_count + 1
  }
}

# ==========================================================================
# 5. 完成总结
# ==========================================================================
message("\n>>> 数据拆分分发完毕！")
message(sprintf("  -> 成功输出 ISH1 目标文件: %d 个", ish1_file_count))