#!/bin/bash

# ==============================================================================
# CRISPR-Cas9 Target Sequencing Analysis Pipeline
# ==============================================================================
# 使用说明：
# 运行该脚本时，请在下方【用户自定义参数区域】修改对应的参考序列文件名和靶点序列。
# 运行命令示例：bash 02_alignment_and_analysis.sh lane1_R1.fq.gz lane1_R2.fq.gz
# ==============================================================================

# --- 1. 用户自定义参数区域 (根据不同靶点修改此处) ---
REF_FASTA="reference.fa"       # 参考基因组/质粒的 fasta 文件名
REF_INDEX="reference.idx"      # 生成的索引文件名
CUT_REF_SEQ="[+sequence]"      # 当前靶点对应的剪切参考序列 (例如: ATCG...)

# --- 2. 自动获取输入文件 ---
# $1 和 $2 代表运行脚本时输入的 R1 和 R2 原始测序文件
INPUT_R1=$1
INPUT_R2=$2

# --- 3. 核心分析流程 ---

# Step 1: 使用 fastp 进行质控、去接头并合并双端测序数据
echo "Starting fastp QC and merging..."
fastp -i "$INPUT_R1" -I "$INPUT_R2" \
      -o drop_1.fq.gz -O drop_2.fq.gz \
      --merge --merged_out merged.fq

# Step 2: 使用 novoalign 建立参考序列索引
echo "Building novoalign index for ${REF_FASTA}..."
novoindex "$REF_INDEX" "$REF_FASTA"

# Step 3: 将合并后的序列比对到参考序列上
echo "Aligning reads with novoalign..."
novoalign -d "$REF_INDEX" -f merged.fq -o sam > map.sam

# Step 4: 使用 samtools 过滤未比对上的 reads (保留 mapped 成功的 reads)
echo "Filtering mapped reads with samtools..."
samtools view -F4 map.sam -h -o mapped.sam

# Step 5: 运行自定义工具分析剪切活性与统计结果
echo "Running castool analysis for target sequence..."
java -jar ~/Biotools.jar castool -table -cutref "$CUT_REF_SEQ" \
     -I mapped.sam \
     -O result.txt \
     -O2 result_stat.txt

echo "Pipeline finished successfully!"
