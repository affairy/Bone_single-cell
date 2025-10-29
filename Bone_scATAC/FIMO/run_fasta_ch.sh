#!/bin/bash

# 设置输入和输出目录
BED_DIR="/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo_Ch"
OUTPUT_DIR="/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo_Ch/getfasta"
FASTA="/home/public/hg38/fasta/hg38.fa"

# 确保输出目录存在
mkdir -p "$OUTPUT_DIR"

# 遍历所有 BED 文件
for bed_file in "$BED_DIR"/*.bed; do
    # 获取文件名（不带路径）
    bed_name=$(basename "$bed_file" .bed)
    
    # 运行 bedtools getfasta
    bedtools getfasta -name+ -fi "$FASTA" -bed "$bed_file" -fo "$OUTPUT_DIR/$bed_name.fa"
    
    echo "Processed: $bed_name.bed -> $bed_name.fa"
done

