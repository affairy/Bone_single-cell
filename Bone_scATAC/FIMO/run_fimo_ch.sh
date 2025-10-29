#!/bin/bash

# 设置相关目录
FASTA_DIR="/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo_Ch/getfasta"
OUTPUT_DIR="/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo_Ch"
FIMO_DB="/mnt/3w/data2/yiyuan/ctDNA/result/SCLC/FIMO/JASPAR/20240403045904_JASPAR2024_combined_matrices_445470_meme.meme"

# 确保输出目录存在
mkdir -p "$OUTPUT_DIR"

# 遍历所有 .fa 文件并运行 fimo
for fa_file in "$FASTA_DIR"/*.fa; do
    # 获取文件名（不带路径）
    fa_name=$(basename "$fa_file" .fa)
    
    # 运行 fimo 命令（不会后台运行）
    fimo --thresh 1e-4 --o "$OUTPUT_DIR/$fa_name" "$FIMO_DB" "$fa_file"
    
    # 输出当前处理的文件名
    echo "Processed: $fa_name.fa"
done


