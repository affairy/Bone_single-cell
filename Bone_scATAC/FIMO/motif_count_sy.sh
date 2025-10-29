#!/bin/bash

# 设置相关目录
FIMO_DIR="/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo"
OUTPUT_FILE="/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo/all_up_motif.tsv"
COUNT_FILE="/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo/motif_count.txt"

# 初始化输出文件（添加表头）
echo -e "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence" > "$OUTPUT_FILE"

# 遍历所有 fimo.tsv 文件
find "$FIMO_DIR" -mindepth 1 -maxdepth 1 -type d | while read gene_dir; do
    fimo_file="$gene_dir/fimo.tsv"
    if [[ -f "$fimo_file" ]]; then
        awk 'NR > 1 && $8 < 1e-4 {print}' "$fimo_file" >> "$OUTPUT_FILE"
    fi
done

echo "所有上调 peak 相关 motif 提取完成，结果保存在 $OUTPUT_FILE"
 
# 统计 motif 频率
cut -f2 "$OUTPUT_FILE" | tail -n +2 | sort | uniq -c | sort -nr > "$COUNT_FILE"

echo "motif 频率统计完成，结果保存在 $COUNT_FILE"
