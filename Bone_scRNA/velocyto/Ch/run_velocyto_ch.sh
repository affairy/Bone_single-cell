#!/bin/bash

# 定义输入路径和输出路径
BAM_DIR="/mnt/data4/yiyuan/Bone/scRNA/velocyto"
BARCODE_DIR="/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/velocyto"
OUTPUT_DIR="/mnt/data4/yiyuan/Bone/scRNA/velocyto"
MASK_GTF="/mnt/data4/yiyuan/Bone/scRNA/velocyto/hg38_rmsk.gtf"
GENCODE_GTF="/home/yiyuan/reference/hg38/gencode/gencode.v38.annotation.gtf"

# 定义样本名列表
SAMPLES=("HAC_0922" "HAC_0928" "HAC1_1130" "HAC2_1226" "HAC3_0228" "OA1_1125" "OA2_1127" "OA3_1130")

# 循环处理每个样本
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing $SAMPLE..."
    
    velocyto run -b "${BARCODE_DIR}/barcode_${SAMPLE}.tsv" \
                 -o "${OUTPUT_DIR}/${SAMPLE}" \
                 -m "${MASK_GTF}" \
                 "${BAM_DIR}/${SAMPLE}/bam/${SAMPLE}_Aligned.sortedByCoord.out.bam" \
                 "${GENCODE_GTF}"
    
    if [ $? -eq 0 ]; then
        echo "$SAMPLE completed successfully."
    else
        echo "Error processing $SAMPLE."
    fi
done

echo "All jobs finished!"
