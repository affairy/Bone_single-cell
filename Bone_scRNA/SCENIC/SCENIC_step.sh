# 软骨
micromamba activate pyscenic
# step1 grn
nohup pyscenic grn \
  --num_workers 15 \
  --output /mnt/data4/yiyuan/Bone/scRNA/scenic/Ch/output/grn.tsv \
  --method grnboost2 \
  /mnt/data4/yiyuan/Bone/scRNA/scenic/Ch/scenic.loom \
  /home/yiyuan/reference/index_genome/cisTarget_databases/hs_hgnc_tfs.txt &


# step2 ctx
nohup pyscenic ctx \
  /mnt/data4/yiyuan/Bone/scRNA/scenic/Ch/output/grn.tsv /home/yiyuan/reference/index_genome/cisTarget_databases/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
  --annotations_fname /home/yiyuan/reference/index_genome/cisTarget_databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname /mnt/data4/yiyuan/Bone/scRNA/scenic/Ch/scenic.loom \
  --mode "dask_multiprocessing" \
  --output /mnt/data4/yiyuan/Bone/scRNA/scenic/Ch/output/ctx.csv \
  --num_workers 20 \
  --mask_dropouts &

micromamba activate pyscenic_env
# step3 AUCell
nohup pyscenic aucell \
  /mnt/data4/yiyuan/Bone/scRNA/scenic/Ch/scenic.loom \
  /mnt/data4/yiyuan/Bone/scRNA/scenic/Ch/output/ctx.csv \
  --output /mnt/data4/yiyuan/Bone/scRNA/scenic/Ch/output/aucell.loom \
  --num_workers 15 &


#####滑膜成纤维
  # step1 grn
nohup pyscenic grn \
  --num_workers 15 \
  --output /mnt/data4/yiyuan/Bone/scRNA/scenic/Fb/output/grn.tsv \
  --method grnboost2 \
  /mnt/data4/yiyuan/Bone/scRNA/scenic/Fb/scenic.loom \
  /home/yiyuan/reference/index_genome/cisTarget_databases/hs_hgnc_tfs.txt &


# step2 ctx
nohup pyscenic ctx \
  /mnt/data4/yiyuan/Bone/scRNA/scenic/Fb/output/grn.tsv /home/yiyuan/reference/index_genome/cisTarget_databases/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
  --annotations_fname /home/yiyuan/reference/index_genome/cisTarget_databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname /mnt/data4/yiyuan/Bone/scRNA/scenic/Fb/scenic.loom \
  --mode "dask_multiprocessing" \
  --output /mnt/data4/yiyuan/Bone/scRNA/scenic/Fb/output/ctx.csv \
  --num_workers 20 \
  --mask_dropouts &

micromamba activate pyscenic_env
# step3 AUCell
nohup pyscenic aucell \
  /mnt/data4/yiyuan/Bone/scRNA/scenic/Fb/scenic.loom \
  /mnt/data4/yiyuan/Bone/scRNA/scenic/Fb/output/ctx.csv \
  --output /mnt/data4/yiyuan/Bone/scRNA/scenic/Fb/output/aucell.loom \
  --num_workers 15 &