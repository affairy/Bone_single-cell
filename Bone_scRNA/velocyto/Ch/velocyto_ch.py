import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import re
import loompy
import matplotlib.pyplot as plt
import seaborn as sns

# 整合多个loom文件
import loompy
files=['/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/HAC1_1130/HAC1_1130_Aligned_APCMJ.loom',
       '/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/HAC2_1226/HAC2_1226_Aligned_OVM6V.loom',
       '/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/HAC3_0228/HAC3_0228_Aligned_6721U.loom',
       '/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/HAC_0922/HAC_0922_Aligned_LOZEE.loom',
       '/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/HAC_0928/HAC_0928_Aligned_VVQ1N.loom',
       '/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/OA1_1125/OA1_1125_Aligned_0BP16.loom',
       '/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/OA2_1127/OA2_1127_Aligned_X0GE0.loom',
       '/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/OA3_1130/OA3_1130_Aligned_1KJX2.loom']
output_filename='/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Ch/combined_Ch.loom'
loompy.combine(files, output_filename, key="Accession")


loom_Ch = sc.read('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Ch/combined_Ch.loom', cache=False)
loom_Ch.obs.head()
loom_Ch.obs = loom_Ch.obs.rename(index=lambda x: re.sub(r'Aligned_.*?:', '', x).replace(':', '_').replace('x', ''))
loom_Ch.obs.head()

loom_Ch.write('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Ch/combined_Ch_modified.loom')

adata = sc.read_h5ad("/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Ch/adata_Ch.h5ad")

####需要修改一下adata的barcode，HAC_0922.HAC_0928的barcode不一致


adata1 = scv.utils.merge(adata, loom_Ch)

import matplotlib.pyplot as plt
import scvelo as scv
import seaborn as sns

# 设置颜色列表
allcolour = ["#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d","#4dae47","#5c9e43" ]

scv.pl.proportions(adata1, groupby='celltype')

# 保存图片
plt.savefig('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Ch/proportions_plot.png', dpi=300, bbox_inches='tight')


scv.pp.filter_and_normalize(adata1)
scv.pp.moments(adata1)

scv.tl.velocity(adata1, mode='stochastic')
scv.tl.velocity_graph(adata1)
adata1
#adata1.obsm['X_UMAP'] = adata1.obsm['X_UMAP'].values  # 转换为 numpy 数组


# 绘制图像，直接将颜色列表传递给 palette
scv.pl.velocity_embedding_grid(adata1, basis='UMAP', color='celltype', 
                               scale=0.25, palette=allcolour)
# 保存图像
plt.savefig('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Ch/umap_velocyto_plot.png', dpi=300, bbox_inches='tight')

scv.tl.velocity_pseudotime(adata1)  # 由0到1
scv.pl.scatter(adata1, color='velocity_pseudotime', cmap='gnuplot', basis='UMAP')
plt.savefig('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Ch/umap_velocity_pseudotime_plot.png', dpi=300, bbox_inches='tight')

adata1.write('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Ch/adata1_Sy.h5ad')

adata1 = sc.read_h5ad("/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Ch/adata1_Sy.h5ad")
scv.tl.paga(adata1, groups='celltype')
scv.pl.paga(adata1, basis='UMAP', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, palette=allcolour)
plt.savefig('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Ch/umap_velocity_plot2.png', dpi=300, bbox_inches='tight')






