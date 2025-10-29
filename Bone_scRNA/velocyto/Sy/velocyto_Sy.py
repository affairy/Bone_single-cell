# 整合多个loom文件
import loompy
files=['/mnt/data4/yiyuan/Bone/scRNA/velocyto/HAS1_1130/HAS1_1130_Aligned_5KPKT.loom',
       '/mnt/data4/yiyuan/Bone/scRNA/velocyto/HAS2_1226/HAS2_1226_Aligned_9Z9S1.loom',
       '/mnt/data4/yiyuan/Bone/scRNA/velocyto/HAS3_0228/HAS3_0228_Aligned_SWVDT.loom',
       '/mnt/data4/yiyuan/Bone/scRNA/velocyto/OAS1_1209/OAS1_1209_Aligned_PJHYD.loom',
       '/mnt/data4/yiyuan/Bone/scRNA/velocyto/OAS2_1209/OAS2_1209_Aligned_6425S.loom',
       '/mnt/data4/yiyuan/Bone/scRNA/velocyto/OAS3_0108/OAS3_0108_Aligned_P0C00.loom']
output_filename='/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Sy/combined_Sy.loom'
loompy.combine(files, output_filename, key="Accession")

adata = sc.read_h5ad("/mnt/data4/yiyuan/Bone/scRNA/velocyto/Sy/adata_Sy.h5ad")

adata1 = scv.utils.merge(adata, loom_Sy)

import matplotlib.pyplot as plt
import scvelo as scv

# 定义你想要的颜色列表
allcolour = ['#DBC9B3','#EED0E0','#EBAEA9','#CB95BB','#EBCC96','#AED0DF','#CBE5DE']

# 绘制 UMAP 图，使用自定义的颜色
sc.pl.embedding(adata, basis='UMAP', color='celltype', 
                 legend_loc='on data', show=False, 
                 palette=allcolour)  # 这里使用 palette 参数来设置颜色

# 保存图像
plt.savefig('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Sy/umap_plot2.png', dpi=300, bbox_inches='tight')

scv.pl.proportions(adata1, groupby='celltype')

# 保存图片
plt.savefig('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Sy/proportions_plot.png', dpi=300, bbox_inches='tight')


scv.pp.filter_and_normalize(adata1)
scv.pp.moments(adata1)

scv.tl.velocity(adata1, mode='stochastic')
scv.tl.velocity_graph(adata1)
adata1
#adata1.obsm['X_UMAP'] = adata1.obsm['X_UMAP'].values  # 转换为 numpy 数组

import matplotlib.pyplot as plt
import scvelo as scv
import seaborn as sns

# 设置颜色列表
allcolour = ['#DBC9B3','#EED0E0','#EBAEA9','#CB95BB','#EBCC96','#AED0DF','#CBE5DE']

# 绘制图像，直接将颜色列表传递给 palette
scv.pl.velocity_embedding_grid(adata1, basis='UMAP', color='celltype', 
                               scale=0.25, palette=allcolour)
# 保存图像
plt.savefig('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Sy/umap_velocyto_plot.png', dpi=300, bbox_inches='tight')

scv.tl.velocity_pseudotime(adata1)  # 由0到1
scv.pl.scatter(adata1, color='velocity_pseudotime', cmap='gnuplot', basis='UMAP')
plt.savefig('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Sy/umap_velocity_pseudotime_plot.png', dpi=300, bbox_inches='tight')

adata1.write('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Sy/adata1_Sy.h5ad')

adata1 = sc.read_h5ad("/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Sy/adata1_Sy.h5ad")
scv.tl.paga(adata1, groups='celltype')
scv.pl.paga(adata1, basis='UMAP', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, palette=allcolour)
plt.savefig('/mnt/8w/data4/yiyuan/Bone/scRNA/velocyto/Sy/umap_velocity_plot2.png', dpi=300, bbox_inches='tight')

