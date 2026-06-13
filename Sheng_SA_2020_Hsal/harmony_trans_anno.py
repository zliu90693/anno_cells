# %%
import scanpy as sc
import harmonypy as hm
import numpy as np
# %%
Hsal_50 = sc.read_h5ad("./Hsal50.h5ad")
Hsal_me = sc.read_h5ad("./cluster_unannotated.h5ad")
# %%
Hsal_50.X = Hsal_50.layers["counts"]
Hsal_me.X = Hsal_me.layers["counts"]
# %%
del Hsal_me.uns
del Hsal_me.obsm
del Hsal_me.varm
del Hsal_me.layers
del Hsal_me.obsp

del Hsal_50.layers
del Hsal_50.obsm
# %%
Hsal_50
# %%
Hsal_me
# %%
# 检测二者共有多少个基因: 
print(
    "Number of genes found in query dataset:",
    Hsal_me.var.index.isin(Hsal_50.var.index).sum(),
)
# 9731
"""
注意, 现在我担心一个问题, 
Hsal_me中有9814个基因, Hsal_50中有11926个基因, 二者共有9731个基因,
但是, 其中有部分基因并非以LOC开头, 如下:
"""
Hsal_50.var.index[~Hsal_50.var.index.str.startswith("LOC")]
"""
我担心的问题是, 可能Hsal_50和Hsal_me对应的参考基因组版本不一致, 而Hsal_50中部分非LOC开头的基因, 在Hal_me中对应的同源基因可能正是以LOC开头的,
因此, 基于名称的比对很可能导致'二者各自独有的基因中, 很可能有部分是同一基因, 只是因为名称不匹配没出现在那9731中'
验证, Hsal_50中全部不以LOC开头的基因:
"""
Hsal_50.var.index[~Hsal_50.var.index.str.startswith("LOC")]
# Index(['Anamorsin', 'Centrin-1', 'Fut8', 'Glomulin', 'Ickin1-1', 'Importin-4',
#        'Traa', 'Trab', 'Trehalase', 'Trnag-ccc', 'Trypsin-1'],
#       dtype='object')
# 11个
"""
9731个基因中不以LOC开头的基因:
"""
9731 - sum([gene.startswith("LOC") for gene in list(set(Hsal_me.var_names) & set(Hsal_50.var_names))])
# 9723 个
"""
也就是说, 9731个基因中有8个不以LOC开头的基因, 剩下3个不以LOC开头的基因想必是Hsal_50特有的,
因此, 这种规模(3个)的不匹配或许没什么问题.
"""
# %%
#! 仅保留二者共有的基因:
# 取共有基因
common_genes = Hsal_50.var_names.intersection(Hsal_me.var_names)
Hsal_50 = Hsal_50[:, common_genes].copy()
Hsal_me = Hsal_me[:, common_genes].copy()
# %%
#! 处理两个数据集中.var和.obs的各列
Hsal_me.var = Hsal_me.var[[]]
Hsal_50.var = Hsal_50.var[[]]
# %%
Hsal_me.obs = Hsal_me.obs[["leiden_res0.25", "leiden_res0.50", "leiden_res1.00"]]
Hsal_50.obs = Hsal_50.obs[["cluster_name"]]
# %%
Hsal_me.obs["batch_key"] = "Hsal_me"
Hsal_50.obs["batch_key"] = "Hsal_50"
# %%
Hsal_me.obs["cluster_name"] = np.nan
Hsal_50.obs[["leiden_res0.25", "leiden_res0.50", "leiden_res1.00"]] = np.nan
# %%
Hsal_me.obs = Hsal_me.obs[["batch_key", "leiden_res0.25", "leiden_res0.50", "leiden_res1.00", "cluster_name"]]
Hsal_50.obs = Hsal_50.obs[["batch_key", "leiden_res0.25", "leiden_res0.50", "leiden_res1.00", "cluster_name"]]
# %%
#! 合并
Hsal_combined = sc.concat(
    [Hsal_50, Hsal_me],
    label='dataset',         # 自动在obs中添加'dataset'列，值为'Hsal_50'/'Hsal_me'
    keys=['Hsal_50', 'Hsal_me'],
    join='inner',            # 只保留共有基因（对应前面取交集的操作）
    merge='same'             # var中相同的列才保留
)
# %%
#! downstream analyse
Hsal_combined.layers["counts"] = Hsal_combined.X.copy()
# %%
Hsal_combined.write_h5ad("./Hsal_combined.h5ad")
# %%
sc.pp.normalize_total(Hsal_combined)
sc.pp.log1p(Hsal_combined)
Hsal_combined.layers['logcounts'] = Hsal_combined.X.copy()
# %%
sc.pp.highly_variable_genes(Hsal_combined, n_top_genes=4000, batch_key='batch_key')
# %%
sc.tl.pca(Hsal_combined, use_highly_variable=True)
# %%
# From https://github.com/slowkow/harmonypy/issues/49
pcs = Hsal_combined.obsm["X_pca"]
harmony_out =  hm.run_harmony(pcs, Hsal_combined.obs, "batch_key")
Hsal_combined.obsm["X_pca_harmony"] = harmony_out.Z_corr
# %%
sc.pp.neighbors(Hsal_combined, use_rep='X_pca_harmony')
sc.tl.umap(Hsal_combined)
# sc.tl.leiden(Hsal_combined)
# %%
for reso in [0.25, 0.5, 1.0]:
    # run_leiden(combined_h5ad, neo_key="harmony_scran", neighbors_key="neighbors_harmony_scran", resolution=reso, flavor="igraph", n_iterations=2)
    sc.tl.leiden(
    Hsal_combined, 
    key_added=f"neo_leiden_res_{reso:.2f}",  # 防止浮点数转字符串时可能出现 res0.6000000000000001
    resolution=reso
)
# %%
#! visualize
sc.pl.umap(Hsal_combined, color=["batch_key", "leiden_res1.00", "cluster_name", "neo_leiden_res_1.00"])
# %%
sc.pl.umap(Hsal_combined, 
        color=["cluster_name"],
        legend_loc="on data"
        )
# %%
sc.pl.umap(Hsal_combined, 
        color=["leiden_res1.00"],
        legend_loc="on data")
# %%
sc.pl.umap(Hsal_combined, 
        color=["neo_leiden_res_1.00"],
        legend_loc="on data")
# %%
# sc.pl.pca(Hsal_combined, 
#         color=["batch_key"])
sc.pl.embedding(
        Hsal_combined,
        basis="X_pca",
        color="batch_key"
    )
sc.pl.embedding(
        Hsal_combined,
        basis="X_pca_harmony",
        color="batch_key"
    )
# %%
Hsal_combined.write_h5ad("./Hsal_me_50_harmony.h5ad")
# %%
Hsal_combined = sc.read_h5ad("./Hsal_me_50_harmony.h5ad")
# %%
