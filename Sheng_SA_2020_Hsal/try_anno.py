# %%
import scanpy as sc
# %%
combined_h5ad = sc.read_h5ad("./cluster_unannotated.h5ad")
# %%
combined_h5ad.X = combined_h5ad.layers['logcounts'].copy() # sc.pl.umap 建议使用 logcounts 而非 raw 数据进行上色
# %%
sc.pl.umap(
    combined_h5ad, 
    color=['MT-nad4l', 'MT-nad4', 'MT-cox1', 'MT-cox2'], 
    # layers='logcounts',  # 【核心】指定读取 logcounts 层的数据进行上色
    cmap='Reds',         # 推荐使用 Reds 或 Viridis 颜色映射
    wspace=0.4
)
# %%
sc.pl.umap(
    combined_h5ad, 
    color=['LOC105186768', 'LOC105180849', 'LOC105186803', 'LOC105188139', "LOC105186752"], 
    cmap='Reds', 
    wspace=0.4
)
# %%
sc.pl.umap(
    combined_h5ad, 
    legend_loc="on data", 
    color=f'leiden_res0.50')
# %%
sc.tl.rank_genes_groups(
    combined_h5ad, 
    groupby="leiden_res0.50", 
    method="wilcoxon", 
    key_added="dea_leiden_res0.50"
)

# %%
sc.tl.dendrogram(
    combined_h5ad,
    groupby="leiden_res0.50",
)
# %%
sc.pl.rank_genes_groups_dotplot(
    combined_h5ad, 
    groupby="leiden_res0.50", 
    standard_scale="var", 
    n_genes=10, 
    key="dea_leiden_res0.50"
)
# %%
sc.tl.filter_rank_genes_groups(
    combined_h5ad,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="dea_leiden_res0.50",
    key_added="dea_leiden_res0.50_filtered",
)
# %%
sc.pl.rank_genes_groups_dotplot(
    combined_h5ad,
    groupby="leiden_res0.50",
    standard_scale="var",
    n_genes=5,
    key="dea_leiden_res0.50_filtered",
)
# %%
sc.pl.umap(
    combined_h5ad, 
    color=['LOC105182601', 'LOC105184919', 'LOC105185107', 'LOC105188896', "LOC105184491"], 
    cmap='Reds',
    wspace=0.4
)
# %%
