# %%
import scanpy as sc
from anndata import AnnData
# %%
Hsal_me = sc.read_h5ad("./cluster_unannotated.h5ad")
# %%
Hsal_me.X = Hsal_me.layers["logcounts"].copy()
# %%
sc.tl.rank_genes_groups(
    Hsal_me, 
    groupby="leiden_res1.00", 
    method="wilcoxon", 
    key_added="leiden_res1.00_DEG"
)
# %%
sc.tl.dendrogram( # 纯粹为了美观和生物学逻辑。它的作用是让后续画图时，相似的 cluster 能够挨在一起，并在图上方画出一个树状分支
    Hsal_me,
    groupby="leiden_res1.00",
)
# %%
sc.pl.rank_genes_groups_dotplot(
    Hsal_me, 
    groupby="leiden_res1.00", 
    standard_scale="var", 
    n_genes=5, 
    key="leiden_res1.00_DEG"
)
# %%
sc.tl.filter_rank_genes_groups(
    Hsal_me,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="leiden_res1.00_DEG",
    key_added="leiden_res1.00_DEG_filtered",
)
# %%
sc.pl.rank_genes_groups_dotplot(
    Hsal_me,
    groupby="leiden_res1.00",
    standard_scale="var",
    n_genes=5,
    key="leiden_res1.00_DEG_filtered",
)
# %%
# Hsal_me_deg_filtered_df = sc.get.rank_genes_groups_df(
#     Hsal_me,
#     group=None, 
#     key="leiden_res1.00_DEG_filtered"
# )

# Hsal_50_filtered_0 = get_specific_cluster_markers(Hsal_50_deg_filtered_df, "0-lKCB")

def get_specific_cluster_markers(
        adata: AnnData, 
        key: str,
        cluster: str
    ):
    df = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    not_na = df[~df["names"].isna()]
    specific_cluster = not_na[not_na["group"] == cluster]
    return specific_cluster
# %%
cluster1_markers = get_specific_cluster_markers(Hsal_me, "leiden_res1.00_DEG_filtered", "1")
# %%
cluster1_markers
# %%
cluster1_markers.to_csv("./markers/Hsal_me_filtered_1.csv")
# %%
sc.pl.embedding(
        Hsal_me,
        basis="X_umap",
        cmap="Reds",
        color=["LOC105184491", "LOC105190891", "LOC105189409", "LOC105184450", "LOC105186139", "leiden_res1.00"],
        legend_loc="on data"
    )
# %%
markers = cluster1_markers["names"].to_list()
# %%
sc.pl.embedding(
        Hsal_me,
        basis="X_umap",
        cmap="Reds",
        color=markers,
        legend_loc="on data"
    )
# %%
cluster1_markers_unfiltered = get_specific_cluster_markers(Hsal_me, "leiden_res1.00_DEG", "1")
# %%
cluster1_markers_unfiltered.to_csv("./markers/Hsal_me_1.csv")
# %%
unfilt_markers = cluster1_markers_unfiltered["names"].to_list()
# %%
unfilt_markers
# %%
sc.pl.embedding(
        Hsal_me,
        basis="X_umap",
        color=["log1p_n_genes_by_counts", "total_counts"],
        # legend_loc="on data"
    )
# %%
Hsal_50 = sc.read_h5ad("./Hsal50_umap_leiden.h5ad")
# %%
sc.pl.embedding(
        Hsal_50,
        basis="X_umap",
        cmap="Reds",
        color=markers,
        legend_loc="on data"
    )
# %%
#? 如果按照原文的“基因数超过200且UMI总数超过500”对该簇进行质控，可以滤去该簇多少细胞?
cluster1 = Hsal_me.obs[Hsal_me.obs["leiden_res1.00"] == "1"]
cluster1 # 1871个
# %%
cluster1[(cluster1["n_genes_by_counts"] > 200) & (cluster1["total_counts"] > 500)]
# 159个
# %%
# test markers:
sc.pl.embedding(
        Hsal_me,
        basis="X_umap",
        cmap="Reds",
        color=["LOC105187517", "LOC105184531", "LOC105187469"]
    )
# %%
