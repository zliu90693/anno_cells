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
get_specific_cluster_markers(Hsal_me, "leiden_res1.00_DEG_filtered", "1")