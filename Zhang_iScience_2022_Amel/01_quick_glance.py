# 尝试直接迁移HY提供的注释，大概有5000个细胞匹配不到
# %%
import scanpy as sc
# %%
Amel_adata = sc.read_h5ad("./matrix/cluster_unannotated.h5ad")
Amel_adataHY = sc.read_h5ad("./matrix/Amel_annobyHY.h5ad")
# %%
sc.pl.umap(
    Amel_adata,
    color="leiden_res1.00",
    legend_loc="on data"
)
# %%
sc.pl.umap(
    Amel_adataHY,
    color=["leiden", "cell_type"],
    legend_loc="on data"
)
# %%
# "F1 -> forager_rep1"
# "F2 -> forager_rep2"
# "N1 -> nurse_rep1"
# "N2 -> nurse_rep2"
# "Q1 -> queen_rep1"
# "Q2 -> queen_rep2"
# %%
Amel_adata.obs["HY_barcode"] = (
    Amel_adata.obs_names
    .to_series()
    .replace({
        "F1": "forager_rep1",
        "F2": "forager_rep2",
        "N1": "nurse_rep1",
        "N2": "nurse_rep2",
        "Q1": "queen_rep1",
        "Q2": "queen_rep2",
    }, regex=True)
)
# %%
Amel_adata_set = set(Amel_adata.obs["HY_barcode"].to_list())
Amel_adataHY_set = set(Amel_adataHY.obs_names.to_list())
# %%
print(len(Amel_adata_set))
print(len(Amel_adataHY_set))
print(len(Amel_adata_set & Amel_adataHY_set))
# %%
index_n_celltype = Amel_adataHY.obs["cell_type"]
index_n_celltype
# %%
Amel_adata.obs["celltype"] = (
    Amel_adata.obs["HY_barcode"]
    .map(index_n_celltype)
    .cat.add_categories("unknown")
    .fillna("unknown")
)
# %%
sc.pl.umap(
    Amel_adata,
    color="celltype",
    # legend_loc="on data"
)
# %%
sc.pl.umap(
    Amel_adata,
    color="celltype",
    groups=["unknown"]
)
# %%
Amel_adata.write_h5ad("./matrix/01_Amel_transfromHY.h5ad")
# %%
