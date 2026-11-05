# %%
import scanpy as sc
import anndata as ad
ad.settings.allow_write_nullable_strings = True
# %%
Acer_clustered = sc.read_h5ad("./matrix/01_Acer_clustered.h5ad")
Acer_HY = sc.read_h5ad("./matrix/Acer_annobyHY.h5ad")
# %%
sc.pl.umap(
    Acer_clustered,
    color="leiden_neo_res_0.50",
    legend_loc="on data"
)
# %%
sc.pl.umap(
    Acer_HY,
    color=["leiden", "cell_type"],
    legend_loc="on data"
)
# %%
Acer_clustered.obs_names
# %%
Acer_HY.obs_names
# %%
Acer_clustered.obs["HY_barcode"] = (
    Acer_clustered.obs_names
    .to_series()
    .replace({
        "_singlet": ""
    }, regex=True)
)
# %%
Acer_set = set(Acer_clustered.obs["HY_barcode"].to_list())
Acer_HY_set = set(Acer_HY.obs_names.to_list())
# %%
print(len(Acer_set)) # 58166
print(len(Acer_HY_set)) # 26805
print(len(Acer_set & Acer_HY_set)) # 26564
# %%
index_n_celltype = Acer_HY.obs["cell_type"]
index_n_celltype
# %%
Acer_clustered.obs["celltype_transfrom_HY"] = (
    Acer_clustered.obs["HY_barcode"]
    .map(index_n_celltype)
    .cat.add_categories("unknown")
    .fillna("unknown")
)
# %%
sc.pl.umap(
    Acer_clustered,
    color="celltype_transfrom_HY",
    # legend_loc="on data"
)
# %%
Acer_clustered.write_h5ad("./matrix/02_Acer_transfromHY.h5ad")
# %%
