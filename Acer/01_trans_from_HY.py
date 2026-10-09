# %%
import scanpy as sc
# %%
Acer = sc.read_h5ad("./matrix/cluster_unannotated.h5ad")
Acer_HY = sc.read_h5ad("./matrix/Acer_annobyHY.h5ad")
# %%
sc.pl.umap(
    Acer,
    color="leiden_res1.00",
    legend_loc="on data"
)
# %%
sc.pl.umap(
    Acer_HY,
    color=["leiden", "cell_type"],
    legend_loc="on data"
)
# %%
Acer.obs_names
# %%
Acer_HY.obs_names
# %%
Acer.obs["HY_barcode"] = (
    Acer.obs_names
    .to_series()
    .replace({
        "_singlet": ""
    }, regex=True)
)
# %%
Acer_set = set(Acer.obs["HY_barcode"].to_list())
Acer_HY_set = set(Acer_HY.obs_names.to_list())
# %%
print(len(Acer_set)) # 58166
print(len(Acer_HY_set)) # 26805
print(len(Acer_set & Acer_HY_set)) # 26564
# %%
index_n_celltype = Acer_HY.obs["cell_type"]
index_n_celltype
# %%
Acer.obs["celltype_transfrom_HY"] = (
    Acer.obs["HY_barcode"]
    .map(index_n_celltype)
    .cat.add_categories("unknown")
    .fillna("unknown")
)
# %%
sc.pl.umap(
    Acer,
    color="celltype_transfrom_HY",
    # legend_loc="on data"
)
# %%
sc.pl.umap(
    Acer,
    color="celltype_transfrom_HY",
    groups=["unknown"]
)
# %%
Acer.write_h5ad("./matrix/01_Acer_transfromHY.h5ad")
# %%
