# %%
import scanpy as sc
from anndata._io.specs.registry import _REGISTRY, IOSpec
from anndata._io.specs.methods import H5Array
import anndata as ad

ad.settings.allow_write_nullable_strings = True
# %%
_REGISTRY.register_read(H5Array, IOSpec("null", "0.1.0"))(lambda elem, _reader: None)
# %%
Acer_intergrated = sc.read_h5ad("./matrix/acer_all_integration.h5ad")
# %%
Acer_intergrated
# %%
#! 降维
sc.pp.neighbors(Acer_intergrated, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(Acer_intergrated)
# %%
sc.set_figure_params(figsize=[6, 6], dpi_save=300)
fig = sc.pl.umap(Acer_intergrated, color=["batch"], wspace=1, return_fig=True)
# %%
#! 聚类
resolutions = [0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20]
for resolution in resolutions:
    sc.tl.leiden(
        Acer_intergrated, 
        key_added=f"leiden_neo_res_{resolution:.2f}", 
        resolution=resolution, 
        flavor="leidenalg", 
        n_iterations=-1,
        random_state=0
    )
# %%
Acer_intergrated.write_h5ad("./matrix/01_Acer_clustered.h5ad", compression="gzip")
# %%
