# %%
import pandas as pd
import scanpy as sc
# %%
Hsal_50 = sc.read_h5ad("./Hsal50.h5ad")
# %%
Hsal_50
# %%
sc.pl.embedding(
        Hsal_50,
        basis="tsne",
        color="cluster_name",
        legend_loc="on data"
    )
# %%
Hsal_50.X = Hsal_50.layers["counts"].copy()
# %%
# sc.pp.normalize_total(Hsal_50)
# sc.pp.log1p(Hsal_50)
# Hsal_50.layers['logcounts'] = Hsal_50.X.copy()
# # %%
# sc.tl.pca(Hsal_50)
# # %%
# sc.pp.neighbors(Hsal_50, use_rep='X_pca')
# sc.tl.umap(Hsal_50)
# # %%
# sc.pl.embedding(
#         Hsal_50,
#         basis="X_umap",
#         color="cluster_name",
#         # legend_loc="on data"
#     )
# %%
sc.pp.neighbors(Hsal_50, use_rep='pca')
sc.tl.umap(Hsal_50)
# %%
sc.pl.embedding(
        Hsal_50,
        basis="X_umap",
        color="cluster_name",
        legend_loc="on data"
    )
# %%
sc.pl.embedding(
        Hsal_50,
        basis="X_umap",
        color="cluster_name",
        # legend_loc="on data"
    )
# %%
for reso in [0.25, 0.50, 1.00]:
    sc.tl.leiden(
        Hsal_50,
        key_added=f"leiden_res{reso:.2f}",
        resolution=reso,
        flavor="leidenalg",
        n_iterations=2
    )
    fig = sc.pl.umap(Hsal_50, legend_loc="on data", color=f'leiden_res{reso:.2f}')
# %%
Hsal_50.write_h5ad("./Hsal50_umap_leiden.h5ad")


# %%
Hsal_me = sc.read_h5ad("./Hsal_me_with_5051_anno.h5ad")
Hsal_me
# %%
sc.pl.embedding(
        Hsal_me,
        basis="X_umap",
        color=["cell_type_50", "leiden_res1.00"],
        legend_loc="on data"
    )
# %%
