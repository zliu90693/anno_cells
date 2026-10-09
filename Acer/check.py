# %%
import scanpy as sc
# %%
Acer_unanno = sc.read_h5ad("./matrix/cluster_unannotated.h5ad")
Acer_unanno
# %%
print(Acer_unanno.obs["batch"].value_counts())
# drone_rep1_singlet     13544
# drone_rep2_singlet     12619
# queen_rep1_singlet     11272
# worker_rep2_singlet    10612
# worker_rep1_singlet    10119
# %%
Acer_transfromHY = sc.read_h5ad("./matrix/01_Acer_transfromHY.h5ad")
print(Acer_transfromHY.obs["batch"].value_counts())
# drone_rep1_singlet     13544
# drone_rep2_singlet     12619
# queen_rep1_singlet     11272
# worker_rep2_singlet    10612
# worker_rep1_singlet    10119
# %%
sc.pl.umap(
    Acer_transfromHY,
    color="celltype_transfrom_HY",
    # legend_loc="on data"
)
# %%
sc.pl.umap(
    Acer_transfromHY,
    color="leiden_res0.25",
    # legend_loc="on data"
)
# %%
