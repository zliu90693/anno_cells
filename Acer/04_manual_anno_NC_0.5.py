# %%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
from functions.anno_utils import *
# %%
Acer_transfromHY = sc.read_h5ad("./matrix/01_Acer_transfromHY.h5ad")
# Acer_transfromHY.X = Acer_transfromHY.layers["counts"].copy()
Acer_transfromHY.X = Acer_transfromHY.layers["logcounts"].copy()
# sc.pp.normalize_total(Acer_transfromHY)
# Acer_transfromHY.layers["norm_counts"] = Acer_transfromHY.X.copy()
# %%
marker_genes = {
    "KC": ["LOC107995382", "LOC107997697", "LOC107998479", "LOC107997070", "LOC108000569", "LOC107999048"],
    "glia": ["LOC107993519", "LOC108002812", "LOC108000234", "LOC107994035", "LOC107996595", "LOC108003667"],
    "OLC": ["LOC107996886", "LOC107999967", "LOC108002029", "LOC108004344", "LOC107999009", "LOC107996386", "LOC107998009", "LOC108000705", "LOC107992410"],
    "OPN": ["LOC108000694"]
}
# %%
sc.pl.dotplot(
    Acer_transfromHY,
    groupby="leiden_res0.50",
    var_names=marker_genes,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
# %%
sc.pl.umap(
    Acer_transfromHY,
    color="leiden_res0.50",
    legend_loc="on data"
)
# %%
KC_clusters_050_dot = {"KCs": ["0", "3", "8", "19", "22", "23", "26"], 
                   "glia": ["40"], 
                   "OLC": ["10", "11", "12", "14", "15", "17", "32", "39"],
                   "maybe_OPN": ["21"]}
# %%
cluster_to_celltype_050_dot = {
    cluster: celltype
    for celltype, clusters in KC_clusters_050_dot.items()
    for cluster in clusters
}
Acer_transfromHY.obs["celltype_050_NC_dot"] = (
    Acer_transfromHY.obs["leiden_res0.50"]
    .map(cluster_to_celltype_050_dot)
    .fillna("Unknown")
)
Acer_transfromHY.obs["celltype_050_NC_dot"].to_csv("./metadata/anno/celltype_050_NC_dot.csv")
# %%
