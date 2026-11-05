# %%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
# %%
Acer_transfromHY = sc.read_h5ad("./matrix/02_Acer_transfromHY.h5ad")
Acer_transfromHY.X = Acer_transfromHY.layers["logcounts"].copy()
# %%
# 见 .archive/03_manual_anno_NC_0.25.py
marker_genes = { 
    "KC": ["LOC107995382", "LOC107997697", "LOC107998479", "LOC107997070", "LOC108000569", "LOC107999048", "LOC107996358", "LOC107999991"],
    "glia": ["LOC107993519", "LOC108002812", "LOC108000234", "LOC107994035", "LOC107996595", "LOC108003667"],
    "OLC": ["LOC107996886", "LOC107999967", "LOC108002029", "LOC108004344", "LOC107999009", "LOC107996386", "LOC107998009", "LOC108000705", "LOC107992410"],
    "OPN": ["LOC108000694"]
}
# %%
sc.pl.dotplot(
    Acer_transfromHY,
    groupby="leiden_neo_res_0.50",
    var_names=marker_genes,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
# %%
