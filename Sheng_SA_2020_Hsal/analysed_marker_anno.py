# %%
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
# %%
Dmel_Hsal_markers = pd.read_csv("./markers/Hsal_Dmel_allcluster_markers.csv")
Dmel_Hsal_markers
# %%
#! 预计将进行双线注释, 分别对 Hsal_me 和 Hsal_50 分别进行可视化
Hsal_me = sc.read_h5ad("./cluster_unannotated.h5ad")
Hsal_50 = sc.read_h5ad("./Hsal50_umap_leiden.h5ad")
# %%
"""
基于标记基因的注释方法有两种。
一种方法是，从包含所有预期细胞类型的标记基因表入手，检查这些基因在哪些细胞簇中表达。
另一种方法是，检查哪些基因在已定义的细胞簇中高表达，然后检查它们是否与已知的细胞类型或状态相关。
如有必要，可以在这两种方法之间来回切换。
"""
# %%
marker_genes = {
    "neuron": ["LOC105189534", "LOC105190174", "LOC105183410"],
    "glia": ["LOC105181500", "LOC105186302"],
    "lKC": ["LOC105187517", "LOC105185095"],
    "DAN": ["LOC105182098", "LOC105187915", "LOC105189735"],
    "Astro glia": ["LOC105186780", 
            "LOC105188823", 
            "LOC105190201", 
            "LOC105184994", "LOC105184266", "LOC105184043"],
    "Ensheathing glia": ["LOC105186058", 
                    "LOC105186184", "LOC105192665",
                    "LOC105182483",
                    "LOC105183203"],
    "Perineurial glia": ["LOC105186476", "LOC105191060", "LOC105183245"],
    "Cortex glia": ["LOC105186657"],
    "Insulin-producing cell": ["LOC105188195"],
    "Hemocytes": ["LOC105182140", "LOC105188417"]
}
# %%
# sum = 0
# for key, value in marker_genes.items():
#     sum+=len(value)
# print(sum)
# %%
#! Hsal_me
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in Hsal_me.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
# %% 
# https://www.sc-best-practices.org/cellular_structure/annotation.html
Hsal_me.X = Hsal_me.layers["logcounts"]
sc.tl.pca(Hsal_me, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(Hsal_me)
sc.tl.umap(Hsal_me)
# %%
#! --------------------------------- Manual 1: neuron or glia? ---------------------------------
# * leiden 1.00
sc.pl.umap(
    Hsal_me,
    color="leiden_res1.00",
    vmin=0,
    vmax="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
    sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
    frameon=False,
    cmap="Reds",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    legend_loc="on data"
)
# %%
for ct in ["neuron", "glia"]:
    sc.pl.umap(
        Hsal_me,
        color=marker_genes_in_data[ct],
        vmin=0,
        vmax="p99", 
        sort_order=False,
        frameon=False,
        cmap="Reds"
    )
# %%
#! --------------------------------- Manual 2: neurons ---------------------------------
for ct in ["lKC", "DAN"]:
    sc.pl.umap(
        Hsal_me,
        color=marker_genes_in_data[ct],
        vmin=0,
        vmax="p99", 
        sort_order=False,
        frameon=False,
        cmap="Reds"
    )
# %%
#! --------------------------------- Manual 3: glias ---------------------------------
for ct in ["Astro glia", "Ensheathing glia", "Perineurial glia", "Cortex glia",]:
    sc.pl.umap(
        Hsal_me,
        color=marker_genes_in_data[ct],
        vmin=0,
        vmax="p99", 
        sort_order=False,
        frameon=False,
        cmap="Reds"
    )
# %%
