# %%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
# %%
Acer_transfromHY = sc.read_h5ad("./matrix/02_Acer_transfromHY.h5ad")
# %%
#! 由于先前发现使用条形图鉴定的准确度欠佳，故本次统一使用0.5分辨率，通过iScience中的marker，以点图形式鉴定
# %%
Acer_transfromHY.X = Acer_transfromHY.layers["logcounts"].copy()
# %%
#! 全脑级别
# 见 .archive/02_manual_anno_iScience.py
marker_genes = {
    "neuron": ["LOC108004334"],
    "KC": ["LOC108002506", "LOC107999048"],
    "glia": ["LOC107993336", "LOC108001026"],
    "OPN": ["LOC108002036", "LOC108001977", "LOC107999745", "LOC108000694"], # 其中 LOC107999745 效应最明显
    "hemo": ["LOC108004139", "LOC108001681"]
}
sc.pl.dotplot(
    Acer_transfromHY,
    groupby="leiden_neo_res_0.50",
    var_names=marker_genes,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
# %%
#! KC亚型
KC_marker_genes = {
    "global": ["LOC107995382"],
    "lKC": ["LOC108002263"],
    "sKC": ["LOC107997067"],
    "FoxP": ["LOC108000201"]
}
sc.pl.dotplot(
    Acer_transfromHY,
    groupby="leiden_neo_res_0.50",
    var_names=KC_marker_genes,
    standard_scale="var",
)
# %%
#! OL亚型
# 见 .archive/02_manual_anno_iScience.py
OLC_marker_genes = {
    "Photoreceptor": ["LOC107993080", "LOC108003337"],
    "Mi1": ["LOC107999009"],
    "T1": ["LOC107995911"],
    "Lawf2": ["LOC107994307", "LOC107994946", "LOC107996247"],
    "PM": ["LOC107993239", "LOC108004330", "LOC108003229"],
    "Tm5c": ["LOC108004323"]
}
sc.pl.dotplot(
    Acer_transfromHY,
    groupby="leiden_neo_res_0.50",
    var_names=OLC_marker_genes,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
# %%
#! 神经递质亚型
# 见 .archive/02_manual_anno_iScience.py
marker_genes_neurotransmitter = {
    "Cholinergic": ["LOC107993140", "LOC107993179"],
    "Glutamatergic": ["LOC107994928"],
    "GABAergic": ["LOC107996223", "LOC107995375"],
    "Monoaminergic": ["LOC108000696"],
    "Dopaminergic": ["LOC107999422"],
    "Octopaminergic": ["LOC107997295"],
    "Serotonergic": ["LOC107994910"],
    "Histaminergic": ["LOC108003440"]
}
sc.pl.dotplot(
    Acer_transfromHY,
    groupby="leiden_neo_res_0.50",
    var_names=marker_genes_neurotransmitter,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
# %%
