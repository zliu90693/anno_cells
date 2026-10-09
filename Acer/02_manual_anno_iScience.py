# %%
#! 将根据leiden_res0.25列的细胞簇，基于Zhang et al., 2022中提供的marker进行注释
# %%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
from functions.anno_utils import *
# %%
Acer_transfromHY = sc.read_h5ad("./matrix/01_Acer_transfromHY.h5ad")
# %%
Acer_transfromHY
# %%
#! 接下来将仿照Zhang et al., 2022 Figure 1中手动注释的部分，使用同源基因查看中蜂不同簇中的表达量
#! Figure 1的那种条形图，使用的数据集并非raw counts，也不是经过log1p normalize后的数据集，而是“经过library-size normalized，但尚未 log1p 的线性表达矩阵”
Acer_transfromHY.X = Acer_transfromHY.layers["counts"].copy()
sc.pp.normalize_total(Acer_transfromHY)
# %%
Acer_transfromHY.layers["norm_counts"] = Acer_transfromHY.X.copy()
Acer_transfromHY
# %%
# Acer_transfromHY.layers["counts"].max()
# %%
marker_genes = [
    "LOC108004334", # LOC410689，即elav同源基因，神经元普遍表达的基因
    # 图D
    "LOC108002506", # LOC408372，即mub同源基因
    "LOC107999048", # LOC408804，即plc同源基因
    "LOC107993336", # LOC410151，即repo同源基因
    "LOC108001026", # LOC408480，即ced6同源基因
    "LOC108002036", # LOC724282，即opt同源基因
    "LOC108001977", # LOC724148，即c15同源基因
    "LOC107999745", # LOC410657，即acj6同源基因
    "LOC108000694", # LOC413466，即oaz同源基因
    "LOC108004139", # LOC551684，即Fer2LCH同源基因
    "LOC108001681", # LOC411597，即Hml同源基因
]

marker_genes_OL = [
    # 图E
    "LOC107993080", # LOC409937，即ninaC同源基因
    "LOC108003337", # LOC724608，即trp同源基因
    "LOC107999009", # LOC724801，即bsh同源基因
    "LOC107995911", # LOC409919，即Eaat-1同源基因
    "LOC107994307", # LOC552079，即hth同源基因
    "LOC107994946", # LOC412528，即eya同源基因
    "LOC107996247", # LOC100577751，即Lim1同源基因
    "LOC107993239", # LOC410658，即Lim3同源基因
    "LOC108004330", # LOC411315，即vsx1同源基因
    "LOC108003229", # LOC408872，即svp同源基因
    "LOC108004323" # 应该是Hiscl1的同源基因，原文没提，在Flybase中确定Hiscl1的gid随后查询同源基因获得
]
# %%

# %%
#! 注释掉的内容已整理为函数形式
# ---------------------------------- 提取表达矩阵 ---------------------------------- 
# Acer_transfromHY_df = sc.get.obs_df(
#     Acer_transfromHY,
#     keys=[cluster_key] + marker_genes,
#     layer=layer
# )
# Acer_transfromHY_df
# # %%
# # ---------------------------------- 宽表转长表 ---------------------------------- 
# Acer_transfromHY_df_long = Acer_transfromHY_df.melt(
#     id_vars=cluster_key,
#     value_vars=marker_genes,
#     var_name="gene",
#     value_name="expression"
# )
# Acer_transfromHY_df_long
# # %%
# # ---------------------------------- 计算每个 cluster 的 mean / SD / N / SEM ---------------------------------- 
# stats = (
#     Acer_transfromHY_df_long
#     .groupby(
#         ["gene", cluster_key],
#         observed=True
#     )["expression"]
#     .agg(
#         mean="mean",
#         sd="std",
#         n="count"
#     )
#     .reset_index()
# )
# # %%
# stats["sem"] = (
#     stats["sd"] /
#     np.sqrt(stats["n"])
# )
# # %%
# # ---------------------------------- cluster 顺序 ---------------------------------- 
# if isinstance(
#     Acer_transfromHY.obs[cluster_key].dtype,
#     pd.CategoricalDtype
# ):
#     cluster_order = list(
#         Acer_transfromHY.obs[cluster_key].cat.categories
#     )
# else:
#     cluster_order = sorted(
#         Acer_transfromHY.obs[cluster_key].unique()
# )
# # %%
# # ---------------------------------- 获取 Scanpy cluster 颜色 ---------------------------------- 
# color_key = f"{cluster_key}_colors"

# if color_key in Acer_transfromHY.uns:
#     colors = Acer_transfromHY.uns[color_key]

#     cluster_colors = dict(
#         zip(
#             Acer_transfromHY.obs[cluster_key].cat.categories,
#             colors
#         )
#     )
# else:
#     cmap = plt.get_cmap("tab20")

#     cluster_colors = {
#         cluster: cmap(i % 20)
#         for i, cluster in enumerate(cluster_order)
#     }
# # %%
# # ---------------------------------- 绘图 ---------------------------------- 
# fig, axes = plt.subplots(
#     nrows=len(marker_genes),
#     ncols=1,
#     figsize=(12, 12),
#     sharex=True
# )
# x = np.arange(len(cluster_order))
# for ax, gene in zip(axes, marker_genes):
#     temp = stats[
#         stats["gene"] == gene
#     ].copy()
#     temp[cluster_key] = pd.Categorical(
#         temp[cluster_key],
#         categories=cluster_order,
#         ordered=True
#     )
#     temp = (
#         temp
#         .sort_values(cluster_key)
#         .set_index(cluster_key)
#         .reindex(cluster_order)
#     )
#     bar_colors = [
#         cluster_colors[c]
#         for c in cluster_order
#     ]
#     ax.bar(
#         x,
#         temp["mean"],
#         yerr=temp["sem"],
#         color=bar_colors,
#         width=0.75,
#         capsize=2,
#         edgecolor="none"
#     )
#     # gene 名称
#     ax.text(
#         0.02,
#         0.80,
#         gene,
#         transform=ax.transAxes,
#         fontsize=10,
#         fontstyle="italic",
#         fontweight="bold"
#     )
#     ax.set_ylabel("")
#     # 去掉上、右边框
#     ax.spines["top"].set_visible(False)
#     ax.spines["right"].set_visible(False)

# # 设置x轴
# axes[-1].set_xticks(x)
# axes[-1].set_xticklabels(
#     cluster_order,
#     fontsize=12
# )
# axes[-1].set_xlabel("Cell cluster")
# fig.supylabel(
#     "Average expression",
#     fontsize=13,
#     fontweight="bold"
# )
# plt.tight_layout()
# plt.show()
# %%
# def get_gene_cluster_table( # 从表达矩阵中提取一个仅包括指定gene和指定簇的子矩阵，格式为长表
#         adata: anndata.AnnData, 
#         cluster_key: str, # “细胞簇划分”在Anndata对象中对应的列名
#         gene_list: list, # marker gene列表
#         layer: str # 提取过程使用的layer
#     ) -> pd.DataFrame:
#     df = sc.get.obs_df(
#         adata,
#         keys=[cluster_key] + gene_list,
#         layer=layer
#     )
#     df_long = df.melt( # 宽表转长表
#         id_vars=cluster_key,
#         value_vars=gene_list,
#         var_name="gene",
#         value_name="expression"
#     )
#     return df_long

# def cal_mean_etal( # 计算每个 cluster 的 mean / SD / N / SEM 
#         df_long: pd.DataFrame,
#         cluster_key: str,
#     ) -> pd.DataFrame:
#     stats = (
#         df_long
#         .groupby(
#             ["gene", cluster_key],
#             observed=True
#         )["expression"]
#         .agg(
#             mean="mean",
#             sd="std",
#             n="count"
#         )
#         .reset_index()
#     )
#     stats["sem"] = (
#         stats["sd"] /
#         np.sqrt(stats["n"])
#     )
#     return stats

# def get_cluster_order( # cluster 顺序
#         adata: anndata.AnnData, 
#         cluster_key: str,
#     ) -> list:
#     if isinstance(
#         adata.obs[cluster_key].dtype,
#         pd.CategoricalDtype
#     ):
#         cluster_order = list(
#             adata.obs[cluster_key].cat.categories
#         )
#     else:
#         cluster_order = sorted(
#             adata.obs[cluster_key].unique()
#     )
#     return cluster_order

# def get_cluster_color(
#         adata: anndata.AnnData, 
#         cluster_key: str,
#         cluster_order: list
#     ) -> list:
#     color_key = f"{cluster_key}_colors"
#     if color_key in adata.uns:
#         colors = adata.uns[color_key]
#         cluster_colors = dict(
#             zip(
#                 adata.obs[cluster_key].cat.categories,
#                 colors
#             )
#         )
#     else:
#         cmap = plt.get_cmap("tab20")
#         cluster_colors = {
#             cluster: cmap(i % 20)
#             for i, cluster in enumerate(cluster_order)
#         }
#     return cluster_colors

# def visualize(
#         gene_list: list, # marker gene列表
#         cluster_key: str,
#         cluster_order: list,
#         stats: pd.DataFrame,
#         cluster_colors: list
#     ):
#     fig, axes = plt.subplots(
#         nrows=len(gene_list),
#         ncols=1,
#         figsize=(12, 12),
#         sharex=True
#     )
#     x = np.arange(len(cluster_order))
#     for ax, gene in zip(axes, gene_list):
#         temp = stats[
#             stats["gene"] == gene
#         ].copy()
#         temp[cluster_key] = pd.Categorical(
#             temp[cluster_key],
#             categories=cluster_order,
#             ordered=True
#         )
#         temp = (
#             temp
#             .sort_values(cluster_key)
#             .set_index(cluster_key)
#             .reindex(cluster_order)
#         )
#         bar_colors = [
#             cluster_colors[c]
#             for c in cluster_order
#         ]
#         ax.bar(
#             x,
#             temp["mean"],
#             yerr=temp["sem"],
#             color=bar_colors,
#             width=0.75,
#             capsize=2,
#             edgecolor="none"
#         )
#         # gene 名称
#         ax.text(
#             0.02,
#             0.80,
#             gene,
#             transform=ax.transAxes,
#             fontsize=10,
#             fontstyle="italic",
#             fontweight="bold"
#         )
#         ax.set_ylabel("")
#         # 去掉上、右边框
#         ax.spines["top"].set_visible(False)
#         ax.spines["right"].set_visible(False)

#     # 设置x轴
#     axes[-1].set_xticks(x)
#     axes[-1].set_xticklabels(
#         cluster_order,
#         fontsize=12
#     )
#     axes[-1].set_xlabel("Cell cluster")
#     fig.supylabel(
#         "Average expression",
#         fontsize=13,
#         fontweight="bold"
#     )
#     plt.tight_layout()
#     plt.show()
# %%
# -----------------------------------------------------------------------
# --------------------- 查看的是全部类型细胞marker表达情况 ------------------
# -----------------------------------------------------------------------
cluster_key="leiden_res0.25"
layer="norm_counts"
# %%
Acer_transfromHY_df_long = get_gene_cluster_table(
    adata=Acer_transfromHY,
    cluster_key=cluster_key,
    gene_list=marker_genes,
    layer=layer
)
stats = cal_mean_etal(
    df_long=Acer_transfromHY_df_long,
    cluster_key=cluster_key
)
cluster_order = get_cluster_order(
    adata=Acer_transfromHY,
    cluster_key=cluster_key
)
cluster_color = get_cluster_color(
    adata=Acer_transfromHY,
    cluster_key=cluster_key,
    cluster_order=cluster_order
)
visualize(
    gene_list=marker_genes,
    cluster_key=cluster_key,
    cluster_order=cluster_order,
    stats=stats,
    cluster_colors=cluster_color
)
# %%
sc.pl.umap(
    Acer_transfromHY,
    color="leiden_res0.25",
    legend_loc="on data"
)
sc.pl.umap(
    Acer_transfromHY,
    color="celltype_transfrom_HY",
    # legend_loc="on data"
)
# %%
#! 经过手动比较markers_expression.png，leiden0.25.png，和celltype_transfrom_HY.png，确定：
#! KC: 0, 3, 8, 20, 21

#! 其余均不确定
#! glia: 可能是34
#! OPN: 我的结果显示可能是 1, 4, 5, 6, 13 但HY的结果仅包括13的一小部分，待进一步确认
#! OLC暂未确定
#! Hemo同样无法确定，可能是4，也可能根本没有Hemo
# %%
# -----------------------------------------------------------------------
# --------------------- 查看的是 OLC 各亚型 marker表达情况 ------------------
# -----------------------------------------------------------------------
Acer_transfromHY_df_long = get_gene_cluster_table(
    adata=Acer_transfromHY,
    cluster_key=cluster_key,
    gene_list=marker_genes_OL,
    layer=layer
)
stats = cal_mean_etal(
    df_long=Acer_transfromHY_df_long,
    cluster_key=cluster_key
)
cluster_order = get_cluster_order(
    adata=Acer_transfromHY,
    cluster_key=cluster_key
)
cluster_color = get_cluster_color(
    adata=Acer_transfromHY,
    cluster_key=cluster_key,
    cluster_order=cluster_order
)
visualize(
    gene_list=marker_genes_OL,
    cluster_key=cluster_key,
    cluster_order=cluster_order,
    stats=stats,
    cluster_colors=cluster_color
)
# %%
sc.pl.umap(
    Acer_transfromHY,
    color="leiden_res0.25",
    legend_loc="on data"
)
sc.pl.umap(
    Acer_transfromHY,
    color="celltype_transfrom_HY",
    # legend_loc="on data"
)
#! 至20260926，发现iScience中提供的OLCs marker效果并不是特别好
#! 这批数据集的glia数似乎极少
#! 有必要参考Mu et al的西蜂marker，似乎她们的明显更全面
# %%
# 将已鉴定的cluster保存为csv格式，避免重复保存h5ad造成巨大写入量
KC_clusters_025 = {"KCs": ["0", "3", "8", "20", "21"], "glia": ["34"]}
cluster_to_celltype_025 = {
    cluster: celltype
    for celltype, clusters in KC_clusters_025.items()
    for cluster in clusters
}
Acer_transfromHY.obs["celltype_025_iSci"] = (
    Acer_transfromHY.obs["leiden_res0.25"]
    .map(cluster_to_celltype_025)
    .fillna("Unknown")
)
# %%
Acer_transfromHY.obs["celltype_025_iSci"].to_csv("./metadata/anno/celltype_025_iSci.csv")
# %%
# %%
# %%
# %%
# %%
# %%
# ------------------------------------------------------------------------------------------------------------------------------------
#! 另：使用best practice中的点图重新验证，本次使用的layer是logcounts，因为best practive中使用的是“scran-normalized, log-transformed counts”
# -----------------------------------------------------------------------------------------------------------------------------------
Acer_transfromHY.X = Acer_transfromHY.layers["logcounts"].copy()
# %%
marker_genes = {
    "neuron": ["LOC108004334"],
    "KC": ["LOC108002506", "LOC107999048"],
    "glia": ["LOC107993336", "LOC108001026"],
    "OPN": ["LOC108002036", "LOC108001977", "LOC107999745", "LOC108000694"],
    "hemo": ["LOC108004139", "LOC108001681"]
}

#* 颜色 = 每个 group 中该 gene 的 mean expression
#* 点大小 = 表达该 gene 的细胞比例
#* standard_scale="var" = 对每个 gene 在 groups 之间缩放到 0–1。

sc.pl.dotplot(
    Acer_transfromHY,
    groupby="leiden_res0.25",
    var_names=marker_genes,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
#! 结果与条形类似.
# KC: 0，3，8，20，21
# glia: 34
#! OPN：与条形不一致，5，6，13
# Hemo：似乎没有
# %%
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
    groupby="leiden_res0.25",
    var_names=OLC_marker_genes,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
#! OLC: 11, 17, 25, 26有可能
#! 有可能其余簇也有OLC，但比例较低（5，6），可能需要更高分辨率下重新查看
# %%
