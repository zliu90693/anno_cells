import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata

def get_gene_cluster_table( # 从表达矩阵中提取一个仅包括指定gene和指定簇的子矩阵，格式为长表
        adata: anndata.AnnData, 
        cluster_key: str, # “细胞簇划分”在Anndata对象中对应的列名
        gene_list: list, # marker gene列表
        layer: str # 提取过程使用的layer
    ) -> pd.DataFrame:
    df = sc.get.obs_df(
        adata,
        keys=[cluster_key] + gene_list,
        layer=layer
    )
    df_long = df.melt( # 宽表转长表
        id_vars=cluster_key,
        value_vars=gene_list,
        var_name="gene",
        value_name="expression"
    )
    return df_long

def cal_mean_etal( # 计算每个 cluster 的 mean / SD / N / SEM 
        df_long: pd.DataFrame,
        cluster_key: str,
    ) -> pd.DataFrame:
    stats = (
        df_long
        .groupby(
            ["gene", cluster_key],
            observed=True
        )["expression"]
        .agg(
            mean="mean",
            sd="std",
            n="count"
        )
        .reset_index()
    )
    stats["sem"] = (
        stats["sd"] /
        np.sqrt(stats["n"])
    )
    return stats

def get_cluster_order( # cluster 顺序
        adata: anndata.AnnData, 
        cluster_key: str,
    ) -> list:
    if isinstance(
        adata.obs[cluster_key].dtype,
        pd.CategoricalDtype
    ):
        cluster_order = list(
            adata.obs[cluster_key].cat.categories
        )
    else:
        cluster_order = sorted(
            adata.obs[cluster_key].unique()
    )
    return cluster_order

def get_cluster_color(
        adata: anndata.AnnData, 
        cluster_key: str,
        cluster_order: list
    ) -> list:
    color_key = f"{cluster_key}_colors"
    if color_key in adata.uns:
        colors = adata.uns[color_key]
        cluster_colors = dict(
            zip(
                adata.obs[cluster_key].cat.categories,
                colors
            )
        )
    else:
        cmap = plt.get_cmap("tab20")
        cluster_colors = {
            cluster: cmap(i % 20)
            for i, cluster in enumerate(cluster_order)
        }
    return cluster_colors

def visualize(
        gene_list: list, # marker gene列表
        cluster_key: str,
        cluster_order: list,
        stats: pd.DataFrame,
        cluster_colors: list
    ):
    fig, axes = plt.subplots(
        nrows=len(gene_list),
        ncols=1,
        figsize=(12, 12),
        sharex=True
    )
    x = np.arange(len(cluster_order))
    for ax, gene in zip(axes, gene_list):
        temp = stats[
            stats["gene"] == gene
        ].copy()
        temp[cluster_key] = pd.Categorical(
            temp[cluster_key],
            categories=cluster_order,
            ordered=True
        )
        temp = (
            temp
            .sort_values(cluster_key)
            .set_index(cluster_key)
            .reindex(cluster_order)
        )
        bar_colors = [
            cluster_colors[c]
            for c in cluster_order
        ]
        ax.bar(
            x,
            temp["mean"],
            yerr=temp["sem"],
            color=bar_colors,
            width=0.75,
            capsize=2,
            edgecolor="none"
        )
        # gene 名称
        ax.text(
            0.02,
            0.80,
            gene,
            transform=ax.transAxes,
            fontsize=10,
            fontstyle="italic",
            fontweight="bold"
        )
        ax.set_ylabel("")
        # 去掉上、右边框
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # 设置x轴
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(
        cluster_order,
        fontsize=12
    )
    axes[-1].set_xlabel("Cell cluster")
    fig.supylabel(
        "Average expression",
        fontsize=13,
        fontweight="bold"
    )
    plt.tight_layout()
    plt.show()