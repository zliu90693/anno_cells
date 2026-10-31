# from pathlib import Path
# import numpy as np
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# from matplotlib.colors import LinearSegmentedColormap
# from matplotlib import cm

# # ---------------------------------------------------------------------------
# # 1. 从长表 res 构建跨物种 AUROC 矩阵
# # ---------------------------------------------------------------------------
# def build_cross_species_matrix(
#     res,
#     species_a="Acer",
#     species_b="Amel",
#     anno_a=None,
#     anno_b=None,
#     value="AUROCs",
#     include_species=False,
# ):
#     """
#     Parameters
#     ----------
#     res : DataFrame
#         02 的输出，至少包含 Var1, Var2, AUROCs（以及 AUROCp）。
#     species_a, species_b : str
#         行、列对应的物种。默认 A 为行、B 为列。
#     anno_a, anno_b : dict or None
#         {cluster: type}。None 时直接用 cluster 作为 type。
#     value : str
#         用 AUROCs 还是 AUROCp 画热图。
#     include_species : bool
#         标签是否带物种前缀。False 时接近 R 的 `type(cluster)` 写法。

#     Returns
#     -------
#     df_cross : DataFrame
#         跨物种注释后的长表。
#     mat : DataFrame
#         行 = species_a 的 cluster，列 = species_b 的 cluster。
#     """
#     df = res.copy()

#     p1 = df["Var1"].str.split(".", n=1, expand=True)
#     p2 = df["Var2"].str.split(".", n=1, expand=True)
#     df["sp1"], df["cl1"] = p1[0], p1[1]
#     df["sp2"], df["cl2"] = p2[0], p2[1]

#     # 关键：只保留跨物种配对，自动排除同物种对角线
#     df = df[df.sp1 != df.sp2].copy()
#     if len(df) == 0:
#         raise ValueError("res 中没有跨物种配对")

#     # 统一方向：行 = species_a，列 = species_b
#     m1 = (df.sp1 == species_a) & (df.sp2 == species_b)
#     m2 = (df.sp1 == species_b) & (df.sp2 == species_a)

#     if m2.any():
#         # R 的方向：Var1 = B, Var2 = A -> 行 A, 列 B
#         df = df[m2].copy()
#         df["cluster_a"], df["cluster_b"] = df.cl2, df.cl1
#     elif m1.any():
#         # 容错：只有 A vs B 时，交换方向
#         df = df[m1].copy()
#         df["cluster_a"], df["cluster_b"] = df.cl1, df.cl2
#     else:
#         raise ValueError("找不到 requested species 的跨物种配对")

#     # 再确认一次：只剩跨物种
#     assert (df.sp1 != df.sp2).all()

#     def map_type(cluster, anno):
#         if anno is None:
#             return cluster
#         return anno.get(cluster, cluster)

#     # df["type_a"] = [map_type(c, anno_a) for c in df.cluster_a]
#     # df["type_b"] = [map_type(c, anno_b) for c in df.cluster_b]
#     df["type_a"] = [anno_a.get(c, c) if anno_a else c for c in df.cluster_a]
#     df["type_b"] = [anno_b.get(c, c) if anno_b else c for c in df.cluster_b]

#     if include_species:
#         df["label_a"] = [f"{species_a}.{t}({c})" for t, c in zip(df.type_a, df.cluster_a)]
#         df["label_b"] = [f"{species_b}.{t}({c})" for t, c in zip(df.type_b, df.cluster_b)]
#     else:
#         # 对应 R 的 paste0(type, '(', cluster, ')')
#         df["label_a"] = [f"{t}({c})" for t, c in zip(df.type_a, df.cluster_a)]
#         df["label_b"] = [f"{t}({c})" for t, c in zip(df.type_b, df.cluster_b)]

#     mat = df.pivot(index="label_a", columns="label_b", values=value)
#     return df, mat

# # ---------------------------------------------------------------------------
# # 2. 蓝-白-红配色，对应 RColorBrewer Blues/Reds
# # ---------------------------------------------------------------------------
# def blue_white_red_cmap():
#     blues = cm.Blues(np.linspace(0.15, 1.0, 6))
#     reds = cm.Reds(np.linspace(0.15, 1.0, 6))
#     return LinearSegmentedColormap.from_list(
#         "blue_white_red",
#         list(blues[::-1]) + ["white"] + list(reds),
#         N=50,
#     )


# # ---------------------------------------------------------------------------
# # 3. 复现 pheatmap
# # ---------------------------------------------------------------------------
# def plot_auroc_heatmap(
#     mat,
#     out_file,
#     title="",
#     annot_cutoff=0.5,
#     figsize=(14, 11),
#     dpi=300,
# ):
#     def fmt(x):
#         # R: if (auroc[i,j] < 0.5) '' else signif(auroc[i,j], 2)
#         if pd.isna(x) or x < annot_cutoff:
#             return ""
#         return f"{x:.2g}"

#     annot = pd.DataFrame(
#         [[fmt(x) for x in row] for row in mat.values],
#         index=mat.index,
#         columns=mat.columns,
#     )

#     g = sns.clustermap(
#         mat,
#         cmap=blue_white_red_cmap(),
#         method="complete",       # R hclust 默认 complete
#         metric="euclidean",      # R dist 默认 euclidean
#         annot=annot,
#         fmt="",                  # annot 已经是字符串
#         annot_kws={"size": 5},
#         xticklabels=True,
#         yticklabels=True,
#         figsize=figsize,
#         linewidths=0.4,
#         linecolor="white",
#     )

#     g.ax_heatmap.set_xticklabels(
#         g.ax_heatmap.get_xticklabels(),
#         rotation=45,
#         ha="right",
#         fontsize=6,
#     )
#     g.ax_heatmap.set_yticklabels(
#         g.ax_heatmap.get_yticklabels(),
#         fontsize=6,
#     )
#     g.ax_heatmap.set_title(title, fontsize=8)

#     Path(out_file).parent.mkdir(parents=True, exist_ok=True)
#     g.savefig(out_file, dpi=dpi, bbox_inches="tight")
#     return g
import re
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap


def infer_type(cluster):
    """OLs-18 -> OLs; 0 -> 0"""
    return cluster.rsplit("-", 1)[0] if "-" in cluster else cluster


def natural_sort_key(s):
    """KCs-2 < KCs-10 < KCs-11 ..."""
    return "".join(
        f"{int(p):010d}" if p.isdigit() else p
        for p in re.split(r"(\d+)", str(s))
    )


def build_cross_species_matrix(
    res,
    species_a="Acer",
    species_b="Amel",
    value="AUROCs",
    include_species=False,
    anno_a=None,
    anno_b=None,
    group_columns_by_type=True,
    type_order=None,
    sort_rows=True,
):
    """
    构建跨物种 AUROC 矩阵。

    参数
    ----
    include_species : bool
        False 时标签只显示 cluster 本身（如 OLs-18），
        对应“只显示原版括号内内容”。
    group_columns_by_type : bool
        是否按 cluster 类型前缀把 Amel 列分组排序。
    type_order : list or None
        自定义类型顺序，例如 ["KCs","OLs","OPNs","UNs","Glial_cells","Hemocyte"]。
        None 时按字符串字母顺序。
    sort_rows : bool
        当 row_cluster=False 时，是否按 cluster 编号自然排序 Acer 行。
    """
    df = res.copy()

    p1 = df["Var1"].str.split(".", n=1, expand=True)
    p2 = df["Var2"].str.split(".", n=1, expand=True)
    df["sp1"], df["cl1"] = p1[0], p1[1]
    df["sp2"], df["cl2"] = p2[0], p2[1]

    # 只保留跨物种配对，自动排除同物种对角线
    df = df[df.sp1 != df.sp2].copy()

    m1 = (df.sp1 == species_a) & (df.sp2 == species_b)
    m2 = (df.sp1 == species_b) & (df.sp2 == species_a)

    if m2.any():
        # R 方向：Var1 = B, Var2 = A -> 行 A, 列 B
        df = df[m2].copy()
        df["cluster_a"], df["cluster_b"] = df.cl2, df.cl1
    elif m1.any():
        df = df[m1].copy()
        df["cluster_a"], df["cluster_b"] = df.cl1, df.cl2
    else:
        raise ValueError("找不到 requested species 的跨物种配对")

    assert (df.sp1 != df.sp2).all()

    # 类型仅用于列排序；不进入标签
    df["type_a"] = [
        anno_a.get(c, infer_type(c)) if anno_a else infer_type(c)
        for c in df.cluster_a
    ]
    df["type_b"] = [
        anno_b.get(c, infer_type(c)) if anno_b else infer_type(c)
        for c in df.cluster_b
    ]

    # 标签：只保留括号内内容 = cluster 本身
    if include_species:
        df["label_a"] = [f"{species_a}.{c}" for c in df.cluster_a]
        df["label_b"] = [f"{species_b}.{c}" for c in df.cluster_b]
    else:
        df["label_a"] = df.cluster_a
        df["label_b"] = df.cluster_b

    mat = df.pivot(index="label_a", columns="label_b", values=value)

    # 按 Amel cluster 的类型前缀分组排序
    if group_columns_by_type:
        info = (
            df.drop_duplicates("label_b")[["label_b", "type_b", "cluster_b"]]
            .copy()
        )
        info["_ckey"] = info["cluster_b"].map(natural_sort_key)

        if type_order is None:
            info["_tkey"] = info["type_b"]
        else:
            rank = {t: i for i, t in enumerate(type_order)}
            info["_tkey"] = info["type_b"].map(
                lambda t: f"{rank.get(t, len(type_order)):04d}"
            )

        info = info.sort_values(["_tkey", "_ckey"])
        mat = mat.loc[:, info["label_b"].tolist()]

    # Acer 行：如果关闭行聚类，则按编号自然排序
    if sort_rows:
        if all(str(x).isdigit() for x in mat.index):
            mat = mat.reindex(sorted(mat.index, key=lambda x: int(str(x))))
        else:
            mat = mat.reindex(sorted(mat.index, key=natural_sort_key))

    return df, mat


def blue_white_red_cmap():
    blues = cm.Blues(np.linspace(0.15, 1.0, 6))
    reds = cm.Reds(np.linspace(0.15, 1.0, 6))
    return LinearSegmentedColormap.from_list(
        "blue_white_red",
        list(blues[::-1]) + ["white"] + list(reds),
        N=50,
    )


def plot_auroc_heatmap(
    mat,
    out_file,
    title="",
    annot_cutoff=0.5,
    figsize=(14, 11),
    dpi=300,
    row_cluster=True,
    col_cluster=False,
):
    def fmt(x):
        if pd.isna(x) or x < annot_cutoff:
            return ""
        return f"{x:.2g}"

    annot = pd.DataFrame(
        [[fmt(x) for x in row] for row in mat.values],
        index=mat.index,
        columns=mat.columns,
    )

    g = sns.clustermap(
        mat,
        cmap=blue_white_red_cmap(),
        method="complete",
        metric="euclidean",
        annot=annot,
        fmt="",
        annot_kws={"size": 5},
        row_cluster=row_cluster,
        col_cluster=col_cluster,   # 关键：关闭列聚类，保留按类型排序
        xticklabels=True,
        yticklabels=True,
        figsize=figsize,
        linewidths=0.4,
        linecolor="white",
    )

    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(),
        rotation=45,
        ha="right",
        fontsize=6,
    )
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(),
        fontsize=6,
    )
    g.ax_heatmap.set_title(title, fontsize=8)

    Path(out_file).parent.mkdir(parents=True, exist_ok=True)
    g.savefig(out_file, dpi=dpi, bbox_inches="tight")
    return g

# df_cross, mat = build_cross_species_matrix(
#     res,
#     species_a="Acer",
#     species_b="Amel",
#     value="AUROCs",
#     include_species=False,      # 只显示 OLs-18 这种括号内内容
#     group_columns_by_type=True, # 把 OLs 放一起、OPNs 放一起
#     type_order=None,            # None: 字母顺序；也可自定义
#     sort_rows=True,
# )

# # 若要自定义类型顺序：
# # type_order=["KCs", "OLs", "OPNs", "UNs", "Glial_cells", "Hemocyte"]

# plot_auroc_heatmap(
#     mat,
#     out_file="metadata/MN_metadata/auroc_heatmap_grouped.pdf",
#     title="Acer (row) vs Amel (col)",
#     row_cluster=True,
#     col_cluster=False,          # 保持按类别排序
# )