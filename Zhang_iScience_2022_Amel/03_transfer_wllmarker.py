# 20260801晚，陈大华的学生发来了完整版的marker，本次将把这一marker迁移到我的数据中
# %%
import scanpy as sc
import pandas as pd
import numpy as np
# %%
Amel = sc.read_h5ad("./matrix/01_Amel_transfromHY.h5ad")
# %%
sc.pl.umap(
    Amel,
    color=["leiden_res1.00", "celltype"],
    # legend_loc="on data"
)
# %%
wll_marker = pd.read_csv("./metadata/metadata.csv")
wll_marker
# %%
Amel.obs["wll_barcode"] = (
    Amel.obs_names.str.split("_").str[1] + "_sc_" + Amel.obs_names.str.split("_").str[0]
)
# %%
#? Amel.obs["wll_barcode"]和wll_marker["Unnamed: 0"]有多少重合？
# %%
print(len(set(Amel.obs["wll_barcode"].to_list()))) # 99200
print(len(set(wll_marker["Unnamed: 0"].to_list()))) # 115169
print(len(set(Amel.obs["wll_barcode"].to_list()) & set(wll_marker["Unnamed: 0"].to_list()))) # 86283
#! 注：HY的数据集与WLL共有细胞数：85158，代码如下：
# Amel_adataHY = sc.read_h5ad("./matrix/Amel_annobyHY.h5ad")
# Amel_adataHY.obs["HYkey"] = (
#     Amel_adataHY.obs_names
#     .to_series()
#     .replace({
#         "forager_rep1": "F1",
#         "forager_rep2": "F2",
#         "nurse_rep1": "N1",
#         "nurse_rep2": "N2",
#         "queen_rep1": "Q1",
#         "queen_rep2": "Q2",
#     }, regex=True)
# )
# Amel_adataHY.obs["wll_barcode"] = (
#     Amel_adataHY.obs["HYkey"].str.split("_").str[1] + "_sc_" + Amel_adataHY.obs_names.str.split("_").str[0]
# )
# print(len(set(Amel_adataHY.obs["wll_barcode"].to_list()) & set(wll_marker["Unnamed: 0"].to_list())))
# %%
# 将wll_marker中的第一列作为index，便于迁移
wll_marker.set_index("Unnamed: 0", inplace=True)
# %%
wll_marker_cell_type = wll_marker["cell_type"]
wll_marker_cell_type
# %%
# 将wll_marker中的cell_type迁移到Amel.obs中
Amel.obs["wll_celltype"] = (
    Amel.obs["wll_barcode"]
    .map(wll_marker_cell_type)
    .fillna("unknown")
)
# %%
sc.pl.umap(
    Amel,
    color=["leiden_res1.00", "celltype", "wll_celltype"],
    # legend_loc="on data"
)
# %%
#? celltype列：HY的标签，和wll_celltype列：WLL的标签，任一列不是unknown的细胞数有多少？
Amel.obs[(Amel.obs["celltype"] != "unknown") | (Amel.obs["wll_celltype"] != "unknown")].shape[0] # 97703
# 综上，Amel共99200个细胞，Amel与Amel_HY共有93447个共同细胞，Amel与wll_marker共86283个共同细胞。
# 如果采用投票的方式综合两组标签，将会注释成功97703个细胞。注释率高达0.98。
# %%
#! 为投票做准备：先将细胞标签统一
Amel.obs["celltype"] # 'UK', 'unknown', 'OPN', 'OLC', 'KC', 'Glia', 'Hemocyte'
Amel.obs["wll_celltype"] # 'unknown', 'OPNs', 'KCs', 'OLs', 'UNs', 'Glial cells', 'Hemocyte'
# Glial cells -> Glia
# UNs -> UK
# OLs -> OLC
# KCs -> KC
# OPNs -> OPN
Amel.obs["wll_celltype"] = (
    Amel.obs["wll_celltype"]
    .replace({
        "Glial cells": "Glia",
        "UNs": "UK",
        "OLs": "OLC",
        "KCs": "KC",
        "OPNs": "OPN",
        # "unknown" 保持不变
    })
)
# %%
#! 投票规则：针对同一细胞，如果HY和WLL注释结果相同，则最终注释采用共同注释结果。
#! 如果注释结果不同，分情况讨论:（1）如果其中一组是unknown，最终注释采用另一组的结果。（2）如果二者都不是unknown，待议。我需要先确定有多少组细胞符合情况（2）
#? 符合情况（2）：“都不是unknown，但注释结果不一样”的细胞有多少？
Amel_neitherunknown = Amel.obs[(Amel.obs["celltype"] != "unknown") & (Amel.obs["wll_celltype"] != "unknown")]
Amel_neitherunknown # 99200 -> 82027
Amel_neitherunknown[Amel_neitherunknown["celltype"] != Amel_neitherunknown["wll_celltype"]]
# 12165
#! 可知：共有12165个细胞符合情况（2）
# %%
#? 这12165个细胞大概是什么情况？
Amel_neitherunknown[Amel_neitherunknown["celltype"] != Amel_neitherunknown["wll_celltype"]][["celltype", "wll_celltype"]].value_counts()
# celltype  wll_celltype
# UK        OPN             5652
#           OLC             3388
# OLC       OPN             1579
#           UK               848
# UK        Glia             335
#           KC               141
# KC        Glia              99
# Glia      UK                58
# KC        UK                22
# UK        Hemocyte           8
# OPN       UK                 8
# KC        Hemocyte           5
# OLC       KC                 5
# KC        OLC                4
# OLC       Glia               3
# Hemocyte  UK                 3
# OLC       Hemocyte           1
# OPN       Glia               1
# Glia      OLC                1
# KC        OPN                1
# Hemocyte  Glia               1
# Glia      OPN                1
#           Hemocyte           1
# Name: count, dtype: int64
#! 可知：这12165个注释不匹配的细胞中，存在很多”HY中是UK，WLL中细胞身份已确定“的情况。而UK是原文中身份不确定的神经细胞，因此，在情况（2）的12165个细胞中，将以wll_celltype列为准
# %%
vote = np.select(
    [
        Amel.obs["celltype"] == "unknown",                          # HY缺失 -> wll
        Amel.obs["wll_celltype"] == "unknown",                      # WLL缺失 -> hy
        Amel.obs["celltype"] == Amel.obs["wll_celltype"],          # 确定且相同 -> 共同
    ],
    [
        Amel.obs["wll_celltype"],
        Amel.obs["celltype"],
        Amel.obs["celltype"],
    ],
    default=Amel.obs["wll_celltype"],   # 规则3：确定但不同 -> wll
)
Amel.obs["vote_celltype"] = vote
# %%
Amel.write_h5ad("./matrix/03_Amel_HYWLLvote.h5ad")
# %%
Amel.obs[Amel.obs["vote_celltype"] != "unknown"].shape[0] # 97703
# %%
