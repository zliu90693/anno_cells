# %%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
from functions.anno_utils import *
# %%
Acer_transfromHY = sc.read_h5ad("./matrix/01_Acer_transfromHY.h5ad")
Acer_transfromHY.X = Acer_transfromHY.layers["counts"].copy()
sc.pp.normalize_total(Acer_transfromHY)
Acer_transfromHY.layers["norm_counts"] = Acer_transfromHY.X.copy()
# %%
# Mu et al的附件中只提供了西方蜜蜂 Marker 的 gene name，我需要根据西方蜜蜂的 GTF 获取其 gene id，随后获取东方蜜蜂的 gene id
Amel_gtf = pd.read_csv("./metadata/Apis_mellifera.tsv", sep="\t")
Amel_gid_gname = Amel_gtf[["gene_name", "gene_id"]]
# %%
Amel_gid_gname.to_csv("./metadata/Amel_gid_gname.tsv", sep="\t")
# %%
#! --------------------------------------- KC Marker --------------------------------------- 
Mu_KC_Marker = {"Pkc", "LOC408976", "E75", "LOC408804", "Cac", "Tk"} #! Mu et al附件中提及的西蜂Marker (gene name格式)
print(Amel_gid_gname[Amel_gid_gname["gene_name"].isin(Mu_KC_Marker)].drop_duplicates())
#   gene_name        gene_id
#   Pkc  GeneID_406106
#   E75  GeneID_410309
#   Tk  GeneID_406083
#   Cac  GeneID_408764
#   LOC408976   LOC408976
#   LOC408804   LOC408804
# %%
KC_Marker = [ #! 东方蜜蜂的同源Marker
    "LOC107995382", # Pkc
    "LOC107997697", # E75
    "LOC107998479", # Tk
    "LOC107997070", # Cac
    "LOC108000569", # LOC408976
    "LOC107999048", # LOC408804
]
# %%
cluster_key="leiden_res0.25"
layer="norm_counts"
Acer_transfromHY_df_long = get_gene_cluster_table(
    adata=Acer_transfromHY,
    cluster_key=cluster_key,
    gene_list=KC_Marker,
    layer=layer
)

cluster_order = get_cluster_order(
    adata=Acer_transfromHY,
    cluster_key=cluster_key
)

visualize(
    gene_list=KC_Marker,
    cluster_key=cluster_key,
    cluster_order=cluster_order,
    stats=cal_mean_etal(
        df_long=Acer_transfromHY_df_long,
        cluster_key=cluster_key
    ),
    cluster_colors=get_cluster_color(
        adata=Acer_transfromHY,
        cluster_key=cluster_key,
        cluster_order=cluster_order
    )
)
#! 与iScience结论相同，0, 3, 8, 20, 21
# %%
#! --------------------------------------- Glia Marker --------------------------------------- 
Mu_Glia_Marker = {"Gat-1B", "Eaat-2", "LOC725117", "LOC411976", "LOC412382", "LOC411704"} #! Mu et al附件中提及的西蜂Marker (gene name格式)
print(Amel_gid_gname[Amel_gid_gname["gene_name"].isin(Mu_Glia_Marker)].drop_duplicates())
#   gene_name   gene_id
#   Gat-1B  GeneID_494511
#   Eaat-2  GeneID_406113
#   LOC725117   LOC725117
#   LOC411976   LOC411976
#   LOC412382   LOC412382
#   LOC411704   LOC411704
Glia_Marker = [
    "LOC107993519", # Gat-1B
    "LOC108002812", # Eaat-2
    "LOC108000234", # LOC725117
    "LOC107994035", # LOC411976
    "LOC107996595", # LOC412382
    "LOC108003667" # LOC411704
]

cluster_key="leiden_res0.25"
layer="norm_counts"
Acer_transfromHY_df_long = get_gene_cluster_table(
    adata=Acer_transfromHY,
    cluster_key=cluster_key,
    gene_list=Glia_Marker,
    layer=layer
)

cluster_order = get_cluster_order(
    adata=Acer_transfromHY,
    cluster_key=cluster_key
)

visualize(
    gene_list=Glia_Marker,
    cluster_key=cluster_key,
    cluster_order=cluster_order,
    stats=cal_mean_etal(
        df_long=Acer_transfromHY_df_long,
        cluster_key=cluster_key
    ),
    cluster_colors=get_cluster_color(
        adata=Acer_transfromHY,
        cluster_key=cluster_key,
        cluster_order=cluster_order
    )
)
#! 只有34是glia，本组东方蜜蜂的glia占比确实极低
# %%
#! --------------------------------------- OLC --------------------------------------- 
Mu_OLC_Marker = {"LOC413215", "LOC724146", "LOC412767", "LOC724801", "Hiscl1", "LOC726953", "LOC724847", "LOC726153", "LOC100577959"}
print(Amel_gid_gname[Amel_gid_gname["gene_name"].isin(Mu_OLC_Marker)].drop_duplicates())
#   gene_name   gene_id
#   Hiscl1  GeneID_409040
#   LOC413215   LOC413215
#   LOC724146   LOC724146
#   LOC412767   LOC412767
#   LOC724801   LOC724801
#   LOC726953   LOC726953
#   LOC724847   LOC724847
#   LOC726153   LOC726153
#   LOC100577959    LOC100577959
OLC_Marker = [
    "LOC107996886", # GeneID_409040
    "LOC107999967", # LOC413215
    "LOC108002029", # LOC724146
    "LOC108004344", # LOC412767
    "LOC107999009", # LOC724801
    "LOC107996386", # LOC726953
    "LOC107998009", # LOC724847
    "LOC108000705", # LOC726153
    "LOC107992410"# LOC100577959
]

cluster_key="leiden_res0.25"
layer="norm_counts"
Acer_transfromHY_df_long = get_gene_cluster_table(
    adata=Acer_transfromHY,
    cluster_key=cluster_key,
    gene_list=OLC_Marker,
    layer=layer
)

cluster_order = get_cluster_order(
    adata=Acer_transfromHY,
    cluster_key=cluster_key
)

visualize(
    gene_list=OLC_Marker,
    cluster_key=cluster_key,
    cluster_order=cluster_order,
    stats=cal_mean_etal(
        df_long=Acer_transfromHY_df_long,
        cluster_key=cluster_key
    ),
    cluster_colors=get_cluster_color(
        adata=Acer_transfromHY,
        cluster_key=cluster_key,
        cluster_order=cluster_order
    )
)
#! OLC Marker效果比iScience好很多
#! OLC: 5, 9, 10, 11, 13, 14, 15, 17, 27, 29, 31, 33
# %%
#! --------------------------------------- OPN --------------------------------------- 
OPN_Marker = [ # 原文只提供了一个：LOC413466
    "LOC108000694", # LOC413466
    "LOC108002036" # iScience的marker，用于防报错
]
cluster_key="leiden_res0.25"
layer="norm_counts"
Acer_transfromHY_df_long = get_gene_cluster_table(
    adata=Acer_transfromHY,
    cluster_key=cluster_key,
    gene_list=OPN_Marker,
    layer=layer
)

cluster_order = get_cluster_order(
    adata=Acer_transfromHY,
    cluster_key=cluster_key
)

visualize(
    gene_list=OPN_Marker,
    cluster_key=cluster_key,
    cluster_order=cluster_order,
    stats=cal_mean_etal(
        df_long=Acer_transfromHY_df_long,
        cluster_key=cluster_key
    ),
    cluster_colors=get_cluster_color(
        adata=Acer_transfromHY,
        cluster_key=cluster_key,
        cluster_order=cluster_order
    )
)
#! OPN: 5, 6, 13
# %%
# 将已鉴定的cluster保存为csv格式，避免重复保存h5ad造成巨大写入量
# 来自NC marker，使用条形图的鉴定结果（可能相对点图不太可信，忽略了表达特定marker的细胞比例）
KC_clusters_025_bar = {"KCs": ["0", "3", "8", "20", "21"], 
                   "glia": ["34"], 
                   "OLC": ["9", "10", "11", "14", "15", "17", "27", "29", "31", "33"],
                   "OPN": ["6"],
                   "OLC?OPN": ["5", "13"]}
cluster_to_celltype_025_bar = {
    cluster: celltype
    for celltype, clusters in KC_clusters_025_bar.items()
    for cluster in clusters
}
Acer_transfromHY.obs["celltype_025_NC_bar"] = (
    Acer_transfromHY.obs["leiden_res0.25"]
    .map(cluster_to_celltype_025_bar)
    .fillna("Unknown")
)
Acer_transfromHY.obs["celltype_025_NC_bar"].to_csv("./metadata/anno/celltype_025_NC_bar.csv")
# %%
#! 另：使用best practice中的点图重新验证，本次使用的layer是logcounts，因为best practive中使用的是“scran-normalized, log-transformed counts”
# -----------------------------------------------------------------------------------------------------------------------------------
# %%
Acer_transfromHY.X = Acer_transfromHY.layers["logcounts"].copy()
marker_genes = {
    "KC": ["LOC107995382", "LOC107997697", "LOC107998479", "LOC107997070", "LOC108000569", "LOC107999048"],
    "glia": ["LOC107993519", "LOC108002812", "LOC108000234", "LOC107994035", "LOC107996595", "LOC108003667"],
    "OLC": ["LOC107996886", "LOC107999967", "LOC108002029", "LOC108004344", "LOC107999009", "LOC107996386", "LOC107998009", "LOC108000705", "LOC107992410"],
    "OPN": ["LOC108000694"]
}

sc.pl.dotplot(
    Acer_transfromHY,
    groupby="leiden_res0.25",
    var_names=marker_genes,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
# KC: 同条形，0，3，8，20，21
# glia: 同，34
# OLC: 同，9, 10, 11, 14, 15, 17, 27, 29, 31, 33
#! OPN: 无法确定，5/6/13表达量高但是比例都极小，需增大分辨率后重新鉴定

# 来自NC marker，使用点图的鉴定结果，相对更准确，因为反映了细胞比例
KC_clusters_025_dot = {"KCs": ["0", "3", "8", "20", "21"], 
                   "glia": ["34"], 
                   "OLC": ["9", "10", "11", "14", "15", "17", "27", "29", "31", "33"],
                   "maybe_OPN": ["5", "6", "13"]}
cluster_to_celltype_025_dot = {
    cluster: celltype
    for celltype, clusters in KC_clusters_025_dot.items()
    for cluster in clusters
}
Acer_transfromHY.obs["celltype_025_NC_dot"] = (
    Acer_transfromHY.obs["leiden_res0.25"]
    .map(cluster_to_celltype_025_dot)
    .fillna("Unknown")
)
Acer_transfromHY.obs["celltype_025_NC_dot"].to_csv("./metadata/anno/celltype_025_NC_dot.csv")
# %%
# %%
# %%
# %%
# %%
# %%
# %%
#? 解决部分簇同时表达两类细胞marker(OLC/OPN)的问题: 查看两类marker在数据集中的表达量
# 簇5和簇13，同时高表达LOC108000694（OPN）和LOC107996386（OLC），需要确定是否是因为leiden分辨率过低导致把两类细胞簇融合到了一起
sc.pl.umap(
    Acer_transfromHY,
    color=["LOC108000694", "LOC107996386"], # OPN OLC
    vmin=0,
    vmax="p99", 
    sort_order=False,
    frameon=False,
    cmap="Reds"
)
#! 结论是，13可能包含多个细胞类型，但无法确定具体身份，5似乎仅一种细胞类型，身份不确定
# %%
# ------------------------------------------------------------------------------------------------------------------------------------# %%
