# %%
from functions.heatmap_pheat import *
# %%
res = pd.read_csv("metadata/MN_metadata/auroc_Acer_vs_Amel.tsv", sep="\t")
# %%
# anno_acer = (
#     pd.read_csv(
#         "archive/results/annotation/Acer_cluster_annotation_res0.50.tsv",
#         sep="\t",
#         dtype={"cluster": str},
#     )
#     .set_index("cluster")["celltype_exploratory"]
#     .to_dict()
# )
# %%
# Amel 的 cluster 名本身包含类型前缀，例如 Glial_cells-23 -> Glial_cells
amel_clusters = (
    res.loc[res["Var1"].str.startswith("Amel."), "Var1"]
    .str.split(".", n=1).str[1]
    .unique()
)
anno_acer = None
anno_amel = {c: c.rsplit("-", 1)[0] for c in amel_clusters}
# %%
df_cross, mat = build_cross_species_matrix(
    res,
    species_a="Acer",
    species_b="Amel",
    value="AUROCs",
    include_species=False,      # 只显示 OLs-18 这种括号内内容
    group_columns_by_type=True, # 把 OLs 放一起、OPNs 放一起
    type_order=None,            # None: 字母顺序；也可自定义
    sort_rows=True,
)

# 若要自定义类型顺序：
# type_order=["KCs", "OLs", "OPNs", "UNs", "Glial_cells", "Hemocyte"]
plot_auroc_heatmap(
    mat,
    out_file="metadata/MN_metadata/auroc_heatmap_grouped.pdf",
    title="Acer (row) vs Amel (col)",
    row_cluster=True,
    col_cluster=False,          # 保持按类别排序
)
# %%
print(df_cross.shape)   # 42 * 49 = 2058
print(mat.shape)        # (42, 49)
# %%
# 可选：保存带 type 的跨物种长表，相当于 R 的 auroc.anno.tsv
df_cross.to_csv(
    "metadata/MN_metadata/auroc.anno.tsv",
    sep="\t",
    index=False,
)
# %%