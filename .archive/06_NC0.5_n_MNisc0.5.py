# %%
import scanpy as sc
import pandas as pd
# %%
Acer_transfromHY = sc.read_h5ad("./matrix/01_Acer_transfromHY.h5ad")
# %%
df = pd.read_csv("./metadata/anno/dot_and_mn_0.5.csv")
# %%
df["leiden_res0.50"] = df["leiden_res0.50"].astype(str)
# # 可选：排除未鉴定簇
# df = df[df["dot_mn"].isin(label_map.keys())]
# %%
# 按目标标签分组，收集簇号列表
agg = (
    df.groupby("dot_mn")["leiden_res0.50"]
    .apply(list)
    .to_dict()
)
# %%
print(agg)
# %%
celltype_050_dotMN = {
    cluster: celltype
    for celltype, clusters in agg.items()
    for cluster in clusters
}
# %%
Acer_transfromHY.obs["celltype_050_dotMN"] = (
    Acer_transfromHY.obs["leiden_res0.50"]
    .map(celltype_050_dotMN)
)
# %%
Acer_transfromHY.obs["celltype_050_dotMN"].to_csv("./metadata/anno/celltype_050_dotMN.csv")
# %%
sc.pl.umap(
    Acer_transfromHY,
    color="leiden_res0.50",
    legend_loc="on data"
)
# %%
