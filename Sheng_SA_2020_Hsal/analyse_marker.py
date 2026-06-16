# %%
import pandas as pd
import scanpy as sc
# %%
Dmel_markers = pd.read_csv("./metadata/Dmel_Markers_GeneID.csv")
Dmel_markers
# %%
Dmel_markers_available = Dmel_markers[~Dmel_markers["GeneID"].isna()]
Dmel_markers_available
# %%
Dmel_Hsal = pd.read_csv("/home/liuzhiyu/Projects/neo_caste/Find_Ortholog/primary_transcripts/OrthoFinder/Results_Jun15/Orthologues/Orthologues_Harpegnathos_saltator/Harpegnathos_saltator__v__Drosophila_melanogaster.tsv", sep="\t")
Dmel_Hsal
# Dmel_Hsal["Drosophila_melanogaster_list"] = Dmel_Hsal["Drosophila_melanogaster"].str.split(",")
# %%
Dmel_Hsal_exploded = Dmel_Hsal.assign( #! explode!!!
    Drosophila_melanogaster=Dmel_Hsal['Drosophila_melanogaster'].str.split(r'\s*,\s*')
).explode('Drosophila_melanogaster').reset_index(drop=True)
# # %%
# Dmel_Hsal_exploded_array = Dmel_Hsal_exploded.set_index("Drosophila_melanogaster")["Harpegnathos_saltator"]
# # %%
# Dmel_markers_available["Hsal_marker"] = Dmel_markers_available["GeneID"].map(Dmel_Hsal_exploded_array)
# Dmel_markers_available
# %%
Dmel_Hsal_exploded = Dmel_Hsal_exploded[["Harpegnathos_saltator", "Drosophila_melanogaster"]]
Dmel_Hsal_exploded.rename(columns={
    'Drosophila_melanogaster': 'GeneID', 
    "Harpegnathos_saltator": "Hsal_marker"
}, inplace=True)
# %%
Dmel_markers_available = Dmel_markers_available.merge(
    Dmel_Hsal_exploded,
    on = "GeneID",
    how = "left"
)
# %%
Dmel_markers_available
#! 注: 仍然有几个果蝇基因没使用Orthofinder结果匹配到同源的Hsal基因
#! 接下来, 我将试图使用OrthoDB找到这几个果蝇基因在跳镰猛蚁中的同源基因
# %%
markers_unobtainable = Dmel_markers_available[Dmel_markers_available["Hsal_marker"].isna()][["Cell_Type", "MarkerName", "GeneID"]]
markers_unobtainable
# %%
# 首先, 我在 NCBI 的 gene database 中查询 GeneID 列的基因名称. 然后在 Homology 栏选择 NCBI Ortholog, 然后搜 harpegnathos saltator
# FBgn0033799 没找到同源基因
# FBgn0026439 没找到同源基因
# FBgn0265767 有, 名字叫zyd, 但未匹配到的原因可能是当前的参考基因组是 Hsal v8.5, 而 zyd的注释出现在 v8.6 中.

# %%
Dmel_markers_available
# %%
Dmel_markers_available.to_csv("./markers/Hsal_Dmel_allcluster_markers.csv")
# %%
