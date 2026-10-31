# %%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
# %%
import scanpy as sc
# %%
Acer = sc.read_h5ad("./matrix/01_Acer_transfromHY.h5ad")
Acer
# %%
#! 基于NC的marker的表达量条形图确定的细胞身份(分辨率0.25)
celltype_025_NC_bar = pd.read_csv("../../Acer/metadata/anno/celltype_025_NC_bar.csv")
celltype_025_NC_bar["Unnamed: 0"].to_list() == Acer.obs.index.to_list() # True
# %%
celltype_map_bar = celltype_025_NC_bar.set_index("Unnamed: 0")["celltype_025_NC_bar"]
Acer.obs["celltype_025_NC_bar"] = Acer.obs.index.map(celltype_map_bar)
# %%
sc.pl.umap(
    Acer, 
    color=["celltype_025_NC_bar"]
)
# %%
#! 基于NC的marker的表达量点图确定的细胞身份(分辨率0.25)
celltype_025_NC_dot = pd.read_csv("../../Acer/metadata/anno/celltype_025_NC_dot.csv")
celltype_025_NC_dot["Unnamed: 0"].to_list() == Acer.obs.index.to_list() # True
# %%
celltype_map_dot = celltype_025_NC_dot.set_index("Unnamed: 0")["celltype_025_NC_dot"]
Acer.obs["celltype_025_NC_dot"] = Acer.obs.index.map(celltype_map_dot)
# %%
sc.pl.umap(
    Acer, 
    color=["celltype_025_NC_dot"]
)
# %%
#! 基于NC的marker的表达量点图确定的细胞身份(分辨率0.50)
celltype_050_NC_dot = pd.read_csv("../../Acer/metadata/anno/celltype_050_NC_dot.csv")
celltype_050_NC_dot["Unnamed: 0"].to_list() == Acer.obs.index.to_list() # True
# %%
celltype_map_dot = celltype_050_NC_dot.set_index("Unnamed: 0")["celltype_050_NC_dot"]
Acer.obs["celltype_050_NC_dot"] = Acer.obs.index.map(celltype_map_dot)
# %%
sc.pl.umap(
    Acer, 
    color=["celltype_050_NC_dot"]
)
# %%
#! 基于NC的marker的表达量点图确定的细胞身份(分辨率1.00)
celltype_100_NC_dot = pd.read_csv("../../Acer/metadata/anno/celltype_100_NC_dot.csv")
celltype_100_NC_dot["Unnamed: 0"].to_list() == Acer.obs.index.to_list() # True
# %%
celltype_map_dot = celltype_100_NC_dot.set_index("Unnamed: 0")["celltype_100_NC_dot"]
Acer.obs["celltype_100_NC_dot"] = Acer.obs.index.map(celltype_map_dot)
# %%
sc.pl.umap(
    Acer, 
    color=["celltype_100_NC_dot"]
)
# %%
#! 使用iScience的数据进行MetaNeighbor(py)(Acer分辨率0.5），其结果与NC marker综合后的结果：
celltype_050_dotMN = pd.read_csv("../../Acer/metadata/anno/celltype_050_dotMN.csv")
celltype_050_dotMN["Unnamed: 0"].to_list() == Acer.obs.index.to_list() # True
# %%
celltype_map_dot = celltype_050_dotMN.set_index("Unnamed: 0")["celltype_050_dotMN"]
Acer.obs["celltype_050_dotMN"] = Acer.obs.index.map(celltype_map_dot)
sc.pl.umap(
    Acer, 
    color=["celltype_050_dotMN"]
)
# %%
