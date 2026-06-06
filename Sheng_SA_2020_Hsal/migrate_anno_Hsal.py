# %%
import scanpy as sc
import pandas as pd
# %%
Hsal_50 = sc.read_h5ad("./Hsal50.h5ad")
Hsal_51 = sc.read_h5ad("./Hsal51.h5ad")
# %%
Hsal_50
# %%
Hsal_51
# %%
# ------------------------------------- 检测为什么原文保留的基因数与我生成的差异那么大 ------------------------------------- 
Hsal_original = sc.read_10x_h5("/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013383.h5")
Hsal_original # 13684个基因, 排除我用的参考基因组本身包含基因过少的原因
# %%
Hsal_50.X = Hsal_50.layers['counts'].copy()
sc.pp.calculate_qc_metrics(
    Hsal_50, 
    inplace=True, 
    percent_top=[20], 
    log1p=True 
)
# %%
sc.pp.filter_genes(Hsal_50, min_cells=20)
Hsal_50
# %%
total_cells = Hsal_50.n_obs
Hsal_50.var['n_cells_not_detected'] = total_cells - Hsal_50.var['n_cells_by_counts']
mask = Hsal_50.var['n_cells_not_detected'] >= 20
num_genes = mask.sum()
# %%
# ------------------------------------- 开始迁移注释 -------------------------------------

# %%
#! 注意! Hsal原文与我产生的数据在 adata.obs_names 的格式上有差异, 需要统一!!!
#! 如果不进行格式统一, 须知不同样本间的细胞条形码有可能是重复的, 仅使用条形码无法作为分辨细胞的唯一要素!!!

# %%
# %%
Hsal_51 = sc.read_h5ad("./Hsal51.h5ad")
# %%
# Hsal_51.obs_names 格式: gam.nb.e13b1-AAACGGGAGTACCGGA
Hsal_51_cluster = Hsal_51.obs[["cell_type"]]
Hsal_51_cluster
# %%
Hsal_me = sc.read_h5ad("./cluster_unannotated.h5ad")
Hsal_me.obs
# %%
# Hsal_me.obs_names 格式: GGACAGATCTATCCCG-1_GSM4013383
Hsal_me.obs["barcode"] = Hsal_me.obs.index.str.split("-").str[0]
Hsal_me.obs["GSM"] = Hsal_me.obs.index.str.split("_").str[-1]
Hsal_me.obs
# %%
# 手动通过 ncbi 整理 GSM 与 internal_name 的对应关系, 得到 GSM_vs_internalname.csv
GSM_vs_internalname = pd.read_csv("./GSM_vs_internalname.csv", index_col="GSM")["internal_name"]
GSM_vs_internalname # GSM_vs_internalname 是一个 Array
# %%
Hsal_me.obs["internal_name"] = Hsal_me.obs["GSM"].map(GSM_vs_internalname)
Hsal_me.obs
# %%
# 生成 gam.nb.e13b1-AAACGGGAGTACCGGA 格式
Hsal_me.obs["key"] = Hsal_me.obs["internal_name"] + '-' + Hsal_me.obs["barcode"]
Hsal_me.obs
# %%
# 检测 Hsal_me.obs["key"] 与 Hsal_51_cluster.index 的重合率
Hsal_me_set = set(Hsal_me.obs["key"])
Hsal_51_set = set(Hsal_51_cluster.index)
print(len(Hsal_me_set))
print(len(Hsal_51_set))
print(len(Hsal_me_set & Hsal_51_set))
print(len(Hsal_me_set) - len(Hsal_me_set & Hsal_51_set)) # 1772
# %%
Hsal_51_celltype = Hsal_51_cluster["cell_type"]
Hsal_me.obs["cell_type_51"] =  Hsal_me.obs["key"].map(Hsal_51_celltype).fillna('Unknown')
Hsal_me.obs
# %%
# 统计分类变量出现频次, Unknown 正好是 1772
Hsal_me.obs.groupby('cell_type_51').size()
# %%
# %%
sc.pl.embedding(
        Hsal_me,
        basis="X_umap",
        color="leiden_res0.50",
        # legend_loc="on data"
    )

sc.pl.embedding(
        Hsal_me,
        basis="X_umap",
        color="cell_type_51",
        # legend_loc="on data"
    )
# %%
# ------------------------------------- Hsal 50 -------------------------------------
# 
Hsal_50_cluster = Hsal_50.obs[["cluster_name"]]
Hsal_50_celltype = Hsal_50_cluster["cluster_name"]
Hsal_me.obs["key_50"] = Hsal_me.obs["key"].str[:-4] # 注意, Hsal_50.obs_names 中 barcode 去除了后四个碱基
Hsal_me.obs["cell_type_50"] =  Hsal_me.obs["key_50"].map(Hsal_50_celltype).fillna('Unknown')
Hsal_me.obs
# %%
Hsal_me.obs.groupby('cell_type_50').size()
# %%
sc.pl.embedding(
        Hsal_me,
        basis="X_umap",
        color="cell_type_50",
        legend_loc="on data"
    )

sc.pl.embedding(
        Hsal_me,
        basis="X_umap",
        color="cell_type_51",
        legend_loc="on data"
    )

sc.pl.embedding(
        Hsal_me,
        basis="X_umap",
        color="leiden_res1.00",
        legend_loc="on data"
    )
# %%
Hsal_me.write_h5ad("./Hsal_me_with_5051_anno.h5ad")
# %%
