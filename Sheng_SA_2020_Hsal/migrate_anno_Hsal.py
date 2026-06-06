# %%
import scanpy as sc
# %%
Hsal_50 = sc.read_h5ad("./Hsal50.h5ad")
Hsal_51 = sc.read_h5ad("./Hsal51.h5ad")
# %%
Hsal_50
# %%
Hsal_51
# %%
# 检测为什么原文保留的基因数与我生成的差异那么大
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
# 开始迁移注释
Hsal_51.obs
# %%
Hsal_me = sc.read_h5ad("./cluster_unannotated.h5ad")
# %%
Hsal_51_cluster = Hsal_51.obs[["cell_type"]]
Hsal_51_cluster["cell_key"] = Hsal_51_cluster.index.str.split("-").str[-1]
Hsal_51_cluster
# %%
mycellset = set(Hsal_me.obs.index.str.split("-").str[0])
cellset51 = set(Hsal_51_cluster["cell_key"])
print(f"mycellset: {len(mycellset)}, cellset51: {len(cellset51)}")
print(len(mycellset.intersection(cellset51)))
# %%
