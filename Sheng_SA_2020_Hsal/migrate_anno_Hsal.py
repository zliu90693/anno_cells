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
Hsal_50_morethan20 = sc.pp.filter_genes(h5ad_concat, min_cells=20)