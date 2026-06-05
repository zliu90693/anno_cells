# !/home/liuzhiyu/Software/miniconda3/envs/rm_ambient_doublet/bin/Rscript
# %%
# .libPaths(c("/home/liuzhiyu/Software/miniconda3/envs/rm_ambient_doublet/lib/R/library", .libPaths()))
# %%
library(anndataR)
# %%
Hsal50 <- readRDS("seurat.HSAL50.2020.rds")
Hsal51 <- readRDS("seurat.HSAL51.2021.rds")
# %%
# https://anndatar.scverse.org/articles/usage_seurat.html
Hsal50_adata <- as_AnnData(Hsal50)
Hsal51_adata <- as_AnnData(Hsal51)
# %%
# dir()
Hsal50_adata$write_h5ad("Hsal50.h5ad")
Hsal51_adata$write_h5ad("Hsal51.h5ad")
# 