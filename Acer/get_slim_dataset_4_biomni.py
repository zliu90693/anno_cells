# %%
import scanpy as sc
import numpy as np
# %%
Acer = sc.read_h5ad("./matrix/02_Acer_transfromHY.h5ad")
# %%
# X 恢复为原始计数（原 X 是 logcounts，counts 存放在 layers 中）
Acer.X = Acer.layers["counts"].copy().astype(np.int32)
Acer.layers.clear()

# 只保留分辨率 0.5 的聚类和 batch
Acer.obs = Acer.obs[["leiden_neo_res_0.50", "batch"]].copy()
# 只保留基因 ID；不需要的话可改成 Acer.var[[]].copy()
Acer.var = Acer.var[["gene_ids"]].copy()

# 丢弃注释用不到的降维/邻接/高变基因等信息
Acer.obsm.clear()
Acer.varm.clear()
Acer.obsp.clear()
Acer.varp.clear()
Acer.raw = None

keep_uns = {
    "leiden_neo_res_0.50",
    "leiden_neo_res_0.50_colors",
    "batch_colors",
}
Acer.uns = {k: Acer.uns[k] for k in keep_uns if k in Acer.uns}

# 把 nullable-string-array 转成普通字符串，方便其它 anndata 版本读取
Acer.obs_names = Acer.obs_names.astype(str)
Acer.var_names = Acer.var_names.astype(str)
Acer.var["gene_ids"] = Acer.var["gene_ids"].astype(object)
# %%
Acer
# %%
Acer.write_h5ad(
    "./matrix/02_Acer_slim.h5ad",
    compression="gzip",
    compression_opts=4,
)
# %%
