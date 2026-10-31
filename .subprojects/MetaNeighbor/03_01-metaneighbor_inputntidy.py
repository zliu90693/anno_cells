# ----------------------------------------------------------------------------
"""
R 的实际数据流是：
pb(pseudo-cell 原始计数) ──saveRDS──> pb.RDS
        │
        └─ CP10K(pb.CP10K) ──> z-score(pb.norm) ──> FindAllMarkers ──> marker.tsv

python 的数据流：
pb_X(pseudo-cell 原始计数)
    ├─ CP10K(1000) ──> z-score(Z) ──> 未被使用
    └─ CP10K(1e4) + log1p ──> Wilcoxon ──> markers
"""
# %%
from functions.tidyAnndata import *
# %%
# 外层只负责 I/O：
# %%
# -------------------------------------------------------------------------------------------
# -------------------------------------- Acer and Amel -------------------------------------- 
# -------------------------------------------------------------------------------------------
# %%
Acer = sc.read_h5ad("./matrix/Acer_4_MetaNeighbor_py.h5ad")
Amel = sc.read_h5ad("./matrix/Amel_wll_4_MetaNeighbor_py.h5ad")
Acer_Amel_121 = pd.read_csv("./orthologs/Acer_vs_Amel_121_4mn.tsv", sep="\t")
# %%
# %%
"""
到这里, 我发现之前先前在notebook中保存 ***_4_MetaNeighbor_py.h5ad 是一个错误的决定, 因为这个文件还需要对基因进一步取子集, 
原因是 Acer 中的基因不会全部出现在那个 121 对照表中, 
这在工程学水平是一个非常丑的实现
"""
# %%
Amel_gid_set = set(Amel.var.index) # 9154 
Acer_Amel_121_sub = Acer_Amel_121[Acer_Amel_121["Apis_mellifera"].isin(Amel_gid_set)] # 9154
# %%
Acer_gid_set = set(Acer.var.index) # 11817
Acer_inter = Acer_gid_set & set(Acer_Amel_121_sub["Apis_cerana"].to_list()) # 8910
Amel_2besave = set(Acer_Amel_121_sub[Acer_Amel_121_sub["Apis_cerana"].isin(Acer_inter)]["Apis_mellifera"].to_list()) # 8910
# %%
Amel4MN = Amel[:, Amel.var.index.isin(Amel_2besave)].copy() # 115169 × 8910
Acer4MN = Acer[:, Acer.var.index.isin(Acer_inter)].copy() # 58166 × 8910
# %%
Amel4MN.write_h5ad("./matrix/MN_input/Amel4MN_8910.h5ad")
Acer4MN.write_h5ad("./matrix/MN_input/Acer4MN_8910.h5ad")
# %%
# -------------------------------------- Amel -------------------------------------- 
pb_Amel, markers_Amel = run_tidy_sce(
    Amel4MN, label_col="cell_type_cluster",
    n_per=10, seed=0, keep_leftover=True,
    target_sum=1e4, min_logfc=1.25, min_pct=0.2,
)
from pathlib import Path
out_pb_Amel = Path("./matrix/MN_dataset/Amel_pseudocell_raw.h5ad")
out_mk_Amel = Path("./metadata/MN_metadata/Amel_markers.tsv")
out_pb_Amel.parent.mkdir(parents=True, exist_ok=True)
out_mk_Amel.parent.mkdir(parents=True, exist_ok=True)
pb_Amel.write_h5ad(out_pb_Amel)
markers_Amel.to_csv(out_mk_Amel, sep="\t", index=False)
# %%
# -------------------------------------- Acer -------------------------------------- 
pb_Acer, markers_Acer = run_tidy_sce(
    Acer4MN, label_col="leiden_res0.50",
    n_per=10, seed=0, keep_leftover=True,
    target_sum=1e4, min_logfc=1.25, min_pct=0.2,
)
out_pb_Acer = Path("./matrix/MN_dataset/Acer_pseudocell_raw(with Amel).h5ad")
out_mk_Acer = Path("./metadata/MN_metadata/Acer_markers(with Amel).tsv")
pb_Acer.write_h5ad(out_pb_Acer)
markers_Acer.to_csv(out_mk_Acer, sep="\t", index=False)
# %%
