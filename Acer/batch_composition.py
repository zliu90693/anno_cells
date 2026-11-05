#! 统计Acer.obs["leiden_neo_res_0.50"]中，各个簇中的三个batch占比，并进行可视化
# %%
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
# %%
Acer = sc.read_h5ad("./matrix/02_Acer_transfromHY.h5ad")
# %%
cluster_key = "leiden_neo_res_0.50"

# 原始 batch 实际有 5 个水平（drone/worker 各两个重复 + queen），
# 按处理组合并成 3 组：drone / worker / queen
Acer.obs["batch_group"] = pd.Categorical(
    Acer.obs["batch"].astype(str).str.replace(r"_rep\d+$", "", regex=True),
    categories=["drone", "worker", "queen"],
    ordered=True,
)

batch_key = "batch_group"  # 若想看原始 5 个 batch，改为 "batch"

cluster_order = sorted(Acer.obs[cluster_key].unique(), key=int)

# 每个 cluster 内各 batch 的细胞数和占比（按行归一化，行和为 1）
counts = pd.crosstab(Acer.obs[cluster_key], Acer.obs[batch_key]).loc[cluster_order]
proportions = counts.div(counts.sum(axis=1), axis=0)

counts.to_csv("./metadata/batch_composition_leiden_neo_res_0.50_counts.csv")
proportions.to_csv("./metadata/batch_composition_leiden_neo_res_0.50_proportions.csv")
proportions.round(3)
# %%
# 堆叠柱状图：每个 cluster 内各 batch 的占比
fig, ax = plt.subplots(figsize=(16, 6))
proportions.plot(
    kind="bar",
    stacked=True,
    ax=ax,
    width=0.9,
    color=["#4C72B0", "#DD8452", "#55A868"],
)
ax.set_xlabel("Leiden cluster (resolution 0.50)")
ax.set_ylabel("Proportion within cluster")
ax.set_title("Batch composition per Leiden cluster\n(drone vs worker vs queen)")
ax.legend(title="Batch group", bbox_to_anchor=(1.01, 1), loc="upper left")
ax.set_ylim(0, 1)
plt.xticks(rotation=90, fontsize=8)
plt.tight_layout()
plt.savefig(
    "./figures/composition/batch_composition_leiden_neo_res_0.50.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
# %%
# 热图：更直观地看每个 cluster 的组成
fig, ax = plt.subplots(figsize=(6, 14))
im = ax.imshow(proportions.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
ax.set_xticks(
    range(proportions.shape[1]), proportions.columns, rotation=45, ha="right"
)
ax.set_yticks(range(proportions.shape[0]), proportions.index, fontsize=7)
ax.set_xlabel("Batch group")
ax.set_ylabel("Leiden cluster (resolution 0.50)")
fig.colorbar(im, ax=ax, label="Proportion within cluster")
plt.tight_layout()
plt.savefig(
    "./figures/composition/batch_composition_leiden_neo_res_0.50_heatmap.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
# %%
