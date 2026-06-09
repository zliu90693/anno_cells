# %%
#! 问题1, 在 NCBI 中, GSM4013392 和 GSM4013393 的 internal name 都是 gam.nb.e13b1, 我将后者作为 gam.nb.e13b2 处理了, 可能是我处理反了导致的细胞簇不匹配吗?
# GEO 直接提供的表达矩阵
import scanpy as sc
# %%
Hsal_GEO_txt = sc.read_text('../../Previous_studies/Sheng_SA_2020_Hsal/data/GSE135513_d30_DGE.txt.gz', delimiter='\t')
# %%
Hsal_GEO_txt = Hsal_GEO_txt.T
# %%
Hsal_GEO_txt.obs["internalname"] = Hsal_GEO_txt.obs.index.str.split("-").str[0]
Hsal_GEO_txt.obs
# %%
(Hsal_GEO_txt.obs["internalname"] == "gam.nb.e13b1").sum() # 8765
# %%
(Hsal_GEO_txt.obs["internalname"] == "gam.nb.e13b2").sum() # 7052
# %%
GSM4013392 = sc.read_10x_h5("../../fastq2matrix/Sheng_SA_2020_Hsal/cellranger-count-out/GSM4013392/outs/filtered_feature_bc_matrix.h5")
GSM4013392 # 1813
# %%
GSM4013393 = sc.read_10x_h5("../../fastq2matrix/Sheng_SA_2020_Hsal/cellranger-count-out/GSM4013393/outs/filtered_feature_bc_matrix.h5")
GSM4013393 # 2110
# %%
Hsal_51 = sc.read_h5ad("./Hsal51.h5ad")
Hsal_51.obs["internalname"] = Hsal_51.obs.index.str.split("-").str[0]
(Hsal_51.obs["internalname"] == "gam.nb.e13b1").sum() # 1906
# %%
(Hsal_51.obs["internalname"] == "gam.nb.e13b2").sum() # 2247
# %%
# ↑ 只看细胞数似乎看不出来GSM和internal name的对应关系
GSM4013392.obs["barcode"] = GSM4013392.obs.index.str.split("-").str[0]
GSM4013393.obs["barcode"] = GSM4013393.obs.index.str.split("-").str[0]
# %%
Hsal_51.obs["barcode"] = Hsal_51.obs.index.str.split("-").str[1]
# %%
Hsal_51_gam_nb_e13b1_barcodeset = set(Hsal_51.obs["barcode"][Hsal_51.obs["internalname"] == "gam.nb.e13b1"])
Hsal_51_gam_nb_e13b1_barcodeset # 1906
# %%
Hsal_51_gam_nb_e13b2_barcodeset = set(Hsal_51.obs["barcode"][Hsal_51.obs["internalname"] == "gam.nb.e13b2"])
Hsal_51_gam_nb_e13b2_barcodeset # 2247
# %%
GSM4013392_barcodeset = set(GSM4013392.obs["barcode"]) # 1813
GSM4013393_barcodeset = set(GSM4013393.obs["barcode"]) # 2110
# %%
print(len(Hsal_51_gam_nb_e13b1_barcodeset & GSM4013392_barcodeset)) # 1767
print(len(Hsal_51_gam_nb_e13b1_barcodeset & GSM4013393_barcodeset)) # 4
print(len(Hsal_51_gam_nb_e13b2_barcodeset & GSM4013392_barcodeset)) # 4
print(len(Hsal_51_gam_nb_e13b2_barcodeset & GSM4013393_barcodeset)) # 2069
# %%
#* 确定了, e13b1 对应 GSM4013392, e13b2 对应 GSM4013393

# %%
#! 问题2, 可能是Hsal 50/51没去除doublet导致的transfer结果不匹配吗?
# %%
import scrublet as scr
# %%
Hsal_50 = sc.read_h5ad("./Hsal50.h5ad")
Hsal_50.X = Hsal_50.layers["counts"].copy()
# %%
# sc.pl.embedding(
#         Hsal_50,
#         basis="tsne",
#         color="cluster_name",
#         legend_loc="on data"
#     )
# %%
sc.pp.scrublet(Hsal_50)
# %%
Hsal_me = sc.read_h5ad("./cluster_unannotated.h5ad")
# %%
Hsal_me.X = Hsal_me.layers["counts"].copy()
# %%
sc.pp.scrublet(Hsal_me)
# %%
#! 问题3, 可视化Marker基因在本次和Hsal50/51中的表达情况
# %%
Hsal_50.X = Hsal_50.layers["counts"].copy()
sc.pp.normalize_total(Hsal_50)
sc.pp.log1p(Hsal_50)
Hsal_50.layers['logcounts'] = Hsal_50.X.copy()
# %%
sc.tl.rank_genes_groups(
    Hsal_50, 
    groupby="cluster_name", 
    method="wilcoxon", 
    key_added="50_cluster_name_DEG"
)
# %%
sc.tl.dendrogram( # 纯粹为了美观和生物学逻辑。它的作用是让后续画图时，相似的 cluster 能够挨在一起，并在图上方画出一个树状分支
    Hsal_50,
    groupby="cluster_name",
)
# %%
Hsal_50.write_h5ad("./Hsal50_with_wilcoxon_DEG.h5ad")
# %%
sc.pl.rank_genes_groups_dotplot(
    Hsal_50, 
    groupby="cluster_name", 
    standard_scale="var", 
    n_genes=5, 
    key="50_cluster_name_DEG"
)
# plt.show()
# %%
sc.tl.filter_rank_genes_groups(
    Hsal_50,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="50_cluster_name_DEG",
    key_added="50_cluster_name_DEG_filtered",
)
# %%
sc.pl.rank_genes_groups_dotplot(
    Hsal_50,
    groupby="cluster_name",
    standard_scale="var",
    n_genes=5,
    key="50_cluster_name_DEG_filtered",
)
# %%
for marker in ["nSyb", "fne", "Syt1", "mub", "PLCε", "ple", "Vmat", "repo", "bdl", "GLaz", "Eaat1", "Gs2", "Rh50", "Gat", "Giα", "egr", "Tsf1", "Idgf", "vkg", "SPARC", "wrapper", "zyd", "Ilp1", "Hml", "Fer2LCH"]:
    if marker in list(Hsal_50.var.index):
        print(f"Marker: {marker} is in the list!")
# 无任何输出.
# %%
for gene in list(Hsal_50.var.index):
    if gene.startswith("LOC"):
        pass
    else:
        print(gene)
# Anamorsin
# Centrin-1
# Fut8
# Glomulin
# Ickin1-1
# Importin-4
# Traa
# Trab
# Trehalase
# Trnag-ccc
# Trypsin-1

# %%
#! 问题4, 可视化本次和Hsal50/51中同源细胞簇识别到的Marker是否一致
Hsal_me = sc.read_h5ad("./cluster_unannotated.h5ad")
Hsal_me
# %%
Hsal_me.X = Hsal_me.layers["logcounts"].copy()
# %%
sc.tl.rank_genes_groups(
    Hsal_me, 
    groupby="leiden_res1.00", 
    method="wilcoxon", 
    key_added="leiden_res1.00_DEG"
)
# %%
sc.tl.dendrogram( # 纯粹为了美观和生物学逻辑。它的作用是让后续画图时，相似的 cluster 能够挨在一起，并在图上方画出一个树状分支
    Hsal_me,
    groupby="leiden_res1.00",
)
# %%
Hsal_me.write_h5ad("./Hsal_me_with_wilcoxon_DEG.h5ad")
# %%
sc.pl.rank_genes_groups_dotplot(
    Hsal_me, 
    groupby="leiden_res1.00", 
    standard_scale="var", 
    n_genes=5, 
    key="leiden_res1.00_DEG"
)
# %%
sc.tl.filter_rank_genes_groups(
    Hsal_me,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="leiden_res1.00_DEG",
    key_added="leiden_res1.00_DEG_filtered",
)
# %%
sc.pl.rank_genes_groups_dotplot(
    Hsal_me, 
    groupby="leiden_res1.00", 
    standard_scale="var", 
    n_genes=5, 
    key="leiden_res1.00_DEG_filtered"
)
# %%
