#? 问题1, 在 NCBI 中, GSM4013392 和 GSM4013393 的 internal name 都是 gam.nb.e13b1, 我将后者作为 gam.nb.e13b2 处理了, 可能是我处理反了导致的细胞簇不匹配吗?
# GEO 直接提供的表达矩阵
# %%
import scanpy as sc

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