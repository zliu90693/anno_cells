# %%
import pandas as pd
# %%
s3_full_cluster_marker = pd.read_excel("./metadata/aba9869_table_s3.xlsx", skiprows=1)[["gene ID", "fly homolog", "human homolog", "manual annotation"]]
s3_full_cluster_marker
# %%
s4_neuron_marker = pd.read_excel("./metadata/aba9869_table_s4.xlsx", skiprows=1)[["gene ID", "fly homolog", "human homolog", "manual annotation"]]
s4_neuron_marker
# %%
s5_glia_marker = pd.read_excel("./metadata/aba9869_table_s5.xlsx", skiprows=1)[["gene ID", "fly homolog", "human homolog", "manual annotation"]]
s5_glia_marker
# %%
s6_astr_marker = pd.read_excel("./metadata/aba9869_table_s6.xlsx", skiprows=1)[["gene ID", "fly homolog", "human homolog", "manual annotation"]]
s6_astr_marker
# %%
# 根据s2得到
KC_markers = ["mub", "Pka-C1", "PLCe", "E75", "Mblk-1", "IP3R"]
DAN_markers = ["ple", "Vmat", "DAT"]
ASTg_markers = ["Eaat1", "Gs2", "Rh50", "Galphai", "Gat"]
Ehg_markers = ["egr", "Tsf1", "Idgf4", "Slc5eg"]
Pg_markers = ["vkg", "Tret", "SPARC"]
Ctxg_markers = ["wrapper", "zyd"]
IPC_markers = ["Ilp1"]
HM_markers = ["Hml", "Fer2LCH"]
# %%
def check_marker(gene_name_list, df):
    cols_to_check = ["fly homolog", "human homolog", "manual annotation"]
    mask_exact = df[cols_to_check].isin(gene_name_list).any(axis=1)
    kc_results_exact = df[mask_exact]
    return kc_results_exact
# %%
print(check_marker(KC_markers, s3_full_cluster_marker))
print(check_marker(DAN_markers, s3_full_cluster_marker))
print(check_marker(ASTg_markers, s3_full_cluster_marker))
print(check_marker(Ehg_markers, s3_full_cluster_marker))
print(check_marker(Pg_markers, s3_full_cluster_marker))
print(check_marker(Ctxg_markers, s3_full_cluster_marker))
print(check_marker(IPC_markers, s3_full_cluster_marker))
print(check_marker(HM_markers, s3_full_cluster_marker))
# %%
print(check_marker(KC_markers, s4_neuron_marker))
print(check_marker(DAN_markers, s4_neuron_marker))
print(check_marker(ASTg_markers, s4_neuron_marker))
print(check_marker(Ehg_markers, s4_neuron_marker))
print(check_marker(Pg_markers, s4_neuron_marker))
print(check_marker(Ctxg_markers, s4_neuron_marker))
print(check_marker(IPC_markers, s4_neuron_marker))
print(check_marker(HM_markers, s4_neuron_marker))
# %%
print(check_marker(KC_markers, s5_glia_marker))
print(check_marker(DAN_markers, s5_glia_marker))
print(check_marker(ASTg_markers, s5_glia_marker))
print(check_marker(Ehg_markers, s5_glia_marker))
print(check_marker(Pg_markers, s5_glia_marker))
print(check_marker(Ctxg_markers, s5_glia_marker))
print(check_marker(IPC_markers, s5_glia_marker))
print(check_marker(HM_markers, s5_glia_marker))
# %%
print(check_marker(KC_markers, s6_astr_marker))
print(check_marker(DAN_markers, s6_astr_marker))
print(check_marker(ASTg_markers, s6_astr_marker))
print(check_marker(Ehg_markers, s6_astr_marker))
print(check_marker(Pg_markers, s6_astr_marker))
print(check_marker(Ctxg_markers, s6_astr_marker))
print(check_marker(IPC_markers, s6_astr_marker))
print(check_marker(HM_markers, s6_astr_marker))
# %%
# import scanpy as sc
# Hsal_50 = sc.read_h5ad("./Hsal50_umap_leiden.h5ad")
# Hsal_50
# %%
