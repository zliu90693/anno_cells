# %%
from functions.auroc_cal import *
# %%

# %%
res, _, _ = run_auroc_pipeline(
    acer_path="matrix/MN_dataset/Acer_pseudocell_raw(with Amel).h5ad",
    refsp_path="matrix/MN_dataset/Amel_pseudocell_raw.h5ad",
    acer_marker_path="metadata/MN_metadata/Acer_markers(with Amel).tsv",
    refsp_marker_path="metadata/MN_metadata/Amel_markers.tsv",
    ortho_path="orthologs/Acer_vs_Amel_121_4mn.tsv",
    out_path="metadata/MN_metadata/auroc_Acer_vs_Amel.tsv",
    name_a="Acer", 
    name_b="Amel",
    ortho_col_a="Apis_cerana",
    ortho_col_b="Apis_mellifera"
)
print(res.shape)
print(res.head())
# %%
# read_orthologs("./orthologs/Acer_vs_Amel_121_4mn.tsv")
# %%
