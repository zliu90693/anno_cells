# 尝试参考best practice中的注释，但发现很难评估不同marker的表达量
# %%
import scanpy as sc
# %%
Amel = sc.read_h5ad("./matrix/01_Amel_transfromHY.h5ad")
# %%
marker_genes = {
    "hemo": ["LOC551684", "LOC411597"],
    "glia": ["LOC410151", "LOC408480"],
    "glia_cortex": ["LOC409235"], # wrapper -(flybase)-> Dmel: FBgn0025878 -(Orthofinder)-> Amel: LOC409235; zyd没找到蜜蜂的同源基因
    "glia_surface": [
        "Tret1", # "Tret1" -(flybase)-> "FBgn0050035" -(Orthofinder)-> "GeneID_551643/Tret1"
        "LOC408551" # "vkg" -(flybase)-> "FBgn0016075" -(Orthofinder)-> "LOC408551"
    ], 
    "glia_astro": [
        "LOC409109", # "ebony" -(flybase)-> "FBgn0000527" -(Orthofinder)-> "LOC409109"
        "GlnS" # "GlnS" in Amel.var["gene_ids"] -> True
    ],
    "glia_ensheathing": [
        "LOC413324", # "Idgf4" -(flybase)-> "FBgn0026415" -(Orthofinder)-> "LOC413324"
        "Tsf1", # "Tsf1" -(flybase)-> "FBgn0022355" -(Orthofinder)-> "GeneID_406078/Tsf1"
    ],
    "neuron": ["LOC410689"],
    "KC": ["LOC408804", "LOC408372"],
    "KC_l": ["Mblk-1"],
    "KC_s": ["E74"],
    "KC_FoxP": ["Foxp"],
    "OPNs": ["LOC724282", "LOC724148", "LOC410657", "LOC413466"],
    "OLCs_photoreceptors": ["LOC409937", "LOC724608"],
    "OLCs_Mi1": ["LOC724801"],
    "OLCs_T1": ["LOC409919"],
    "OLCs_PM": ["LOC410658", "LOC411315", "LOC408872"],
    "OLCs_Lawf2": ["LOC552079", "LOC412528", "LOC100577751"],
    "OLCs_Tm5c": ["Hiscl1"]
}
# %%
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in Amel.var.index: # 注意，这里的marker名称是以.var.index为标准的，而非.var["gene_ids "]
            markers_found.append(marker)
        else:
            print(f"{marker} not in Amel.var.index!!")
    marker_genes_in_data[ct] = markers_found
# %%
Amel.X = Amel.layers["logcounts"]
# %%
sc.tl.pca(Amel, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(Amel)
sc.tl.umap(Amel)
# %%
sc.pl.umap(
    Amel,
    color=["celltype", "leiden_res0.25"],
)
# %%
glia_cts = ["glia", "glia_cortex", "glia_surface", "glia_astro", "glia_ensheathing"]
for ct in glia_cts:
    print(f"{ct.upper()}:")
    sc.pl.umap(
        Amel,
        color=marker_genes_in_data[ct],
        vmin=0,
        vmax="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
        sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
        frameon=False,
        cmap="Reds",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    )
    print("\n\n\n")
# %%
neuron_cts = ["neuron", "KC", "KC_l", "KC_s", "KC_FoxP", "OPNs", "OLCs_photoreceptors", "OLCs_Mi1", "OLCs_T1", "OLCs_PM", "OLCs_Lawf2", "OLCs_Tm5c"]
for ct in neuron_cts:
    print(f"{ct.upper()}:")
    sc.pl.umap(
        Amel,
        color=marker_genes_in_data[ct],
        vmin=0,
        vmax="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
        sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
        frameon=False,
        cmap="Reds",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    )
    print("\n\n\n")
# %%
