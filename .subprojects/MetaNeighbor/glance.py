# %%
# ---------------------------------------------------------------------------
# glance.py
# Quick visualisation of the MetaNeighbor AUROC results written by 02.py.
# The output of 02.py lives in:
#   ./results/metaneighbor/<tag>/auroc.matrix.tsv
# where <tag> = "Acer-<ref_label>_vs_<other_name>-<other_label>".
# ---------------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

# %%
RES_DIR = "./results/metaneighbor"

# Acer (rows) vs Amel (cols) AUROC, one file per Acer leiden resolution.
AMEL_TAGS = [
    "Acer-leiden_res0.25_vs_Amel-cell_type",
    "Acer-leiden_res0.50_vs_Amel-cell_type",
    "Acer-leiden_res1.00_vs_Amel-cell_type",
]

HSAL_TAGS = [
    "Acer-leiden_res0.25_vs_Hsal-cluster_name",
    "Acer-leiden_res0.50_vs_Hsal-cluster_name",
    "Acer-leiden_res1.00_vs_Hsal-cluster_name"
]

DMEL_TAGS = [
    "Acer-leiden_res0.25_vs_Dmel-cell_type",
    "Acer-leiden_res0.50_vs_Dmel-cell_type",
    "Acer-leiden_res1.00_vs_Dmel-cell_type",
]

# %%
def load_auroc(tag, res_dir=RES_DIR):
    """Load one AUROC matrix (index/columns = study|cell_type labels)."""
    path = os.path.join(res_dir, tag, "auroc.matrix.tsv")
    return pd.read_csv(path, sep="\t", index_col=0)


def _cluster_order(sub):
    """Leaf order of a square AUROC sub-matrix (distance = 1 - AUROC)."""
    if sub.shape[0] < 2:
        return list(sub.index)
    d = 1.0 - sub.values
    d = (d + d.T) / 2.0
    np.fill_diagonal(d, 0.0)
    d = np.clip(d, 0.0, None)
    z = linkage(squareform(d, checks=False), method="average")
    return list(sub.index[leaves_list(z)])


def _auto_fontsize(n_rows, n_cols, figsize, factor=0.7, min_pt=2.0, max_pt=10.0):
    """Pick a label size so text roughly fits one cell, capped to a readable range."""
    row_pt = figsize[1] * 72 / max(n_rows, 1) * factor
    col_pt = figsize[0] * 72 / max(n_cols, 1) * factor
    return float(np.clip(min(row_pt, col_pt), min_pt, max_pt))


def _species_orders(auroc, ref_species="Acer"):
    """Cluster each species' own sub-matrix separately, returning both orders."""
    labels = list(auroc.index)
    ref = [l for l in labels if l.split("|", 1)[0] == ref_species]
    oth = [l for l in labels if l.split("|", 1)[0] != ref_species]
    ref_order = _cluster_order(auroc.loc[ref, ref])
    oth_order = _cluster_order(auroc.loc[oth, oth])
    return ref_order, oth_order


def plot_auroc_heatmap(tag, res_dir=RES_DIR, mode="blocks", threshold=None,
                       ref_species="Acer", figsize=(20, 18), cmap="RdBu_r",
                       label_fontsize=None, cbar_label="AUROC"):
    """Draw and save the AUROC heatmap for one run.

    mode="blocks" -> full matrix, rows/cols ordered as [ref_species block,
                     other species block]; the two species never interleave.
    mode="cross"  -> cross-species block only: rows = ref_species clusters,
                     columns = other-species clusters (like 03heatmap_pheat.R).
    mode="clustermap" -> original pymn global clustering (species interleave).

    threshold=None  -> continuous AUROC values
    threshold=0.5   -> binary high/low similarity view
    """
    auroc = load_auroc(tag, res_dir)
    ref_order, oth_order = _species_orders(auroc, ref_species)

    if mode == "blocks":
        order = ref_order + oth_order
        data = auroc.loc[order, order]
    elif mode == "cross":
        data = auroc.loc[ref_order, oth_order]
    elif mode == "clustermap":
        # Original pymn behaviour: global clustering, species interleave.
        from pymn.plotting import plotMetaNeighborUS
        cm_figsize = (9, 9)
        fs = label_fontsize if label_fontsize is not None else _auto_fontsize(
            auroc.shape[0], auroc.shape[1], cm_figsize)
        cm = plotMetaNeighborUS(
            auroc, threshold=threshold, show=False,
            figsize=cm_figsize, cmap=cmap,
            xticklabels=list(auroc.columns), yticklabels=list(auroc.index),
        )
        cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xticklabels(),
                                      rotation=90, fontsize=fs)
        cm.ax_heatmap.set_yticklabels(cm.ax_heatmap.get_yticklabels(),
                                      rotation=0, fontsize=fs)
        cm.fig.suptitle(tag, y=1.02)
        suffix = "continuous" if threshold is None else f"thr{threshold}"
        out = os.path.join(res_dir, tag, f"auroc.heatmap.{suffix}.pdf")
        cm.fig.savefig(out, bbox_inches="tight")
        print(f"saved {out}")
        return cm

    else:
        raise ValueError("mode must be 'blocks', 'cross' or 'clustermap'")

    plot_data = data if threshold is None else (data >= threshold)
    if label_fontsize is None:
        label_fontsize = _auto_fontsize(data.shape[0], data.shape[1], figsize)

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        plot_data,
        ax=ax,
        cmap=cmap,
        vmin=0, vmax=1,
        square=(mode == "cross"),
        xticklabels=list(plot_data.columns),
        yticklabels=list(plot_data.index),
        cbar_kws={"label": cbar_label},
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90,
                       fontsize=label_fontsize)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0,
                       fontsize=label_fontsize)

    if mode == "blocks":
        # separator between the two species blocks
        n_ref = len(ref_order)
        ax.axhline(n_ref, color="black", lw=1)
        ax.axvline(n_ref, color="black", lw=1)

    ax.set_title(tag, fontsize=10)
    fig.tight_layout()

    suffix = "continuous" if threshold is None else f"thr{threshold}"
    out = os.path.join(res_dir, tag, f"auroc.heatmap.{mode}.{suffix}.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")
    return data

# %%
# Continuous AUROC heatmap for every run (species blocks kept together).
# ALL_TAGS = AMEL_TAGS + HSAL_TAGS + DMEL_TAGS
for tag in HSAL_TAGS:
    if os.path.isdir(os.path.join(RES_DIR, tag)):
        plot_auroc_heatmap(tag, mode="blocks")
    else:
        print(f"missing: {tag} (run 02.py first)")

# %%
# Optional: cross-species block only (rows = Acer, columns = other species),
# mirroring 03heatmap_pheat.R.
for tag in HSAL_TAGS:
    if os.path.isdir(os.path.join(RES_DIR, tag)):
        plot_auroc_heatmap(tag, mode="cross")

# %%
# Optional: threshold at 0.5 to highlight strong cross-species matches.
for tag in HSAL_TAGS:
    if os.path.isdir(os.path.join(RES_DIR, tag)):
        plot_auroc_heatmap(tag, mode="cross", threshold=0.5)

# %%
