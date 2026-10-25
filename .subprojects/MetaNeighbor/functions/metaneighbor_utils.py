"""Helpers for reproducing the cross-species MetaNeighbor pipeline with pymn.

The flow mirrors ../../.reference/snRNA_NF/05_cross_species:
  - 01tidySCE.R   : unify raw counts -> build pseudo-cells -> select marker genes
  - 02auroc_cal.R : merge two species on shared orthologs -> run MetaNeighbor (AUROC)

pymn.MetaNeighborUS implements the same algorithm as the reference
`run_MetaNeighbor_s` (Spearman rank correlation + neighbor voting), so the
only pieces we need to re-create are the pseudo-cell construction and the
orthologous-marker-gene selection.
"""

from __future__ import annotations

import os
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import pymn


# --------------------------------------------------------------------------- #
# 1. Unify the expression matrix to raw counts
# --------------------------------------------------------------------------- #
def ensure_counts(adata, count_layer: str = "counts"):
    """Make sure `adata.X` holds raw UMI counts.

    The four datasets come from different sources, so their `.X` is not
    consistent: Acer stores log-normalised counts in `.X` (raw counts in
    `layers["counts"]`), Hsal has no `.X` at all (raw counts in
    `layers["counts"]`), while Amel/Dmel already store raw counts in `.X`.
    We always prefer an explicit `counts` layer when it exists.
    """
    if count_layer in adata.layers:
        adata.X = adata.layers[count_layer].copy()
    elif adata.X is None:
        raise ValueError(
            f"adata has neither an explicit '{count_layer}' layer nor a .X matrix"
        )
    return adata


# --------------------------------------------------------------------------- #
# 2. Harmonise two species onto a shared set of one-to-one orthologs
# --------------------------------------------------------------------------- #
def drop_duplicate_genes(adata):
    """Keep the first occurrence of each `var_name` (some gene ids repeat)."""
    names = pd.Index(adata.var_names)
    if names.is_unique:
        return adata
    return adata[:, ~names.duplicated()].copy()


def harmonize_orthologs(ref, other, ortho, ref_col, other_col):
    """Subset two AnnData objects to shared 1:1 orthologs and rename genes.

    Genes in `other` are renamed to their `ref` ortholog id, so that both
    objects afterwards share the exact same `var_names` (the Acer gene ids).

    Parameters
    ----------
    ref, other : AnnData
        Count matrices of the two species (ref = Acer).
    ortho : pd.DataFrame
        One-to-one ortholog table (columns include `ref_col` and `other_col`).
    ref_col, other_col : str
        Column names of the two species inside `ortho`.

    Returns
    -------
    (ref_sub, other_sub) : tuple of AnnData, with identical ordered `var_names`.
    """
    ref = drop_duplicate_genes(ref)
    other = drop_duplicate_genes(other)

    pair = ortho[[ref_col, other_col]].dropna()
    pair = pair.drop_duplicates(subset=ref_col).drop_duplicates(subset=other_col)

    ref_keep = pair[pair[ref_col].isin(ref.var_names)]
    other_keep = pair[pair[other_col].isin(other.var_names)]
    pair = ref_keep[ref_keep[other_col].isin(other_keep[other_col])]

    ref_sub = ref[:, ref.var_names.isin(pair[ref_col])].copy()
    other_sub = other[:, other.var_names.isin(pair[other_col])].copy()

    other_map = dict(zip(pair[other_col], pair[ref_col]))
    other_sub.var_names = [other_map[g] for g in other_sub.var_names]

    genes = sorted(set(ref_sub.var_names) & set(other_sub.var_names))
    ref_sub = ref_sub[:, genes].copy()
    other_sub = other_sub[:, genes].copy()
    return ref_sub, other_sub


# --------------------------------------------------------------------------- #
# 3. Pseudo-cell construction (10 randomly chosen cells summed per pseudo-cell)
# --------------------------------------------------------------------------- #
def make_pseudocells(adata, label_col: str, n_cells: int = 10, seed: int = 0):
    """Aggregate cells into pseudo-cells by summing UMI counts.

    For a cluster with `cnt` cells we create ``cnt // n_cells`` pseudo-cells,
    each built from `n_cells` randomly chosen cells (without replacement);
    the remaining `cnt % n_cells` cells are discarded. Clusters with <=
    `n_cells` cells collapse into a single pseudo-cell. This reproduces
    ``sample(c(rep(seq(1, floor(cnt/10)), 10), rep(0, cnt %% 10)))`` in
    ``01tidySCE.R``.

    Returns an AnnData whose `.X` holds pseudo-cell summed counts, `obs["cell_type"]`
    holds the label and `var_names` are inherited from `adata`.
    """
    rng = np.random.default_rng(seed)
    labels = adata.obs[label_col].astype(str).values
    uniq, inv = np.unique(labels, return_inverse=True)

    group_id = np.full(labels.shape[0], -1, dtype=np.int64)
    for gi in range(len(uniq)):
        idx = np.where(inv == gi)[0]
        cnt = idx.shape[0]
        if cnt <= n_cells:
            group_id[idx] = 0
        else:
            k = cnt // n_cells
            assign = np.concatenate(
                [np.tile(np.arange(k), n_cells), np.full(cnt % n_cells, -1)]
            )
            group_id[idx] = rng.permutation(assign)

    keep = group_id >= 0
    n_cells_kept = int(keep.sum())
    lab_idx = inv[keep]
    grp_idx = group_id[keep]

    n_groups = grp_idx.max() + 1
    key = lab_idx.astype(np.int64) * n_groups + grp_idx
    uniq_key, codes = np.unique(key, return_inverse=True)
    n_pseudo = uniq_key.shape[0]
    ct_per = uniq[uniq_key // n_groups]
    grp_per = uniq_key % n_groups

    membership = sp.csr_matrix(
        (np.ones(n_cells_kept), (codes, np.arange(n_cells_kept))),
        shape=(n_pseudo, n_cells_kept),
    )
    pseudo_X = membership @ adata[keep].X

    pb = ad.AnnData(X=pseudo_X.tocsr())
    pb.obs_names = [f"{ct}_{g}" for ct, g in zip(ct_per, grp_per)]
    pb.obs["cell_type"] = np.asarray(ct_per, dtype=str)
    pb.obs["pseudo_group"] = grp_per
    pb.var_names = adata.var_names
    return pb


# --------------------------------------------------------------------------- #
# 4. Orthologous marker genes (Python analogue of Seurat::FindAllMarkers)
# --------------------------------------------------------------------------- #
def select_marker_genes(
    pb,
    label_col: str = "cell_type",
    n_top: int = 30,
    min_logfc: float = 1.25,
    max_padj: float = 0.05,
):
    """Select cluster marker genes from a pseudo-cell AnnData.

    Pseudo-cell counts are CP10K + log1p normalised, then a Wilcoxon rank-sum
    test (``scanpy.tl.rank_genes_groups``) is run per cluster. Genes with
    adjusted p-value < `max_padj` and log2 fold change >= `min_logfc` are kept,
    and the top `n_top` genes per cluster (by score) are returned.
    """
    labels = pb.obs[label_col].astype(str)
    counts = labels.value_counts()
    valid = counts[counts >= 2].index
    if len(valid) < 2:
        return set()

    tmp = pb[labels.isin(valid)].copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
        sc.tl.rank_genes_groups(
            tmp, groupby=label_col, method="wilcoxon", n_genes=tmp.n_vars
        )
    df = sc.get.rank_genes_groups_df(tmp, group=None)
    df = df[df["pvals_adj"] < max_padj]
    df = df[df["logfoldchanges"] >= min_logfc]
    df = df.sort_values("scores", ascending=False).groupby("group", observed=True)
    df = df.head(n_top)
    return set(df["names"])


# --------------------------------------------------------------------------- #
# 5. End-to-end single pair (Acer vs another species)
# --------------------------------------------------------------------------- #
def run_metaneighbor_pair(
    ref,
    ref_label: str,
    other,
    other_label: str,
    ortho,
    ref_col: str,
    other_col: str,
    ref_name: str,
    other_name: str,
    out_dir: str,
    count_layer: str = "counts",
    n_cells: int = 10,
    seed: int = 0,
    n_top: int = 30,
    min_logfc: float = 1.25,
    max_padj: float = 0.05,
    fast_version: bool = False,
    save_inputs: bool = False,
):
    """Run the full pre-processing + AUROC step for one species pair.

    Returns the cell-type x cell-type AUROC DataFrame and writes
    ``auroc.tsv`` (long format) plus ``auroc.matrix.tsv`` to `out_dir`.
    """
    os.makedirs(out_dir, exist_ok=True)

    ref = ensure_counts(ref, count_layer=count_layer)
    other = ensure_counts(other, count_layer=count_layer)

    ref_h, other_h = harmonize_orthologs(ref, other, ortho, ref_col, other_col)
    print(f"[{ref_name} vs {other_name}] shared ortholog genes: {ref_h.n_vars}")

    ref_pb = make_pseudocells(ref_h, ref_label, n_cells=n_cells, seed=seed)
    other_pb = make_pseudocells(other_h, other_label, n_cells=n_cells, seed=seed)
    print(
        f"[{ref_name} vs {other_name}] pseudo-cells: "
        f"{ref_name}={ref_pb.n_obs}, {other_name}={other_pb.n_obs}"
    )

    markers = select_marker_genes(
        ref_pb, n_top=n_top, min_logfc=min_logfc, max_padj=max_padj
    )
    markers |= select_marker_genes(
        other_pb, n_top=n_top, min_logfc=min_logfc, max_padj=max_padj
    )

    merged = ad.concat(
        [ref_pb, other_pb],
        label="study",
        keys=[ref_name, other_name],
        join="inner",
        index_unique="-",
    )
    merged.obs["study"] = merged.obs["study"].astype(str)
    merged.obs["cell_type"] = merged.obs["cell_type"].astype(str)

    var_genes = [g for g in merged.var_names if g in markers]
    print(f"[{ref_name} vs {other_name}] marker genes used: {len(var_genes)}")

    pymn.MetaNeighborUS(
        merged,
        study_col="study",
        ct_col="cell_type",
        var_genes=var_genes,
        symmetric_output=True,
        node_degree_normalization=True,
        fast_version=fast_version,
        mn_key="MetaNeighborUS",
    )
    auroc = merged.uns["MetaNeighborUS"]

    if save_inputs:
        ref_pb.write_h5ad(os.path.join(out_dir, f"pseudocells_{ref_name}.h5ad"))
        other_pb.write_h5ad(os.path.join(out_dir, f"pseudocells_{other_name}.h5ad"))

    pd.Series(sorted(var_genes)).to_csv(
        os.path.join(out_dir, "marker_genes.txt"), index=False, header=False
    )
    auroc.to_csv(os.path.join(out_dir, "auroc.matrix.tsv"), sep="\t")
    auroc.stack().rename("AUROC").reset_index().rename(
        columns={"level_0": "train", "level_1": "test"}
    ).to_csv(os.path.join(out_dir, "auroc.tsv"), sep="\t", index=False)

    return auroc
