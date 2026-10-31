from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc


def load_counts(adata, label_col, count_layer=None, copy=True):
    """返回 (adata, counts_csr, labels, var_names)。"""
    a = adata.copy() if copy else adata
    X = a.layers[count_layer] if count_layer is not None else a.X
    if X is None:
        raise ValueError("no expression matrix found")
    labels = a.obs[label_col].astype(str).values
    return a, X.tocsr(), labels, a.var_names


def assign_pseudo_groups(labels, n_per=10, seed=0):
    """每 n_per 个细胞一组；余数细胞标记为 0。"""
    rng = np.random.default_rng(seed)
    out = np.full(labels.size, -1, dtype=np.int64)
    for ct in pd.unique(labels):
        idx = np.flatnonzero(labels == ct)
        cnt = idx.size
        if cnt <= n_per:
            out[idx] = 1
        else:
            k = cnt // n_per
            assign = np.concatenate([
                np.repeat(np.arange(1, k + 1), n_per),
                np.zeros(cnt % n_per, dtype=np.int64),
            ])
            out[idx] = rng.permutation(assign)
    return out


def build_pseudocells(counts, labels, group_ids, var_names,
                    keep_leftover=False):
    """按 (label, group) 求和，返回 raw-count AnnData。"""
    mask = group_ids >= 0 if keep_leftover else group_ids > 0 # 如果显式指定 keep_leftover=True，将余数中凑不满 10 个的真细胞作为伪细胞处理，如果默认，则丢弃凑不满 10 个的真细胞
    keep = np.flatnonzero(mask)

    keys = np.array([f"{labels[i]}\x1f{group_ids[i]}" for i in keep])
    uniq_keys, key_codes = np.unique(keys, return_inverse=True)
    membership = sp.csr_matrix(
        (np.ones(keep.size), (key_codes, np.arange(keep.size))),
        shape=(uniq_keys.size, keep.size),
    )
    pb_X = (membership @ counts[keep]).tocsr()

    ct_per = np.array([k.split("\x1f")[0] for k in uniq_keys])
    grp_per = np.array([int(k.split("\x1f")[1]) for k in uniq_keys])

    pb = ad.AnnData(X=pb_X.astype(np.float32))
    pb.obs_names = [f"{ct}_{g}" for ct, g in zip(ct_per, grp_per)]
    pb.obs["cell_type"] = ct_per
    pb.obs["pseudo_group"] = grp_per
    pb.var_names = var_names
    return pb


def compute_markers(pb, label_col="cell_type", target_sum=1e4,
                    min_logfc=1.25, max_padj=0.05, min_pct=0.2,
                    p_adj_method="bh"):
    """在 pseudo-cell 上计算 marker；不修改输入 pb，不写文件。"""
    a = pb.copy()
    sc.pp.normalize_total(a, target_sum=target_sum)
    sc.pp.log1p(a)

    sc.tl.rank_genes_groups(
        a, groupby=label_col, method="wilcoxon",
        use_raw=False, pts=True,
    )
    m = sc.get.rank_genes_groups_df(a, group=None).rename(columns={
        "names": "gene",
        "group": "cluster",
        "pvals_adj": "p_val_adj",
        "pvals": "p_val",
        "logfoldchanges": "avg_log2FC",
    })
    if p_adj_method == "bonferroni":
        m["p_val_adj"] = np.minimum(m["p_val"] * a.n_vars, 1.0)

    m["min_pct"] = m[["pct_nz_group", "pct_nz_reference"]].max(axis=1)
    m = m[
        (m["avg_log2FC"] >= min_logfc)
        & (m["p_val_adj"] < max_padj)
        & (m["min_pct"] >= min_pct)
    ]
    return m[["cluster", "gene", "p_val_adj", "p_val", "avg_log2FC"]]


def run_tidy_sce(adata, label_col, n_per=10, seed=0,
                 keep_leftover=False, **marker_kwargs):
    """01tidySCE 的完整流程；返回 (pb_raw, markers)。"""
    a, counts, labels, var_names = load_counts(adata, label_col)
    gid = assign_pseudo_groups(labels, n_per=n_per, seed=seed)
    pb = build_pseudocells(
        counts, labels, gid, var_names, keep_leftover=keep_leftover
    )
    markers = compute_markers(pb, **marker_kwargs)
    return pb, markers