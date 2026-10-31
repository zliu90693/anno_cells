from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from scipy.stats import rankdata

# ---------------------------------------------------------------------------
# 模块 2：同源基因映射
# ---------------------------------------------------------------------------
def read_orthologs(path, col_a="Apis_cerana", col_b="Apis_mellifera"):
    """返回 {col_a 基因 ID: col_b 基因 ID}，要求 1:1。"""
    ortho = (
        pd.read_csv(path, sep="\t")
        [[col_a, col_b]]
        .astype(str)
        .drop_duplicates()
    )
    if ortho[col_a].duplicated().any() or ortho[col_b].duplicated().any():
        raise ValueError("ortholog table is not 1:1")
    return dict(zip(ortho[col_a], ortho[col_b]))


def harmonize_genes(acer, amel, a2m, label_col="cell_type"):
    """把 Acer 基因名映射成 Amel ID，并按共同基因对齐。"""
    acer = acer.copy()
    acer.var_names = [a2m.get(g, g) for g in acer.var_names]
    if not acer.var_names.is_unique:
        raise ValueError("mapped Acer gene names are not unique")

    genes = [g for g in amel.var_names if g in set(acer.var_names)]
    if len(genes) == 0:
        raise ValueError("no common genes after ortholog mapping")

    acer = acer[:, genes].copy()
    amel = amel[:, genes].copy()
    return acer, amel, genes

# ---------------------------------------------------------------------------
# 模块 3：合并物种 + pheno
# ---------------------------------------------------------------------------
def merge_species(acer, amel, name_a="Acer", name_b="Amel",
                label_col="cell_type"):
    X = sp.vstack([acer.X, amel.X], format="csr").astype(np.float64)

    names = (
        [f"{name_a}.{x}" for x in acer.obs_names]
        + [f"{name_b}.{x}" for x in amel.obs_names]
    )
    study_id = np.array(
        [name_a] * acer.n_obs + [name_b] * amel.n_obs
    )
    cell_type = np.array(
        [f"{name_a}.{x}" for x in acer.obs[label_col].astype(str)]
        + [f"{name_b}.{x}" for x in amel.obs[label_col].astype(str)]
    )

    pheno = pd.DataFrame({
        "Sample_ID": names,
        "Study_ID": study_id,
        "Celltype": cell_type,
    }).set_index("Sample_ID")

    return X, pheno

# ---------------------------------------------------------------------------
# 模块 4：联合 CP10K + 按基因 z-score
# ---------------------------------------------------------------------------
def joint_normalize_zscore(X, target_sum=1000.0):
    """X: cells x genes，返回 cells x genes 的 z-score。"""
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0

    Xn = np.asarray(
        X.multiply(target_sum / lib[:, None]).todense(),
        dtype=np.float64,
    )
    mu = Xn.mean(axis=0, keepdims=True)
    sd = Xn.std(axis=0, ddof=1, keepdims=True)

    Z = np.zeros_like(Xn)
    np.divide(Xn - mu, sd, out=Z, where=sd != 0)
    return Z

# ---------------------------------------------------------------------------
# 模块 5：marker 选择
# ---------------------------------------------------------------------------
def select_markers(marker_df, mapping=None, top_n=30,
                   min_logfc=1.25, max_padj=0.05):
    d = marker_df.copy()
    if mapping is not None:
        d["gene"] = d["gene"].map(mapping)
    d = d.dropna(subset=["gene"])
    d = d[(d["p_val_adj"] < max_padj) & (d["avg_log2FC"] >= min_logfc)]

    # 模拟 dplyr::top_n(30, wt=avg_log2FC)，并列时保留
    d["_rank"] = (
        d.groupby("cluster")["avg_log2FC"]
        .rank(method="min", ascending=False)
    )
    d = d[d["_rank"] <= top_n]
    return set(d["gene"])


def build_marker_set(mk_a, mk_b, a2m, genes,
                     top_n=30, min_logfc=1.25, max_padj=0.05):
    a = select_markers(mk_a, mapping=a2m, top_n=top_n,
                       min_logfc=min_logfc, max_padj=max_padj)
    b = select_markers(mk_b, mapping=None, top_n=top_n,
                       min_logfc=min_logfc, max_padj=max_padj)
    return (a | b) & set(genes)


# ---------------------------------------------------------------------------
# 模块 6：MetaNeighbor AUROC（对应 run_MetaNeighbor_s / _p）
# ---------------------------------------------------------------------------
def _rank_global(a):
    """R: rank(cor.dat, ties.method='average', na.last='keep')。"""
    flat = a.ravel()
    mask = ~np.isnan(flat)
    r = np.full(flat.shape, np.nan, dtype=np.float64)
    r[mask] = rankdata(flat[mask], method="average")
    r = r.reshape(a.shape)
    r[np.isnan(r)] = 0.0
    m = r.max()
    return r / m if m > 0 else r


def _rank_axis0(a):
    """R: apply(abs(predicts.temp), 2, rank, na.last='keep')。"""
    out = np.full(a.shape, np.nan, dtype=np.float64)
    for j in range(a.shape[1]):
        col = a[:, j]
        mask = ~np.isnan(col)
        out[mask, j] = rankdata(col[mask], method="average")
    return out


def run_metaneighbor(marker_genes, expr, gene_names, celltypes, pheno,
                     method="spearman"):
    """
    expr: genes x cells（与 R 的 data.zscore 方向一致）
    pheno: index 与 expr 的列一一对应，含 Study_ID / Celltype
    """
    assert expr.shape[0] == len(gene_names)
    assert expr.shape[1] == len(pheno)
    assert method in ("spearman", "pearson")

    ct_values = pheno["Celltype"].values
    n_cells = expr.shape[1]
    n_types = len(celltypes)

    # cell.labels
    labels = np.zeros((n_cells, n_types), dtype=np.float64)
    for i, ct in enumerate(celltypes):
        labels[:, i] = (ct_values == ct)

    # 只保留 marker genes
    Xm = expr[np.isin(gene_names, list(marker_genes)), :]
    if Xm.shape[0] < 2:
        raise ValueError("fewer than 2 marker genes overlap the matrix")

    # Spearman: 先按细胞（列）在基因间排秩，再算 Pearson
    if method == "spearman":
        Xm = rankdata(Xm, axis=0, method="average")

    # R: cor(data[marker, ], method=...)
    cor_dat = np.corrcoef(Xm, rowvar=False)  # cells x cells

    # R: rank.dat
    rank_dat = _rank_global(cor_dat)

    # R: sumin / sumall
    sumin = rank_dat @ labels
    colsum = rank_dat.sum(axis=0)
    sumall = np.tile(colsum[:, None], (1, n_types))
    predicts = np.divide(
        sumin, sumall, out=np.zeros_like(sumin), where=sumall != 0
    )

    cell_nv = np.zeros((n_types, n_types), dtype=np.float64)

    for i, ct in enumerate(celltypes):
        is_ct = (ct_values == ct)
        studies = pd.unique(pheno.loc[is_ct, "Study_ID"])
        in_study = pheno["Study_ID"].isin(studies).values
        target = is_ct[in_study]

        # R: apply(abs(predicts.temp), 2, rank)
        pred_rank = _rank_axis0(np.abs(predicts[in_study, :]))
        pred_rank[~target, :] = 0.0

        np_ = int(target.sum())          # 目标 cell type 的细胞数
        nn_ = int((~target).sum())       # 其他 cell type 的细胞数
        if np_ == 0 or nn_ == 0:
            cell_nv[i, :] = np.nan
            continue

        p = np.nansum(pred_rank, axis=0)
        cell_nv[i, :] = (p / np_ - (np_ + 1) / 2.0) / nn_

    # R: cell.NV = (cell.NV + t(cell.NV)) / 2
    cell_nv = (cell_nv + cell_nv.T) / 2.0
    return pd.DataFrame(cell_nv, index=celltypes, columns=celltypes)


# ---------------------------------------------------------------------------
# 模块 7：整理与输出
# ---------------------------------------------------------------------------
def melt_auroc(auroc_p, auroc_s):
    rows, cols = np.meshgrid(
        auroc_p.index, auroc_p.columns, indexing="ij"
    )
    return pd.DataFrame({
        "Var1": rows.ravel(),
        "Var2": cols.ravel(),
        "AUROCp": auroc_p.values.ravel(),
        "AUROCs": auroc_s.values.ravel(),
    })

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

def run_auroc_pipeline(
    acer_path, refsp_path,
    acer_marker_path, refsp_marker_path,
    ortho_path, out_path=None,
    name_a=None, name_b=None,
    label_col="cell_type",
    ortho_col_a=None,
    ortho_col_b=None,
    top_n=30, min_logfc=1.25, max_padj=0.05,
    target_sum=1000.0,
):
    # 1. read
    acer = ad.read_h5ad(acer_path)
    refsp = ad.read_h5ad(refsp_path)
    mk_a = pd.read_csv(acer_marker_path, sep="\t")
    mk_b = pd.read_csv(refsp_marker_path, sep="\t")

    # 2. ortholog + gene alignment
    a2m = read_orthologs(ortho_path, ortho_col_a, ortho_col_b)
    acer, refsp, genes = harmonize_genes(acer, refsp, a2m, label_col=label_col)

    # 3. merge + pheno
    X, pheno = merge_species(acer, refsp, name_a, name_b, label_col=label_col)

    # 4. joint CP10K + z-score
    Z = joint_normalize_zscore(X, target_sum=target_sum)
    expr = Z.T                    # genes x cells
    del Z, X

    # 5. markers
    markers = build_marker_set(
        mk_a, mk_b, a2m, genes,
        top_n=top_n, min_logfc=min_logfc, max_padj=max_padj,
    )
    celltypes = list(pd.unique(pheno["Celltype"]))

    # 6. AUROC
    auroc_s = run_metaneighbor(
        markers, expr, genes, celltypes, pheno, method="spearman"
    )
    auroc_p = run_metaneighbor(
        markers, expr, genes, celltypes, pheno, method="pearson"
    )

    # 7. melt + save
    res = melt_auroc(auroc_p, auroc_s)
    if out_path is not None:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        res.to_csv(out_path, sep="\t", index=False)

    return res, auroc_s, auroc_p