# %%
"""Alluvial plot of Leiden cluster correspondence across resolutions.

Dots = clusters at each resolution.
Lines/ribbons = cells shared between clusters, colored by the 0.25 cluster
so that the same origin can be followed across all three columns.
"""
# %%
from pathlib import Path as FilePath

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import anndata as ad

HERE = FilePath(__file__).resolve().parent
H5AD = HERE / "matrix/01_Acer_transfromHY.h5ad"
RES_COLS = ["leiden_res0.25", "leiden_res0.50", "leiden_res1.00"]
RES_LABELS = ["0.25", "0.50", "1.00"]
OUT = HERE / "figures/leiden_resolution_alluvial.png"
# %%
# ---------------------------------------------------------------- read obs
adata = ad.read_h5ad(H5AD, backed="r")
obs = adata.obs[RES_COLS].astype(str)
adata.file.close()
# %%
N = len(obs)


def cat_order(series):
    cats = list(pd.Categorical(series).categories)
    try:
        return sorted(cats, key=float)
    except ValueError:
        return sorted(cats)


orders = [cat_order(obs[c]) for c in RES_COLS]
counts = [obs[c].value_counts() for c in RES_COLS]

# ---------------------------------------------------------------- layout
GAP = 0.004
avail = [1.0 - GAP * (len(o) - 1) for o in orders]

# node y-extents: {cluster: (y0, y1)}, y grows downwards (axis inverted later)
pos = []
for o, cnt, av in zip(orders, counts, avail):
    p, y = {}, 0.0
    for c in o:
        h = cnt.get(c, 0) / N * av
        p[c] = (y, y + h)
        y += h + GAP
    pos.append(p)

# ---------------------------------------------------------------- colors
palette = (
    list(plt.get_cmap("tab20").colors)
    + list(plt.get_cmap("tab20b").colors)
    + list(plt.get_cmap("tab20c").colors)
)
color = {c: palette[i % len(palette)] for i, c in enumerate(orders[0])}

# ---------------------------------------------------------------- ribbons
def ribbon(ax, x0, y0a, y0b, x1, y1a, y1b, col, alpha=0.45):
    xm = (x0 + x1) / 2
    verts = [
        (x0, y0a),
        (xm, y0a), (xm, y1a), (x1, y1a),
        (x1, y1b),
        (xm, y1b), (xm, y0b), (x0, y0b),
        (x0, y0a),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    ax.add_patch(
        PathPatch(Path(verts, codes), facecolor=col, edgecolor="none",
                  alpha=alpha, zorder=1)
    )


fig, ax = plt.subplots(figsize=(10, 12))
XS = [0.0, 1.4, 2.8]
PAD = 0.20  # gap between dot and ribbon start

df = obs.copy()

# transition 0 -> 1
ct01 = pd.crosstab(df[RES_COLS[0]], df[RES_COLS[1]]).reindex(
    index=orders[0], columns=orders[1], fill_value=0
)
tgt = {c: pos[1][c][0] for c in orders[1]}
block01 = {}  # (c0, c1) -> (y0, y1) on the right side of column 1
for c0 in orders[0]:
    src = pos[0][c0][0]
    for c1 in orders[1]:
        n = int(ct01.loc[c0, c1])
        if n == 0:
            continue
        h0 = n / N * avail[0]
        h1 = n / N * avail[1]
        y0a, y0b = src, src + h0
        y1a, y1b = tgt[c1], tgt[c1] + h1
        if (c0, c1) not in block01:
            block01[(c0, c1)] = [y1a, y1b]
        else:
            block01[(c0, c1)][1] = y1b
        ribbon(ax, XS[0] + PAD, y0a, y0b, XS[1] - PAD, y1a, y1b, color[c0])
        src += h0
        tgt[c1] += h1

# transition 1 -> 2 (split each column-1 block by column-2 cluster,
# keeping the same origin color)
ct12 = df.groupby(RES_COLS).size()
tgt = {c: pos[2][c][0] for c in orders[2]}
for c1 in orders[1]:
    for c0 in orders[0]:
        if (c0, c1) not in block01:
            continue
        src = block01[(c0, c1)][0]
        for c2 in orders[2]:
            n = int(ct12.get((c0, c1, c2), 0))
            if n == 0:
                continue
            h0 = n / N * avail[1]
            h1 = n / N * avail[2]
            y1a, y1b = src, src + h0
            y2a, y2b = tgt[c2], tgt[c2] + h1
            ribbon(ax, XS[1] + PAD, y1a, y1b, XS[2] - PAD, y2a, y2b, color[c0])
            src += h0
            tgt[c2] += h1

# ---------------------------------------------------------------- dots
max_count = max(c.max() for c in counts)
for j, (o, cnt, p) in enumerate(zip(orders, counts, pos)):
    for c in o:
        n = int(cnt.get(c, 0))
        if n == 0:
            continue
        y = (p[c][0] + p[c][1]) / 2
        s = 25 + 220 * (n / max_count) ** 0.6
        ax.scatter(XS[j], y, s=s, color="white", edgecolor="black",
                   linewidth=0.8, zorder=3)
        ax.scatter(XS[j], y, s=s * 0.72, color=color[c] if j == 0 else "0.35",
                   edgecolor="none", zorder=4)
        ax.text(XS[j] + 0.04, y, f"{c} ({n})", ha="left", va="center",
                fontsize=6.2, zorder=5,
                bbox=dict(facecolor="white", edgecolor="none",
                          alpha=0.65, pad=0.3))

# ---------------------------------------------------------------- headers
for x, lab in zip(XS, RES_LABELS):
    ax.text(x - 0.02, -0.015, f"leiden {lab}", ha="center", va="bottom",
            fontsize=11, fontweight="bold")
ax.text(XS[0], 1.035, "dot = cluster (n cells); ribbon color follows the leiden 0.25 cluster",
        ha="left", va="top", fontsize=7, color="0.35")

ax.set_xlim(-0.12, 3.15)
ax.set_ylim(-0.05, 1.02)
ax.invert_yaxis()
ax.axis("off")
fig.savefig(OUT, dpi=300, bbox_inches="tight")
print(f"saved: {OUT}")

# %%
