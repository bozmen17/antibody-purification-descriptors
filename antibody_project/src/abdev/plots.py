"""Figures for the report.

House rules applied throughout: one y-scale per panel, categorical hues assigned in a
fixed order and never cycled (a ninth group folds into "Other"), a legend whenever two
or more series share a panel, signed quantities on a diverging blue/orange pair with a
neutral midpoint, magnitudes on a single-hue ramp, and recessive grid and axes.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

# Categorical hues, in fixed assignment order.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e4e3df"

SEQUENTIAL = LinearSegmentedColormap.from_list(
    "blue_ramp", ["#f2f7fe", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#184f95"])
DIVERGING_NEG, DIVERGING_MID, DIVERGING_POS = "#2a78d6", "#f0efec", "#eb6834"

MAX_SERIES = 8


def use_house_style() -> None:
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "axes.labelcolor": INK_SOFT,
        "axes.titlecolor": INK,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "xtick.color": INK_SOFT,
        "ytick.color": INK_SOFT,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "lines.linewidth": 2.0,
        "font.size": 10,
        "figure.dpi": 140,
        "savefig.bbox": "tight",
    })


def top_groups(labels: pd.Series, max_groups: int = MAX_SERIES) -> list[str]:
    """The groups that get their own hue; everything else is folded into "Other"."""
    counts = labels.value_counts()
    if len(counts) <= max_groups:
        return list(counts.index)
    return list(counts.index[: max_groups - 1]) + ["Other"]


def fold_groups(labels: pd.Series, keep: list[str]) -> pd.Series:
    return labels.where(labels.isin(keep), "Other")


def _save(fig, out_dir: Path, name: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_scree(explained: pd.Series, out_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    pct = explained * 100
    ax.bar(range(len(pct)), pct.to_numpy(), width=0.6, color=SERIES[0])
    for i, v in enumerate(pct.to_numpy()):
        ax.text(i, v + 0.6, f"{v:.0f}%", ha="center", va="bottom", fontsize=8.5, color=INK_SOFT)
    cumulative = pct.cumsum().iloc[min(2, len(pct) - 1)]
    ax.set_xticks(range(len(pct)), list(pct.index))
    ax.set_ylabel("variance explained (%)")
    ax.set_ylim(0, pct.max() * 1.22)
    ax.set_title(f"First three components carry {cumulative:.0f}% of descriptor variance")
    ax.xaxis.grid(False)
    return _save(fig, out_dir, "01_pca_scree")


def plot_pc_scatter(scores: pd.DataFrame, labels: pd.Series, explained: pd.Series,
                    out_dir: Path, x: str = "PC1", y: str = "PC2",
                    name: str = "02_pca_by_format") -> Path:
    keep = top_groups(labels)
    folded = fold_groups(labels, keep)
    order = [g for g in keep if (folded == g).any()]

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    for i, group in enumerate(order):
        m = folded == group
        ax.scatter(scores.loc[m, x], scores.loc[m, y],
                   s=26 if group != order[0] else 18,
                   c=SERIES[i % len(SERIES)], marker=MARKERS[i % len(MARKERS)],
                   alpha=0.55 if group == order[0] else 0.9,
                   linewidths=0.6, edgecolors=SURFACE,
                   label=f"{group} (n={int(m.sum())})", zorder=2 + i)
    ax.axhline(0, color=GRID, linewidth=0.8, zorder=1)
    ax.axvline(0, color=GRID, linewidth=0.8, zorder=1)
    ax.set_xlabel(f"{x} ({explained[x] * 100:.0f}% of variance)")
    ax.set_ylabel(f"{y} ({explained[y] * 100:.0f}% of variance)")
    ax.set_title("Descriptor space by format family")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), handletextpad=0.4)
    return _save(fig, out_dir, name)


def plot_loadings(loadings: pd.DataFrame, out_dir: Path, component: str = "PC1",
                  n: int = 10) -> Path:
    s = loadings[component]
    s = s.reindex(s.abs().sort_values(ascending=False).index)[:n].sort_values()
    colors = [DIVERGING_POS if v > 0 else DIVERGING_NEG for v in s]

    fig, ax = plt.subplots(figsize=(6.6, 0.38 * len(s) + 1.6))
    ax.barh(range(len(s)), s.to_numpy(), color=colors, height=0.62)
    ax.axvline(0, color=INK_SOFT, linewidth=0.9)
    ax.set_yticks(range(len(s)), [i.replace("_", " ") for i in s.index])
    ax.set_xlabel(f"loading on {component}")
    ax.set_title(f"What {component} is made of")
    ax.yaxis.grid(False)
    handles = [Line2D([], [], marker="s", linestyle="none", color=DIVERGING_POS,
                      label="positive loading"),
               Line2D([], [], marker="s", linestyle="none", color=DIVERGING_NEG,
                      label="negative loading")]
    ax.legend(handles=handles, loc="lower right")
    return _save(fig, out_dir, f"03_loadings_{component.lower()}")


def plot_k_sweep(sweep: pd.DataFrame, out_dir: Path) -> Path:
    """Silhouette and ARI share one axis: both are unitless indices on the same scale."""
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for i, (col, label) in enumerate([("silhouette", "silhouette (cluster separation)"),
                                      ("ari_vs_format", "ARI vs format family")]):
        ax.plot(sweep["k"], sweep[col], color=SERIES[i], marker=MARKERS[i],
                markersize=5, markeredgecolor=SURFACE, markeredgewidth=0.8, label=label)
        last = sweep.iloc[-1]
        ax.annotate(f"{last[col]:.2f}", (last["k"], last[col]),
                    textcoords="offset points", xytext=(6, 0),
                    fontsize=8.5, color=INK_SOFT, va="center")
    ax.axhline(0, color=GRID, linewidth=0.8)
    ax.set_xlabel("number of clusters (k)")
    ax.set_ylabel("index value")
    ax.set_xticks(sweep["k"].tolist())
    ax.set_title("Cluster separation and agreement with format, across k")
    ax.legend(loc="upper right")
    return _save(fig, out_dir, "04_cluster_k_sweep")


def plot_descriptor_by_group(df: pd.DataFrame, column: str, group: str, out_dir: Path,
                             name: str | None = None, min_n: int = 8) -> Path:
    counts = df[group].value_counts()
    candidates = [g for g in counts.index if counts[g] >= min_n]
    pairs = [(g, df.loc[df[group] == g, column].dropna().to_numpy()) for g in candidates]
    pairs = [(g, d) for g, d in pairs if len(d) >= min_n]
    if not pairs:
        raise ValueError(f"no group has at least {min_n} values of {column}")
    groups = [g for g, _ in pairs]
    data = [d for _, d in pairs]
    medians = [float(np.median(d)) for d in data]
    order = np.argsort(medians)
    data = [data[i] for i in order]
    groups = [groups[i] for i in order]

    fig, ax = plt.subplots(figsize=(6.8, 0.44 * len(groups) + 1.8))
    bp = ax.boxplot(data, vert=False, widths=0.56, patch_artist=True,
                    positions=list(range(len(data))), showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor("#cde2fb")
        patch.set_edgecolor(SERIES[0])
        patch.set_linewidth(1.2)
    for element in ("whiskers", "caps"):
        for artist in bp[element]:
            artist.set_color(SERIES[0])
            artist.set_linewidth(1.0)
    for median in bp["medians"]:
        median.set_color(INK)
        median.set_linewidth(1.6)
    for i, d in enumerate(data):
        rng = np.random.default_rng(0)
        ax.scatter(d, i + rng.uniform(-0.14, 0.14, size=len(d)), s=7,
                   color=SERIES[0], alpha=0.35, linewidths=0, zorder=3)
    ax.set_yticks(range(len(groups)), [f"{g}  (n={int(counts[g])})" for g in groups])
    ax.set_xlabel(column.replace("_", " "))
    ax.set_title(f"{column.replace('_', ' ')} by {group.replace('_', ' ')}")
    ax.yaxis.grid(False)
    return _save(fig, out_dir, name or f"05_{column}_by_{group}")


def plot_composition(table: pd.DataFrame, out_dir: Path,
                     name: str = "06_cluster_composition") -> Path:
    """Row-normalised cluster x format table on a single-hue ramp."""
    values = table.to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(0.9 * table.shape[1] + 3.2, 0.5 * table.shape[0] + 2.2))
    im = ax.imshow(values, cmap=SEQUENTIAL, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(table.shape[1]),
                  [str(c).replace(" (", "\n(") for c in table.columns],
                  rotation=35, ha="right")
    ax.set_yticks(range(table.shape[0]), [f"cluster {i}" for i in table.index])
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            v = values[i, j]
            if v >= 0.05:
                ax.text(j, i, f"{v:.0%}", ha="center", va="center", fontsize=8,
                        color="#f0efec" if v > 0.55 else INK)
    ax.grid(False)
    ax.set_title("Share of each cluster by format family")
    fig.colorbar(im, ax=ax, shrink=0.8, label="share of cluster")
    return _save(fig, out_dir, name)


def plot_arm_pi(df: pd.DataFrame, out_dir: Path, band: float,
                name: str = "07_bispecific_arm_pI") -> Path:
    """Fv pI of arm A against arm B; points near the diagonal are hard to separate."""
    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    lo = float(min(df["pI_arm_a"].min(), df["pI_arm_b"].min())) - 0.3
    hi = float(max(df["pI_arm_a"].max(), df["pI_arm_b"].max())) + 0.3
    grid = np.array([lo, hi])
    ax.fill_between(grid, grid - 2 * band, grid + 2 * band, color=GRID, alpha=0.7,
                    linewidth=0, label=f"arm ΔpI < {2 * band:g}", zorder=1)
    ax.plot(grid, grid, color=INK_SOFT, linewidth=0.9, zorder=2)
    for i, family in enumerate(sorted(df["format_family"].unique())):
        m = df["format_family"] == family
        ax.scatter(df.loc[m, "pI_arm_a"], df.loc[m, "pI_arm_b"], s=24,
                   c=SERIES[i % len(SERIES)], marker=MARKERS[i % len(MARKERS)],
                   alpha=0.85, linewidths=0.6, edgecolors=SURFACE,
                   label=f"{family} (n={int(m.sum())})", zorder=3 + i)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel("Fv pI, arm A")
    ax.set_ylabel("Fv pI, arm B")
    ax.set_title("How far apart are a bispecific's two arms?")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    return _save(fig, out_dir, name)


def plot_charge_gap(df: pd.DataFrame, out_dir: Path, cutoff: float,
                    name: str = "08_bispecific_homodimer_charge_gap") -> Path:
    """Distribution of the Fv charge difference between product and homodimers at pH 5.5."""
    values = df["delta_charge_homodimer_pH5p5"].to_numpy()
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    bins = np.arange(0, max(values.max(), cutoff) + 1.0, 0.5)
    ax.hist(values, bins=bins, color=SERIES[0], edgecolor=SURFACE, linewidth=0.8)
    ax.axvline(cutoff, color=DIVERGING_POS, linewidth=1.6)
    share = float((values < cutoff).mean())
    ax.text(0.98, 0.92, f"{share:.0%} of bispecifics below {cutoff:g} e (orange line)",
            transform=ax.transAxes, ha="right", color=INK, fontsize=9)
    ax.set_xlabel("|charge(product) − charge(homodimer)| at pH 5.5, Fv level (e)")
    ax.set_ylabel("bispecifics")
    ax.set_title("Charge handle for removing heavy-chain homodimers")
    ax.xaxis.grid(False)
    return _save(fig, out_dir, name)
