"""Step 2: PCA, clustering and figures over the descriptor table.

    python scripts/02_analyse.py
    python scripts/02_analyse.py --k 6 --dims 5

Writes figures to `figures/`, result tables to `results/`, and prints the numbers
quoted in the README so they can never drift from the code.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abdev import analysis, plots  # noqa: E402
from abdev.dataset import feature_columns  # noqa: E402

IN = ROOT / "data" / "processed" / "fv_descriptors.csv"
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"

HEADLINE_DESCRIPTORS = [
    "pI_fv", "charge_pH5p5", "gravy_fv", "max_hydrophobic_patch",
    "hydrophobicity_index_fv", "aromatic_fraction", "instability_index",
    "n_glycosylation_sequon_fv",
]


def section(title: str) -> None:
    print(f"\n{title}\n{'-' * len(title)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5, help="clusters for the reported solution")
    ap.add_argument("--dims", type=int, default=5, help="PCA dimensions used for clustering")
    args = ap.parse_args()

    if not IN.exists():
        raise SystemExit("run scripts/01_build_dataset.py first")

    plots.use_house_style()
    RESULTS.mkdir(exist_ok=True)
    df = pd.read_csv(IN)

    core, loop_dependent = feature_columns(df)
    features, dropped = analysis.prune_collinear(df, core)

    section("descriptor set")
    print(f"{len(core)} core descriptors, {len(features)} kept after pruning |r| > 0.95")
    for col, against, r in dropped:
        print(f"  dropped {col:<28} (r = {r} with {against})")
    print(f"{len(loop_dependent)} loop- or light-chain-dependent descriptors held back; "
          f"{int(df['loops_resolved'].sum())} arms have resolved loops")

    pca = analysis.run_pca(df, features, n_components=6)
    section("PCA")
    for name, frac in pca.explained.items():
        print(f"  {name}: {frac:.1%}")
    print(f"  first three components: {pca.explained[:3].sum():.1%}")
    for component in ("PC1", "PC2", "PC3"):
        print(f"\n  top loadings on {component}:")
        for descriptor, value in pca.top_loadings(component, 6).items():
            print(f"    {value:+.3f}  {descriptor}")

    sweep = analysis.sweep_k(pca.scores, df["format_family"], dims=args.dims)
    section("k-means sweep")
    print(sweep.to_string(index=False))
    best = sweep.loc[sweep["silhouette"].idxmax()]
    print(f"\n  best silhouette at k = {int(best['k'])} ({best['silhouette']:.3f}), "
          f"where ARI vs format = {best['ari_vs_format']:.3f}")

    df["cluster_kmeans"] = analysis.cluster(pca.scores, args.k, args.dims, "kmeans")
    df["cluster_ward"] = analysis.cluster(pca.scores, args.k, args.dims, "ward")
    section(f"agreement with format labels at k = {args.k}")
    for method in ("kmeans", "ward"):
        scores = analysis.agreement(df["format_family"], df[f"cluster_{method}"])
        print(f"  {method:<8} ARI = {scores['ari']:+.3f}   AMI = {scores['ami']:+.3f}")
    baseline = analysis.permutation_baseline(df["format_family"], df["cluster_kmeans"])
    print(f"  label-shuffled baseline: ARI = {baseline['mean']:+.3f} "
          f"(sd {baseline['sd']:.3f})")
    print("\n  same test against isotype and against molecular origin:")
    placeholders = {"", "na", "nan", "none", "tbc"}
    for column in ("isotype", "light_type", "origin"):
        known = ~df[column].astype(str).str.strip().str.lower().isin(placeholders)
        subset = df[known]
        scores = analysis.agreement(subset[column], subset["cluster_kmeans"])
        print(f"    {column:<11} ARI = {scores['ari']:+.3f}   AMI = {scores['ami']:+.3f}")

    composition = analysis.cluster_composition(df, "cluster_kmeans", "format_family")
    summary = analysis.group_summary(df, "format_family", HEADLINE_DESCRIPTORS)
    cluster_profile = analysis.group_summary(df, "cluster_kmeans", HEADLINE_DESCRIPTORS)

    section("median descriptors by format family")
    print(summary.to_string())
    section("median descriptors by cluster")
    print(cluster_profile.to_string())

    sweep.to_csv(RESULTS / "kmeans_sweep.csv", index=False)
    pca.loadings.to_csv(RESULTS / "pca_loadings.csv")
    pca.explained.rename("fraction").to_csv(RESULTS / "pca_explained_variance.csv")
    composition.to_csv(RESULTS / "cluster_composition.csv")
    summary.to_csv(RESULTS / "median_descriptors_by_format.csv")
    cluster_profile.to_csv(RESULTS / "median_descriptors_by_cluster.csv")
    df.drop(columns=["vh", "vl"]).to_csv(RESULTS / "arms_with_clusters.csv", index=False)

    written = [
        plots.plot_scree(pca.explained, FIGURES),
        plots.plot_pc_scatter(pca.scores, df["format_family"], pca.explained, FIGURES),
        plots.plot_loadings(pca.loadings, FIGURES, "PC1"),
        plots.plot_loadings(pca.loadings, FIGURES, "PC2"),
        plots.plot_k_sweep(sweep, FIGURES),
        plots.plot_descriptor_by_group(df, "pI_fv", "format_family", FIGURES,
                                       name="05_pI_by_format"),
        plots.plot_descriptor_by_group(df, "max_hydrophobic_patch", "format_family", FIGURES,
                                       name="05_hydrophobic_patch_by_format"),
        plots.plot_composition(composition, FIGURES),
    ]

    section("written")
    for path in written:
        print(f"  {path.relative_to(ROOT)}")
    print(f"  {len(list(RESULTS.glob('*.csv')))} tables in results/")


if __name__ == "__main__":
    main()
