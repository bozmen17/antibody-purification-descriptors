"""PCA and clustering over the descriptor table, with the checks that make it honest.

Two questions are kept apart:

1. *What varies?* PCA on standardised descriptors, reporting explained variance and
   the loadings that define each component. Highly collinear descriptors are pruned
   first, otherwise a block of near-duplicates dictates PC1 by weight of numbers.
2. *Does the variation follow format?* k-means and hierarchical clusters are compared
   with the engineering format family by adjusted Rand index and by adjusted mutual
   information. A silhouette score says whether the clusters are real at all. Both
   are reported even when they are unflattering: the interesting result here is that
   format is only weakly recoverable from Fv sequence descriptors, because the
   properties the descriptors capture live mostly in the variable domain, while the
   format is largely a constant-region decision.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    silhouette_score,
)
from sklearn.preprocessing import StandardScaler

RANDOM_STATE = 0


def prune_collinear(df: pd.DataFrame, columns: list[str], threshold: float = 0.95
                    ) -> tuple[list[str], list[tuple[str, str, float]]]:
    """Drop one of each pair of descriptors correlating above `threshold`.

    The survivor is the earlier column, so the kept name is the more primary descriptor
    (e.g. charge at a pH is kept over charge per kDa at the same pH).
    """
    corr = df[columns].corr().abs()
    kept: list[str] = []
    dropped: list[tuple[str, str, float]] = []
    for col in columns:
        hit = next(((k, corr.loc[col, k]) for k in kept if corr.loc[col, k] > threshold), None)
        if hit is None:
            kept.append(col)
        else:
            dropped.append((col, hit[0], round(float(hit[1]), 3)))
    return kept, dropped


@dataclass
class PCAResult:
    scores: pd.DataFrame          # one row per arm, columns PC1..PCn
    loadings: pd.DataFrame        # descriptors x components
    explained: pd.Series          # fraction of variance per component
    features: list[str]
    model: PCA
    scaler: StandardScaler

    def top_loadings(self, component: str = "PC1", n: int = 8) -> pd.Series:
        return self.loadings[component].reindex(
            self.loadings[component].abs().sort_values(ascending=False).index)[:n]


def run_pca(df: pd.DataFrame, features: list[str], n_components: int = 6) -> PCAResult:
    scaler = StandardScaler()
    x = scaler.fit_transform(df[features].to_numpy(dtype=float))
    n = min(n_components, x.shape[1], x.shape[0])
    model = PCA(n_components=n, random_state=RANDOM_STATE)
    scores = model.fit_transform(x)
    names = [f"PC{i + 1}" for i in range(n)]
    return PCAResult(
        scores=pd.DataFrame(scores, columns=names, index=df.index),
        loadings=pd.DataFrame(model.components_.T, index=features, columns=names),
        explained=pd.Series(model.explained_variance_ratio_, index=names),
        features=features,
        model=model,
        scaler=scaler,
    )


def sweep_k(scores: pd.DataFrame, labels: pd.Series, dims: int = 5,
            k_values=range(2, 11)) -> pd.DataFrame:
    """Silhouette and agreement-with-format across a range of k, for k-means.

    Clustering runs in the PCA space (the leading `dims` components), which denoises
    the descriptor set without changing what is being measured.
    """
    x = scores.iloc[:, :dims].to_numpy()
    rows = []
    for k in k_values:
        km = KMeans(n_clusters=k, n_init=10, random_state=RANDOM_STATE).fit(x)
        rows.append({
            "k": k,
            "silhouette": round(float(silhouette_score(x, km.labels_)), 4),
            "inertia": round(float(km.inertia_), 1),
            "ari_vs_format": round(float(adjusted_rand_score(labels, km.labels_)), 4),
            "ami_vs_format": round(float(adjusted_mutual_info_score(labels, km.labels_)), 4),
        })
    return pd.DataFrame(rows)


def cluster(scores: pd.DataFrame, k: int, dims: int = 5, method: str = "kmeans") -> np.ndarray:
    x = scores.iloc[:, :dims].to_numpy()
    if method == "kmeans":
        return KMeans(n_clusters=k, n_init=10, random_state=RANDOM_STATE).fit_predict(x)
    if method == "ward":
        return AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(x)
    raise ValueError(f"unknown clustering method: {method}")


def agreement(labels_true, labels_pred) -> dict[str, float]:
    return {
        "ari": round(float(adjusted_rand_score(labels_true, labels_pred)), 4),
        "ami": round(float(adjusted_mutual_info_score(labels_true, labels_pred)), 4),
    }


def permutation_baseline(labels_true, labels_pred, n: int = 500) -> dict[str, float]:
    """What ARI would look like if the clusters carried no information about format.

    Adjusted indices are already chance-corrected, so this is a sanity check on the
    correction rather than the headline number: it should sit at roughly zero.
    """
    rng = np.random.default_rng(RANDOM_STATE)
    pred = np.asarray(labels_pred)
    scores = [adjusted_rand_score(labels_true, rng.permutation(pred)) for _ in range(n)]
    return {"mean": round(float(np.mean(scores)), 4), "sd": round(float(np.std(scores)), 4)}


def group_summary(df: pd.DataFrame, by: str, columns: list[str]) -> pd.DataFrame:
    """Median descriptor value per group, with group sizes, for the write-up table."""
    out = df.groupby(by)[columns].median(numeric_only=True)
    out.insert(0, "n", df.groupby(by).size())
    return out.sort_values("n", ascending=False)


def cluster_composition(df: pd.DataFrame, cluster_col: str, label_col: str) -> pd.DataFrame:
    """Row-normalised contingency table: what each cluster is made of."""
    table = pd.crosstab(df[cluster_col], df[label_col])
    return table.div(table.sum(axis=1), axis=0).round(3)
