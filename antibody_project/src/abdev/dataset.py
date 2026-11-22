"""Load Thera-SAbDab and reduce it to one row per variable-domain arm.

Thera-SAbDab gives one row per therapeutic, with up to two Fv arms (a second pair of
sequence columns is filled for bispecifics). Descriptors are a property of an arm, so
the table is exploded: one row per arm, carrying the parent molecule's format and
metadata. The free-text `Format` field (a hundred-odd distinct strings) is collapsed
into a small number of families, which is what the clustering is asked to recover.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .descriptors import clean_sequence, fv_descriptors

THERA_SABDAB_URL = (
    "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/"
    "TheraSAbDab_SeqStruc_OnlineDownload.csv"
)

ARMS = (
    ("arm1", "HeavySequence", "LightSequence"),
    ("arm2", "HeavySequence(ifbispec)", "LightSequence(ifbispec)"),
)


def classify_format(raw: str) -> str:
    """Collapse the free-text Format field into one of a few engineering families.

    Order matters: the most specific test wins, so a bispecific built from scFvs is
    counted as a bispecific fragment rather than as an scFv.
    """
    if not isinstance(raw, str) or not raw.strip():
        return "Unclassified"
    f = raw.lower()
    multispecific = any(k in f for k in ("bispec", "trispec", "tetraspec", "pentavalent"))
    fragment_only = not any(k in f for k in ("mab", "igg", "ch2", "ch3", "fc", "igm"))

    if multispecific:
        return "Multispecific (fragment-based)" if fragment_only else "Multispecific (IgG-like)"
    if "single domain" in f or "nanobody" in f or "vhh" in f or re.search(r"\bvh-", f):
        return "Single domain"
    if "scfv" in f:
        return "scFv / scFv fusion"
    if "fab" in f:
        return "Fab / Fab fusion"
    if "adc" in f or "radiolabel" in f or "conjugat" in f or "toxin" in f:
        return "IgG conjugate"
    if "fusion" in f:
        return "IgG fusion"
    if "mab" in f or "igg" in f or "igm" in f:
        return "Whole mAb"
    return "Other"


def _find_column(df: pd.DataFrame, prefix: str) -> str | None:
    """Thera-SAbDab stamps a date into some column names; match on the stable prefix."""
    for col in df.columns:
        if col.lower().startswith(prefix.lower()):
            return col
    return None


def load_raw(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = {c for _, h, l in ARMS for c in (h, l)} - set(df.columns)
    if missing:
        raise ValueError(f"Thera-SAbDab download is missing expected columns: {sorted(missing)}")
    return df


def explode_arms(df: pd.DataFrame) -> pd.DataFrame:
    """One row per Fv arm, with usable sequences only."""
    trial_col = _find_column(df, "Highest_Clin_Trial")
    genetics_col = _find_column(df, "Genetics")
    records = []

    for _, row in df.iterrows():
        for arm, h_col, l_col in ARMS:
            vh = clean_sequence(row[h_col])
            if vh is None:
                continue  # no heavy variable domain means nothing to describe
            vl = clean_sequence(row[l_col])
            records.append({
                "therapeutic": row.get("Therapeutic"),
                "arm": arm,
                "format_raw": row.get("Format"),
                "format_family": classify_format(row.get("Format")),
                "isotype": str(row.get("CH1 Isotype", "")).split(";")[0].strip(),
                "light_type": str(row.get("VD LC", "")).split(";")[0].strip(),
                "origin": (str(row.get(genetics_col, "")).split(";")[0].strip()
                           if genetics_col else ""),
                "highest_trial": row.get(trial_col) if trial_col else None,
                "status": row.get("Est. Status"),
                "target": row.get("Target"),
                "year_proposed": pd.to_numeric(row.get("Year Proposed"), errors="coerce"),
                "vh": vh,
                "vl": vl,
            })

    out = pd.DataFrame.from_records(records)
    # A handful of entries appear twice (one row per specificity of the same molecule).
    out = out.drop_duplicates(subset=["vh", "vl"], keep="first").reset_index(drop=True)
    return out


def add_descriptors(arms: pd.DataFrame) -> pd.DataFrame:
    """Append every descriptor column to the exploded arm table."""
    # pandas stores a missing light chain as NaN, which is truthy; map it back to None
    # so single-domain arms are described as unpaired rather than crashing.
    rows = [fv_descriptors(r.vh, r.vl if isinstance(r.vl, str) else None)
            for r in arms.itertuples()]
    desc = pd.DataFrame(rows, index=arms.index)
    overlap = set(desc.columns) & set(arms.columns)
    if overlap:
        raise ValueError(f"descriptor names collide with metadata columns: {sorted(overlap)}")
    return pd.concat([arms, desc], axis=1)


def build(path: str | Path) -> pd.DataFrame:
    return add_descriptors(explode_arms(load_raw(path)))


DESCRIPTOR_START = "paired"  # first descriptor column produced by fv_descriptors

#: Flags and bookkeeping columns that are descriptors but not features.
NOT_FEATURES = {"year_proposed", "loops_resolved", "paired", "vh3_framework", "extinction_280"}


def feature_columns(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    """Split numeric descriptors into two feature sets.

    Returns `(core, loop_dependent)`. Core descriptors are defined for every arm,
    including single-domain ones, so the whole table can be projected at once.
    Loop-dependent descriptors (and the VL-specific ones) are only defined for a
    subset, so they get their own analysis on that subset rather than being imputed.
    """
    start = df.columns.get_loc(DESCRIPTOR_START)
    numeric = [c for c in df.iloc[:, start:].select_dtypes("number").columns
               if c not in NOT_FEATURES]
    core = [c for c in numeric if df[c].notna().all()]
    loop_dependent = [c for c in numeric if c not in core]
    return core, loop_dependent
