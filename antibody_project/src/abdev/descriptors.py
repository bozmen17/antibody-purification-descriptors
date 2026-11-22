"""Sequence-only descriptors for an Fv, chosen for relevance to downstream purification.

Why these descriptors: the chromatography steps that dominate an antibody platform
process separate on charge (pI and charge at the pH of the step), on hydrophobicity
(HIC, and the aggregation risk that follows from exposed apolar surface), and on the
framework family the capture resin recognises. Structure-based surface patches are
approximated here by sliding-window sums, which is what you can do from sequence alone.
"""
from __future__ import annotations

import re

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .charge import charge_patch, extreme_window, isoelectric_point, net_charge
from .regions import heavy_loops, light_loops
from .scales import EISENBERG, KYTE_DOOLITTLE, NEGATIVE, POSITIVE, STANDARD_AA

#: pH values of interest: a typical elution pH, a wash/load pH, and formulation pH.
PH_POINTS = (3.6, 5.5, 7.4)

#: Sequence motifs that tend to cause heterogeneity during processing and storage.
MOTIFS = {
    "glycosylation_sequon": r"N[^P][ST]",
    "deamidation_site": r"N[GSNTH]",
    "isomerisation_site": r"D[GSTDH]",
    "cysteine": r"C",
}


def clean_sequence(value) -> str | None:
    """Normalise one sequence field, returning None when it is missing or non-standard."""
    if not isinstance(value, str):
        return None
    seq = re.sub(r"[^A-Za-z]", "", value).upper()
    if len(seq) < 60 or seq in {"NA", "NONE"} or set(seq) - STANDARD_AA:
        return None
    return seq


def hydrophobicity_index(seqs: list[str]) -> float:
    """Ratio of hydrophobic to hydrophilic Eisenberg contributions over the given regions."""
    joined = "".join(seqs)
    if not joined:
        return float("nan")
    apolar = sum(v for v in (EISENBERG[a] for a in joined) if v > 0)
    polar = sum(v for v in (EISENBERG[a] for a in joined) if v < 0)
    return float(-apolar / polar) if polar else float("nan")


def count_motif(seqs: list[str], pattern: str) -> int:
    """Count a motif across regions, allowing overlaps."""
    return sum(len(re.findall(f"(?=({pattern}))", s)) for s in seqs)


def aromatic_fraction(seq: str) -> float:
    return sum(seq.count(a) for a in "FWY") / len(seq)


def vh3_framework(vh: str) -> bool:
    """Flag a VH3-type framework, whose affinity-capture behaviour differs from other families."""
    return bool(re.search(r"[LVQ]ESGGG", vh[:20]))


def fv_descriptors(vh: str, vl: str | None = None) -> dict:
    """Compute every descriptor for one Fv. `vl` is None for single-domain formats."""
    chains = [vh] + ([vl] if vl else [])
    fv = "".join(chains)
    protein = ProteinAnalysis(fv)

    hl = heavy_loops(vh)
    ll = light_loops(vl) if vl else None
    loops = hl.as_list() + (ll.as_list() if ll else [])
    loops_ok = hl.ok and (ll.ok if ll else True)

    d: dict[str, object] = {
        "paired": vl is not None,
        "fv_length": len(fv),
        "fv_mass_kda": round(protein.molecular_weight() / 1000.0, 3),
        "loops_resolved": loops_ok,
        "loop_note": "; ".join(n for n in (hl.note, ll.note if ll else "") if n),
    }

    # --- charge
    d["pI_fv"] = isoelectric_point(*chains)
    d["pI_vh"] = isoelectric_point(vh)
    d["pI_vl"] = isoelectric_point(vl) if vl else np.nan
    d["pI_chain_gap"] = abs(d["pI_vh"] - d["pI_vl"]) if vl else np.nan
    for ph in PH_POINTS:
        key = f"{ph:g}".replace(".", "p")
        d[f"charge_pH{key}"] = round(sum(net_charge(c, ph) for c in chains), 3)
        d[f"charge_per_kda_pH{key}"] = round(d[f"charge_pH{key}"] / d["fv_mass_kda"], 3)
    d["pos_patch_pH5p5"] = round(max(charge_patch(c, 5.5, mode="max") for c in chains), 3)
    d["neg_patch_pH5p5"] = round(min(charge_patch(c, 5.5, mode="min") for c in chains), 3)
    d["frac_positive"] = round(sum(fv.count(a) for a in POSITIVE) / len(fv), 4)
    d["frac_negative"] = round(sum(fv.count(a) for a in NEGATIVE) / len(fv), 4)
    d["charge_asymmetry"] = round(d["frac_positive"] - d["frac_negative"], 4)

    # --- hydrophobicity
    d["gravy_fv"] = round(protein.gravy(), 4)
    d["max_hydrophobic_patch"] = round(
        max(extreme_window(c, KYTE_DOOLITTLE, 7, "max") for c in chains), 3)
    d["aromatic_fraction"] = round(aromatic_fraction(fv), 4)
    d["hydrophobicity_index_fv"] = round(hydrophobicity_index(chains), 4)
    d["hydrophobicity_index_loops"] = (
        round(hydrophobicity_index(loops), 4) if loops_ok else np.nan)

    # --- loops
    d["loop_total_length"] = sum(len(c) for c in loops) if loops_ok else np.nan
    d["cdr3_h_length"] = len(hl.cdr3) if hl.ok else np.nan
    d["cdr3_h_charge_pH7p4"] = round(net_charge(hl.cdr3, 7.4, False), 3) if hl.ok else np.nan
    d["loop_frac_positive"] = (
        round(sum("".join(loops).count(a) for a in POSITIVE) / max(1, len("".join(loops))), 4)
        if loops_ok else np.nan)
    d["loop_aromatic_fraction"] = (
        round(aromatic_fraction("".join(loops)), 4) if loops_ok and loops else np.nan)

    # --- secondary structure propensity and stability proxies
    helix, turn, sheet = protein.secondary_structure_fraction()
    d["frac_sheet"] = round(sheet, 4)
    d["frac_turn"] = round(turn, 4)
    d["instability_index"] = round(protein.instability_index(), 3)
    d["extinction_280"] = protein.molar_extinction_coefficient()[1]

    # --- processing liabilities
    for name, pattern in MOTIFS.items():
        d[f"n_{name}_fv"] = count_motif(chains, pattern)
        d[f"n_{name}_loops"] = count_motif(loops, pattern) if loops_ok else np.nan
    d["n_oxidation_sites_loops"] = (
        sum("".join(loops).count(a) for a in "MW") if loops_ok else np.nan)
    d["vh3_framework"] = vh3_framework(vh)

    return d
