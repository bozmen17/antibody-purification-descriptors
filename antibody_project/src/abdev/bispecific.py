"""Charge separation between a bispecific and its mispaired by-products.

A two-arm bispecific (arms A and B) is expressed alongside the species its chains can
mispair into. Two kinds matter for purification:

* **Heavy-chain homodimers** AA and BB, next to the wanted heterodimer AB. With identical
  constant regions, the Fv-level charge of AA is 2*q_A and of AB is q_A + q_B, so the
  charge difference between the product and either homodimer is exactly |q_A - q_B|. That
  is the handle an ion-exchange step has for removing them.
* **Light-chain mispairs**, where an arm carries the other arm's light chain. Their charge
  differs from the product by |q(VL_A) - q(VL_B)|. With a common light chain this is zero
  by construction, which is why common-LC formats only have to solve the homodimer problem.

Everything is computed on the variable domains Thera-SAbDab provides. Constant-region
engineering (heterodimerisation or pI-shifting mutations) is invisible here and can add to
or replace the Fv-level difference; see the README.
"""
from __future__ import annotations

import pandas as pd

from .charge import isoelectric_point, net_charge
from .dataset import classify_format
from .descriptors import clean_sequence

#: Load/wash pH for cation exchange and formulation pH, matching descriptors.PH_POINTS.
PH_CEX = 5.5
PH_NEUTRAL = 7.4

#: Illustrative cut-offs for "hard to resolve by charge alone" on the Fv contribution.
#: They are not process limits; they make the distribution easy to read.
CHARGE_GAP_CUTOFF = 1.0     # elementary charges at pH 5.5
PI_GAP_CUTOFF = 0.2         # pH units between product and nearest homodimer


def _chains(vh: str, vl: str | None) -> list[str]:
    return [vh] + ([vl] if vl else [])


def _charge(chains: list[str], ph: float) -> float:
    return sum(net_charge(c, ph) for c in chains)


def bispecific_pairs(raw: pd.DataFrame) -> pd.DataFrame:
    """One row per therapeutic with two sequenced arms."""
    rows = []
    for _, r in raw.iterrows():
        vh1 = clean_sequence(r["HeavySequence"])
        vh2 = clean_sequence(r["HeavySequence(ifbispec)"])
        if vh1 is None or vh2 is None:
            continue
        vl1 = clean_sequence(r["LightSequence"])
        vl2 = clean_sequence(r["LightSequence(ifbispec)"])
        rows.append({
            "therapeutic": r.get("Therapeutic"),
            "format_raw": r.get("Format"),
            "format_family": classify_format(r.get("Format")),
            "vh1": vh1, "vl1": vl1, "vh2": vh2, "vl2": vl2,
        })
    return pd.DataFrame.from_records(rows)


def pair_descriptors(vh1: str, vl1: str | None, vh2: str, vl2: str | None) -> dict:
    """Charge separation between product AB and its mispaired species, Fv level."""
    a, b = _chains(vh1, vl1), _chains(vh2, vl2)
    d: dict[str, object] = {
        "common_light_chain": vl1 is not None and vl1 == vl2,
        "both_light_chains": vl1 is not None and vl2 is not None,
        "pI_arm_a": isoelectric_point(*a),
        "pI_arm_b": isoelectric_point(*b),
    }
    d["delta_pI_arms"] = round(abs(d["pI_arm_a"] - d["pI_arm_b"]), 3)

    # Heavy-chain homodimers vs product: Fv charge difference is exactly |q_A - q_B|.
    for ph, key in ((PH_CEX, "pH5p5"), (PH_NEUTRAL, "pH7p4")):
        d[f"delta_charge_homodimer_{key}"] = round(abs(_charge(a, ph) - _charge(b, ph)), 3)

    # pI of each assembled species (two Fvs each); chains passed separately for termini.
    pi_aa = isoelectric_point(*a, *a)
    pi_ab = isoelectric_point(*a, *b)
    pi_bb = isoelectric_point(*b, *b)
    d.update(pI_homodimer_aa=pi_aa, pI_heterodimer=pi_ab, pI_homodimer_bb=pi_bb)
    d["pI_gap_to_nearest_homodimer"] = round(min(abs(pi_ab - pi_aa), abs(pi_ab - pi_bb)), 3)

    # Light-chain mispair vs product (only meaningful with two distinct light chains).
    if d["both_light_chains"] and not d["common_light_chain"]:
        d["delta_charge_lc_mispair_pH5p5"] = round(
            abs(net_charge(vl1, PH_CEX) - net_charge(vl2, PH_CEX)), 3)
    else:
        d["delta_charge_lc_mispair_pH5p5"] = 0.0 if d["common_light_chain"] else float("nan")

    d["hard_by_charge_homodimer"] = d["delta_charge_homodimer_pH5p5"] < CHARGE_GAP_CUTOFF
    return d


def build_pairs(raw: pd.DataFrame) -> pd.DataFrame:
    pairs = bispecific_pairs(raw)
    def seq(v):  # pandas turns a missing light chain into NaN
        return v if isinstance(v, str) else None

    desc = pd.DataFrame([pair_descriptors(r.vh1, seq(r.vl1), r.vh2, seq(r.vl2))
                         for r in pairs.itertuples()], index=pairs.index)
    return pd.concat([pairs, desc], axis=1)
