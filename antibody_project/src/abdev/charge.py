"""Charge, isoelectric point and sliding-window patch descriptors."""
from __future__ import annotations

import numpy as np

from .scales import (
    PKA_ACIDIC,
    PKA_BASIC,
    PKA_C_TERMINUS,
    PKA_N_TERMINUS,
)


def residue_charge(aa: str, ph: float) -> float:
    """Fractional charge of one side chain at a given pH."""
    if aa in PKA_BASIC:
        return 1.0 / (1.0 + 10 ** (ph - PKA_BASIC[aa]))
    if aa in PKA_ACIDIC:
        return -1.0 / (1.0 + 10 ** (PKA_ACIDIC[aa] - ph))
    return 0.0


def net_charge(chain: str, ph: float, include_termini: bool = True) -> float:
    """Net charge of a polypeptide chain at a given pH (Henderson-Hasselbalch)."""
    q = sum(residue_charge(aa, ph) for aa in chain)
    if include_termini:
        q += 1.0 / (1.0 + 10 ** (ph - PKA_N_TERMINUS))
        q -= 1.0 / (1.0 + 10 ** (PKA_C_TERMINUS - ph))
    return q


def isoelectric_point(*chains: str, tol: float = 1e-4) -> float:
    """pI of one or more chains treated as a single molecule, by bisection.

    Passing chains separately matters: each one contributes its own pair of termini.
    """
    low, high = 0.0, 14.0
    while high - low > tol:
        mid = 0.5 * (low + high)
        if sum(net_charge(c, mid) for c in chains) > 0:
            low = mid
        else:
            high = mid
    return round(0.5 * (low + high), 3)


def charge_profile(chain: str, ph_values) -> list[float]:
    """Net charge across a pH range, e.g. to compare binding windows."""
    return [net_charge(chain, ph) for ph in ph_values]


def extreme_window(chain: str, scale: dict[str, float], width: int = 7,
                   mode: str = "max") -> float:
    """Most extreme sliding-window sum of a per-residue scale.

    A crude sequence-only stand-in for a surface patch: `mode="max"` with a
    hydrophobicity scale approximates the largest hydrophobic patch, and
    `mode="min"`/`"max"` on per-residue charge approximates negative/positive patches.
    """
    if len(chain) < width:
        return float("nan")
    values = np.array([scale.get(aa, 0.0) for aa in chain], dtype=float)
    sums = np.convolve(values, np.ones(width), mode="valid")
    return float(sums.max() if mode == "max" else sums.min())


def charge_patch(chain: str, ph: float, width: int = 7, mode: str = "max") -> float:
    """Most extreme windowed charge sum along a chain at a given pH."""
    if len(chain) < width:
        return float("nan")
    values = np.array([residue_charge(aa, ph) for aa in chain], dtype=float)
    sums = np.convolve(values, np.ones(width), mode="valid")
    return float(sums.max() if mode == "max" else sums.min())
