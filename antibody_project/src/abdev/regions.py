"""Locate the three CDR loops in a variable domain from conserved sequence anchors.

This uses the well-known anchor rules (conserved Cys, Trp, and the WGxG / FGxG
framework-4 motifs) rather than full numbering. It is quick and dependency-free;
every call reports whether the loops it found have plausible lengths, so entries
that fail can be dropped or inspected instead of silently producing nonsense.

ANARCI-based IMGT numbering is the obvious upgrade; see the README.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

# Framework-4 motifs that close the third loop.
FR4_HEAVY = re.compile(r"W[GA][QRKCSG][GAQ]")
FR4_LIGHT = re.compile(r"F[GA][GAQSCV][GA]")
# Motif that closes the second heavy loop (framework 3 entry).
FR3_HEAVY = re.compile(r"[KRQ][LIVFTAM][TSIAN]")

LENGTH_LIMITS = {
    "H": {"cdr1": (8, 14), "cdr2": (14, 21), "cdr3": (1, 40)},
    "L": {"cdr1": (8, 18), "cdr2": (5, 9), "cdr3": (4, 16)},
}


@dataclass
class Loops:
    """The three loops of one variable domain, plus a flag saying whether to trust them."""

    chain_type: str
    cdr1: str = ""
    cdr2: str = ""
    cdr3: str = ""
    ok: bool = False
    note: str = ""

    def as_list(self) -> list[str]:
        return [c for c in (self.cdr1, self.cdr2, self.cdr3) if c]

    @property
    def concatenated(self) -> str:
        return "".join(self.as_list())


def _plausible(loops: Loops) -> bool:
    limits = LENGTH_LIMITS[loops.chain_type]
    for name, (lo, hi) in limits.items():
        if not lo <= len(getattr(loops, name)) <= hi:
            loops.note = f"{name} length {len(getattr(loops, name))} outside {lo}-{hi}"
            return False
    return True


def heavy_loops(seq: str) -> Loops:
    """CDR-H1 (AbM-like), CDR-H2 and CDR-H3 of a VH domain."""
    out = Loops("H")
    try:
        cys1 = seq.index("C", 17, 30)
        trp = seq.index("W", cys1 + 11, cys1 + 21)
        out.cdr1 = seq[cys1 + 4: trp]

        # CDR-H2 is taken as a fixed 16-residue window from the Trp anchor (Kabat 50-65).
        # The framework-3 motif is only a cross-check: it is absent in a number of
        # murine frameworks, so making the boundary depend on it loses those entries.
        start2 = trp + 14
        out.cdr2 = seq[start2: start2 + 16]
        fr3 = FR3_HEAVY.search(seq, start2 + 15, start2 + 22)
        if fr3 is not None and fr3.start() != start2 + 17:
            out.note = "framework-3 anchor offset from the expected position"

        fr4 = FR4_HEAVY.search(seq, 90)
        if fr4 is None:
            raise ValueError("no WGxG framework-4 motif")
        cys2 = seq.rindex("C", 80, fr4.start())
        out.cdr3 = seq[cys2 + 3: fr4.start()]
    except ValueError as exc:
        out.note = str(exc)
        return out
    out.ok = _plausible(out)
    return out


def light_loops(seq: str) -> Loops:
    """CDR-L1, CDR-L2 and CDR-L3 of a VL (kappa or lambda) domain."""
    out = Loops("L")
    try:
        cys1 = seq.index("C", 17, 30)
        trp = seq.index("W", cys1 + 9, cys1 + 21)
        out.cdr1 = seq[cys1 + 1: trp]
        out.cdr2 = seq[trp + 15: trp + 22]

        fr4 = FR4_LIGHT.search(seq, 88)
        if fr4 is None:
            raise ValueError("no FGxG framework-4 motif")
        cys2 = seq.rindex("C", 75, fr4.start())
        out.cdr3 = seq[cys2 + 1: fr4.start()]
    except ValueError as exc:
        out.note = str(exc)
        return out
    out.ok = _plausible(out)
    return out
