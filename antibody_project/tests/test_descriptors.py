"""Checks on the descriptor code, including loop boundaries against known antibodies.

The CDR assertions use published Kabat/AbM loops for trastuzumab and adalimumab, so a
change to the anchor rules that breaks a real antibody fails here rather than silently
shifting every descriptor downstream.
"""
from __future__ import annotations

import math

import pytest

from abdev.charge import isoelectric_point, net_charge, residue_charge
from abdev.dataset import classify_format
from abdev.descriptors import clean_sequence, fv_descriptors, hydrophobicity_index
from abdev.regions import heavy_loops, light_loops

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
ADALIMUMAB_VH = (
    "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEG"
    "RFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS"
)


# --------------------------------------------------------------------------- charge
def test_residue_charge_signs():
    assert residue_charge("K", 7.4) > 0.9
    assert residue_charge("E", 7.4) < -0.9
    assert residue_charge("A", 7.4) == 0.0
    # Histidine is half-protonated near its pKa, which is why it dominates pH 6-7 behaviour.
    assert 0.3 < residue_charge("H", 6.5) < 0.7


def test_net_charge_falls_monotonically_with_ph():
    peptide = "KRHDEYACGGKK"
    charges = [net_charge(peptide, ph) for ph in (3.0, 5.0, 7.0, 9.0, 11.0)]
    assert charges == sorted(charges, reverse=True)


def test_pi_matches_composition():
    assert isoelectric_point("K" * 20) > 10.0
    assert isoelectric_point("E" * 20) < 4.5
    # At its own pI a molecule is uncharged, by definition.
    for seq in ("KRHDEYACGG", TRASTUZUMAB_VH, TRASTUZUMAB_VL):
        assert abs(net_charge(seq, isoelectric_point(seq))) < 0.02


def test_pi_of_paired_fv_accounts_for_both_chains():
    """Two chains carry two sets of termini, so they are passed separately, not concatenated."""
    paired = isoelectric_point(TRASTUZUMAB_VH, TRASTUZUMAB_VL)
    fused = isoelectric_point(TRASTUZUMAB_VH + TRASTUZUMAB_VL)
    assert abs(sum(net_charge(c, paired) for c in (TRASTUZUMAB_VH, TRASTUZUMAB_VL))) < 0.02
    # The two routes agree closely: the difference is one extra pair of termini.
    assert abs(paired - fused) < 0.5


# --------------------------------------------------------------------------- loops
def test_trastuzumab_heavy_loops():
    loops = heavy_loops(TRASTUZUMAB_VH)
    assert loops.ok, loops.note
    assert loops.cdr1 == "GFNIKDTYIH"          # AbM H1, Kabat 26-35
    assert loops.cdr2 == "RIYPTNGYTRYADSVK"    # Kabat H2, 50-65
    assert loops.cdr3 == "WGGDGFYAMDY"         # Kabat H3, 95-102


def test_trastuzumab_light_loops():
    loops = light_loops(TRASTUZUMAB_VL)
    assert loops.ok, loops.note
    assert loops.cdr1 == "RASQDVNTAVA"         # Kabat L1, 24-34
    assert loops.cdr2 == "SASFLYS"             # Kabat L2, 50-56
    assert loops.cdr3 == "QQHYTTPPT"           # Kabat L3, 89-97


def test_adalimumab_heavy_cdr3():
    loops = heavy_loops(ADALIMUMAB_VH)
    assert loops.ok, loops.note
    assert loops.cdr3 == "VSYLSTASSLDY"


def test_loops_are_substrings_in_order():
    for seq, finder in ((TRASTUZUMAB_VH, heavy_loops), (TRASTUZUMAB_VL, light_loops)):
        loops = finder(seq)
        positions = [seq.index(c) for c in loops.as_list()]
        assert positions == sorted(positions)


def test_nonsense_sequence_is_flagged_not_guessed():
    loops = heavy_loops("A" * 120)
    assert not loops.ok
    assert loops.note


# --------------------------------------------------------------------------- descriptors
def test_clean_sequence_rejects_placeholders_and_junk():
    assert clean_sequence("na") is None
    assert clean_sequence(None) is None
    assert clean_sequence(float("nan")) is None
    assert clean_sequence("ACDEF") is None                  # too short to be a domain
    assert clean_sequence("B" * 120) is None                # non-standard residue
    assert clean_sequence(TRASTUZUMAB_VH.lower()) == TRASTUZUMAB_VH


def test_hydrophobicity_index_ordering():
    apolar_rich = hydrophobicity_index(["IIIILLLLVVVV", "KE"])
    polar_rich = hydrophobicity_index(["IL", "KKKKEEEERRRR"])
    assert apolar_rich > polar_rich > 0
    assert math.isnan(hydrophobicity_index(["IIII"]))       # no hydrophilic term to divide by
    assert math.isnan(hydrophobicity_index([]))


def test_fv_descriptors_paired_and_unpaired():
    paired = fv_descriptors(TRASTUZUMAB_VH, TRASTUZUMAB_VL)
    single = fv_descriptors(TRASTUZUMAB_VH)

    assert paired["paired"] and not single["paired"]
    assert paired["fv_length"] == len(TRASTUZUMAB_VH) + len(TRASTUZUMAB_VL)
    assert paired["loops_resolved"]
    assert 4 < paired["pI_fv"] < 10
    assert paired["charge_pH3p6"] > paired["charge_pH7p4"]  # protonated at low pH
    assert paired["n_cysteine_fv"] >= 4                     # two intradomain disulfides
    assert math.isnan(single["pI_vl"])
    assert paired["cdr3_h_length"] == 11
    # Descriptor set is the same shape whether or not a light chain is present.
    assert set(paired) == set(single)


def test_descriptor_values_are_finite_where_defined():
    d = fv_descriptors(TRASTUZUMAB_VH, TRASTUZUMAB_VL)
    for key, value in d.items():
        if isinstance(value, float):
            assert math.isfinite(value), key


def test_missing_light_chain_in_table_is_treated_as_unpaired():
    """pandas stores a missing VL as NaN; it must reach fv_descriptors as None."""
    import pandas as pd
    from abdev.dataset import add_descriptors

    arms = pd.DataFrame({"vh": [TRASTUZUMAB_VH, TRASTUZUMAB_VH],
                         "vl": [TRASTUZUMAB_VL, None]})
    out = add_descriptors(arms)
    assert list(out["paired"]) == [True, False]


# --------------------------------------------------------------------------- format labels
@pytest.mark.parametrize("raw,expected", [
    ("Whole mAb", "Whole mAb"),
    ("Whole mAb ADC", "IgG conjugate"),
    ("Whole mAb Radiolabelled", "IgG conjugate"),
    ("Whole mAb Fusion", "IgG fusion"),
    ("Bispecific mAb", "Multispecific (IgG-like)"),
    ("Bispecific scFv", "Multispecific (fragment-based)"),
    ("Bispecific Mixed mAb and scFv", "Multispecific (IgG-like)"),
    ("Fab", "Fab / Fab fusion"),
    ("scFv", "scFv / scFv fusion"),
    ("Nanobody", "Single domain"),
    ("Single Domain Variable Fragment;H", "Single domain"),
    ("", "Unclassified"),
])
def test_classify_format(raw, expected):
    assert classify_format(raw) == expected


# --------------------------------------------------------------------------- bispecific pairs
def test_identical_arms_have_no_charge_handle():
    from abdev.bispecific import pair_descriptors

    d = pair_descriptors(TRASTUZUMAB_VH, TRASTUZUMAB_VL, TRASTUZUMAB_VH, TRASTUZUMAB_VL)
    assert d["common_light_chain"]
    assert d["delta_pI_arms"] == 0
    assert d["delta_charge_homodimer_pH5p5"] == 0
    assert d["pI_gap_to_nearest_homodimer"] == 0
    assert d["delta_charge_lc_mispair_pH5p5"] == 0
    assert d["hard_by_charge_homodimer"]


def test_homodimer_charge_gap_equals_arm_charge_difference():
    """Charge(AB) - charge(AA) = q_B - q_A exactly, so the gap is the arm difference."""
    from abdev.bispecific import pair_descriptors

    d = pair_descriptors(TRASTUZUMAB_VH, TRASTUZUMAB_VL, ADALIMUMAB_VH, None)
    q_a = net_charge(TRASTUZUMAB_VH, 5.5) + net_charge(TRASTUZUMAB_VL, 5.5)
    q_b = net_charge(ADALIMUMAB_VH, 5.5)
    assert abs(d["delta_charge_homodimer_pH5p5"] - abs(q_a - q_b)) < 1e-3
    # Heterodimer pI lies between the two homodimers.
    lo, hi = sorted((d["pI_homodimer_aa"], d["pI_homodimer_bb"]))
    assert lo <= d["pI_heterodimer"] <= hi
    assert math.isnan(d["delta_charge_lc_mispair_pH5p5"])   # arm B has no light chain
