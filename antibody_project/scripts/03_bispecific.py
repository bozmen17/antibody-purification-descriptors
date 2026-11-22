"""Step 3: charge separation between bispecifics and their mispaired species.

    python scripts/03_bispecific.py

Reads the raw Thera-SAbDab download (run step 1 first), writes
`results/bispecific_charge_separation.csv` and two figures, and prints the numbers
quoted in the README.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abdev import plots  # noqa: E402
from abdev.bispecific import CHARGE_GAP_CUTOFF, PI_GAP_CUTOFF, build_pairs  # noqa: E402
from abdev.dataset import load_raw  # noqa: E402

RAW = ROOT / "data" / "raw" / "TheraSAbDab_SeqStruc_OnlineDownload.csv"
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"


def main() -> None:
    if not RAW.exists():
        raise SystemExit("run scripts/01_build_dataset.py first")
    plots.use_house_style()
    RESULTS.mkdir(exist_ok=True)

    df = build_pairs(load_raw(RAW))
    n = len(df)
    print(f"{n} bispecifics with both arms sequenced")
    print(df["format_family"].value_counts().to_string())
    print(f"  common light chain: {int(df['common_light_chain'].sum())}   "
          f"two distinct light chains: "
          f"{int((df['both_light_chains'] & ~df['common_light_chain']).sum())}   "
          f"at least one heavy-only arm: {int((~df['both_light_chains']).sum())}")

    def describe(col: str, unit: str) -> None:
        q = df[col].dropna().quantile([0.25, 0.5, 0.75])
        print(f"  {col:<34} median {q[0.5]:.2f} {unit} (IQR {q[0.25]:.2f}-{q[0.75]:.2f})")

    print("\nseparation between product and mispaired species (Fv level)")
    describe("delta_pI_arms", "pH")
    describe("pI_gap_to_nearest_homodimer", "pH")
    describe("delta_charge_homodimer_pH5p5", "e")
    describe("delta_charge_homodimer_pH7p4", "e")
    describe("delta_charge_lc_mispair_pH5p5", "e")

    hard_q = df["delta_charge_homodimer_pH5p5"] < CHARGE_GAP_CUTOFF
    hard_pi = df["pI_gap_to_nearest_homodimer"] < PI_GAP_CUTOFF
    print(f"\n  homodimer charge gap at pH 5.5 < {CHARGE_GAP_CUTOFF:g} e: "
          f"{int(hard_q.sum())}/{n} ({hard_q.mean():.0%})")
    print(f"  pI gap to nearest homodimer < {PI_GAP_CUTOFF:g}: "
          f"{int(hard_pi.sum())}/{n} ({hard_pi.mean():.0%})")
    print(f"  both: {int((hard_q & hard_pi).sum())}/{n}")

    print("\nby format family (median homodimer charge gap at pH 5.5, share < cut-off):")
    by = df.groupby("format_family").agg(
        n=("therapeutic", "size"),
        median_gap=("delta_charge_homodimer_pH5p5", "median"),
        share_hard=("hard_by_charge_homodimer", "mean"))
    print(by.round(2).to_string())
    cl = df.groupby("common_light_chain")["delta_charge_homodimer_pH5p5"].median()
    print(f"\n  median gap, common LC: {cl.get(True, float('nan')):.2f} e   "
          f"distinct LCs: {cl.get(False, float('nan')):.2f} e")

    out = df.drop(columns=["vh1", "vl1", "vh2", "vl2"])
    out.to_csv(RESULTS / "bispecific_charge_separation.csv", index=False)
    written = [plots.plot_arm_pi(df, FIGURES, PI_GAP_CUTOFF),
               plots.plot_charge_gap(df, FIGURES, CHARGE_GAP_CUTOFF)]
    print("\nwritten")
    print("  results/bispecific_charge_separation.csv")
    for p in written:
        print(f"  {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
