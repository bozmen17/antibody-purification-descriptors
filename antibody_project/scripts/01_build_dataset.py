"""Step 1: download Thera-SAbDab if needed, compute descriptors, write the table.

    python scripts/01_build_dataset.py
    python scripts/01_build_dataset.py --refresh     # re-download the source CSV
"""
from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abdev.dataset import THERA_SABDAB_URL, build, feature_columns  # noqa: E402

RAW = ROOT / "data" / "raw" / "TheraSAbDab_SeqStruc_OnlineDownload.csv"
OUT = ROOT / "data" / "processed" / "fv_descriptors.csv"


def download(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {THERA_SABDAB_URL}")
    with urllib.request.urlopen(THERA_SABDAB_URL, timeout=120) as response:
        path.write_bytes(response.read())
    print(f"  -> {path.relative_to(ROOT)} ({path.stat().st_size / 1e6:.1f} MB)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true", help="re-download the source CSV")
    args = ap.parse_args()

    if args.refresh or not RAW.exists():
        download(RAW)

    df = build(RAW)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)

    core, loop_dependent = feature_columns(df)
    print(f"\n{len(df)} variable-domain arms from {df['therapeutic'].nunique()} therapeutics")
    print(f"  paired Fv: {int(df['paired'].sum())}   "
          f"heavy-only: {int((~df['paired']).sum())}")
    print(f"  loops resolved by the anchor rules: {int(df['loops_resolved'].sum())} "
          f"({df['loops_resolved'].mean():.1%})")
    print(f"  descriptors: {len(core)} defined for every arm, "
          f"{len(loop_dependent)} needing resolved loops or a light chain")
    print("\nformat families:")
    print(df["format_family"].value_counts().to_string())
    print(f"\nwrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
