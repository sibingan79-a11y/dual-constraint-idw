"""Run all reproducibility scripts."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    "reproduce_table2.py",
    "reproduce_fig5_sensitivity.py",
    "reproduce_fig6_maps.py",
    "reproduce_fig7_profile.py",
]


def main() -> None:
    for script in SCRIPTS:
        path = ROOT / "scripts" / script
        print(f"\n[RUN] {script}")
        subprocess.run([sys.executable, str(path)], cwd=str(ROOT), check=True)
    print("\nAll reproducibility scripts completed.")
    print(f"Outputs are available in: {ROOT / 'outputs'}")


if __name__ == "__main__":
    main()

