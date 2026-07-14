"""One-off generator for e2e test fixtures under example_pdb_files/fixed/.

Runs OpenBabel and PDBFixer once against each raw structure used by the e2e
suite and commits the result as a static fixture, so tests read a pre-fixed
file with fix_pdb_enabled=False instead of re-running (stochastic) fixing at
test time. Re-run manually if a source PDB in example_pdb_files/ changes.

Usage: PYENV_VERSION=hbat python scripts/generate_fixed_pdb_fixtures.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hbat.core.pdb_fixer import PDBFixer
from hbat.core.analysis import AnalysisParameters
from hbat.core.analyzer import MolecularInteractionAnalyzer

SOURCE_DIR = "example_pdb_files"
OUTPUT_DIR = os.path.join(SOURCE_DIR, "fixed")

STEMS = ["6rsa", "1ubi", "4hhb", "4laz", "7nwd"]

INTERACTION_ATTRS = [
    "hydrogen_bonds",
    "halogen_bonds",
    "pi_interactions",
    "pi_pi_interactions",
    "carbonyl_interactions",
    "n_pi_interactions",
    "water_bridges",
]


def generate() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fixer = PDBFixer()

    for stem in STEMS:
        input_path = os.path.join(SOURCE_DIR, f"{stem}.pdb")
        if not os.path.exists(input_path):
            print(f"SKIP {stem}: {input_path} not found")
            continue

        for method in ("openbabel", "pdbfixer"):
            output_path = os.path.join(OUTPUT_DIR, f"{stem}_{method}.pdb")
            fixer.fix_pdb_file_to_file(
                input_pdb_path=input_path,
                output_pdb_path=output_path,
                method=method,
                add_hydrogens=True,
                add_heavy_atoms=(method == "pdbfixer"),
            )
            print(f"wrote {output_path}")


def report() -> None:
    print("\nInteraction counts on generated fixtures (fix_pdb_enabled=False):")
    for stem in STEMS:
        for method in ("openbabel", "pdbfixer"):
            path = os.path.join(OUTPUT_DIR, f"{stem}_{method}.pdb")
            if not os.path.exists(path):
                continue
            params = AnalysisParameters(fix_pdb_enabled=False)
            analyzer = MolecularInteractionAnalyzer(params)
            success = analyzer.analyze_file(path)
            counts = {
                attr: len(getattr(analyzer, attr) or []) for attr in INTERACTION_ATTRS
            }
            lig = (
                len(analyzer.ligand_interactions.interactions)
                if analyzer.ligand_interactions
                else 0
            )
            counts["ligand_interactions"] = lig
            print(f"{stem}_{method}.pdb success={success} {counts}")


if __name__ == "__main__":
    generate()
    report()
