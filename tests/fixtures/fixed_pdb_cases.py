"""Expected results for the pre-fixed PDB regression structures."""

from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]

PRIMARY_INTERACTION_ATTRIBUTES = (
    "hydrogen_bonds",
    "halogen_bonds",
    "pi_interactions",
    "pi_pi_interactions",
    "carbonyl_interactions",
    "n_pi_interactions",
    "water_bridges",
)

COUNT_ATTRIBUTES = PRIMARY_INTERACTION_ATTRIBUTES + (
    "ligand_interactions",
    "cooperativity_chains",
)

FIXED_PDB_CASES = [
    {
        "name": "6rsa.pdb",
        "file": str(REPOSITORY_ROOT / "example_pdb_files/fixed/6rsa_openbabel.pdb"),
        "expected_interactions": [
            "hydrogen_bonds",
            "pi_interactions",
            "carbonyl_interactions",
            "n_pi_interactions",
            "water_bridges",
        ],
        "expected_counts": {
            "hydrogen_bonds": 212,
            "halogen_bonds": 0,
            "pi_interactions": 19,
            "pi_pi_interactions": 0,
            "carbonyl_interactions": 35,
            "n_pi_interactions": 1,
            "water_bridges": 65,
            "ligand_interactions": 21,
            "cooperativity_chains": 40,
        },
        "expected_ligand_interactions": True,
        "expected_ligand_interactions_with_water_bridges": True,
    },
    {
        "name": "7nwd.pdb",
        "file": str(REPOSITORY_ROOT / "example_pdb_files/fixed/7nwd_openbabel.pdb"),
        "expected_interactions": [
            "hydrogen_bonds",
            "pi_pi_interactions",
            "pi_interactions",
        ],
        "expected_counts": {
            "hydrogen_bonds": 27,
            "halogen_bonds": 0,
            "pi_interactions": 3,
            "pi_pi_interactions": 2,
            "carbonyl_interactions": 0,
            "n_pi_interactions": 0,
            "water_bridges": 0,
            "ligand_interactions": 0,
            "cooperativity_chains": 1,
        },
        "expected_ligand_interactions": False,
        "expected_ligand_interactions_with_water_bridges": False,
    },
    {
        "name": "1ubi.pdb",
        "file": str(REPOSITORY_ROOT / "example_pdb_files/fixed/1ubi_openbabel.pdb"),
        "expected_interactions": [
            "hydrogen_bonds",
            "pi_interactions",
            "carbonyl_interactions",
            "water_bridges",
        ],
        "expected_counts": {
            "hydrogen_bonds": 125,
            "halogen_bonds": 0,
            "pi_interactions": 2,
            "pi_pi_interactions": 0,
            "carbonyl_interactions": 21,
            "n_pi_interactions": 0,
            "water_bridges": 15,
            "ligand_interactions": 0,
            "cooperativity_chains": 22,
        },
        "expected_ligand_interactions": False,
        "expected_ligand_interactions_with_water_bridges": False,
    },
    {
        "name": "4laz.pdb",
        "file": str(REPOSITORY_ROOT / "example_pdb_files/fixed/4laz_openbabel.pdb"),
        "expected_interactions": [
            "hydrogen_bonds",
            "halogen_bonds",
            "pi_pi_interactions",
            "pi_interactions",
            "carbonyl_interactions",
            "n_pi_interactions",
            "water_bridges",
        ],
        "expected_counts": {
            "hydrogen_bonds": 863,
            "halogen_bonds": 1,
            "pi_interactions": 65,
            "pi_pi_interactions": 1,
            "carbonyl_interactions": 157,
            "n_pi_interactions": 1,
            "water_bridges": 211,
            "ligand_interactions": 34,
            "cooperativity_chains": 148,
        },
        "expected_ligand_interactions": True,
        "expected_ligand_interactions_with_water_bridges": True,
    },
    {
        "name": "4hhb.pdb",
        "file": str(REPOSITORY_ROOT / "example_pdb_files/fixed/4hhb_openbabel.pdb"),
        "expected_interactions": [
            "hydrogen_bonds",
            "pi_interactions",
            "carbonyl_interactions",
            "water_bridges",
        ],
        "expected_counts": {
            "hydrogen_bonds": 781,
            "halogen_bonds": 0,
            "pi_interactions": 114,
            "pi_pi_interactions": 0,
            "carbonyl_interactions": 242,
            "n_pi_interactions": 0,
            "water_bridges": 79,
            "ligand_interactions": 17,
            "cooperativity_chains": 137,
        },
        "expected_ligand_interactions": True,
        "expected_ligand_interactions_with_water_bridges": True,
    },
]
