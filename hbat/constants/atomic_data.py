"""
Atomic Data Constants

This module contains atomic properties and constants for all elements
commonly found in protein, DNA, RNA, and water molecules in PDB structures.
"""


class AtomicData:
    """Atomic properties and constants.

    This class contains atomic data for all elements commonly found in
    protein, DNA, RNA, and water molecules in PDB structures.
    """

    # Covalent radii in Angstroms
    COVALENT_RADII = {
        # Main biomolecule elements
        "H": 0.31,  # Hydrogen
        "D": 0.31,  # Deuterium (same as H)
        "C": 0.76,  # Carbon
        "N": 0.71,  # Nitrogen
        "O": 0.66,  # Oxygen
        "P": 1.07,  # Phosphorus (DNA/RNA backbone)
        "S": 1.05,  # Sulfur (Cys, Met)
        # Halogens
        "F": 0.57,  # Fluorine
        "CL": 0.99,  # Chlorine
        "BR": 1.14,  # Bromine
        "I": 1.33,  # Iodine
        # Common metals in proteins
        "NA": 1.66,  # Sodium
        "MG": 1.41,  # Magnesium
        "K": 2.03,  # Potassium
        "CA": 1.76,  # Calcium
        "MN": 1.39,  # Manganese
        "FE": 1.32,  # Iron
        "CO": 1.26,  # Cobalt
        "NI": 1.24,  # Nickel
        "CU": 1.32,  # Copper
        "ZN": 1.22,  # Zinc
    }

    # Van der Waals radii in Angstroms
    VDW_RADII = {
        # Main biomolecule elements
        "H": 1.09,  # Hydrogen
        "D": 1.09,  # Deuterium (same as H)
        "C": 1.70,  # Carbon
        "N": 1.55,  # Nitrogen
        "O": 1.52,  # Oxygen
        "P": 1.80,  # Phosphorus (DNA/RNA backbone)
        "S": 1.80,  # Sulfur (Cys, Met)
        # Halogens
        "F": 1.47,  # Fluorine
        "CL": 1.75,  # Chlorine
        "BR": 1.85,  # Bromine
        "I": 1.98,  # Iodine
        # Common metals in proteins
        "NA": 2.27,  # Sodium
        "MG": 1.73,  # Magnesium
        "K": 2.75,  # Potassium
        "CA": 2.31,  # Calcium
        "MN": 2.05,  # Manganese
        "FE": 2.04,  # Iron
        "CO": 2.00,  # Cobalt
        "NI": 1.63,  # Nickel
        "CU": 1.40,  # Copper
        "ZN": 1.39,  # Zinc
    }

    # Electronegativity values (Pauling scale)
    ELECTRONEGATIVITY = {
        # Main biomolecule elements
        "H": 2.20,  # Hydrogen
        "D": 2.20,  # Deuterium (same as H)
        "C": 2.55,  # Carbon
        "N": 3.04,  # Nitrogen
        "O": 3.44,  # Oxygen
        "P": 2.19,  # Phosphorus (DNA/RNA backbone)
        "S": 2.58,  # Sulfur (Cys, Met)
        # Halogens
        "F": 3.98,  # Fluorine
        "CL": 3.16,  # Chlorine
        "BR": 2.96,  # Bromine
        "I": 2.66,  # Iodine
        # Common metals in proteins (less relevant for covalent interactions)
        "NA": 0.93,  # Sodium
        "MG": 1.31,  # Magnesium
        "K": 0.82,  # Potassium
        "CA": 1.00,  # Calcium
        "MN": 1.55,  # Manganese
        "FE": 1.83,  # Iron
        "CO": 1.88,  # Cobalt
        "NI": 1.91,  # Nickel
        "CU": 1.90,  # Copper
        "ZN": 1.65,  # Zinc
    }

    # Atomic masses in amu (atomic mass units)
    ATOMIC_MASSES = {
        # Main biomolecule elements
        "H": 1.008,  # Hydrogen
        "D": 2.014,  # Deuterium (heavy hydrogen)
        "C": 12.011,  # Carbon
        "N": 14.007,  # Nitrogen
        "O": 15.999,  # Oxygen
        "P": 30.974,  # Phosphorus (DNA/RNA backbone)
        "S": 32.065,  # Sulfur (Cys, Met)
        # Halogens
        "F": 18.998,  # Fluorine
        "CL": 35.453,  # Chlorine
        "BR": 79.904,  # Bromine
        "I": 126.904,  # Iodine
        # Common metals in proteins
        "NA": 22.990,  # Sodium
        "MG": 24.305,  # Magnesium
        "K": 39.098,  # Potassium
        "CA": 40.078,  # Calcium
        "MN": 54.938,  # Manganese
        "FE": 55.845,  # Iron
        "CO": 58.933,  # Cobalt
        "NI": 58.693,  # Nickel
        "CU": 63.546,  # Copper
        "ZN": 65.38,  # Zinc
    }

    # Default atomic mass for unknown elements
    DEFAULT_ATOMIC_MASS = 12.011  # Carbon mass

    # Hydrogen detection threshold
    MIN_HYDROGEN_RATIO = 0.25  # 25% of atoms must be hydrogen
