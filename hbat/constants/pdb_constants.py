"""
PDB Structure Constants

This module contains constants specifically related to PDB file processing,
including residue mappings, atom classifications, and molecular recognition
patterns used throughout HBAT's structure analysis components.
"""

from typing import Dict, List

# Non-standard to standard protein residue substitutions
PROTEIN_SUBSTITUTIONS: Dict[str, str] = {
    "2AS": "ASP",
    "3AH": "HIS",
    "5HP": "GLU",
    "5OW": "LYS",
    "ACL": "ARG",
    "AGM": "ARG",
    "AIB": "ALA",
    "ALM": "ALA",
    "ALO": "THR",
    "ALY": "LYS",
    "ARM": "ARG",
    "ASA": "ASP",
    "ASB": "ASP",
    "ASK": "ASP",
    "ASL": "ASP",
    "ASQ": "ASP",
    "AYA": "ALA",
    "BCS": "CYS",
    "BHD": "ASP",
    "BMT": "THR",
    "BNN": "ALA",
    "BUC": "CYS",
    "BUG": "LEU",
    "C5C": "CYS",
    "C6C": "CYS",
    "CAS": "CYS",
    "CCS": "CYS",
    "CEA": "CYS",
    "CGU": "GLU",
    "CHG": "ALA",
    "CLE": "LEU",
    "CME": "CYS",
    "CSD": "ALA",
    "CSO": "CYS",
    "CSP": "CYS",
    "CSS": "CYS",
    "CSW": "CYS",
    "CSX": "CYS",
    "CXM": "MET",
    "CY1": "CYS",
    "CY3": "CYS",
    "CYG": "CYS",
    "CYM": "CYS",
    "CYQ": "CYS",
    "DAH": "PHE",
    "DAL": "ALA",
    "DAR": "ARG",
    "DAS": "ASP",
    "DCY": "CYS",
    "DGL": "GLU",
    "DGN": "GLN",
    "DHA": "ALA",
    "DHI": "HIS",
    "DIL": "ILE",
    "DIV": "VAL",
    "DLE": "LEU",
    "DLY": "LYS",
    "DNP": "ALA",
    "DPN": "PHE",
    "DPR": "PRO",
    "DSN": "SER",
    "DSP": "ASP",
    "DTH": "THR",
    "DTR": "TRP",
    "DTY": "TYR",
    "DVA": "VAL",
    "EFC": "CYS",
    "FLA": "ALA",
    "FME": "MET",
    "GGL": "GLU",
    "GL3": "GLY",
    "GLZ": "GLY",
    "GMA": "GLU",
    "GSC": "GLY",
    "HAC": "ALA",
    "HAR": "ARG",
    "HIC": "HIS",
    "HIP": "HIS",
    "HMR": "ARG",
    "HPQ": "PHE",
    "HTR": "TRP",
    "HYP": "PRO",
    "IAS": "ASP",
    "IIL": "ILE",
    "IYR": "TYR",
    "KCX": "LYS",
    "LLP": "LYS",
    "LLY": "LYS",
    "LTR": "TRP",
    "LYM": "LYS",
    "LYZ": "LYS",
    "MAA": "ALA",
    "MEN": "ASN",
    "MHS": "HIS",
    "MIS": "SER",
    "MK8": "LEU",
    "MLE": "LEU",
    "MPQ": "GLY",
    "MSA": "GLY",
    "MSE": "MET",
    "MVA": "VAL",
    "NEM": "HIS",
    "NEP": "HIS",
    "NLE": "LEU",
    "NLN": "LEU",
    "NLP": "LEU",
    "NMC": "GLY",
    "OAS": "SER",
    "OCS": "CYS",
    "OMT": "MET",
    "PAQ": "TYR",
    "PCA": "GLU",
    "PEC": "CYS",
    "PHI": "PHE",
    "PHL": "PHE",
    "PR3": "CYS",
    "PRR": "ALA",
    "PTR": "TYR",
    "PYX": "CYS",
    "SAC": "SER",
    "SAR": "GLY",
    "SCH": "CYS",
    "SCS": "CYS",
    "SCY": "CYS",
    "SEL": "SER",
    "SEP": "SER",
    "SET": "SER",
    "SHC": "CYS",
    "SHR": "LYS",
    "SMC": "CYS",
    "SOC": "CYS",
    "STY": "TYR",
    "SVA": "SER",
    "TIH": "ALA",
    "TPL": "TRP",
    "TPO": "THR",
    "TPQ": "ALA",
    "TRG": "LYS",
    "TRO": "TRP",
    "TYB": "TYR",
    "TYI": "TYR",
    "TYQ": "TYR",
    "TYS": "TYR",
    "TYY": "TYR",
}
"""Dict[str, str]: Mapping of non-standard protein residue codes to their standard amino acid equivalents.

This comprehensive dictionary provides substitutions for modified, methylated, phosphorylated,
and other chemically altered amino acid residues commonly found in PDB structures. Used by
PDB fixing operations to standardize protein residue names for consistent analysis.

Examples:
    - MSE (selenomethionine) → MET (methionine)
    - CSO (cysteine sulfenic acid) → CYS (cysteine)
    - HYP (hydroxyproline) → PRO (proline)
    - PCA (pyroglutamic acid) → GLU (glutamic acid)

Note: This dictionary contains only protein residue substitutions. Nucleotide
modifications are handled separately.
"""

PROTEIN_RESIDUES: List[str] = [
    "ALA",
    "ASN",
    "CYS",
    "GLU",
    "HIS",
    "LEU",
    "MET",
    "PRO",
    "THR",
    "TYR",
    "ARG",
    "ASP",
    "GLN",
    "GLY",
    "ILE",
    "LYS",
    "PHE",
    "SER",
    "TRP",
    "VAL",
]
"""List[str]: Standard three-letter codes for the 20 canonical amino acid residues.

This list contains all naturally occurring protein amino acids in their standard
three-letter abbreviation format as used in PDB files. Used for residue type
validation, protein chain identification, and analysis scope determination.

The 20 amino acids are:
    - Alanine (ALA), Arginine (ARG), Asparagine (ASN), Aspartic acid (ASP)
    - Cysteine (CYS), Glutamic acid (GLU), Glutamine (GLN), Glycine (GLY)
    - Histidine (HIS), Isoleucine (ILE), Leucine (LEU), Lysine (LYS)
    - Methionine (MET), Phenylalanine (PHE), Proline (PRO), Serine (SER)
    - Threonine (THR), Tryptophan (TRP), Tyrosine (TYR), Valine (VAL)
"""

RNA_RESIDUES: List[str] = ["A", "G", "C", "U", "I"]
"""List[str]: Standard single-letter codes for RNA nucleotide residues.

Contains the five RNA nucleotides commonly found in PDB structures:
    - A (Adenine): Purine base forming A-U base pairs
    - G (Guanine): Purine base forming G-C base pairs  
    - C (Cytosine): Pyrimidine base forming C-G base pairs
    - U (Uracil): Pyrimidine base forming U-A base pairs
    - I (Inosine): Modified nucleotide, wobble base pairing

Used for nucleic acid chain identification and RNA structure analysis.
"""

DNA_RESIDUES: List[str] = ["DA", "DG", "DC", "DT", "DI"]
"""List[str]: Standard two-letter codes for DNA nucleotide residues.

Contains the five DNA nucleotides commonly found in PDB structures:
    - DA (Deoxyadenosine): Purine base forming A-T base pairs
    - DG (Deoxyguanosine): Purine base forming G-C base pairs
    - DC (Deoxycytidine): Pyrimidine base forming C-G base pairs
    - DT (Deoxythymidine): Pyrimidine base forming T-A base pairs
    - DI (Deoxyinosine): Modified nucleotide, wobble base pairing

Used for nucleic acid chain identification and DNA structure analysis.
The 'D' prefix distinguishes DNA nucleotides from RNA nucleotides.
"""


# Pre-computed mapping for common PDB atoms (for performance)
_COMMON_PDB_ATOMS: Dict[str, str] = {
    # Protein backbone
    "N": "N",
    "CA": "C",
    "C": "C",
    "O": "O",
    # Common side chain atoms
    "CB": "C",
    "CG": "C",
    "CD": "C",
    "CE": "C",
    "CZ": "C",
    "ND1": "N",
    "ND2": "N",
    "NE": "N",
    "NE1": "N",
    "NE2": "N",
    "NH1": "N",
    "NH2": "N",
    "NZ": "N",
    "OD1": "O",
    "OD2": "O",
    "OE1": "O",
    "OE2": "O",
    "OG": "O",
    "OG1": "O",
    "OH": "O",
    "SG": "S",
    "SD": "S",
    # DNA/RNA backbone
    "P": "P",
    "OP1": "O",
    "OP2": "O",
    "O5'": "O",
    "C5'": "C",
    "C4'": "C",
    "O4'": "O",
    "C3'": "C",
    "O3'": "O",
    "C2'": "C",
    "O2'": "O",
    "C1'": "C",
    # Nucleotide bases
    "N1": "N",
    "C2": "C",
    "N2": "N",
    "O2": "O",
    "N3": "N",
    "C4": "C",
    "N4": "N",
    "O4": "O",
    "C5": "C",
    "C5M": "C",
    "C6": "C",
    "N6": "N",
    "O6": "O",
    "N7": "N",
    "C8": "C",
    "N9": "N",
    # Common hydrogens
    "H": "H",
    "HA": "H",
    "HB": "H",
    "HG": "H",
    "HD": "H",
    "HE": "H",
    "HZ": "H",
    "HH": "H",
    "HN": "H",
    "HO": "H",
    "HS": "H",
    "HOH": "H",
    # Water
    "OH2": "O",
    # Common heteroatoms and metals
    "F": "F",
    "CL": "CL",
    "BR": "BR",
    "I": "I",
    "D": "D",
}


# Maintain backward compatibility
PDB_ATOM_TO_ELEMENT: Dict[str, str] = _COMMON_PDB_ATOMS
"""Dict[str, str]: Pre-computed mapping of common PDB atom names to their element types.

This dictionary provides fast lookup for the most frequently encountered PDB atoms.
For comprehensive coverage including unusual atoms, use pdb_atom_to_element() function
which uses regex-based pattern matching.

Coverage includes:
    - Protein backbone and common side chain atoms
    - DNA/RNA backbone and nucleotide base atoms  
    - Standard hydrogen atoms
    - Water molecules

For full pattern-based mapping that handles:
    - Greek letter remoteness indicators (CA, CB, CG, CD, CE, CZ, CH)
    - Numbered variants (C1', H2'', OP1, etc.)
    - Ion charges (CA2+, MG2+, etc.)
    - IUPAC hydrogen naming conventions
    - Uncommon PDB atom names

Use pdb_atom_to_element() function instead.

Used for:
    - Looking up atomic properties (radius, mass, electronegativity)
    - Covalent bond detection
    - Van der Waals calculations
    - Molecular mass calculations
"""

PROTEIN_BACKBONE_ATOMS: List[str] = ["N", "CA", "C", "O"]
"""List[str]: Standard protein backbone atom names in PDB format.

Defines the four atoms that form the protein backbone (main chain):
    - N: Amino nitrogen atom
    - CA: Alpha carbon atom (central carbon)
    - C: Carbonyl carbon atom
    - O: Carbonyl oxygen atom

These atoms are present in all amino acid residues (except proline's modified N)
and form the peptide bonds that connect residues.
"""

DNA_RNA_BACKBONE_ATOMS: List[str] = [
    "P",
    "OP1",
    "OP2",
    "O5'",
    "C5'",
    "C4'",
    "O4'",
    "C3'",
    "O3'",
    "C2'",
    "O2'",
    "C1'",
]
"""List[str]: Standard DNA/RNA backbone atom names in PDB format.

Sugar-phosphate backbone atoms:
    - P: Phosphorus atom
    - OP1, OP2: Non-bridging phosphate oxygens
    - O5': 5' phosphate oxygen (bridging)
    - C5': 5' carbon of ribose/deoxyribose
    - C4': 4' carbon of ribose/deoxyribose
    - O4': 4' oxygen of ribose/deoxyribose (ring oxygen)
    - C3': 3' carbon of ribose/deoxyribose
    - O3': 3' phosphate oxygen (bridging)
    - C2': 2' carbon of ribose/deoxyribose
    - O2': 2' hydroxyl oxygen (RNA only, absent in DNA)
    - C1': 1' carbon of ribose/deoxyribose (anomeric carbon)

Note: O2' is present in RNA but absent in DNA (deoxyribose lacks 2' hydroxyl).
"""

BACKBONE_ATOMS: List[str] = PROTEIN_BACKBONE_ATOMS + DNA_RNA_BACKBONE_ATOMS
"""List[str]: Combined backbone atom names for proteins, DNA, and RNA in PDB format.

This list is the combination of PROTEIN_BACKBONE_ATOMS and DNA_RNA_BACKBONE_ATOMS,
providing a comprehensive set of backbone atoms for all major biomolecule types.

Used for:
    - Backbone hydrogen bond identification across all molecule types
    - Secondary structure analysis
    - Main chain vs side chain/base classification
    - Nucleic acid backbone conformation analysis
"""

PROTEIN_SIDECHAIN_ATOMS: List[str] = [
    "CB",
    "CG",
    "CD",
    "NE",
    "CZ",
    "NH1",
    "NH2",
    "OD1",
    "ND2",
    "OD2",
    "SG",
    "OE1",
    "NE2",
    "OE2",
    "CD2",
    "ND1",
    "CE1",
    "CG1",
    "CG2",
    "CD1",
    "CE",
    "NZ",
    "SD",
    "CE2",
    "OG",
    "OG1",
    "NE1",
    "CE3",
    "CZ2",
    "CZ3",
    "CH2",
    "OH",
]
"""List[str]: Common protein side chain atom names in PDB format.

Comprehensive list of side chain (R-group) atoms found in the 20 standard amino acids:
    - Aliphatic carbons: CB, CG, CD, CE, CZ (branching from CA)
    - Aromatic carbons: CD1/CD2, CE1/CE2/CE3, CZ2/CZ3, CH2 (ring systems)
    - Nitrogen atoms: NE, NH1, NH2, ND1, ND2, NE1, NE2, NZ (basic groups)
    - Oxygen atoms: OD1, OD2, OE1, OE2, OG, OG1, OH (acidic/hydroxyl groups)
    - Sulfur atoms: SG, SD (cysteine, methionine)

Used for:
    - Side chain interaction analysis
    - Functional group identification
    - Hydrogen bond donor/acceptor classification
"""

DNA_RNA_BASE_ATOMS: List[str] = [
    # Purine base atoms (Adenine, Guanine)
    "N1",
    "C2",
    "N3",
    "C4",
    "C5",
    "C6",
    "N6",
    "N7",
    "C8",
    "N9",
    "O6",
    "N2",
    # Pyrimidine base atoms (Cytosine, Thymine, Uracil)
    "O2",
    "N4",
    "O4",
    "C5M",
]
"""List[str]: Common DNA/RNA base atom names in PDB format.

Base atoms found in nucleotides:

Purine bases (Adenine, Guanine):
  - N1, C2, N3, C4, C5, C6: Six-membered ring atoms
  - N7, C8, N9: Five-membered ring atoms
  - N6: Amino group on adenine
  - O6, N2: Functional groups on guanine

Pyrimidine bases (Cytosine, Thymine, Uracil):
  - N1, C2, N3, C4, C5, C6: Six-membered ring atoms
  - O2: Carbonyl oxygen at position 2
  - N4: Amino group on cytosine
  - O4: Carbonyl oxygen at position 4 (thymine/uracil)
  - C5M: Methyl group on thymine (also called C7)

Used for:
  - Base-base interactions (hydrogen bonding, stacking)
  - Protein-nucleic acid recognition
  - Base functional group identification
"""

SIDECHAIN_ATOMS: List[str] = PROTEIN_SIDECHAIN_ATOMS + DNA_RNA_BASE_ATOMS
"""List[str]: Combined side chain and base atoms for proteins and nucleic acids.

This list is the combination of PROTEIN_SIDECHAIN_ATOMS and DNA_RNA_BASE_ATOMS,
providing a comprehensive set of non-backbone atoms for all major biomolecule types.

Used for:
    - Side chain/base interaction analysis
    - Distinguishing backbone from functional groups
    - Molecular recognition studies
"""


WATER_MOLECULES: List[str] = ["HOH", "WAT", "DOD", "TIP3", "TIP4", "TIP5", "W"]
"""List[str]: Standard water molecule residue names in PDB files.

Recognition patterns for different water representations:
    - HOH: Standard PDB water molecule designation
    - WAT: Alternative water molecule name
    - DOD: Deuterated water (heavy water)
    - TIP3: TIP3P water model (3-point)
    - TIP4: TIP4P water model (4-point)
    - TIP5: TIP5P water model (5-point)
    - W: Abbreviated water designation

Used for:
    - Water molecule identification in PDB structures
    - Solvent exclusion during analysis
    - Water-mediated interaction detection
    - Hydration shell analysis
"""

RESIDUES: List[str] = PROTEIN_RESIDUES + DNA_RESIDUES + RNA_RESIDUES + WATER_MOLECULES
"""List[str]: Combined list of all standard residue codes for proteins, DNA, and RNA.

This list is the combination of PROTEIN_RESIDUES, DNA_RESIDUES, WATER_MOLECULES, and RNA_RESIDUES,
providing a comprehensive set of standard residues found in biomolecular structures.

Used for:
    - General residue type validation
    - Distinguishing standard residues from heterogens
    - Biomolecule type identification
"""

RESIDUES_WITH_AROMATIC_RINGS: List[str] = [
    # Protein residues
    "PHE",  # Phenylalanine
    "TYR",  # Tyrosine
    "TRP",  # Tryptophan
    "HIS",  # Histidine
    "HID",  # Histidine (deprotonated)
    "HIE",  # Histidine (protonated)
    "HIP",  # Histidine (protonated)
    "TYI",  # Tyrosine (ionized)
    "TYQ",  # Tyrosine (quinone)
    "TYB",  # Tyrosine (brominated)
    # DNA nucleotides
    "DA",  # Deoxyadenosine (purine)
    "DG",  # Deoxyguanosine (purine)
    "DC",  # Deoxycytidine (pyrimidine)
    "DT",  # Deoxythymidine (pyrimidine)
    # RNA nucleotides
    "A",  # Adenine (purine)
    "G",  # Guanine (purine)
    "C",  # Cytosine (pyrimidine)
    "U",  # Uracil (pyrimidine)
]

"""List[str]: Residues containing aromatic rings in their structures.

This list includes:

Protein residues:
  - PHE: Phenylalanine (benzene ring)
  - TYR: Tyrosine (phenolic ring)
  - TRP: Tryptophan (indole ring)
  - HIS: Histidine (imidazole ring)
  - HID, HIE, HIP: Different protonation states of histidine
  - TYI, TYQ, TYB: Variants of tyrosine with modifications

DNA nucleotides:
  - DA: Deoxyadenosine (purine ring: adenine)
  - DG: Deoxyguanosine (purine ring: guanine)
  - DC: Deoxycytidine (pyrimidine ring: cytosine)
  - DT: Deoxythymidine (pyrimidine ring: thymine)

RNA nucleotides:
  - A: Adenine (purine ring)
  - G: Guanine (purine ring)
  - C: Cytosine (pyrimidine ring)
  - U: Uracil (pyrimidine ring)

Used for:
  - Aromatic interaction analysis
  - π-π stacking detection between proteins and nucleic acids
  - DNA/RNA-protein interface studies
"""

# Molecular interaction element lists
HYDROGEN_ELEMENTS: List[str] = ["H", "D"]
"""List[str]: Hydrogen element types including isotopes.

Contains the hydrogen element symbols commonly found in PDB structures:
- H: Standard hydrogen (protium)
- D: Deuterium (heavy hydrogen isotope)

Used for:
    - Hydrogen bond donor/acceptor detection
    - Identifying hydrogen atoms in molecular interactions
    - Mass calculations and isotope effects
    - NMR-related structural analysis
"""

HALOGEN_ELEMENTS: List[str] = ["F", "CL", "BR", "I"]
"""List[str]: Elements that can participate in halogen bonding as donors.

These halogens can act as electron acceptors in halogen bonds when covalently
bonded to carbon (C-X...Y geometry). The halogen forms a σ-hole that can interact
with electron-rich regions on acceptor atoms.

- F: Fluorine (weakest halogen bond donor due to high electronegativity)
- CL: Chlorine (common in drug design, moderate halogen bonding)
- BR: Bromine (strong halogen bond donor, commonly studied)
- I: Iodine (strongest halogen bond donor due to large, polarizable electron cloud)
"""

HYDROGEN_BOND_DONOR_ELEMENTS: List[str] = ["N", "O", "S", "F", "C"]
"""List[str]: Elements that can act as hydrogen bond donors.

These elements can form hydrogen bonds when covalently bonded to hydrogen atoms
(D-H...A geometry). They are electronegative enough to polarize the D-H bond,
creating a partial positive charge on the hydrogen that can interact with
electron-rich acceptor atoms.

- N: Nitrogen (amino groups, ring nitrogens, strong donors)
- O: Oxygen (hydroxyl groups, moderate to strong donors)
- S: Sulfur (thiol groups, weak donors due to lower electronegativity)
"""

HYDROGEN_BOND_ACCEPTOR_ELEMENTS: List[str] = ["N", "O", "S", "F", "CL"]
"""List[str]: Elements that can act as hydrogen bond acceptors.

These electronegative elements have lone pairs of electrons that can accept
hydrogen bonds from donor atoms (D-H...A geometry). They can form favorable
electrostatic interactions with the partial positive charge on hydrogen.

- N: Nitrogen (lone pairs on amino groups, ring nitrogens)
- O: Oxygen (lone pairs on carbonyl, hydroxyl, ether groups - strongest acceptors)
- S: Sulfur (lone pairs on thiol, sulfide groups - weaker acceptors)
- F: Fluorine (strongest electronegativity, excellent acceptor but rare in proteins)
- CL: Chlorine (moderate acceptor, sometimes found in modified residues)
"""

HALOGEN_BOND_ACCEPTOR_ELEMENTS: List[str] = ["N", "O", "S", "P", "SE"]
"""List[str]: Elements that can act as halogen bond acceptors.

These electronegative atoms can donate electron density to the σ-hole of
halogen atoms in halogen bonds. They typically have lone pairs of electrons
that can interact with the positive electrostatic potential of the halogen.

- N: Nitrogen (lone pairs on amino groups, ring nitrogens)
- O: Oxygen (lone pairs on carbonyl, hydroxyl, ether groups)
- S: Sulfur (lone pairs on thiol, sulfide groups, weaker than N/O)
- P: Phosphorus (lone pairs in phosphate groups, moderate acceptor)
- SE: Selenium (lone pairs in selenocysteine, rare but possible acceptor)
"""

PI_INTERACTION_DONOR: List[str] = ["C", "N", "O", "S"]
"""List[str]: Elements that can act as π-interaction donors.

These atoms can participate in π-interactions when bonded to interaction atoms:
- C: Carbon atoms (C-Cl...π, C-Br...π, C-I...π, C-H...π)
- N: Nitrogen atoms (N-H...π)
- O: Oxygen atoms (O-H...π)
- S: Sulfur atoms (S-H...π)
"""

PI_INTERACTION_ATOMS: List[str] = ["H", "F", "CL", "BR", "I"]
"""List[str]: Elements that can participate in π-interactions.

Contains elements that can form π interactions when bonded to donor atoms:
- H: Hydrogen (C-H...π, N-H...π, O-H...π, S-H...π)
- F: Fluorine (C-F...π, rare but possible)
- CL: Chlorine (C-Cl...π halogen-π interactions)
- BR: Bromine (C-Br...π halogen-π interactions)
- I: Iodine (C-I...π halogen-π interactions)
"""

RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS: Dict[str, List[str]] = {
    # Protein residues
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HID": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HIE": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HIP": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "TYI": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYQ": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYB": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    # DNA nucleotides (purine bases)
    "DA": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],  # Adenine
    "DG": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],  # Guanine
    # DNA nucleotides (pyrimidine bases)
    "DC": ["N1", "C2", "N3", "C4", "C5", "C6"],  # Cytosine
    "DT": ["N1", "C2", "N3", "C4", "C5", "C6"],  # Thymine
    # RNA nucleotides (purine bases)
    "A": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],  # Adenine
    "G": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],  # Guanine
    # RNA nucleotides (pyrimidine bases)
    "C": ["N1", "C2", "N3", "C4", "C5", "C6"],  # Cytosine
    "U": ["N1", "C2", "N3", "C4", "C5", "C6"],  # Uracil
}
"""Dict[str, List[str]]: Mapping of aromatic residues to their ring atom names.

This dictionary provides the specific atom names that form aromatic ring systems
for each residue type containing aromatic groups:

Protein residues:
  Phenylalanine (PHE) and variants:
    - 6-membered benzene ring: CG-CD1-CE1-CZ-CE2-CD2
    
  Tyrosine (TYR, TYI, TYQ, TYB) and variants:
    - 6-membered phenolic ring: CG-CD1-CE1-CZ-CE2-CD2
    - TYI: Ionized tyrosine (deprotonated hydroxyl)
    - TYQ: Quinone form of tyrosine
    - TYB: Brominated tyrosine
    
  Tryptophan (TRP):
    - 5-membered pyrrole ring: CG-CD1-NE1-CE2-CD2
    - 6-membered benzene ring: CD2-CE2-CZ2-CH2-CZ3-CE3
    - Forms bicyclic indole system
    
  Histidine (HIS, HID, HIE, HIP):
    - 5-membered imidazole ring: CG-ND1-CE1-NE2-CD2
    - HID: Delta protonated (H on ND1)
    - HIE: Epsilon protonated (H on NE2)
    - HIP: Both nitrogens protonated (positive charge)

DNA nucleotides:
  Adenine (DA) and Guanine (DG) - Purine bases:
    - 5-membered ring: N9-C8-N7-C5-C4
    - 6-membered ring: C5-C6-N1-C2-N3-C4
    - Forms bicyclic purine system

  Cytosine (DC) and Thymine (DT) - Pyrimidine bases:
    - 6-membered ring: N1-C2-N3-C4-C5-C6

RNA nucleotides:
  Adenine (A) and Guanine (G) - Purine bases:
    - Same purine ring system as DNA counterparts
    
  Cytosine (C) and Uracil (U) - Pyrimidine bases:
    - Same pyrimidine ring system as DNA counterparts

Used for:
  - Calculating aromatic ring centroids for π interactions
  - Identifying atoms involved in π-π stacking
  - Determining ring plane orientations
  - X-H...π interaction analysis where these atoms form the π system
  - DNA/RNA-protein interface interactions
  - Nucleotide base stacking analysis
"""

HYDROPHOBIC_RESIDUES: List[str] = [
    "VAL",
    "LEU",
    "ILE",
    "MET",
    "PHE",
    "TRP",
    "PRO",
    "ALA",
]
"""List[str]: Hydrophobic amino acid residues with nonpolar side chains.

These amino acids have side chains that are predominantly nonpolar and hydrophobic:
    - VAL (Valine): Branched aliphatic chain
    - LEU (Leucine): Branched aliphatic chain
    - ILE (Isoleucine): Branched aliphatic chain
    - MET (Methionine): Sulfur-containing nonpolar chain
    - PHE (Phenylalanine): Aromatic benzyl group
    - TRP (Tryptophan): Aromatic indole group
    - PRO (Proline): Cyclic imino acid structure
    - ALA (Alanine): Simple methyl group

Used for:
    - Hydrophobic interaction analysis
    - Protein folding studies
    - Membrane protein analysis
    - Hydrophobic patch identification
"""

CHARGED_RESIDUES: List[str] = ["ARG", "LYS", "ASP", "GLU", "HIS"]
"""List[str]: Charged amino acid residues with ionizable side chains.

These amino acids carry formal charges at physiological pH:
    - ARG (Arginine): Positively charged guanidinium group (+1)
    - LYS (Lysine): Positively charged amino group (+1)
    - ASP (Aspartic acid): Negatively charged carboxylate group (-1)
    - GLU (Glutamic acid): Negatively charged carboxylate group (-1)
    - HIS (Histidine): Can be positively charged imidazolium group (pKa ~6)

Used for:
    - Electrostatic interaction analysis
    - Salt bridge identification
    - pH-dependent behavior studies
    - Ion binding site analysis
"""

RESIDUE_TYPES: List[str] = ["DNA", "RNA", "PROTEIN", "LIGAND", "WATER"]
"""List[str]: Standard residue type classifications for molecular analysis.

Classification categories for different types of molecular residues:
    - DNA: Deoxyribonucleotide residues (DA, DG, DC, DT, DI)
    - RNA: Ribonucleotide residues (A, G, C, U, I)
    - PROTEIN: Amino acid residues (20 standard amino acids and variants)
    - LIGAND: Ligands, cofactors, metals, and other heteroatom residues
    - WATER: Water molecules and solvent


Used for:
    - Residue type identification and classification
    - Molecular component analysis
    - Structure validation and processing
    - Interaction type determination
"""

# Single letter codes for atom properties
RESIDUE_TYPE_CODES: Dict[str, str] = {
    "DNA": "D",
    "RNA": "R",
    "PROTEIN": "P",
    "WATER": "W",
    "LIGAND": "L",
}
"""Dict[str, str]: Single letter codes for residue types.

Mapping of full residue type names to compact single letter codes:
    - "DNA" → "D": Deoxyribonucleotide residues
    - "RNA" → "R": Ribonucleotide residues  
    - "PROTEIN" → "P": Amino acid residues
    - "WATER" → "W": Water molecules and solvent
    - "LIGAND" → "L": Ligands, cofactors, metals, and other heteroatom residues

Used for compact representation in hydrogen bond descriptions and atom records.
"""

BACKBONE_SIDECHAIN_CODES: Dict[str, str] = {
    "BACKBONE": "B",
    "SIDECHAIN": "S",
    "NOT_APPLICABLE": "N",
}
"""Dict[str, str]: Single letter codes for backbone vs sidechain classification.

Mapping of atom structural classification to compact single letter codes:
    - "BACKBONE" → "B": Main chain atoms (protein backbone, DNA/RNA sugar-phosphate)
    - "SIDECHAIN" → "S": Side chain atoms (protein R-groups, nucleotide bases)

Used for describing hydrogen bond donor-acceptor relationships (e.g., S-S, S-B, B-B).
"""

AROMATIC_CODES: Dict[str, str] = {"AROMATIC": "A", "NON-AROMATIC": "N"}
"""Dict[str, str]: Single letter codes for aromatic classification.

Mapping of aromatic property classification to compact single letter codes:
    - "AROMATIC" → "A": Atoms that are part of aromatic ring systems
    - "NON-AROMATIC" → "N": Atoms that are not part of aromatic ring systems

Used for identifying atoms involved in π-interactions and aromatic stacking.
"""

# Carbonyl interaction constants
RESIDUES_WITH_BACKBONE_CARBONYLS: List[str] = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    # Common variants
    "HID",
    "HIE",
    "HIP",
    "CYX",
    "MSE",
]
"""List[str]: Residues that have backbone carbonyl groups (C=O).

All standard amino acid residues contain a backbone carbonyl group as part
of the peptide bond. This list includes:
- All 20 standard amino acids
- Common histidine protonation variants (HID, HIE, HIP)
- Cysteine disulfide variant (CYX)
- Selenomethionine (MSE)

Used for:
    - Backbone carbonyl identification
    - n→π* interaction detection between peptide bonds
    - Secondary structure analysis
"""

RESIDUES_WITH_SIDECHAIN_CARBONYLS: Dict[str, List[str]] = {
    "ASN": ["CG", "OD1"],  # Asparagine amide
    "GLN": ["CD", "OE1"],  # Glutamine amide
    "ASP": ["CG", "OD1"],  # Aspartate carboxylate (representative)
    "GLU": ["CD", "OE1"],  # Glutamate carboxylate (representative)
}
"""Dict[str, List[str]]: Residues with sidechain carbonyl groups and their atom names.

Maps residue names to [carbon_atom, oxygen_atom] pairs for sidechain carbonyls:
    - ASN: CG-OD1 (asparagine amide group)
    - GLN: CD-OE1 (glutamine amide group)
    - ASP: CG-OD1 (aspartate carboxylate, representative oxygen)
    - GLU: CD-OE1 (glutamate carboxylate, representative oxygen)

Note: Carboxylates (ASP, GLU) have two oxygen atoms; OD1/OE1 is used as 
representative for interaction detection.

Used for:
    - Sidechain carbonyl identification
    - n→π* interaction detection
    - Electrostatic interaction analysis
"""

BACKBONE_CARBONYL_ATOMS: Dict[str, str] = {"C": "O"}
"""Dict[str, str]: Backbone carbonyl carbon to oxygen atom mapping.

Standard peptide bond atoms:
    - C: Backbone carbonyl carbon
    - O: Backbone carbonyl oxygen

Used for:
    - Backbone carbonyl group identification
    - Peptide bond interaction analysis
"""

CARBONYL_BOND_LENGTH_RANGE: Dict[str, tuple] = {
    "amide": (1.15, 1.35),  # C=O in peptide bonds and amide sidechains
    "carboxylate": (1.15, 1.40),  # C-O in carboxylate groups (longer due to resonance)
}
"""Dict[str, tuple]: Acceptable C-O bond length ranges for different carbonyl types.

Bond length validation ranges in Angstroms:
    - "amide": (1.15, 1.35) Å - Peptide bonds and amide sidechains (ASN, GLN)
    - "carboxylate": (1.15, 1.40) Å - Carboxylate groups (ASP, GLU)

Carboxylate bonds are slightly longer due to resonance delocalization.

Used for:
    - Carbonyl group validation
    - Chemical bond verification
"""
