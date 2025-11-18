"""
Analysis parameters for HBAT molecular interaction analysis.

This module contains the AnalysisParameters class that defines all configurable
parameters for hydrogen bond, halogen bond, and π interaction analysis.
"""

from typing import List


class ParametersDefault:
    """Default values for molecular interaction analysis parameters."""

    # Hydrogen bond parameters
    HB_DISTANCE_CUTOFF = 2.5  # Å - H...A distance cutoff
    HB_ANGLE_CUTOFF = 120.0  # degrees - D-H...A angle cutoff
    HB_DA_DISTANCE = 3.5  # Å - Donor-acceptor distance cutoff

    # Weak hydrogen bond parameters (for carbon donors)
    WHB_DISTANCE_CUTOFF = 3.6  # Å - H...A distance cutoff
    WHB_ANGLE_CUTOFF = 150.0  # degrees - D-H...A angle cutoff
    WHB_DA_DISTANCE = 3.5  # Å - Donor-acceptor distance cutoff

    # Halogen bond parameters (updated defaults)
    XB_DISTANCE_CUTOFF = (
        3.9  # Å - X...A distance cutoff <= vdW sum (Max vdW sum: I with SE → 3.88 Å)
    )
    XB_ANGLE_CUTOFF = (
        150.0  # degrees - C-X...A angle cutoff (updated from 120° to 150°)
    )

    # π interaction parameters
    PI_DISTANCE_CUTOFF = (
        3.5  # Å - H...π distance cutoff (legacy, kept for compatibility)
    )
    PI_ANGLE_CUTOFF = (
        110.0  # degrees - D-H...π angle cutoff (legacy, kept for compatibility)
    )

    # π interaction subtype parameters
    # Halogen-π interactions
    PI_CCL_DISTANCE_CUTOFF = 3.5  # Å - C-Cl...π distance cutoff
    PI_CCL_ANGLE_CUTOFF = 145  # degrees - C-Cl...π angle cutoff
    PI_CBR_DISTANCE_CUTOFF = 3.5  # Å - C-Br...π distance cutoff
    PI_CBR_ANGLE_CUTOFF = 155  # degrees - C-Br...π angle cutoff
    PI_CI_DISTANCE_CUTOFF = 3.6  # Å - C-I...π distance cutoff
    PI_CI_ANGLE_CUTOFF = 165.0  # degrees - C-I...π angle cutoff

    # Hydrogen-π interactions
    PI_CH_DISTANCE_CUTOFF = 3.5  # Å - C-H...π distance cutoff
    PI_CH_ANGLE_CUTOFF = 110.0  # degrees - C-H...π angle cutoff
    PI_NH_DISTANCE_CUTOFF = 3.2  # Å - N-H...π distance cutoff
    PI_NH_ANGLE_CUTOFF = 115.0  # degrees - N-H...π angle cutoff
    PI_OH_DISTANCE_CUTOFF = 3.0  # Å - O-H...π distance cutoff
    PI_OH_ANGLE_CUTOFF = 115.0  # degrees - O-H...π angle cutoff
    PI_SH_DISTANCE_CUTOFF = 3.8  # Å - S-H...π distance cutoff
    PI_SH_ANGLE_CUTOFF = 105.0  # degrees - S-H...π angle cutoff

    # π-π stacking interaction parameters
    PI_PI_DISTANCE_CUTOFF = 3.8  # Å - π-π centroid-to-centroid distance cutoff
    PI_PI_PARALLEL_ANGLE_CUTOFF = 30.0  # degrees - max angle for parallel stacking
    PI_PI_TSHAPED_ANGLE_MIN = 60.0  # degrees - min angle for T-shaped stacking
    PI_PI_TSHAPED_ANGLE_MAX = 90.0  # degrees - max angle for T-shaped stacking
    PI_PI_OFFSET_CUTOFF = 2.0  # Å - max lateral offset for parallel stacking

    # Carbonyl-carbonyl interaction parameters (n→π*)
    CARBONYL_DISTANCE_CUTOFF = (
        3.2  # Å - O···C distance cutoff (Bürgi-Dunitz trajectory)
    )
    CARBONYL_ANGLE_MIN = 95.0  # degrees - min O···C=O angle (Bürgi-Dunitz)
    CARBONYL_ANGLE_MAX = 125.0  # degrees - max O···C=O angle (Bürgi-Dunitz)

    # n→π* interaction parameters (lone pair to π system)
    N_PI_DISTANCE_CUTOFF = 3.6  # Å - lone pair to π center distance cutoff
    N_PI_DISTANCE_MIN = 2.5  # Å - minimum distance to avoid unrealistic close contacts
    N_PI_SULFUR_DISTANCE_CUTOFF = (
        4.0  # Å - sulfur-specific distance cutoff (larger vdW radius)
    )
    N_PI_ANGLE_MIN = 0.0  # degrees - min angle to π plane (perpendicular approach)
    N_PI_ANGLE_MAX = 45.0  # degrees - max angle to π plane (optimal geometry)

    # General analysis parameters
    COVALENT_CUTOFF_FACTOR = 0.85  # Covalent bond detection factor (0.0-1.0)
    ANALYSIS_MODE = "complete"  # Analysis mode: "complete" or "local"

    # Bond distance thresholds
    MAX_BOND_DISTANCE = 2.5  # Reasonable maximum for most covalent bonds (Angstroms)
    MIN_BOND_DISTANCE = 0.5  # Minimum realistic bond distance (Angstroms)

    # PDB structure fixing parameters (updated defaults)
    FIX_PDB_ENABLED = True  # Enable PDB structure fixing (changed from False to True)
    FIX_PDB_METHOD = "pdbfixer"  # Method: "openbabel" or "pdbfixer" (changed from "openbabel" to "pdbfixer")

    # Fixing operations (explicit control)
    FIX_PDB_ADD_HYDROGENS = (
        True  # Add missing hydrogen atoms (both OpenBabel and PDBFixer)
    )
    FIX_PDB_ADD_HEAVY_ATOMS = False  # Add missing heavy atoms (PDBFixer only)
    FIX_PDB_REPLACE_NONSTANDARD = False  # Replace nonstandard residues (PDBFixer only)
    FIX_PDB_REMOVE_HETEROGENS = False  # Remove heterogens (PDBFixer only)
    FIX_PDB_KEEP_WATER = True  # Keep water when removing heterogens (PDBFixer only)


class AnalysisParameters:
    """Parameters for comprehensive molecular interaction analysis.

    This class contains all configurable parameters for detecting and analyzing
    molecular interactions in protein structures. Supports multiple interaction
    types with subtype-specific parameters.

    **Hydrogen Bonds (Classical):**

    - Strong N-H···O, O-H···O, N-H···N interactions
    - Default: 2.5Å H···A distance, 120° D-H···A angle

    **Weak Hydrogen Bonds (C-H···O):**

    - Carbon-hydrogen donor interactions with oxygen acceptors
    - Default: 3.6Å H···A distance, 150° D-H···A angle
    - Important in protein-ligand binding and aromatic interactions

    **Halogen Bonds:**
    
    - C-X···A interactions where X is Cl, Br, I
    - Default: 3.9Å X···A distance, 150° C-X···A angle (updated default)

    **π Interactions (Multiple Subtypes):**

    - Hydrogen-π: C-H···π, N-H···π, O-H···π, S-H···π
    - Halogen-π: C-Cl···π, C-Br···π, C-I···π
    - Each subtype has optimized distance/angle cutoffs

    **π-π Stacking Interactions:**

    - Aromatic ring-ring interactions (parallel, T-shaped, offset)
    - Default: 3.8Å centroid distance, 30° parallel angle cutoff
    - Important for protein folding and stability

    **Carbonyl Interactions (n→π*):**

    - Carbonyl-carbonyl n→π* interactions (Bürgi-Dunitz trajectory)
    - Default: 3.2Å O···C distance, 95°-125° O···C=O angle
    - Found in polyproline II helices and protein backbones

    **n→π* Interactions:**

    - Lone pair to π system interactions
    - Default: 3.6Å distance, 0°-45° angle to π plane
    - Includes O, N, S lone pairs interacting with aromatic systems

    :param hb_distance_cutoff: Maximum H...A distance for hydrogen bonds (Å)
    :type hb_distance_cutoff: float
    :param hb_angle_cutoff: Minimum D-H...A angle for hydrogen bonds (degrees)
    :type hb_angle_cutoff: float
    :param hb_donor_acceptor_cutoff: Maximum D...A distance for hydrogen bonds (Å)
    :type hb_donor_acceptor_cutoff: float
    :param whb_distance_cutoff: Maximum H...A distance for weak hydrogen bonds (Å)
    :type whb_distance_cutoff: float
    :param whb_angle_cutoff: Minimum D-H...A angle for weak hydrogen bonds (degrees)
    :type whb_angle_cutoff: float
    :param whb_donor_acceptor_cutoff: Maximum D...A distance for weak hydrogen bonds (Å)
    :type whb_donor_acceptor_cutoff: float
    :param xb_distance_cutoff: Maximum X...A distance for halogen bonds (Å)
    :type xb_distance_cutoff: float
    :param xb_angle_cutoff: Minimum C-X...A angle for halogen bonds (degrees, default: 150°)
    :type xb_angle_cutoff: float
    :param pi_distance_cutoff: Maximum H...π distance for π interactions (Å, legacy)
    :type pi_distance_cutoff: float
    :param pi_angle_cutoff: Minimum D-H...π angle for π interactions (degrees, legacy)
    :type pi_angle_cutoff: float
    :param pi_ccl_distance_cutoff: Maximum C-Cl...π distance for halogen-π interactions (Å)
    :type pi_ccl_distance_cutoff: float
    :param pi_ccl_angle_cutoff: Minimum C-Cl...π angle for halogen-π interactions (degrees)
    :type pi_ccl_angle_cutoff: float
    :param pi_cbr_distance_cutoff: Maximum C-Br...π distance for halogen-π interactions (Å)
    :type pi_cbr_distance_cutoff: float
    :param pi_cbr_angle_cutoff: Minimum C-Br...π angle for halogen-π interactions (degrees)
    :type pi_cbr_angle_cutoff: float
    :param pi_ci_distance_cutoff: Maximum C-I...π distance for halogen-π interactions (Å)
    :type pi_ci_distance_cutoff: float
    :param pi_ci_angle_cutoff: Minimum C-I...π angle for halogen-π interactions (degrees)
    :type pi_ci_angle_cutoff: float
    :param pi_ch_distance_cutoff: Maximum C-H...π distance for hydrogen-π interactions (Å)
    :type pi_ch_distance_cutoff: float
    :param pi_ch_angle_cutoff: Minimum C-H...π angle for hydrogen-π interactions (degrees)
    :type pi_ch_angle_cutoff: float
    :param pi_nh_distance_cutoff: Maximum N-H...π distance for hydrogen-π interactions (Å)
    :type pi_nh_distance_cutoff: float
    :param pi_nh_angle_cutoff: Minimum N-H...π angle for hydrogen-π interactions (degrees)
    :type pi_nh_angle_cutoff: float
    :param pi_oh_distance_cutoff: Maximum O-H...π distance for hydrogen-π interactions (Å)
    :type pi_oh_distance_cutoff: float
    :param pi_oh_angle_cutoff: Minimum O-H...π angle for hydrogen-π interactions (degrees)
    :type pi_oh_angle_cutoff: float
    :param pi_sh_distance_cutoff: Maximum S-H...π distance for hydrogen-π interactions (Å)
    :type pi_sh_distance_cutoff: float
    :param pi_sh_angle_cutoff: Minimum S-H...π angle for hydrogen-π interactions (degrees)
    :type pi_sh_angle_cutoff: float
    :param pi_pi_distance_cutoff: Maximum centroid-to-centroid distance for π-π stacking (Å)
    :type pi_pi_distance_cutoff: float
    :param pi_pi_parallel_angle_cutoff: Maximum angle for parallel π-π stacking (degrees)
    :type pi_pi_parallel_angle_cutoff: float
    :param pi_pi_tshaped_angle_min: Minimum angle for T-shaped π-π stacking (degrees)
    :type pi_pi_tshaped_angle_min: float
    :param pi_pi_tshaped_angle_max: Maximum angle for T-shaped π-π stacking (degrees)
    :type pi_pi_tshaped_angle_max: float
    :param pi_pi_offset_cutoff: Maximum lateral offset for parallel π-π stacking (Å)
    :type pi_pi_offset_cutoff: float
    :param carbonyl_distance_cutoff: Maximum O···C distance for carbonyl n→π* interactions (Å)
    :type carbonyl_distance_cutoff: float
    :param carbonyl_angle_min: Minimum O···C=O angle for carbonyl interactions (degrees, Bürgi-Dunitz)
    :type carbonyl_angle_min: float
    :param carbonyl_angle_max: Maximum O···C=O angle for carbonyl interactions (degrees, Bürgi-Dunitz)
    :type carbonyl_angle_max: float
    :param n_pi_distance_cutoff: Maximum lone pair to π center distance for n→π* interactions (Å)
    :type n_pi_distance_cutoff: float
    :param n_pi_sulfur_distance_cutoff: Maximum sulfur lone pair to π center distance for n→π* interactions (Å)
    :type n_pi_sulfur_distance_cutoff: float
    :param n_pi_angle_min: Minimum angle to π plane for n→π* interactions (degrees)
    :type n_pi_angle_min: float
    :param n_pi_angle_max: Maximum angle to π plane for n→π* interactions (degrees)
    :type n_pi_angle_max: float
    :param covalent_cutoff_factor: Factor for covalent bond detection
    :type covalent_cutoff_factor: float
    :param analysis_mode: Analysis mode ('local' or 'global')
    :type analysis_mode: str
    :param fix_pdb_enabled: Enable PDB structure fixing
    :type fix_pdb_enabled: bool
    :param fix_pdb_method: PDB fixing method ('openbabel' or 'pdbfixer')
    :type fix_pdb_method: str
    :param fix_pdb_add_hydrogens: Add missing hydrogen atoms (both methods)
    :type fix_pdb_add_hydrogens: bool
    :param fix_pdb_add_heavy_atoms: Add missing heavy atoms (PDBFixer only)
    :type fix_pdb_add_heavy_atoms: bool
    :param fix_pdb_replace_nonstandard: Replace nonstandard residues (PDBFixer only)
    :type fix_pdb_replace_nonstandard: bool
    :param fix_pdb_remove_heterogens: Remove heterogens (PDBFixer only)
    :type fix_pdb_remove_heterogens: bool
    :param fix_pdb_keep_water: Keep water when removing heterogens (PDBFixer only)
    :type fix_pdb_keep_water: bool
    """

    def __init__(
        self,
        # Hydrogen bond parameters
        hb_distance_cutoff: float = ParametersDefault.HB_DISTANCE_CUTOFF,
        hb_angle_cutoff: float = ParametersDefault.HB_ANGLE_CUTOFF,
        hb_donor_acceptor_cutoff: float = ParametersDefault.HB_DA_DISTANCE,
        # Weak hydrogen bond parameters (for carbon donors)
        whb_distance_cutoff: float = ParametersDefault.WHB_DISTANCE_CUTOFF,
        whb_angle_cutoff: float = ParametersDefault.WHB_ANGLE_CUTOFF,
        whb_donor_acceptor_cutoff: float = ParametersDefault.WHB_DA_DISTANCE,
        # Halogen bond parameters
        xb_distance_cutoff: float = ParametersDefault.XB_DISTANCE_CUTOFF,
        xb_angle_cutoff: float = ParametersDefault.XB_ANGLE_CUTOFF,
        # Pi interaction parameters (legacy)
        pi_distance_cutoff: float = ParametersDefault.PI_DISTANCE_CUTOFF,
        pi_angle_cutoff: float = ParametersDefault.PI_ANGLE_CUTOFF,
        # Pi interaction subtype parameters
        pi_ccl_distance_cutoff: float = ParametersDefault.PI_CCL_DISTANCE_CUTOFF,
        pi_ccl_angle_cutoff: float = ParametersDefault.PI_CCL_ANGLE_CUTOFF,
        pi_cbr_distance_cutoff: float = ParametersDefault.PI_CBR_DISTANCE_CUTOFF,
        pi_cbr_angle_cutoff: float = ParametersDefault.PI_CBR_ANGLE_CUTOFF,
        pi_ci_distance_cutoff: float = ParametersDefault.PI_CI_DISTANCE_CUTOFF,
        pi_ci_angle_cutoff: float = ParametersDefault.PI_CI_ANGLE_CUTOFF,
        pi_ch_distance_cutoff: float = ParametersDefault.PI_CH_DISTANCE_CUTOFF,
        pi_ch_angle_cutoff: float = ParametersDefault.PI_CH_ANGLE_CUTOFF,
        pi_nh_distance_cutoff: float = ParametersDefault.PI_NH_DISTANCE_CUTOFF,
        pi_nh_angle_cutoff: float = ParametersDefault.PI_NH_ANGLE_CUTOFF,
        pi_oh_distance_cutoff: float = ParametersDefault.PI_OH_DISTANCE_CUTOFF,
        pi_oh_angle_cutoff: float = ParametersDefault.PI_OH_ANGLE_CUTOFF,
        pi_sh_distance_cutoff: float = ParametersDefault.PI_SH_DISTANCE_CUTOFF,
        pi_sh_angle_cutoff: float = ParametersDefault.PI_SH_ANGLE_CUTOFF,
        # π-π stacking parameters
        pi_pi_distance_cutoff: float = ParametersDefault.PI_PI_DISTANCE_CUTOFF,
        pi_pi_parallel_angle_cutoff: float = ParametersDefault.PI_PI_PARALLEL_ANGLE_CUTOFF,
        pi_pi_tshaped_angle_min: float = ParametersDefault.PI_PI_TSHAPED_ANGLE_MIN,
        pi_pi_tshaped_angle_max: float = ParametersDefault.PI_PI_TSHAPED_ANGLE_MAX,
        pi_pi_offset_cutoff: float = ParametersDefault.PI_PI_OFFSET_CUTOFF,
        # Carbonyl-carbonyl interaction parameters
        carbonyl_distance_cutoff: float = ParametersDefault.CARBONYL_DISTANCE_CUTOFF,
        carbonyl_angle_min: float = ParametersDefault.CARBONYL_ANGLE_MIN,
        carbonyl_angle_max: float = ParametersDefault.CARBONYL_ANGLE_MAX,
        # n→π* interaction parameters
        n_pi_distance_cutoff: float = ParametersDefault.N_PI_DISTANCE_CUTOFF,
        n_pi_sulfur_distance_cutoff: float = ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF,
        n_pi_angle_min: float = ParametersDefault.N_PI_ANGLE_MIN,
        n_pi_angle_max: float = ParametersDefault.N_PI_ANGLE_MAX,
        # General parameters
        covalent_cutoff_factor: float = ParametersDefault.COVALENT_CUTOFF_FACTOR,
        analysis_mode: str = ParametersDefault.ANALYSIS_MODE,
        # PDB structure fixing parameters
        fix_pdb_enabled: bool = ParametersDefault.FIX_PDB_ENABLED,
        fix_pdb_method: str = ParametersDefault.FIX_PDB_METHOD,
        fix_pdb_add_hydrogens: bool = ParametersDefault.FIX_PDB_ADD_HYDROGENS,
        fix_pdb_add_heavy_atoms: bool = ParametersDefault.FIX_PDB_ADD_HEAVY_ATOMS,
        fix_pdb_replace_nonstandard: bool = ParametersDefault.FIX_PDB_REPLACE_NONSTANDARD,
        fix_pdb_remove_heterogens: bool = ParametersDefault.FIX_PDB_REMOVE_HETEROGENS,
        fix_pdb_keep_water: bool = ParametersDefault.FIX_PDB_KEEP_WATER,
        **kwargs,
    ) -> None:
        """Initialize analysis parameters.

        :param hb_distance_cutoff: Maximum H...A distance for hydrogen bonds (Å)
        :type hb_distance_cutoff: float
        :param hb_angle_cutoff: Minimum D-H...A angle for hydrogen bonds (degrees)
        :type hb_angle_cutoff: float
        :param hb_donor_acceptor_cutoff: Maximum D...A distance for hydrogen bonds (Å)
        :type hb_donor_acceptor_cutoff: float
        :param whb_distance_cutoff: Maximum H...A distance for weak hydrogen bonds (Å)
        :type whb_distance_cutoff: float
        :param whb_angle_cutoff: Minimum D-H...A angle for weak hydrogen bonds (degrees)
        :type whb_angle_cutoff: float
        :param whb_donor_acceptor_cutoff: Maximum D...A distance for weak hydrogen bonds (Å)
        :type whb_donor_acceptor_cutoff: float
        :param xb_distance_cutoff: Maximum X...A distance for halogen bonds (Å)
        :type xb_distance_cutoff: float
        :param xb_angle_cutoff: Minimum C-X...A angle for halogen bonds (degrees)
        :type xb_angle_cutoff: float
        :param pi_distance_cutoff: Maximum H...π distance for π interactions (Å)
        :type pi_distance_cutoff: float
        :param pi_angle_cutoff: Minimum D-H...π angle for π interactions (degrees)
        :type pi_angle_cutoff: float
        :param covalent_cutoff_factor: Factor for covalent bond detection
        :type covalent_cutoff_factor: float
        :param analysis_mode: Analysis mode ('local' or 'global')
        :type analysis_mode: str
        :param fix_pdb_enabled: Enable PDB structure fixing
        :type fix_pdb_enabled: bool
        :param fix_pdb_method: PDB fixing method ('openbabel' or 'pdbfixer')
        :type fix_pdb_method: str
        :param fix_pdb_add_hydrogens: Add missing hydrogen atoms (both methods)
        :type fix_pdb_add_hydrogens: bool
        :param fix_pdb_add_heavy_atoms: Add missing heavy atoms (PDBFixer only)
        :type fix_pdb_add_heavy_atoms: bool
        :param fix_pdb_replace_nonstandard: Replace nonstandard residues (PDBFixer only)
        :type fix_pdb_replace_nonstandard: bool
        :param fix_pdb_remove_heterogens: Remove heterogens (PDBFixer only)
        :type fix_pdb_remove_heterogens: bool
        :param fix_pdb_keep_water: Keep water when removing heterogens (PDBFixer only)
        :type fix_pdb_keep_water: bool
        :param kwargs: Additional parameters (for future extensibility)
        :type kwargs: dict
        """
        # Hydrogen bond parameters
        self.hb_distance_cutoff = hb_distance_cutoff
        self.hb_angle_cutoff = hb_angle_cutoff
        self.hb_donor_acceptor_cutoff = hb_donor_acceptor_cutoff

        # Weak hydrogen bond parameters (for carbon donors)
        self.whb_distance_cutoff = whb_distance_cutoff
        self.whb_angle_cutoff = whb_angle_cutoff
        self.whb_donor_acceptor_cutoff = whb_donor_acceptor_cutoff

        # Halogen bond parameters
        self.xb_distance_cutoff = xb_distance_cutoff
        self.xb_angle_cutoff = xb_angle_cutoff

        # π interaction parameters (legacy)
        self.pi_distance_cutoff = pi_distance_cutoff
        self.pi_angle_cutoff = pi_angle_cutoff

        # π interaction subtype parameters
        self.pi_ccl_distance_cutoff = pi_ccl_distance_cutoff
        self.pi_ccl_angle_cutoff = pi_ccl_angle_cutoff
        self.pi_cbr_distance_cutoff = pi_cbr_distance_cutoff
        self.pi_cbr_angle_cutoff = pi_cbr_angle_cutoff
        self.pi_ci_distance_cutoff = pi_ci_distance_cutoff
        self.pi_ci_angle_cutoff = pi_ci_angle_cutoff
        self.pi_ch_distance_cutoff = pi_ch_distance_cutoff
        self.pi_ch_angle_cutoff = pi_ch_angle_cutoff
        self.pi_nh_distance_cutoff = pi_nh_distance_cutoff
        self.pi_nh_angle_cutoff = pi_nh_angle_cutoff
        self.pi_oh_distance_cutoff = pi_oh_distance_cutoff
        self.pi_oh_angle_cutoff = pi_oh_angle_cutoff
        self.pi_sh_distance_cutoff = pi_sh_distance_cutoff
        self.pi_sh_angle_cutoff = pi_sh_angle_cutoff

        # π-π stacking parameters
        self.pi_pi_distance_cutoff = pi_pi_distance_cutoff
        self.pi_pi_parallel_angle_cutoff = pi_pi_parallel_angle_cutoff
        self.pi_pi_tshaped_angle_min = pi_pi_tshaped_angle_min
        self.pi_pi_tshaped_angle_max = pi_pi_tshaped_angle_max
        self.pi_pi_offset_cutoff = pi_pi_offset_cutoff

        # Carbonyl-carbonyl interaction parameters
        self.carbonyl_distance_cutoff = carbonyl_distance_cutoff
        self.carbonyl_angle_min = carbonyl_angle_min
        self.carbonyl_angle_max = carbonyl_angle_max

        # n→π* interaction parameters
        self.n_pi_distance_cutoff = n_pi_distance_cutoff
        self.n_pi_sulfur_distance_cutoff = n_pi_sulfur_distance_cutoff
        self.n_pi_angle_min = n_pi_angle_min
        self.n_pi_angle_max = n_pi_angle_max

        # General parameters
        self.covalent_cutoff_factor = covalent_cutoff_factor
        self.analysis_mode = analysis_mode

        # PDB structure fixing parameters
        self.fix_pdb_enabled = fix_pdb_enabled
        self.fix_pdb_method = fix_pdb_method
        self.fix_pdb_add_hydrogens = fix_pdb_add_hydrogens
        self.fix_pdb_add_heavy_atoms = fix_pdb_add_heavy_atoms
        self.fix_pdb_replace_nonstandard = fix_pdb_replace_nonstandard
        self.fix_pdb_remove_heterogens = fix_pdb_remove_heterogens
        self.fix_pdb_keep_water = fix_pdb_keep_water

        # Store any additional parameters
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self) -> str:
        """Return string representation of the parameters object.

        :returns: String representation showing all parameters
        :rtype: str
        """
        params = []
        # Hydrogen bond parameters
        params.append(f"hb_distance_cutoff={self.hb_distance_cutoff}")
        params.append(f"hb_angle_cutoff={self.hb_angle_cutoff}")
        params.append(f"hb_donor_acceptor_cutoff={self.hb_donor_acceptor_cutoff}")
        # Weak hydrogen bond parameters
        params.append(f"whb_distance_cutoff={self.whb_distance_cutoff}")
        params.append(f"whb_angle_cutoff={self.whb_angle_cutoff}")
        params.append(f"whb_donor_acceptor_cutoff={self.whb_donor_acceptor_cutoff}")
        # Halogen bond parameters
        params.append(f"xb_distance_cutoff={self.xb_distance_cutoff}")
        params.append(f"xb_angle_cutoff={self.xb_angle_cutoff}")
        # π interaction parameters (legacy)
        params.append(f"pi_distance_cutoff={self.pi_distance_cutoff}")
        params.append(f"pi_angle_cutoff={self.pi_angle_cutoff}")
        # π interaction subtype parameters
        params.append(f"pi_ccl_distance_cutoff={self.pi_ccl_distance_cutoff}")
        params.append(f"pi_ccl_angle_cutoff={self.pi_ccl_angle_cutoff}")
        params.append(f"pi_cbr_distance_cutoff={self.pi_cbr_distance_cutoff}")
        params.append(f"pi_cbr_angle_cutoff={self.pi_cbr_angle_cutoff}")
        params.append(f"pi_ci_distance_cutoff={self.pi_ci_distance_cutoff}")
        params.append(f"pi_ci_angle_cutoff={self.pi_ci_angle_cutoff}")
        params.append(f"pi_ch_distance_cutoff={self.pi_ch_distance_cutoff}")
        params.append(f"pi_ch_angle_cutoff={self.pi_ch_angle_cutoff}")
        params.append(f"pi_nh_distance_cutoff={self.pi_nh_distance_cutoff}")
        params.append(f"pi_nh_angle_cutoff={self.pi_nh_angle_cutoff}")
        params.append(f"pi_oh_distance_cutoff={self.pi_oh_distance_cutoff}")
        params.append(f"pi_oh_angle_cutoff={self.pi_oh_angle_cutoff}")
        params.append(f"pi_sh_distance_cutoff={self.pi_sh_distance_cutoff}")
        params.append(f"pi_sh_angle_cutoff={self.pi_sh_angle_cutoff}")
        # π-π stacking parameters
        params.append(f"pi_pi_distance_cutoff={self.pi_pi_distance_cutoff}")
        params.append(f"pi_pi_parallel_angle_cutoff={self.pi_pi_parallel_angle_cutoff}")
        params.append(f"pi_pi_tshaped_angle_min={self.pi_pi_tshaped_angle_min}")
        params.append(f"pi_pi_tshaped_angle_max={self.pi_pi_tshaped_angle_max}")
        params.append(f"pi_pi_offset_cutoff={self.pi_pi_offset_cutoff}")
        # Carbonyl interaction parameters
        params.append(f"carbonyl_distance_cutoff={self.carbonyl_distance_cutoff}")
        params.append(f"carbonyl_angle_min={self.carbonyl_angle_min}")
        params.append(f"carbonyl_angle_max={self.carbonyl_angle_max}")
        # n→π* interaction parameters
        params.append(f"n_pi_distance_cutoff={self.n_pi_distance_cutoff}")
        params.append(f"n_pi_sulfur_distance_cutoff={self.n_pi_sulfur_distance_cutoff}")
        params.append(f"n_pi_angle_min={self.n_pi_angle_min}")
        params.append(f"n_pi_angle_max={self.n_pi_angle_max}")
        # General parameters
        params.append(f"covalent_cutoff_factor={self.covalent_cutoff_factor}")
        params.append(f"analysis_mode='{self.analysis_mode}'")
        # PDB fixing parameters
        params.append(f"fix_pdb_enabled={self.fix_pdb_enabled}")
        params.append(f"fix_pdb_method='{self.fix_pdb_method}'")
        params.append(f"fix_pdb_add_hydrogens={self.fix_pdb_add_hydrogens}")
        params.append(f"fix_pdb_add_heavy_atoms={self.fix_pdb_add_heavy_atoms}")
        params.append(f"fix_pdb_replace_nonstandard={self.fix_pdb_replace_nonstandard}")
        params.append(f"fix_pdb_remove_heterogens={self.fix_pdb_remove_heterogens}")
        params.append(f"fix_pdb_keep_water={self.fix_pdb_keep_water}")

        return f"AnalysisParameters({', '.join(params)})"

    def __eq__(self, other: object) -> bool:
        """Compare two AnalysisParameters objects for equality.

        :param other: Other AnalysisParameters object to compare
        :type other: AnalysisParameters
        :returns: True if all parameters are equal
        :rtype: bool
        """
        if not isinstance(other, AnalysisParameters):
            return False

        return (
            self.hb_distance_cutoff == other.hb_distance_cutoff
            and self.hb_angle_cutoff == other.hb_angle_cutoff
            and self.hb_donor_acceptor_cutoff == other.hb_donor_acceptor_cutoff
            and self.whb_distance_cutoff == other.whb_distance_cutoff
            and self.whb_angle_cutoff == other.whb_angle_cutoff
            and self.whb_donor_acceptor_cutoff == other.whb_donor_acceptor_cutoff
            and self.xb_distance_cutoff == other.xb_distance_cutoff
            and self.xb_angle_cutoff == other.xb_angle_cutoff
            and self.pi_distance_cutoff == other.pi_distance_cutoff
            and self.pi_angle_cutoff == other.pi_angle_cutoff
            and self.pi_ccl_distance_cutoff == other.pi_ccl_distance_cutoff
            and self.pi_ccl_angle_cutoff == other.pi_ccl_angle_cutoff
            and self.pi_cbr_distance_cutoff == other.pi_cbr_distance_cutoff
            and self.pi_cbr_angle_cutoff == other.pi_cbr_angle_cutoff
            and self.pi_ci_distance_cutoff == other.pi_ci_distance_cutoff
            and self.pi_ci_angle_cutoff == other.pi_ci_angle_cutoff
            and self.pi_ch_distance_cutoff == other.pi_ch_distance_cutoff
            and self.pi_ch_angle_cutoff == other.pi_ch_angle_cutoff
            and self.pi_nh_distance_cutoff == other.pi_nh_distance_cutoff
            and self.pi_nh_angle_cutoff == other.pi_nh_angle_cutoff
            and self.pi_oh_distance_cutoff == other.pi_oh_distance_cutoff
            and self.pi_oh_angle_cutoff == other.pi_oh_angle_cutoff
            and self.pi_sh_distance_cutoff == other.pi_sh_distance_cutoff
            and self.pi_sh_angle_cutoff == other.pi_sh_angle_cutoff
            and self.pi_pi_distance_cutoff == other.pi_pi_distance_cutoff
            and self.pi_pi_parallel_angle_cutoff == other.pi_pi_parallel_angle_cutoff
            and self.pi_pi_tshaped_angle_min == other.pi_pi_tshaped_angle_min
            and self.pi_pi_tshaped_angle_max == other.pi_pi_tshaped_angle_max
            and self.pi_pi_offset_cutoff == other.pi_pi_offset_cutoff
            and self.carbonyl_distance_cutoff == other.carbonyl_distance_cutoff
            and self.carbonyl_angle_min == other.carbonyl_angle_min
            and self.carbonyl_angle_max == other.carbonyl_angle_max
            and self.n_pi_distance_cutoff == other.n_pi_distance_cutoff
            and self.n_pi_sulfur_distance_cutoff == other.n_pi_sulfur_distance_cutoff
            and self.n_pi_angle_min == other.n_pi_angle_min
            and self.n_pi_angle_max == other.n_pi_angle_max
            and self.covalent_cutoff_factor == other.covalent_cutoff_factor
            and self.analysis_mode == other.analysis_mode
            and self.fix_pdb_enabled == other.fix_pdb_enabled
            and self.fix_pdb_method == other.fix_pdb_method
            and self.fix_pdb_add_hydrogens == other.fix_pdb_add_hydrogens
            and self.fix_pdb_add_heavy_atoms == other.fix_pdb_add_heavy_atoms
            and self.fix_pdb_replace_nonstandard == other.fix_pdb_replace_nonstandard
            and self.fix_pdb_remove_heterogens == other.fix_pdb_remove_heterogens
            and self.fix_pdb_keep_water == other.fix_pdb_keep_water
        )

    def __hash__(self) -> int:
        """Return hash of the parameters object for use in sets/dicts.

        :returns: Hash value based on all parameters
        :rtype: int
        """
        return hash(
            (
                self.hb_distance_cutoff,
                self.hb_angle_cutoff,
                self.hb_donor_acceptor_cutoff,
                self.whb_distance_cutoff,
                self.whb_angle_cutoff,
                self.whb_donor_acceptor_cutoff,
                self.xb_distance_cutoff,
                self.xb_angle_cutoff,
                self.pi_distance_cutoff,
                self.pi_angle_cutoff,
                self.pi_ccl_distance_cutoff,
                self.pi_ccl_angle_cutoff,
                self.pi_cbr_distance_cutoff,
                self.pi_cbr_angle_cutoff,
                self.pi_ci_distance_cutoff,
                self.pi_ci_angle_cutoff,
                self.pi_ch_distance_cutoff,
                self.pi_ch_angle_cutoff,
                self.pi_nh_distance_cutoff,
                self.pi_nh_angle_cutoff,
                self.pi_oh_distance_cutoff,
                self.pi_oh_angle_cutoff,
                self.pi_sh_distance_cutoff,
                self.pi_sh_angle_cutoff,
                self.covalent_cutoff_factor,
                self.analysis_mode,
                self.fix_pdb_enabled,
                self.fix_pdb_method,
                self.fix_pdb_add_hydrogens,
                self.fix_pdb_add_heavy_atoms,
                self.fix_pdb_replace_nonstandard,
                self.fix_pdb_remove_heterogens,
                self.fix_pdb_keep_water,
            )
        )

    def to_dict(self) -> dict:
        """Convert parameters to dictionary format.

        :returns: Dictionary representation of all parameters
        :rtype: dict
        """
        return {
            # Hydrogen bond parameters
            "hb_distance_cutoff": self.hb_distance_cutoff,
            "hb_angle_cutoff": self.hb_angle_cutoff,
            "hb_donor_acceptor_cutoff": self.hb_donor_acceptor_cutoff,
            # Weak hydrogen bond parameters
            "whb_distance_cutoff": self.whb_distance_cutoff,
            "whb_angle_cutoff": self.whb_angle_cutoff,
            "whb_donor_acceptor_cutoff": self.whb_donor_acceptor_cutoff,
            # Halogen bond parameters
            "xb_distance_cutoff": self.xb_distance_cutoff,
            "xb_angle_cutoff": self.xb_angle_cutoff,
            # π interaction parameters
            "pi_distance_cutoff": self.pi_distance_cutoff,
            "pi_angle_cutoff": self.pi_angle_cutoff,
            # π interaction subtype parameters
            "pi_ccl_distance_cutoff": self.pi_ccl_distance_cutoff,
            "pi_ccl_angle_cutoff": self.pi_ccl_angle_cutoff,
            "pi_cbr_distance_cutoff": self.pi_cbr_distance_cutoff,
            "pi_cbr_angle_cutoff": self.pi_cbr_angle_cutoff,
            "pi_ci_distance_cutoff": self.pi_ci_distance_cutoff,
            "pi_ci_angle_cutoff": self.pi_ci_angle_cutoff,
            "pi_ch_distance_cutoff": self.pi_ch_distance_cutoff,
            "pi_ch_angle_cutoff": self.pi_ch_angle_cutoff,
            "pi_nh_distance_cutoff": self.pi_nh_distance_cutoff,
            "pi_nh_angle_cutoff": self.pi_nh_angle_cutoff,
            "pi_oh_distance_cutoff": self.pi_oh_distance_cutoff,
            "pi_oh_angle_cutoff": self.pi_oh_angle_cutoff,
            "pi_sh_distance_cutoff": self.pi_sh_distance_cutoff,
            "pi_sh_angle_cutoff": self.pi_sh_angle_cutoff,
            # General parameters
            "covalent_cutoff_factor": self.covalent_cutoff_factor,
            "analysis_mode": self.analysis_mode,
            # PDB fixing parameters
            "fix_pdb_enabled": self.fix_pdb_enabled,
            "fix_pdb_method": self.fix_pdb_method,
            "fix_pdb_add_hydrogens": self.fix_pdb_add_hydrogens,
            "fix_pdb_add_heavy_atoms": self.fix_pdb_add_heavy_atoms,
            "fix_pdb_replace_nonstandard": self.fix_pdb_replace_nonstandard,
            "fix_pdb_remove_heterogens": self.fix_pdb_remove_heterogens,
            "fix_pdb_keep_water": self.fix_pdb_keep_water,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "AnalysisParameters":
        """Create AnalysisParameters object from dictionary.

        :param data: Dictionary containing parameter values
        :type data: dict
        :returns: New AnalysisParameters object
        :rtype: AnalysisParameters
        """
        return cls(**data)

    def validate(self) -> List[str]:
        """Validate parameter values and return list of validation errors.

        Checks all parameter values for validity and logical consistency.
        Returns a list of validation error messages, empty list if all valid.

        :returns: List of validation error messages
        :rtype: List[str]
        """

        errors = []

        # Distance parameter validation
        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.hb_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Hydrogen bond distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.hb_donor_acceptor_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Donor-acceptor distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.whb_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Weak hydrogen bond distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.whb_donor_acceptor_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Weak hydrogen bond donor-acceptor distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.xb_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Halogen bond distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.pi_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"π interaction distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        # π interaction subtype distance validation
        pi_subtypes = [
            ("pi_ccl_distance_cutoff", "C-Cl...π"),
            ("pi_cbr_distance_cutoff", "C-Br...π"),
            ("pi_ci_distance_cutoff", "C-I...π"),
            ("pi_ch_distance_cutoff", "C-H...π"),
            ("pi_nh_distance_cutoff", "N-H...π"),
            ("pi_oh_distance_cutoff", "O-H...π"),
            ("pi_sh_distance_cutoff", "S-H...π"),
        ]

        for param_name, desc in pi_subtypes:
            param_value = getattr(self, param_name)
            if not (
                ParameterRanges.MIN_DISTANCE
                <= param_value
                <= ParameterRanges.MAX_DISTANCE
            ):
                errors.append(
                    f"{desc} distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
                )

        # Angle parameter validation
        if not (
            ParameterRanges.MIN_ANGLE
            <= self.hb_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"Hydrogen bond angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.whb_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"Weak hydrogen bond angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.xb_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"Halogen bond angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.pi_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"π interaction angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        # π interaction subtype angle validation
        pi_subtype_angles = [
            ("pi_ccl_angle_cutoff", "C-Cl...π"),
            ("pi_cbr_angle_cutoff", "C-Br...π"),
            ("pi_ci_angle_cutoff", "C-I...π"),
            ("pi_ch_angle_cutoff", "C-H...π"),
            ("pi_nh_angle_cutoff", "N-H...π"),
            ("pi_oh_angle_cutoff", "O-H...π"),
            ("pi_sh_angle_cutoff", "S-H...π"),
        ]

        for param_name, desc in pi_subtype_angles:
            param_value = getattr(self, param_name)
            if not (
                ParameterRanges.MIN_ANGLE <= param_value <= ParameterRanges.MAX_ANGLE
            ):
                errors.append(
                    f"{desc} angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
                )

        # Covalent factor validation
        if not (
            ParameterRanges.MIN_COVALENT_FACTOR
            <= self.covalent_cutoff_factor
            <= ParameterRanges.MAX_COVALENT_FACTOR
        ):
            errors.append(
                f"Covalent bond factor must be between {ParameterRanges.MIN_COVALENT_FACTOR}-{ParameterRanges.MAX_COVALENT_FACTOR}"
            )

        # Analysis mode validation

        if self.analysis_mode not in AnalysisModes.ALL_MODES:
            errors.append(
                f"Analysis mode must be one of: {', '.join(AnalysisModes.ALL_MODES)}"
            )

        # PDB fixing parameter validation
        if self.fix_pdb_method not in PDBFixingModes.ALL_METHODS:
            errors.append(
                f"PDB fixing method must be one of: {', '.join(PDBFixingModes.ALL_METHODS)}"
            )

        # Logical consistency validation for PDB fixing
        if self.fix_pdb_enabled:
            if self.fix_pdb_method == "openbabel":
                # OpenBabel only supports hydrogen addition
                if self.fix_pdb_add_heavy_atoms:
                    errors.append(
                        "OpenBabel does not support adding heavy atoms - only PDBFixer supports this"
                    )
                if self.fix_pdb_replace_nonstandard:
                    errors.append(
                        "OpenBabel does not support replacing nonstandard residues - only PDBFixer supports this"
                    )
                if self.fix_pdb_remove_heterogens:
                    errors.append(
                        "OpenBabel does not support removing heterogens - only PDBFixer supports this"
                    )

            # At least one operation must be selected when fixing is enabled
            operations_selected = [
                self.fix_pdb_add_hydrogens,
                self.fix_pdb_add_heavy_atoms,
                self.fix_pdb_replace_nonstandard,
                self.fix_pdb_remove_heterogens,
            ]
            if not any(operations_selected):
                errors.append(
                    "At least one PDB fixing operation must be selected when PDB fixing is enabled"
                )

        # π-π stacking interaction parameter validation
        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.pi_pi_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"π-π stacking distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.pi_pi_parallel_angle_cutoff
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"π-π parallel angle cutoff must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.pi_pi_tshaped_angle_min
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"π-π T-shaped angle minimum must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.pi_pi_tshaped_angle_max
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"π-π T-shaped angle maximum must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.pi_pi_offset_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"π-π offset cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        # Logical consistency for π-π T-shaped angles
        if self.pi_pi_tshaped_angle_min >= self.pi_pi_tshaped_angle_max:
            errors.append("π-π T-shaped angle minimum must be less than maximum")

        # Carbonyl interaction parameter validation
        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.carbonyl_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"Carbonyl interaction distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.carbonyl_angle_min
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"Carbonyl interaction angle minimum must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.carbonyl_angle_max
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"Carbonyl interaction angle maximum must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        # Logical consistency for carbonyl angles
        if self.carbonyl_angle_min >= self.carbonyl_angle_max:
            errors.append(
                "Carbonyl interaction angle minimum must be less than maximum"
            )

        # n→π* interaction parameter validation
        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.n_pi_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"n→π* interaction distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_DISTANCE
            <= self.n_pi_sulfur_distance_cutoff
            <= ParameterRanges.MAX_DISTANCE
        ):
            errors.append(
                f"n→π* sulfur interaction distance cutoff must be between {ParameterRanges.MIN_DISTANCE}-{ParameterRanges.MAX_DISTANCE}Å"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.n_pi_angle_min
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"n→π* interaction angle minimum must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        if not (
            ParameterRanges.MIN_ANGLE
            <= self.n_pi_angle_max
            <= ParameterRanges.MAX_ANGLE
        ):
            errors.append(
                f"n→π* interaction angle maximum must be between {ParameterRanges.MIN_ANGLE}-{ParameterRanges.MAX_ANGLE}°"
            )

        # Logical consistency for n→π* angles
        if self.n_pi_angle_min >= self.n_pi_angle_max:
            errors.append("n→π* interaction angle minimum must be less than maximum")

        return errors


# Parameter validation ranges
class ParameterRanges:
    """Valid ranges for analysis parameters."""

    # Distance ranges (Angstroms)
    MIN_DISTANCE = 0.5
    MAX_DISTANCE = 6

    # Angle ranges (degrees)
    MIN_ANGLE = 0.0
    MAX_ANGLE = 180.0

    # Factor ranges
    MIN_COVALENT_FACTOR = 0.0
    MAX_COVALENT_FACTOR = 1.0


# PDB fixing mode constants
class PDBFixingModes:
    """Available PDB structure fixing modes."""

    OPENBABEL = "openbabel"
    PDBFIXER = "pdbfixer"

    ALL_METHODS = [OPENBABEL, PDBFIXER]

    # Available operations
    ADD_HYDROGENS = "add_hydrogens"
    ADD_HEAVY_ATOMS = "add_heavy_atoms"  # PDBFixer only
    REPLACE_NONSTANDARD = "replace_nonstandard"  # PDBFixer only
    REMOVE_HETEROGENS = "remove_heterogens"  # PDBFixer only

    # Operations available for each method
    OPENBABEL_OPERATIONS = [ADD_HYDROGENS]
    PDBFIXER_OPERATIONS = [
        ADD_HYDROGENS,
        ADD_HEAVY_ATOMS,
        REPLACE_NONSTANDARD,
        REMOVE_HETEROGENS,
    ]


# Analysis mode constants
class AnalysisModes:
    """Available analysis modes."""

    COMPLETE = "complete"
    LOCAL = "local"

    ALL_MODES = [COMPLETE, LOCAL]


# Bond detection method constants
class BondDetectionMethods:
    """Available bond detection methods."""

    CONECT_RECORDS = "conect_records"
    RESIDUE_LOOKUP = "residue_lookup"
    DISTANCE_BASED = "distance_based"

    ALL_METHODS = [CONECT_RECORDS, RESIDUE_LOOKUP, DISTANCE_BASED]
