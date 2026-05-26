"""
Centralized configuration for HBAT UI, export, and parameter management.

This module defines the single source of truth for:
- Interaction type columns and display properties
- Parameter definitions and validation ranges
- Data extraction patterns for all interaction types

This eliminates duplication across GUI, Server, and Export modules.
"""

import math
from dataclasses import dataclass, field, asdict
from typing import Any, Callable, Dict, List, Optional, Union


@dataclass
class ColumnConfig:
    """Configuration for a single column in results display/export.

    Defines how a field should be displayed, exported, and formatted
    across all interfaces (GUI, Server, Export).
    """

    name: str  # Python attribute name (e.g., "donor_residue")
    label: str  # Display label (e.g., "Donor Residue")
    csv_header: str = ""  # CSV header name (defaults to label with underscores)
    json_key: str = ""  # JSON key (defaults to name)
    width: int = 100  # UI column width (pixels)
    precision: Optional[int] = None  # Decimal places for numbers (None = as-is)
    accessor: Optional[Callable[[Any], Any]] = (
        None  # Function to extract value from object
    )

    def __post_init__(self):
        """Set defaults for optional fields."""
        if not self.csv_header:
            self.csv_header = self.label.replace(" ", "_")
        if not self.json_key:
            self.json_key = self.name


@dataclass
class InteractionConfig:
    """Configuration for a single interaction type.

    Defines all metadata needed to display, export, and format
    a specific type of molecular interaction.
    """

    id: str  # Unique identifier (e.g., "hydrogen_bonds")
    label: str  # Display name (e.g., "Hydrogen Bonds")
    icon: str = "link"  # Icon name for UI
    analyzer_attr: str = ""  # Attribute name on analyzer (defaults to id)
    filename: str = ""  # Base filename for export (defaults to id)
    columns: List[ColumnConfig] = field(default_factory=list)
    summary_keys: List[str] = field(default_factory=list)  # Keys in get_summary()

    def __post_init__(self):
        """Set defaults for optional fields."""
        if not self.analyzer_attr:
            self.analyzer_attr = self.id
        if not self.filename:
            self.filename = self.id.replace("_", "_")


@dataclass
class ParameterConfig:
    """Configuration for a single analysis parameter.

    Defines ranges, defaults, and metadata for parameters in
    the AnalysisParameters class.
    """

    name: str  # Parameter name
    label: str  # Display label
    default: Union[float, bool, str]  # Default value
    min_val: Optional[Union[float, int]] = None  # Minimum value (None = no limit)
    max_val: Optional[Union[float, int]] = None  # Maximum value (None = no limit)
    category: str = "General"  # Parameter category for grouping
    description: str = ""  # Detailed description
    param_type: str = "float"  # "float", "int", "bool", "str"


# ============================================================================
# INTERACTION CONFIGURATIONS
# ============================================================================

HYDROGEN_BONDS_CONFIG = InteractionConfig(
    id="hydrogen_bonds",
    label="Hydrogen Bonds",
    icon="link",
    analyzer_attr="hydrogen_bonds",
    filename="h_bonds",
    columns=[
        ColumnConfig(
            name="donor_residue",
            label="Donor Residue",
            width=120,
            accessor=lambda hb: hb.donor_residue,
        ),
        ColumnConfig(
            name="donor_atom",
            label="Donor Atom",
            width=100,
            accessor=lambda hb: hb.donor.name,
        ),
        ColumnConfig(
            name="hydrogen",
            label="Hydrogen Atom",
            width=120,
            accessor=lambda hb: hb.hydrogen.name,
        ),
        ColumnConfig(
            name="acceptor_residue",
            label="Acceptor Residue",
            width=120,
            accessor=lambda hb: hb.acceptor_residue,
        ),
        ColumnConfig(
            name="acceptor_atom",
            label="Acceptor Atom",
            width=100,
            accessor=lambda hb: hb.acceptor.name,
        ),
        ColumnConfig(
            name="distance",
            label="H...A (Å)",
            precision=3,
            width=90,
            accessor=lambda hb: f"{hb.distance:.3f}",
        ),
        ColumnConfig(
            name="angle",
            label="Angle (°)",
            precision=1,
            width=90,
            accessor=lambda hb: f"{math.degrees(hb.angle):.1f}",
        ),
        ColumnConfig(
            name="da_distance",
            label="D...A (Å)",
            precision=3,
            width=90,
            accessor=lambda hb: f"{hb.donor_acceptor_distance:.3f}",
        ),
        ColumnConfig(
            name="type", label="Type", width=120, accessor=lambda hb: hb.bond_type
        ),
        ColumnConfig(
            name="da_props",
            label="D-A Props",
            width=100,
            csv_header="D-A_Properties",
            accessor=lambda hb: hb.donor_acceptor_properties,
        ),
        ColumnConfig(
            name="bs_int",
            label="B/S",
            width=70,
            csv_header="B/S_Interaction",
            accessor=lambda hb: hb.get_backbone_sidechain_interaction(),
        ),
    ],
    summary_keys=["hydrogen_bonds"],
)

WATER_BRIDGES_CONFIG = InteractionConfig(
    id="water_bridges",
    label="Water Bridges",
    icon="water",
    analyzer_attr="water_bridges",
    filename="water_bridges",
    columns=[
        ColumnConfig(
            name="start_res",
            label="Start Residue",
            width=140,
            accessor=lambda wb: wb.get_donor_residue(),
        ),
        ColumnConfig(
            name="end_res",
            label="End Residue",
            width=140,
            accessor=lambda wb: wb.get_acceptor_residue(),
        ),
        ColumnConfig(
            name="hops", label="Hops", width=60, accessor=lambda wb: wb.bridge_length
        ),
        ColumnConfig(
            name="water_residues",
            label="Water Residues",
            width=300,
            accessor=lambda wb: "; ".join(wb.water_residues),
        ),
        ColumnConfig(
            name="distance",
            label="Distance (Å)",
            precision=2,
            width=100,
            accessor=lambda wb: f"{wb.get_donor_acceptor_distance():.2f}",
        ),
    ],
    summary_keys=["water_bridges"],
)

HALOGEN_BONDS_CONFIG = InteractionConfig(
    id="halogen_bonds",
    label="Halogen Bonds",
    icon="atom",
    analyzer_attr="halogen_bonds",
    filename="x_bonds",
    columns=[
        ColumnConfig(
            name="halogen_res",
            label="Halogen Residue",
            width=140,
            accessor=lambda xb: xb.donor_residue,
        ),
        ColumnConfig(
            name="donor_atom",
            label="Donor Atom",
            width=100,
            accessor=lambda xb: xb.donor.name,
        ),
        ColumnConfig(
            name="halogen_atom",
            label="Halogen Atom",
            width=120,
            accessor=lambda xb: xb.halogen.name,
        ),
        ColumnConfig(
            name="acceptor_res",
            label="Acceptor Residue",
            width=140,
            accessor=lambda xb: xb.acceptor_residue,
        ),
        ColumnConfig(
            name="acceptor_atom",
            label="Acceptor Atom",
            width=120,
            accessor=lambda xb: xb.acceptor.name,
        ),
        ColumnConfig(
            name="distance",
            label="X...A (Å)",
            precision=3,
            width=90,
            accessor=lambda xb: f"{xb.distance:.3f}",
        ),
        ColumnConfig(
            name="angle",
            label="Angle (°)",
            precision=1,
            width=90,
            accessor=lambda xb: f"{math.degrees(xb.angle):.1f}",
        ),
        ColumnConfig(
            name="type", label="Type", width=120, accessor=lambda xb: xb.bond_type
        ),
        ColumnConfig(
            name="bs_interaction",
            label="B/S Interaction",
            width=100,
            accessor=lambda xb: xb.get_backbone_sidechain_interaction(),
        ),
        ColumnConfig(
            name="da_properties",
            label="D-A Properties",
            width=120,
            accessor=lambda xb: xb.donor_acceptor_properties,
        ),
    ],
    summary_keys=["halogen_bonds"],
)

PI_INTERACTIONS_CONFIG = InteractionConfig(
    id="pi_interactions",
    label="π Interactions",
    icon="hexagon",
    analyzer_attr="pi_interactions",
    filename="pi_interactions",
    columns=[
        ColumnConfig(
            name="donor_res",
            label="Donor Residue",
            width=140,
            accessor=lambda pi: pi.donor_residue,
        ),
        ColumnConfig(
            name="donor_atom",
            label="Donor Atom",
            width=120,
            accessor=lambda pi: pi.donor.name,
        ),
        ColumnConfig(
            name="pi_res",
            label="π Residue",
            width=140,
            accessor=lambda pi: pi.acceptor_residue,
        ),
        ColumnConfig(
            name="distance",
            label="H...π (Å)",
            precision=3,
            width=110,
            accessor=lambda pi: f"{pi.distance:.3f}",
        ),
        ColumnConfig(
            name="angle",
            label="Angle (°)",
            precision=1,
            width=110,
            accessor=lambda pi: f"{math.degrees(pi.angle):.1f}",
        ),
        ColumnConfig(
            name="type",
            label="Type",
            width=90,
            accessor=lambda pi: pi.get_interaction_type_display(),
        ),
        ColumnConfig(
            name="da_props",
            label="D-A Props",
            width=100,
            csv_header="D-A_Properties",
            accessor=lambda pi: pi.donor_acceptor_properties,
        ),
        ColumnConfig(
            name="bs_int",
            label="B/S",
            width=70,
            csv_header="B/S_Interaction",
            accessor=lambda pi: pi.get_backbone_sidechain_interaction(),
        ),
    ],
    summary_keys=["pi_interactions"],
)

PI_PI_STACKING_CONFIG = InteractionConfig(
    id="pi_pi_interactions",
    label="π-π Stacking",
    icon="ring",
    analyzer_attr="pi_pi_interactions",
    filename="pi_pi_interactions",
    columns=[
        ColumnConfig(
            name="ring1_residue",
            label="Ring 1 Residue",
            width=140,
            accessor=lambda pi_pi: pi_pi.ring1_residue,
        ),
        ColumnConfig(
            name="ring1_type",
            label="Ring 1 Type",
            width=100,
            accessor=lambda pi_pi: pi_pi.ring1_type,
        ),
        ColumnConfig(
            name="ring2_residue",
            label="Ring 2 Residue",
            width=140,
            accessor=lambda pi_pi: pi_pi.ring2_residue,
        ),
        ColumnConfig(
            name="ring2_type",
            label="Ring 2 Type",
            width=100,
            accessor=lambda pi_pi: pi_pi.ring2_type,
        ),
        ColumnConfig(
            name="distance",
            label="Distance (Å)",
            precision=3,
            width=100,
            accessor=lambda pi_pi: f"{pi_pi.distance:.3f}",
        ),
        ColumnConfig(
            name="plane_angle",
            label="Plane Angle (°)",
            precision=1,
            width=120,
            csv_header="Plane_Angle_Degrees",
            accessor=lambda pi_pi: f"{pi_pi.plane_angle:.1f}",
        ),
        ColumnConfig(
            name="offset",
            label="Offset (Å)",
            precision=3,
            width=100,
            accessor=lambda pi_pi: f"{pi_pi.offset:.3f}",
        ),
        ColumnConfig(
            name="stacking_type",
            label="Stacking Type",
            width=120,
            accessor=lambda pi_pi: pi_pi.stacking_type,
        ),
    ],
    summary_keys=["pi_pi_stacking"],
)

CARBONYL_INTERACTIONS_CONFIG = InteractionConfig(
    id="carbonyl_interactions",
    label="Carbonyl Interactions",
    icon="molecule",
    analyzer_attr="carbonyl_interactions",
    filename="carbonyl_interactions",
    columns=[
        ColumnConfig(
            name="acceptor_res",
            label="Acceptor Residue",
            width=130,
            accessor=lambda c: c.donor_residue,
        ),
        ColumnConfig(
            name="acceptor_atom",
            label="Acceptor Atom",
            width=120,
            accessor=lambda c: c.donor_carbon.name,
        ),
        ColumnConfig(
            name="carbonyl_res",
            label="Carbonyl Residue",
            width=130,
            accessor=lambda c: c.acceptor_residue,
        ),
        ColumnConfig(
            name="carbonyl_atoms",
            label="Carbonyl C=O",
            width=120,
            accessor=lambda c: f"{c.acceptor_carbon.name}={c.acceptor_oxygen.name}",
        ),
        ColumnConfig(
            name="distance",
            label="O···C Distance (Å)",
            precision=3,
            width=140,
            accessor=lambda c: f"{c.distance:.3f}",
        ),
        ColumnConfig(
            name="angle",
            label="Bürgi-Dunitz Angle (°)",
            precision=1,
            width=150,
            csv_header="Burgi_Dunitz_Angle_Degrees",
            accessor=lambda c: f"{c.burgi_dunitz_angle:.1f}",
        ),
        ColumnConfig(
            name="carbonyl_type",
            label="Carbonyl Type",
            width=120,
            csv_header="Interaction_Type",
            accessor=lambda c: c.interaction_classification,
        ),
        ColumnConfig(
            name="bs_int",
            label="B/S",
            width=70,
            csv_header="B/S_Interaction",
            accessor=lambda c: "Backbone" if c.is_backbone else "Sidechain",
        ),
    ],
    summary_keys=["carbonyl_interactions"],
)

N_PI_INTERACTIONS_CONFIG = InteractionConfig(
    id="n_pi_interactions",
    label="n→π* Interactions",
    icon="electric",
    analyzer_attr="n_pi_interactions",
    filename="n_pi_interactions",
    columns=[
        ColumnConfig(
            name="donor_res",
            label="Donor Residue",
            width=120,
            accessor=lambda n_pi: n_pi.donor_residue,
        ),
        ColumnConfig(
            name="donor_atom",
            label="Donor Atom",
            width=120,
            csv_header="Lone_Pair_Atom",
            accessor=lambda n_pi: n_pi.lone_pair_atom.name,
        ),
        ColumnConfig(
            name="pi_res",
            label="π Residue",
            width=120,
            csv_header="Pi_System",
            accessor=lambda n_pi: n_pi.acceptor_residue,
        ),
        ColumnConfig(
            name="distance",
            label="Distance (Å)",
            precision=3,
            width=110,
            accessor=lambda n_pi: f"{n_pi.distance:.3f}",
        ),
        ColumnConfig(
            name="angle",
            label="Angle (°)",
            precision=1,
            width=110,
            csv_header="Angle_To_Plane_Degrees",
            accessor=lambda n_pi: f"{n_pi.angle_to_plane:.1f}",
        ),
        ColumnConfig(
            name="donor_element",
            label="Donor Element",
            width=120,
            accessor=lambda n_pi: n_pi.lone_pair_atom.element,
        ),
        ColumnConfig(
            name="bs_int",
            label="B/S",
            width=70,
            csv_header="B/S_Interaction",
            accessor=lambda n_pi: "N/A",
        ),  # Not currently available
    ],
    summary_keys=["n_pi_interactions"],
)

COOPERATIVITY_CHAINS_CONFIG = InteractionConfig(
    id="cooperativity_chains",
    label="Cooperativity Chains",
    icon="link-2",
    analyzer_attr="cooperativity_chains",
    filename="cooperativity_chains",
    columns=[
        ColumnConfig(
            name="chain_id",
            label="Chain ID",
            width=100,
            accessor=lambda chain, idx=0: idx + 1,
        ),
        ColumnConfig(
            name="chain_length",
            label="Length",
            width=100,
            accessor=lambda chain: chain.chain_length,
        ),
        ColumnConfig(
            name="chain_description",
            label="Chain Description",
            width=1000,
            accessor=lambda chain: " → ".join(
                f"{i.get_donor_residue()}({i.get_donor_atom().name if i.get_donor_atom() else '?'})"
                for i in chain.interactions
            ),
        ),
    ],
    summary_keys=["cooperativity_chains"],
)

# Dictionary mapping interaction type IDs to their configurations
INTERACTION_CONFIGS: Dict[str, InteractionConfig] = {
    "hydrogen_bonds": HYDROGEN_BONDS_CONFIG,
    "water_bridges": WATER_BRIDGES_CONFIG,
    "halogen_bonds": HALOGEN_BONDS_CONFIG,
    "pi_interactions": PI_INTERACTIONS_CONFIG,
    "pi_pi_interactions": PI_PI_STACKING_CONFIG,
    "carbonyl_interactions": CARBONYL_INTERACTIONS_CONFIG,
    "n_pi_interactions": N_PI_INTERACTIONS_CONFIG,
    "cooperativity_chains": COOPERATIVITY_CHAINS_CONFIG,
}


# ============================================================================
# PARAMETER CONFIGURATIONS
# ============================================================================

PARAMETER_CONFIGS: List[ParameterConfig] = [
    # Hydrogen Bond Parameters
    ParameterConfig(
        name="hb_distance_cutoff",
        label="H...A Distance Cutoff",
        default=2.5,
        min_val=1.5,
        max_val=4.0,
        category="Hydrogen Bonds",
        description="Maximum H...A distance for hydrogen bonds (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="hb_angle_cutoff",
        label="D-H...A Angle Cutoff",
        default=120.0,
        min_val=90.0,
        max_val=180.0,
        category="Hydrogen Bonds",
        description="Minimum D-H...A angle for hydrogen bonds (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="hb_donor_acceptor_cutoff",
        label="D...A Distance Cutoff",
        default=3.5,
        min_val=2.5,
        max_val=4.5,
        category="Hydrogen Bonds",
        description="Maximum D...A distance for hydrogen bonds (Å)",
        param_type="float",
    ),
    # Weak Hydrogen Bond Parameters
    ParameterConfig(
        name="whb_distance_cutoff",
        label="C-H...A Distance Cutoff",
        default=3.6,
        min_val=2.5,
        max_val=4.5,
        category="Weak Hydrogen Bonds",
        description="Maximum H...A distance for weak hydrogen bonds (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="whb_angle_cutoff",
        label="C-H...A Angle Cutoff",
        default=150.0,
        min_val=90.0,
        max_val=180.0,
        category="Weak Hydrogen Bonds",
        description="Minimum D-H...A angle for weak hydrogen bonds (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="whb_donor_acceptor_cutoff",
        label="D...A Distance Cutoff",
        default=3.5,
        min_val=2.5,
        max_val=4.5,
        category="Weak Hydrogen Bonds",
        description="Maximum D...A distance for weak hydrogen bonds (Å)",
        param_type="float",
    ),
    # Halogen Bond Parameters
    ParameterConfig(
        name="xb_distance_cutoff",
        label="X...A Distance Cutoff",
        default=3.9,
        min_val=3.0,
        max_val=4.5,
        category="Halogen Bonds",
        description="Maximum X...A distance for halogen bonds (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="xb_angle_cutoff",
        label="C-X...A Angle Cutoff",
        default=150.0,
        min_val=90.0,
        max_val=180.0,
        category="Halogen Bonds",
        description="Minimum C-X...A angle for halogen bonds (degrees)",
        param_type="float",
    ),
    # π Interaction Parameters (Legacy)
    ParameterConfig(
        name="pi_distance_cutoff",
        label="H...π Distance Cutoff (Legacy)",
        default=3.5,
        min_val=2.5,
        max_val=4.5,
        category="π Interactions",
        description="Maximum H...π distance for π interactions (Å, legacy)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_angle_cutoff",
        label="D-H...π Angle Cutoff (Legacy)",
        default=110.0,
        min_val=80.0,
        max_val=180.0,
        category="π Interactions",
        description="Minimum D-H...π angle for π interactions (degrees, legacy)",
        param_type="float",
    ),
    # π Interaction Subtype Parameters - Halogen-π
    ParameterConfig(
        name="pi_ccl_distance_cutoff",
        label="C-Cl...π Distance Cutoff",
        default=3.5,
        min_val=2.5,
        max_val=4.5,
        category="π Interactions - Halogen",
        description="Maximum C-Cl...π distance (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_ccl_angle_cutoff",
        label="C-Cl...π Angle Cutoff",
        default=145.0,
        min_val=100.0,
        max_val=180.0,
        category="π Interactions - Halogen",
        description="Minimum C-Cl...π angle (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_cbr_distance_cutoff",
        label="C-Br...π Distance Cutoff",
        default=3.5,
        min_val=2.5,
        max_val=4.5,
        category="π Interactions - Halogen",
        description="Maximum C-Br...π distance (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_cbr_angle_cutoff",
        label="C-Br...π Angle Cutoff",
        default=155.0,
        min_val=100.0,
        max_val=180.0,
        category="π Interactions - Halogen",
        description="Minimum C-Br...π angle (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_ci_distance_cutoff",
        label="C-I...π Distance Cutoff",
        default=3.6,
        min_val=2.5,
        max_val=4.5,
        category="π Interactions - Halogen",
        description="Maximum C-I...π distance (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_ci_angle_cutoff",
        label="C-I...π Angle Cutoff",
        default=165.0,
        min_val=100.0,
        max_val=180.0,
        category="π Interactions - Halogen",
        description="Minimum C-I...π angle (degrees)",
        param_type="float",
    ),
    # π Interaction Subtype Parameters - Hydrogen-π
    ParameterConfig(
        name="pi_ch_distance_cutoff",
        label="C-H...π Distance Cutoff",
        default=3.5,
        min_val=2.5,
        max_val=4.5,
        category="π Interactions - Hydrogen",
        description="Maximum C-H...π distance (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_ch_angle_cutoff",
        label="C-H...π Angle Cutoff",
        default=110.0,
        min_val=80.0,
        max_val=180.0,
        category="π Interactions - Hydrogen",
        description="Minimum C-H...π angle (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_nh_distance_cutoff",
        label="N-H...π Distance Cutoff",
        default=3.2,
        min_val=2.5,
        max_val=4.0,
        category="π Interactions - Hydrogen",
        description="Maximum N-H...π distance (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_nh_angle_cutoff",
        label="N-H...π Angle Cutoff",
        default=115.0,
        min_val=80.0,
        max_val=180.0,
        category="π Interactions - Hydrogen",
        description="Minimum N-H...π angle (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_oh_distance_cutoff",
        label="O-H...π Distance Cutoff",
        default=3.0,
        min_val=2.5,
        max_val=4.0,
        category="π Interactions - Hydrogen",
        description="Maximum O-H...π distance (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_oh_angle_cutoff",
        label="O-H...π Angle Cutoff",
        default=115.0,
        min_val=80.0,
        max_val=180.0,
        category="π Interactions - Hydrogen",
        description="Minimum O-H...π angle (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_sh_distance_cutoff",
        label="S-H...π Distance Cutoff",
        default=3.8,
        min_val=2.5,
        max_val=4.5,
        category="π Interactions - Hydrogen",
        description="Maximum S-H...π distance (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_sh_angle_cutoff",
        label="S-H...π Angle Cutoff",
        default=105.0,
        min_val=80.0,
        max_val=180.0,
        category="π Interactions - Hydrogen",
        description="Minimum S-H...π angle (degrees)",
        param_type="float",
    ),
    # π-π Stacking Parameters
    ParameterConfig(
        name="pi_pi_distance_cutoff",
        label="Centroid Distance Cutoff",
        default=3.8,
        min_val=3.0,
        max_val=4.5,
        category="π-π Stacking",
        description="Maximum centroid-to-centroid distance for π-π stacking (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_pi_parallel_angle_cutoff",
        label="Parallel Angle Cutoff",
        default=30.0,
        min_val=0.0,
        max_val=90.0,
        category="π-π Stacking",
        description="Maximum angle for parallel π-π stacking (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_pi_tshaped_angle_min",
        label="T-shaped Angle Min",
        default=60.0,
        min_val=0.0,
        max_val=90.0,
        category="π-π Stacking",
        description="Minimum angle for T-shaped π-π stacking (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_pi_tshaped_angle_max",
        label="T-shaped Angle Max",
        default=90.0,
        min_val=0.0,
        max_val=90.0,
        category="π-π Stacking",
        description="Maximum angle for T-shaped π-π stacking (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="pi_pi_offset_cutoff",
        label="Offset Cutoff",
        default=2.0,
        min_val=0.0,
        max_val=3.0,
        category="π-π Stacking",
        description="Maximum lateral offset for parallel π-π stacking (Å)",
        param_type="float",
    ),
    # Carbonyl Interaction Parameters
    ParameterConfig(
        name="carbonyl_distance_cutoff",
        label="O···C Distance Cutoff",
        default=3.2,
        min_val=2.5,
        max_val=4.0,
        category="Carbonyl Interactions",
        description="Maximum O···C distance for carbonyl n→π* interactions (Å, Bürgi-Dunitz)",
        param_type="float",
    ),
    ParameterConfig(
        name="carbonyl_angle_min",
        label="O···C=O Angle Min",
        default=95.0,
        min_val=70.0,
        max_val=120.0,
        category="Carbonyl Interactions",
        description="Minimum O···C=O angle for carbonyl interactions (degrees, Bürgi-Dunitz)",
        param_type="float",
    ),
    ParameterConfig(
        name="carbonyl_angle_max",
        label="O···C=O Angle Max",
        default=125.0,
        min_val=100.0,
        max_val=180.0,
        category="Carbonyl Interactions",
        description="Maximum O···C=O angle for carbonyl interactions (degrees, Bürgi-Dunitz)",
        param_type="float",
    ),
    # n→π* Interaction Parameters
    ParameterConfig(
        name="n_pi_distance_cutoff",
        label="Lone Pair to π Distance Cutoff",
        default=3.6,
        min_val=2.5,
        max_val=4.5,
        category="n→π* Interactions",
        description="Maximum lone pair to π center distance (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="n_pi_sulfur_distance_cutoff",
        label="Sulfur Distance Cutoff",
        default=4.0,
        min_val=3.0,
        max_val=4.5,
        category="n→π* Interactions",
        description="Maximum sulfur lone pair to π center distance (Å)",
        param_type="float",
    ),
    ParameterConfig(
        name="n_pi_angle_min",
        label="Angle to Plane Min",
        default=0.0,
        min_val=0.0,
        max_val=45.0,
        category="n→π* Interactions",
        description="Minimum angle to π plane for n→π* interactions (degrees)",
        param_type="float",
    ),
    ParameterConfig(
        name="n_pi_angle_max",
        label="Angle to Plane Max",
        default=45.0,
        min_val=0.0,
        max_val=90.0,
        category="n→π* Interactions",
        description="Maximum angle to π plane for n→π* interactions (degrees)",
        param_type="float",
    ),
    # General Parameters
    ParameterConfig(
        name="covalent_cutoff_factor",
        label="Covalent Bond Detection Factor",
        default=0.85,
        min_val=0.5,
        max_val=1.0,
        category="General",
        description="Factor for covalent bond detection (0.0-1.0)",
        param_type="float",
    ),
    ParameterConfig(
        name="analysis_mode",
        label="Analysis Mode",
        default="local",
        category="General",
        description="Analysis mode: 'local' or 'complete'",
        param_type="str",
    ),
    # PDB Structure Fixing Parameters
    ParameterConfig(
        name="fix_pdb_enabled",
        label="Enable PDB Fixing",
        default=True,
        category="PDB Structure Fixing",
        description="Enable PDB structure fixing",
        param_type="bool",
    ),
    ParameterConfig(
        name="fix_pdb_method",
        label="PDB Fixing Method",
        default="openbabel",
        category="PDB Structure Fixing",
        description="Method for PDB fixing: 'openbabel' or 'pdbfixer'",
        param_type="str",
    ),
    ParameterConfig(
        name="fix_pdb_add_hydrogens",
        label="Add Hydrogens",
        default=True,
        category="PDB Structure Fixing",
        description="Add missing hydrogen atoms",
        param_type="bool",
    ),
    ParameterConfig(
        name="fix_pdb_add_heavy_atoms",
        label="Add Heavy Atoms",
        default=False,
        category="PDB Structure Fixing",
        description="Add missing heavy atoms (PDBFixer only)",
        param_type="bool",
    ),
    ParameterConfig(
        name="fix_pdb_replace_nonstandard",
        label="Replace Nonstandard Residues",
        default=False,
        category="PDB Structure Fixing",
        description="Replace nonstandard residues (PDBFixer only)",
        param_type="bool",
    ),
    ParameterConfig(
        name="fix_pdb_remove_heterogens",
        label="Remove Heterogens",
        default=False,
        category="PDB Structure Fixing",
        description="Remove heterogens (PDBFixer only)",
        param_type="bool",
    ),
    ParameterConfig(
        name="fix_pdb_keep_water",
        label="Keep Water",
        default=True,
        category="PDB Structure Fixing",
        description="Keep water when removing heterogens (PDBFixer only)",
        param_type="bool",
    ),
]

# Create lookup dictionary for fast parameter access
PARAMETER_LOOKUP: Dict[str, ParameterConfig] = {
    param.name: param for param in PARAMETER_CONFIGS
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def get_interaction_config(interaction_type: str) -> Optional[InteractionConfig]:
    """Get configuration for a specific interaction type.

    Args:
        interaction_type: Interaction type ID (e.g., "hydrogen_bonds")

    Returns:
        InteractionConfig or None if not found
    """
    return INTERACTION_CONFIGS.get(interaction_type)


def get_parameter_config(param_name: str) -> Optional[ParameterConfig]:
    """Get configuration for a specific parameter.

    Args:
        param_name: Parameter name

    Returns:
        ParameterConfig or None if not found
    """
    return PARAMETER_LOOKUP.get(param_name)


def extract_interaction_data(
    interaction: Any, config: InteractionConfig, interaction_index: int = 0
) -> Dict[str, Any]:
    """Extract formatted data from an interaction using config.

    Extracts all column values from an interaction object using the
    accessor functions defined in the config.

    Args:
        interaction: MolecularInteraction object
        config: InteractionConfig for this interaction type
        interaction_index: Index for interactions without unique IDs

    Returns:
        Dictionary mapping column names to formatted values
    """
    data = {}

    for col in config.columns:
        try:
            if col.accessor:
                # Use custom accessor if provided
                if col.name == "chain_id":  # Special case for cooperativity chains
                    value = col.accessor(interaction, interaction_index)
                else:
                    value = col.accessor(interaction)
            else:
                # Default: get attribute directly
                value = getattr(interaction, col.name, None)

            # Apply precision formatting if needed
            if value is not None and col.precision is not None:
                if isinstance(value, (int, float)):
                    value = round(value, col.precision)

            data[col.name] = value
        except Exception as e:
            # Log accessor errors but don't fail
            data[col.name] = None

    return data


def get_parameters_by_category() -> Dict[str, List[ParameterConfig]]:
    """Group parameters by their category.

    Returns:
        Dictionary mapping category names to lists of ParameterConfig
    """
    categories: Dict[str, List[ParameterConfig]] = {}

    for param in PARAMETER_CONFIGS:
        if param.category not in categories:
            categories[param.category] = []
        categories[param.category].append(param)

    return categories
