"""
Command-line interface for HBAT.

This module provides a command-line interface for running HBAT analysis
without the GUI, suitable for batch processing and scripting.
"""

import argparse
import json
import os
import sys
import time
from typing import Any, Dict, List, Optional

from .. import __version__
from ..constants.parameters import ParametersDefault
from ..core.analysis import AnalysisParameters, NPMolecularInteractionAnalyzer
from ..core.pdb_parser import PDBParser
from ..export.results import (
    export_to_csv_files,
    export_to_json_files,
    export_to_json_single_file,
    export_to_txt_single_file,
)


class ProgressBar:
    """Simple CLI progress bar for analysis operations."""

    def __init__(self, width: int = 50):
        """Initialize progress bar.

        :param width: Width of the progress bar in characters
        :type width: int
        """
        self.width = width
        self.current_step = ""
        self.last_progress = -1

    def update(self, message: str, progress: Optional[int] = None) -> None:
        """Update progress bar with new message and optional percentage.

        :param message: Current operation message
        :type message: str
        :param progress: Progress percentage (0-100), optional
        :type progress: Optional[int]
        """
        # Clear previous line and print new status
        print(f"\r\033[K[INFO] {message}", end="", flush=True)

        if progress is not None and progress != self.last_progress:
            # Add progress bar for percentage updates with emoji
            filled_width = int(self.width * progress / 100)
            bar = "●" * filled_width + "○" * (self.width - filled_width)
            print(f" [{bar}] {progress}%", end="", flush=True)
            self.last_progress = progress

    def finish(self, message: str) -> None:
        """Finish progress bar with final message.

        :param message: Final completion message
        :type message: str
        """
        print(f"\r\033[K[INFO] {message}")
        self.last_progress = -1


def create_parser() -> argparse.ArgumentParser:
    """Create command-line argument parser.

    Creates and configures an ArgumentParser with all CLI options for HBAT,
    including input/output options, analysis parameters, and preset management.

    :returns: Configured argument parser
    :rtype: argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description="HBAT - Hydrogen Bond Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.pdb                         # Display results to console
  %(prog)s input.pdb -o results.txt          # Save results to text file
  %(prog)s input.pdb -o results.json         # Save results to JSON file (single file)
  %(prog)s input.pdb --csv results           # Export to multiple CSV files (one per interaction type)
  %(prog)s input.pdb --json results          # Export to multiple JSON files (one per interaction type)
  %(prog)s input.pdb --hb-distance 3.0       # Custom H-bond distance cutoff
  %(prog)s input.pdb --mode local            # Local interactions only
  %(prog)s --list-presets                    # List available presets
  %(prog)s input.pdb --preset high_resolution # Use preset with custom overrides
        """,
    )

    # Version
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    # Input file (optional when listing presets)
    parser.add_argument("input", nargs="?", help="Input PDB file")

    # Output options
    parser.add_argument(
        "-o",
        "--output",
        help="Output file (format auto-detected from extension: .txt, .json)",
    )
    parser.add_argument(
        "--json", help="Export to multiple JSON files (base name for files)"
    )
    parser.add_argument(
        "--csv", help="Export to multiple CSV files (base name for files)"
    )

    # Preset options
    preset_group = parser.add_argument_group("Preset Options")
    preset_group.add_argument(
        "--preset", type=str, help="Load parameters from preset file (.hbat or .json)"
    )
    preset_group.add_argument(
        "--list-presets",
        action="store_true",
        help="List available example presets and exit",
    )

    # Analysis parameters
    param_group = parser.add_argument_group("Analysis Parameters")
    param_group.add_argument(
        "--hb-distance",
        type=float,
        default=ParametersDefault.HB_DISTANCE_CUTOFF,
        help=f"Hydrogen bond H...A distance cutoff in Å (default: {ParametersDefault.HB_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--hb-angle",
        type=float,
        default=ParametersDefault.HB_ANGLE_CUTOFF,
        help=f"Hydrogen bond D-H...A angle cutoff in degrees (default: {ParametersDefault.HB_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--da-distance",
        type=float,
        default=ParametersDefault.HB_DA_DISTANCE,
        help=f"Donor-acceptor distance cutoff in Å (default: {ParametersDefault.HB_DA_DISTANCE})",
    )
    param_group.add_argument(
        "--whb-distance",
        type=float,
        default=ParametersDefault.WHB_DISTANCE_CUTOFF,
        help=f"Weak hydrogen bond H...A distance cutoff in Å for carbon donors (default: {ParametersDefault.WHB_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--whb-angle",
        type=float,
        default=ParametersDefault.WHB_ANGLE_CUTOFF,
        help=f"Weak hydrogen bond D-H...A angle cutoff in degrees for carbon donors (default: {ParametersDefault.WHB_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--whb-da-distance",
        type=float,
        default=ParametersDefault.WHB_DA_DISTANCE,
        help=f"Weak hydrogen bond donor-acceptor distance cutoff in Å for carbon donors (default: {ParametersDefault.WHB_DA_DISTANCE})",
    )
    param_group.add_argument(
        "--xb-distance",
        type=float,
        default=ParametersDefault.XB_DISTANCE_CUTOFF,
        help=f"Halogen bond X...A distance cutoff in Å (default: {ParametersDefault.XB_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--xb-angle",
        type=float,
        default=ParametersDefault.XB_ANGLE_CUTOFF,
        help=f"Halogen bond C-X...A angle cutoff in degrees (default: {ParametersDefault.XB_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-distance",
        type=float,
        default=ParametersDefault.PI_DISTANCE_CUTOFF,
        help=f"π interaction H...π distance cutoff in Å (default: {ParametersDefault.PI_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-angle",
        type=float,
        default=ParametersDefault.PI_ANGLE_CUTOFF,
        help=f"π interaction D-H...π angle cutoff in degrees (default: {ParametersDefault.PI_ANGLE_CUTOFF})",
    )

    # π interaction subtype parameters
    param_group.add_argument(
        "--pi-ccl-distance",
        type=float,
        default=ParametersDefault.PI_CCL_DISTANCE_CUTOFF,
        help=f"C-Cl...π interaction distance cutoff in Å (default: {ParametersDefault.PI_CCL_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-ccl-angle",
        type=float,
        default=ParametersDefault.PI_CCL_ANGLE_CUTOFF,
        help=f"C-Cl...π interaction angle cutoff in degrees (default: {ParametersDefault.PI_CCL_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-cbr-distance",
        type=float,
        default=ParametersDefault.PI_CBR_DISTANCE_CUTOFF,
        help=f"C-Br...π interaction distance cutoff in Å (default: {ParametersDefault.PI_CBR_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-cbr-angle",
        type=float,
        default=ParametersDefault.PI_CBR_ANGLE_CUTOFF,
        help=f"C-Br...π interaction angle cutoff in degrees (default: {ParametersDefault.PI_CBR_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-ci-distance",
        type=float,
        default=ParametersDefault.PI_CI_DISTANCE_CUTOFF,
        help=f"C-I...π interaction distance cutoff in Å (default: {ParametersDefault.PI_CI_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-ci-angle",
        type=float,
        default=ParametersDefault.PI_CI_ANGLE_CUTOFF,
        help=f"C-I...π interaction angle cutoff in degrees (default: {ParametersDefault.PI_CI_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-ch-distance",
        type=float,
        default=ParametersDefault.PI_CH_DISTANCE_CUTOFF,
        help=f"C-H...π interaction distance cutoff in Å (default: {ParametersDefault.PI_CH_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-ch-angle",
        type=float,
        default=ParametersDefault.PI_CH_ANGLE_CUTOFF,
        help=f"C-H...π interaction angle cutoff in degrees (default: {ParametersDefault.PI_CH_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-nh-distance",
        type=float,
        default=ParametersDefault.PI_NH_DISTANCE_CUTOFF,
        help=f"N-H...π interaction distance cutoff in Å (default: {ParametersDefault.PI_NH_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-nh-angle",
        type=float,
        default=ParametersDefault.PI_NH_ANGLE_CUTOFF,
        help=f"N-H...π interaction angle cutoff in degrees (default: {ParametersDefault.PI_NH_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-oh-distance",
        type=float,
        default=ParametersDefault.PI_OH_DISTANCE_CUTOFF,
        help=f"O-H...π interaction distance cutoff in Å (default: {ParametersDefault.PI_OH_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-oh-angle",
        type=float,
        default=ParametersDefault.PI_OH_ANGLE_CUTOFF,
        help=f"O-H...π interaction angle cutoff in degrees (default: {ParametersDefault.PI_OH_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-sh-distance",
        type=float,
        default=ParametersDefault.PI_SH_DISTANCE_CUTOFF,
        help=f"S-H...π interaction distance cutoff in Å (default: {ParametersDefault.PI_SH_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-sh-angle",
        type=float,
        default=ParametersDefault.PI_SH_ANGLE_CUTOFF,
        help=f"S-H...π interaction angle cutoff in degrees (default: {ParametersDefault.PI_SH_ANGLE_CUTOFF})",
    )

    # π-π stacking interaction parameters
    param_group.add_argument(
        "--pi-pi-distance",
        type=float,
        default=ParametersDefault.PI_PI_DISTANCE_CUTOFF,
        help=f"π-π stacking centroid distance cutoff in Å (default: {ParametersDefault.PI_PI_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-pi-parallel-angle",
        type=float,
        default=ParametersDefault.PI_PI_PARALLEL_ANGLE_CUTOFF,
        help=f"π-π parallel stacking max angle cutoff in degrees (default: {ParametersDefault.PI_PI_PARALLEL_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-pi-tshaped-angle-min",
        type=float,
        default=ParametersDefault.PI_PI_TSHAPED_ANGLE_MIN,
        help=f"π-π T-shaped stacking min angle in degrees (default: {ParametersDefault.PI_PI_TSHAPED_ANGLE_MIN})",
    )
    param_group.add_argument(
        "--pi-pi-tshaped-angle-max",
        type=float,
        default=ParametersDefault.PI_PI_TSHAPED_ANGLE_MAX,
        help=f"π-π T-shaped stacking max angle in degrees (default: {ParametersDefault.PI_PI_TSHAPED_ANGLE_MAX})",
    )
    param_group.add_argument(
        "--pi-pi-offset",
        type=float,
        default=ParametersDefault.PI_PI_OFFSET_CUTOFF,
        help=f"π-π parallel stacking max offset in Å (default: {ParametersDefault.PI_PI_OFFSET_CUTOFF})",
    )

    # Carbonyl interaction parameters
    param_group.add_argument(
        "--carbonyl-distance",
        type=float,
        default=ParametersDefault.CARBONYL_DISTANCE_CUTOFF,
        help=f"Carbonyl n→π* O···C distance cutoff in Å (default: {ParametersDefault.CARBONYL_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--carbonyl-angle-min",
        type=float,
        default=ParametersDefault.CARBONYL_ANGLE_MIN,
        help=f"Carbonyl n→π* min O···C=O angle in degrees (default: {ParametersDefault.CARBONYL_ANGLE_MIN})",
    )
    param_group.add_argument(
        "--carbonyl-angle-max",
        type=float,
        default=ParametersDefault.CARBONYL_ANGLE_MAX,
        help=f"Carbonyl n→π* max O···C=O angle in degrees (default: {ParametersDefault.CARBONYL_ANGLE_MAX})",
    )

    # n→π* interaction parameters
    param_group.add_argument(
        "--n-pi-distance",
        type=float,
        default=ParametersDefault.N_PI_DISTANCE_CUTOFF,
        help=f"n→π* interaction distance cutoff in Å (default: {ParametersDefault.N_PI_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--n-pi-sulfur-distance",
        type=float,
        default=ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF,
        help=f"n→π* sulfur-specific distance cutoff in Å (default: {ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--n-pi-angle-min",
        type=float,
        default=ParametersDefault.N_PI_ANGLE_MIN,
        help=f"n→π* min angle to π plane in degrees (default: {ParametersDefault.N_PI_ANGLE_MIN})",
    )
    param_group.add_argument(
        "--n-pi-angle-max",
        type=float,
        default=ParametersDefault.N_PI_ANGLE_MAX,
        help=f"n→π* max angle to π plane in degrees (default: {ParametersDefault.N_PI_ANGLE_MAX})",
    )

    param_group.add_argument(
        "--covalent-factor",
        type=float,
        default=ParametersDefault.COVALENT_CUTOFF_FACTOR,
        help=f"Covalent bond detection factor (default: {ParametersDefault.COVALENT_CUTOFF_FACTOR})",
    )

    # Analysis mode
    param_group.add_argument(
        "--mode",
        choices=["complete", "local"],
        default=ParametersDefault.ANALYSIS_MODE,
        help="Analysis mode: complete (all interactions) or local (intra-residue only)",
    )

    # PDB structure fixing options
    fix_group = parser.add_argument_group("PDB Structure Fixing")
    fix_group.add_argument(
        "--fix-pdb",
        action="store_true",
        help="Enable PDB structure fixing",
    )
    fix_group.add_argument(
        "--fix-method",
        choices=["openbabel", "pdbfixer"],
        default=ParametersDefault.FIX_PDB_METHOD,
        help=f"PDB fixing method: openbabel or pdbfixer (default: {ParametersDefault.FIX_PDB_METHOD})",
    )
    fix_group.add_argument(
        "--fix-add-hydrogens",
        action="store_true",
        default=ParametersDefault.FIX_PDB_ADD_HYDROGENS,
        help="Add missing hydrogen atoms (both OpenBabel and PDBFixer)",
    )
    fix_group.add_argument(
        "--fix-add-heavy-atoms",
        action="store_true",
        help="Add missing heavy atoms (PDBFixer only)",
    )
    fix_group.add_argument(
        "--fix-replace-nonstandard",
        action="store_true",
        help="Replace nonstandard residues (PDBFixer only)",
    )
    fix_group.add_argument(
        "--fix-remove-heterogens",
        action="store_true",
        help="Remove heterogens (PDBFixer only)",
    )
    fix_group.add_argument(
        "--fix-keep-water",
        action="store_true",
        default=ParametersDefault.FIX_PDB_KEEP_WATER,
        help="Keep water when removing heterogens (PDBFixer only)",
    )

    # Output control
    output_group = parser.add_argument_group("Output Control")
    output_group.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Verbose output with detailed progress",
    )
    output_group.add_argument(
        "--quiet", "-q", action="store_true", help="Quiet mode with minimal output"
    )
    output_group.add_argument(
        "--summary-only", action="store_true", help="Output summary statistics only"
    )

    # Analysis filters
    filter_group = parser.add_argument_group("Analysis Filters")
    filter_group.add_argument(
        "--no-hydrogen-bonds", action="store_true", help="Skip hydrogen bond analysis"
    )
    filter_group.add_argument(
        "--no-halogen-bonds", action="store_true", help="Skip halogen bond analysis"
    )
    filter_group.add_argument(
        "--no-pi-interactions", action="store_true", help="Skip π interaction analysis"
    )
    filter_group.add_argument(
        "--no-pi-pi-stacking", action="store_true", help="Skip π-π stacking analysis"
    )
    filter_group.add_argument(
        "--no-carbonyl-interactions",
        action="store_true",
        help="Skip carbonyl n→π* interaction analysis",
    )
    filter_group.add_argument(
        "--no-n-pi-interactions",
        action="store_true",
        help="Skip n→π* interaction analysis",
    )

    return parser


def get_example_presets_directory() -> str:
    """Get the example presets directory relative to the package.

    Locates the example_presets directory that contains predefined
    analysis parameter configurations.

    :returns: Path to the example presets directory
    :rtype: str
    """
    # Get the directory of this file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Go up to the hbat package root, then to example_presets
    package_root = os.path.dirname(os.path.dirname(current_dir))
    example_presets_dir = os.path.join(package_root, "example_presets")

    return example_presets_dir


def list_available_presets() -> None:
    """List all available example presets.

    Displays a formatted list of all available preset files with descriptions
    and icons to help users choose appropriate analysis parameters.

    :returns: None
    :rtype: None
    """
    presets_dir = get_example_presets_directory()

    if not os.path.exists(presets_dir):
        print("No example presets directory found.")
        return

    print("🎯 Available HBAT Presets:")
    print("=" * 50)

    preset_files = [f for f in os.listdir(presets_dir) if f.endswith(".hbat")]
    if not preset_files:
        print("No preset files found in example_presets directory.")
        return

    # Define icons for different preset types
    preset_icons = {
        "high_resolution": "🔬",
        "standard_resolution": "⚙️",
        "low_resolution": "📐",
        "nmr_structures": "🧬",
        "strong_interactions_only": "💪",
        "drug_design_strict": "💊",
        "membrane_proteins": "🧱",
        "weak_interactions_permissive": "🌐",
    }

    for preset_file in sorted(preset_files):
        preset_name = preset_file.replace(".hbat", "")
        icon = preset_icons.get(preset_name, "📄")

        # Try to load and show description
        preset_path = os.path.join(presets_dir, preset_file)
        try:
            with open(preset_path, "r") as f:
                preset_data = json.load(f)
                description = preset_data.get("description", "No description available")
                print(f"{icon} {preset_name}")
                print(f"   {description}")
                print()
        except Exception:
            print(f"{icon} {preset_name}")
            print("   (Unable to load description)")
            print()


def load_preset_file(preset_path: str) -> AnalysisParameters:
    """Load parameters from a preset file.

    Reads and parses a JSON preset file to create AnalysisParameters
    with predefined values for various analysis scenarios.

    :param preset_path: Path to the preset file to load
    :type preset_path: str
    :returns: Analysis parameters loaded from the preset
    :rtype: AnalysisParameters
    :raises ValueError: If preset file format is invalid
    :raises FileNotFoundError: If preset file doesn't exist
    """
    try:
        with open(preset_path, "r", encoding="utf-8") as f:
            preset_data = json.load(f)

        # Validate basic structure
        if "parameters" not in preset_data:
            raise ValueError("Invalid preset file format: missing 'parameters' section")

        params = preset_data["parameters"]

        # Extract parameters with defaults
        hb_params = params.get("hydrogen_bonds", {})
        whb_params = params.get("weak_hydrogen_bonds", {})
        xb_params = params.get("halogen_bonds", {})
        pi_params = params.get("pi_interactions", {})
        pi_pi_params = params.get("pi_pi_stacking", {})
        carbonyl_params = params.get("carbonyl_interactions", {})
        n_pi_params = params.get("n_pi_interactions", {})
        general_params = params.get("general", {})
        fix_params = params.get("pdb_fixing", {})

        return AnalysisParameters(
            hb_distance_cutoff=hb_params.get(
                "h_a_distance_cutoff", ParametersDefault.HB_DISTANCE_CUTOFF
            ),
            hb_angle_cutoff=hb_params.get(
                "dha_angle_cutoff", ParametersDefault.HB_ANGLE_CUTOFF
            ),
            hb_donor_acceptor_cutoff=hb_params.get(
                "d_a_distance_cutoff", ParametersDefault.HB_DA_DISTANCE
            ),
            whb_distance_cutoff=whb_params.get(
                "h_a_distance_cutoff", ParametersDefault.WHB_DISTANCE_CUTOFF
            ),
            whb_angle_cutoff=whb_params.get(
                "dha_angle_cutoff", ParametersDefault.WHB_ANGLE_CUTOFF
            ),
            whb_donor_acceptor_cutoff=whb_params.get(
                "d_a_distance_cutoff", ParametersDefault.WHB_DA_DISTANCE
            ),
            xb_distance_cutoff=xb_params.get(
                "x_a_distance_cutoff", ParametersDefault.XB_DISTANCE_CUTOFF
            ),
            xb_angle_cutoff=xb_params.get(
                "dxa_angle_cutoff", ParametersDefault.XB_ANGLE_CUTOFF
            ),
            pi_distance_cutoff=pi_params.get(
                "h_pi_distance_cutoff", ParametersDefault.PI_DISTANCE_CUTOFF
            ),
            pi_angle_cutoff=pi_params.get(
                "dh_pi_angle_cutoff", ParametersDefault.PI_ANGLE_CUTOFF
            ),
            # π interaction subtype parameters
            pi_ccl_distance_cutoff=pi_params.get(
                "ccl_pi_distance_cutoff", ParametersDefault.PI_CCL_DISTANCE_CUTOFF
            ),
            pi_ccl_angle_cutoff=pi_params.get(
                "ccl_pi_angle_cutoff", ParametersDefault.PI_CCL_ANGLE_CUTOFF
            ),
            pi_cbr_distance_cutoff=pi_params.get(
                "cbr_pi_distance_cutoff", ParametersDefault.PI_CBR_DISTANCE_CUTOFF
            ),
            pi_cbr_angle_cutoff=pi_params.get(
                "cbr_pi_angle_cutoff", ParametersDefault.PI_CBR_ANGLE_CUTOFF
            ),
            pi_ci_distance_cutoff=pi_params.get(
                "ci_pi_distance_cutoff", ParametersDefault.PI_CI_DISTANCE_CUTOFF
            ),
            pi_ci_angle_cutoff=pi_params.get(
                "ci_pi_angle_cutoff", ParametersDefault.PI_CI_ANGLE_CUTOFF
            ),
            pi_ch_distance_cutoff=pi_params.get(
                "ch_pi_distance_cutoff", ParametersDefault.PI_CH_DISTANCE_CUTOFF
            ),
            pi_ch_angle_cutoff=pi_params.get(
                "ch_pi_angle_cutoff", ParametersDefault.PI_CH_ANGLE_CUTOFF
            ),
            pi_nh_distance_cutoff=pi_params.get(
                "nh_pi_distance_cutoff", ParametersDefault.PI_NH_DISTANCE_CUTOFF
            ),
            pi_nh_angle_cutoff=pi_params.get(
                "nh_pi_angle_cutoff", ParametersDefault.PI_NH_ANGLE_CUTOFF
            ),
            pi_oh_distance_cutoff=pi_params.get(
                "oh_pi_distance_cutoff", ParametersDefault.PI_OH_DISTANCE_CUTOFF
            ),
            pi_oh_angle_cutoff=pi_params.get(
                "oh_pi_angle_cutoff", ParametersDefault.PI_OH_ANGLE_CUTOFF
            ),
            pi_sh_distance_cutoff=pi_params.get(
                "sh_pi_distance_cutoff", ParametersDefault.PI_SH_DISTANCE_CUTOFF
            ),
            pi_sh_angle_cutoff=pi_params.get(
                "sh_pi_angle_cutoff", ParametersDefault.PI_SH_ANGLE_CUTOFF
            ),
            covalent_cutoff_factor=general_params.get(
                "covalent_cutoff_factor", ParametersDefault.COVALENT_CUTOFF_FACTOR
            ),
            analysis_mode=general_params.get(
                "analysis_mode", ParametersDefault.ANALYSIS_MODE
            ),
            # PDB fixing parameters
            fix_pdb_enabled=fix_params.get(
                "enabled", ParametersDefault.FIX_PDB_ENABLED
            ),
            fix_pdb_method=fix_params.get("method", ParametersDefault.FIX_PDB_METHOD),
            fix_pdb_add_hydrogens=fix_params.get(
                "add_hydrogens", ParametersDefault.FIX_PDB_ADD_HYDROGENS
            ),
            fix_pdb_add_heavy_atoms=fix_params.get(
                "add_heavy_atoms", ParametersDefault.FIX_PDB_ADD_HEAVY_ATOMS
            ),
            fix_pdb_replace_nonstandard=fix_params.get(
                "replace_nonstandard", ParametersDefault.FIX_PDB_REPLACE_NONSTANDARD
            ),
            fix_pdb_remove_heterogens=fix_params.get(
                "remove_heterogens", ParametersDefault.FIX_PDB_REMOVE_HETEROGENS
            ),
            fix_pdb_keep_water=fix_params.get(
                "keep_water", ParametersDefault.FIX_PDB_KEEP_WATER
            ),
            # π-π stacking parameters
            pi_pi_distance_cutoff=pi_pi_params.get(
                "distance_cutoff", ParametersDefault.PI_PI_DISTANCE_CUTOFF
            ),
            pi_pi_parallel_angle_cutoff=pi_pi_params.get(
                "parallel_angle_cutoff", ParametersDefault.PI_PI_PARALLEL_ANGLE_CUTOFF
            ),
            pi_pi_tshaped_angle_min=pi_pi_params.get(
                "tshaped_angle_min", ParametersDefault.PI_PI_TSHAPED_ANGLE_MIN
            ),
            pi_pi_tshaped_angle_max=pi_pi_params.get(
                "tshaped_angle_max", ParametersDefault.PI_PI_TSHAPED_ANGLE_MAX
            ),
            pi_pi_offset_cutoff=pi_pi_params.get(
                "offset_cutoff", ParametersDefault.PI_PI_OFFSET_CUTOFF
            ),
            # Carbonyl interaction parameters
            carbonyl_distance_cutoff=carbonyl_params.get(
                "distance_cutoff", ParametersDefault.CARBONYL_DISTANCE_CUTOFF
            ),
            carbonyl_angle_min=carbonyl_params.get(
                "angle_min", ParametersDefault.CARBONYL_ANGLE_MIN
            ),
            carbonyl_angle_max=carbonyl_params.get(
                "angle_max", ParametersDefault.CARBONYL_ANGLE_MAX
            ),
            # n→π* interaction parameters
            n_pi_distance_cutoff=n_pi_params.get(
                "distance_cutoff", ParametersDefault.N_PI_DISTANCE_CUTOFF
            ),
            n_pi_sulfur_distance_cutoff=n_pi_params.get(
                "sulfur_distance_cutoff", ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF
            ),
            n_pi_angle_min=n_pi_params.get(
                "angle_min", ParametersDefault.N_PI_ANGLE_MIN
            ),
            n_pi_angle_max=n_pi_params.get(
                "angle_max", ParametersDefault.N_PI_ANGLE_MAX
            ),
        )

    except Exception as e:
        print_error(f"Failed to load preset file '{preset_path}': {str(e)}")
        sys.exit(1)


def resolve_preset_path(preset_name: str) -> str:
    """Resolve preset name to full path.

    Takes a preset name or partial path and resolves it to a full
    path, searching in the example presets directory if needed.

    :param preset_name: Name or path of the preset to resolve
    :type preset_name: str
    :returns: Full path to the preset file
    :rtype: str
    :raises SystemExit: If preset file cannot be found
    """
    # If it's already a full path, use it
    if os.path.isabs(preset_name) and os.path.exists(preset_name):
        return preset_name

    # If it's a relative path and exists, use it
    if os.path.exists(preset_name):
        return preset_name

    # Try to find it in the example presets directory
    presets_dir = get_example_presets_directory()

    # If preset_name doesn't have extension, try adding .hbat
    if not preset_name.endswith((".hbat", ".json")):
        preset_name += ".hbat"

    preset_path = os.path.join(presets_dir, preset_name)
    if os.path.exists(preset_path):
        return preset_path

    # Try without directory (basename only)
    basename = os.path.basename(preset_name)
    preset_path = os.path.join(presets_dir, basename)
    if os.path.exists(preset_path):
        return preset_path

    print_error(f"Preset file not found: {preset_name}")
    print_error(f"Looked in: {presets_dir}")
    print_error("Use --list-presets to see available presets")
    sys.exit(1)


def load_parameters_from_args(args: argparse.Namespace) -> AnalysisParameters:
    """Create AnalysisParameters from command-line arguments.

    Processes command-line arguments to create analysis parameters,
    with support for preset files and parameter overrides.

    :param args: Parsed command-line arguments
    :type args: argparse.Namespace
    :returns: Analysis parameters configured from arguments
    :rtype: AnalysisParameters
    """
    # If preset is specified, load from preset file first
    if hasattr(args, "preset") and args.preset:
        preset_path = resolve_preset_path(args.preset)
        print_progress(
            f"Loading parameters from preset: {preset_path}",
            not args.quiet if hasattr(args, "quiet") else True,
        )
        params = load_preset_file(preset_path)

        # Override preset parameters with any explicitly provided CLI arguments
        # Only override if the argument was explicitly set (not default)
        parser = create_parser()
        defaults = vars(parser.parse_args([]))  # Get default values

        if args.hb_distance != defaults.get("hb_distance"):
            params.hb_distance_cutoff = args.hb_distance
        if args.hb_angle != defaults.get("hb_angle"):
            params.hb_angle_cutoff = args.hb_angle
        if args.da_distance != defaults.get("da_distance"):
            params.hb_donor_acceptor_cutoff = args.da_distance
        if args.whb_distance != defaults.get("whb_distance"):
            params.whb_distance_cutoff = args.whb_distance
        if args.whb_angle != defaults.get("whb_angle"):
            params.whb_angle_cutoff = args.whb_angle
        if args.whb_da_distance != defaults.get("whb_da_distance"):
            params.whb_donor_acceptor_cutoff = args.whb_da_distance
        if args.xb_distance != defaults.get("xb_distance"):
            params.xb_distance_cutoff = args.xb_distance
        if args.xb_angle != defaults.get("xb_angle"):
            params.xb_angle_cutoff = args.xb_angle
        if args.pi_distance != defaults.get("pi_distance"):
            params.pi_distance_cutoff = args.pi_distance
        if args.pi_angle != defaults.get("pi_angle"):
            params.pi_angle_cutoff = args.pi_angle

        # π interaction subtype parameter overrides
        if args.pi_ccl_distance != defaults.get("pi_ccl_distance"):
            params.pi_ccl_distance_cutoff = args.pi_ccl_distance
        if args.pi_ccl_angle != defaults.get("pi_ccl_angle"):
            params.pi_ccl_angle_cutoff = args.pi_ccl_angle
        if args.pi_cbr_distance != defaults.get("pi_cbr_distance"):
            params.pi_cbr_distance_cutoff = args.pi_cbr_distance
        if args.pi_cbr_angle != defaults.get("pi_cbr_angle"):
            params.pi_cbr_angle_cutoff = args.pi_cbr_angle
        if args.pi_ci_distance != defaults.get("pi_ci_distance"):
            params.pi_ci_distance_cutoff = args.pi_ci_distance
        if args.pi_ci_angle != defaults.get("pi_ci_angle"):
            params.pi_ci_angle_cutoff = args.pi_ci_angle
        if args.pi_ch_distance != defaults.get("pi_ch_distance"):
            params.pi_ch_distance_cutoff = args.pi_ch_distance
        if args.pi_ch_angle != defaults.get("pi_ch_angle"):
            params.pi_ch_angle_cutoff = args.pi_ch_angle
        if args.pi_nh_distance != defaults.get("pi_nh_distance"):
            params.pi_nh_distance_cutoff = args.pi_nh_distance
        if args.pi_nh_angle != defaults.get("pi_nh_angle"):
            params.pi_nh_angle_cutoff = args.pi_nh_angle
        if args.pi_oh_distance != defaults.get("pi_oh_distance"):
            params.pi_oh_distance_cutoff = args.pi_oh_distance
        if args.pi_oh_angle != defaults.get("pi_oh_angle"):
            params.pi_oh_angle_cutoff = args.pi_oh_angle
        if args.pi_sh_distance != defaults.get("pi_sh_distance"):
            params.pi_sh_distance_cutoff = args.pi_sh_distance
        if args.pi_sh_angle != defaults.get("pi_sh_angle"):
            params.pi_sh_angle_cutoff = args.pi_sh_angle

        # π-π stacking parameter overrides
        if args.pi_pi_distance != defaults.get("pi_pi_distance"):
            params.pi_pi_distance_cutoff = args.pi_pi_distance
        if args.pi_pi_parallel_angle != defaults.get("pi_pi_parallel_angle"):
            params.pi_pi_parallel_angle_cutoff = args.pi_pi_parallel_angle
        if args.pi_pi_tshaped_angle_min != defaults.get("pi_pi_tshaped_angle_min"):
            params.pi_pi_tshaped_angle_min = args.pi_pi_tshaped_angle_min
        if args.pi_pi_tshaped_angle_max != defaults.get("pi_pi_tshaped_angle_max"):
            params.pi_pi_tshaped_angle_max = args.pi_pi_tshaped_angle_max
        if args.pi_pi_offset != defaults.get("pi_pi_offset"):
            params.pi_pi_offset_cutoff = args.pi_pi_offset

        # Carbonyl interaction parameter overrides
        if args.carbonyl_distance != defaults.get("carbonyl_distance"):
            params.carbonyl_distance_cutoff = args.carbonyl_distance
        if args.carbonyl_angle_min != defaults.get("carbonyl_angle_min"):
            params.carbonyl_angle_min = args.carbonyl_angle_min
        if args.carbonyl_angle_max != defaults.get("carbonyl_angle_max"):
            params.carbonyl_angle_max = args.carbonyl_angle_max

        # n→π* interaction parameter overrides
        if args.n_pi_distance != defaults.get("n_pi_distance"):
            params.n_pi_distance_cutoff = args.n_pi_distance
        if args.n_pi_sulfur_distance != defaults.get("n_pi_sulfur_distance"):
            params.n_pi_sulfur_distance_cutoff = args.n_pi_sulfur_distance
        if args.n_pi_angle_min != defaults.get("n_pi_angle_min"):
            params.n_pi_angle_min = args.n_pi_angle_min
        if args.n_pi_angle_max != defaults.get("n_pi_angle_max"):
            params.n_pi_angle_max = args.n_pi_angle_max

        if args.covalent_factor != defaults.get("covalent_factor"):
            params.covalent_cutoff_factor = args.covalent_factor
        if args.mode != defaults.get("mode"):
            params.analysis_mode = args.mode

        return params
    else:
        # Use CLI arguments only
        return AnalysisParameters(
            hb_distance_cutoff=args.hb_distance,
            hb_angle_cutoff=args.hb_angle,
            hb_donor_acceptor_cutoff=args.da_distance,
            whb_distance_cutoff=args.whb_distance,
            whb_angle_cutoff=args.whb_angle,
            whb_donor_acceptor_cutoff=args.whb_da_distance,
            xb_distance_cutoff=args.xb_distance,
            xb_angle_cutoff=args.xb_angle,
            pi_distance_cutoff=args.pi_distance,
            pi_angle_cutoff=args.pi_angle,
            # π interaction subtype parameters
            pi_ccl_distance_cutoff=args.pi_ccl_distance,
            pi_ccl_angle_cutoff=args.pi_ccl_angle,
            pi_cbr_distance_cutoff=args.pi_cbr_distance,
            pi_cbr_angle_cutoff=args.pi_cbr_angle,
            pi_ci_distance_cutoff=args.pi_ci_distance,
            pi_ci_angle_cutoff=args.pi_ci_angle,
            pi_ch_distance_cutoff=args.pi_ch_distance,
            pi_ch_angle_cutoff=args.pi_ch_angle,
            pi_nh_distance_cutoff=args.pi_nh_distance,
            pi_nh_angle_cutoff=args.pi_nh_angle,
            pi_oh_distance_cutoff=args.pi_oh_distance,
            pi_oh_angle_cutoff=args.pi_oh_angle,
            pi_sh_distance_cutoff=args.pi_sh_distance,
            pi_sh_angle_cutoff=args.pi_sh_angle,
            # π-π stacking parameters
            pi_pi_distance_cutoff=args.pi_pi_distance,
            pi_pi_parallel_angle_cutoff=args.pi_pi_parallel_angle,
            pi_pi_tshaped_angle_min=args.pi_pi_tshaped_angle_min,
            pi_pi_tshaped_angle_max=args.pi_pi_tshaped_angle_max,
            pi_pi_offset_cutoff=args.pi_pi_offset,
            # Carbonyl interaction parameters
            carbonyl_distance_cutoff=args.carbonyl_distance,
            carbonyl_angle_min=args.carbonyl_angle_min,
            carbonyl_angle_max=args.carbonyl_angle_max,
            # n→π* interaction parameters
            n_pi_distance_cutoff=args.n_pi_distance,
            n_pi_sulfur_distance_cutoff=args.n_pi_sulfur_distance,
            n_pi_angle_min=args.n_pi_angle_min,
            n_pi_angle_max=args.n_pi_angle_max,
            covalent_cutoff_factor=args.covalent_factor,
            analysis_mode=args.mode,
            # PDB fixing parameters
            fix_pdb_enabled=args.fix_pdb,
            fix_pdb_method=args.fix_method,
            fix_pdb_add_hydrogens=args.fix_add_hydrogens,
            fix_pdb_add_heavy_atoms=args.fix_add_heavy_atoms,
            fix_pdb_replace_nonstandard=args.fix_replace_nonstandard,
            fix_pdb_remove_heterogens=args.fix_remove_heterogens,
            fix_pdb_keep_water=args.fix_keep_water,
        )


def print_progress(message: str, verbose: bool = True) -> None:
    """Print progress message if verbose mode enabled.

    :param message: Progress message to display
    :type message: str
    :param verbose: Whether to actually print the message
    :type verbose: bool
    :returns: None
    :rtype: None
    """
    if verbose:
        print(f"[INFO] {message}")


def print_error(message: str) -> None:
    """Print error message to stderr.

    :param message: Error message to display
    :type message: str
    :returns: None
    :rtype: None
    """
    print(f"[ERROR] {message}", file=sys.stderr)


def validate_input_file(filename: str) -> bool:
    """Validate input PDB file.

    Checks if the input file exists, is readable, and contains
    valid PDB-format content.

    :param filename: Path to the PDB file to validate
    :type filename: str
    :returns: True if file is valid, False otherwise
    :rtype: bool
    """
    if not os.path.exists(filename):
        print_error(f"Input file '{filename}' not found")
        return False

    if not os.path.isfile(filename):
        print_error(f"'{filename}' is not a regular file")
        return False

    try:
        with open(filename, "r") as f:
            # Check if file contains PDB-like content
            # Read more lines to account for headers
            content = f.read()
            has_atoms = "ATOM" in content or "HETATM" in content
            has_pdb_keywords = any(
                keyword in content
                for keyword in ["HEADER", "TITLE", "COMPND", "ATOM", "HETATM"]
            )
            if not (has_atoms or has_pdb_keywords):
                print_error(f"'{filename}' does not appear to be a valid PDB file")
                return False
    except Exception as e:
        print_error(f"Cannot read input file: {e}")
        return False

    return True


def format_results_text(
    analyzer: NPMolecularInteractionAnalyzer,
    input_file: str,
    summary_only: bool = False,
) -> str:
    """Format analysis results as text.

    Creates a human-readable text report of analysis results,
    with options for summary or detailed output.

    :param analyzer: Analysis results to format
    :type analyzer: MolecularInteractionAnalyzer
    :param input_file: Path to the input file analyzed
    :type input_file: str
    :param summary_only: Whether to include only summary statistics
    :type summary_only: bool
    :returns: Formatted text report
    :rtype: str
    """
    lines = []
    lines.append("HBAT Analysis Results")
    lines.append("=" * 50)
    lines.append(f"Input file: {input_file}")
    lines.append(f"Analysis time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    # Statistics summary
    summary = analyzer.get_summary()
    lines.append("Summary:")
    lines.append(f"  Hydrogen bonds: {summary['hydrogen_bonds']['count']}")
    lines.append(f"  Halogen bonds: {summary['halogen_bonds']['count']}")
    lines.append(f"  π interactions: {summary['pi_interactions']['count']}")

    # Add new interaction types if they exist
    if "pi_pi_stacking" in summary:
        lines.append(f"  π-π stacking: {summary['pi_pi_stacking']['count']}")
    elif hasattr(analyzer, "pi_pi_interactions"):
        lines.append(f"  π-π stacking: {len(analyzer.pi_pi_interactions)}")

    if "carbonyl_interactions" in summary:
        lines.append(
            f"  Carbonyl interactions: {summary['carbonyl_interactions']['count']}"
        )
    elif hasattr(analyzer, "carbonyl_interactions"):
        lines.append(f"  Carbonyl interactions: {len(analyzer.carbonyl_interactions)}")

    if "n_pi_interactions" in summary:
        lines.append(f"  n→π* interactions: {summary['n_pi_interactions']['count']}")
    elif hasattr(analyzer, "n_pi_interactions"):
        lines.append(f"  n→π* interactions: {len(analyzer.n_pi_interactions)}")

    lines.append(f"  Cooperativity chains: {summary['cooperativity_chains']['count']}")
    lines.append(f"  Total interactions: {summary['total_interactions']}")
    lines.append("")

    # Bond detection statistics
    if "bond_detection" in summary:
        bond_stats = summary["bond_detection"]
        lines.append("Bond Detection:")
        lines.append(f"  Total bonds detected: {bond_stats['total_bonds']}")
        if bond_stats["breakdown"]:
            for method, stats in bond_stats["breakdown"].items():
                method_name = method.replace("_", " ").title()
                lines.append(
                    f"    {method_name}: {stats['count']} ({stats['percentage']}%)"
                )
        lines.append("")

    if summary_only:
        return "\n".join(lines)

    # Detailed results
    if analyzer.hydrogen_bonds:
        lines.append("Hydrogen Bonds:")
        lines.append("-" * 30)
        for i, hb in enumerate(analyzer.hydrogen_bonds, 1):
            lines.append(f"{i:3d}. {hb}")
        lines.append("")

    if analyzer.halogen_bonds:
        lines.append("Halogen Bonds:")
        lines.append("-" * 30)
        for i, xb in enumerate(analyzer.halogen_bonds, 1):
            lines.append(f"{i:3d}. {xb}")
        lines.append("")

    if analyzer.pi_interactions:
        lines.append("π Interactions:")
        lines.append("-" * 30)
        for i, pi in enumerate(analyzer.pi_interactions, 1):
            lines.append(f"{i:3d}. {pi}")
        lines.append("")

    # Add new interaction type details
    if hasattr(analyzer, "pi_pi_interactions") and analyzer.pi_pi_interactions:
        lines.append("π-π Stacking Interactions:")
        lines.append("-" * 35)
        for i, pi_pi in enumerate(analyzer.pi_pi_interactions, 1):
            lines.append(f"{i:3d}. {pi_pi}")
        lines.append("")

    if hasattr(analyzer, "carbonyl_interactions") and analyzer.carbonyl_interactions:
        lines.append("Carbonyl-Carbonyl Interactions:")
        lines.append("-" * 35)
        for i, carbonyl in enumerate(analyzer.carbonyl_interactions, 1):
            lines.append(f"{i:3d}. {carbonyl}")
        lines.append("")

    if hasattr(analyzer, "n_pi_interactions") and analyzer.n_pi_interactions:
        lines.append("n→π* Interactions:")
        lines.append("-" * 20)
        for i, n_pi in enumerate(analyzer.n_pi_interactions, 1):
            lines.append(f"{i:3d}. {n_pi}")
        lines.append("")

    if analyzer.cooperativity_chains:
        lines.append("Cooperativity Chains:")
        lines.append("-" * 30)
        for i, chain in enumerate(analyzer.cooperativity_chains, 1):
            lines.append(f"{i:3d}. {chain}")
        lines.append("")

    return "\n".join(lines)


def detect_output_format(filename: str) -> str:
    """Detect output format from file extension.

    :param filename: Output filename
    :type filename: str
    :returns: Format type ('text', 'json')
    :rtype: str
    :raises ValueError: If file extension is not supported
    """
    import os

    _, ext = os.path.splitext(filename)
    ext_lower = ext.lower()

    if ext_lower == ".txt":
        return "text"
    elif ext_lower == ".csv":
        raise ValueError(
            "Single CSV file output is not supported. "
            "Use --csv flag to export multiple CSV files (one per interaction type).\n"
            "Example: hbat input.pdb --csv output"
        )
    elif ext_lower == ".json":
        return "json"
    else:
        raise ValueError(
            f"Unsupported output format '{ext}'. Use .txt or .json for single file output"
        )


def run_analysis(args: argparse.Namespace) -> int:
    """Run the analysis with given arguments.

    Performs the complete analysis workflow including parameter loading,
    analysis execution, and result output based on command-line arguments.

    :param args: Parsed command-line arguments
    :type args: argparse.Namespace
    :returns: Exit code (0 for success, non-zero for failure)
    :rtype: int
    """
    # Validate input
    if not validate_input_file(args.input):
        return 1

    verbose = args.verbose and not args.quiet

    try:
        # Load parameters
        parameters = load_parameters_from_args(args)

        print_progress(f"Starting analysis of {args.input}", verbose)
        print_progress(f"Analysis mode: {parameters.analysis_mode}", verbose)

        # Create analyzer
        analyzer = NPMolecularInteractionAnalyzer(parameters)

        # Set up progress tracking (show unless quiet mode)
        progress_bar = None
        if not args.quiet:
            progress_bar = ProgressBar()

            def cli_progress_callback(message: str) -> None:
                """Progress callback for CLI updates."""
                # Use progress bar for clean display
                if progress_bar:
                    # Extract percentage if present in message
                    if "%" in message:
                        try:
                            # Split on space to get the percentage part
                            parts = message.split()
                            percent_part = None
                            for part in parts:
                                if "%" in part:
                                    percent_part = part.replace("%", "")
                                    break

                            if percent_part and percent_part.isdigit():
                                progress = int(percent_part)
                                # Get the step name (everything before the percentage)
                                step_name = message.split(f" {percent_part}%")[0]
                                progress_bar.update(step_name, progress)
                            else:
                                progress_bar.update(message)
                        except (ValueError, IndexError):
                            progress_bar.update(message)
                    else:
                        progress_bar.update(message)
                else:
                    # Fallback if no progress bar
                    print_progress(message, verbose)

            analyzer.progress_callback = cli_progress_callback

        # Run analysis
        start_time = time.time()
        success = analyzer.analyze_file(args.input)
        analysis_time = time.time() - start_time

        if not success:
            if progress_bar:
                progress_bar.finish("Analysis failed")
            print_error("Analysis failed")
            return 1

        if progress_bar:
            progress_bar.finish(f"Analysis completed in {analysis_time:.2f} seconds")
        else:
            print_progress(
                f"Analysis completed in {analysis_time:.2f} seconds", verbose
            )

        # Apply analysis filters (clear results for disabled interaction types)
        if hasattr(args, "no_hydrogen_bonds") and args.no_hydrogen_bonds:
            analyzer.hydrogen_bonds = []
        if hasattr(args, "no_halogen_bonds") and args.no_halogen_bonds:
            analyzer.halogen_bonds = []
        if hasattr(args, "no_pi_interactions") and args.no_pi_interactions:
            analyzer.pi_interactions = []
        if hasattr(args, "no_pi_pi_stacking") and args.no_pi_pi_stacking:
            analyzer.pi_pi_interactions = []
        if hasattr(args, "no_carbonyl_interactions") and args.no_carbonyl_interactions:
            analyzer.carbonyl_interactions = []
        if hasattr(args, "no_n_pi_interactions") and args.no_n_pi_interactions:
            analyzer.n_pi_interactions = []

        # Get results
        summary = analyzer.get_summary()

        if not args.quiet:
            # Get counts for all interaction types
            hb_count = summary["hydrogen_bonds"]["count"]
            xb_count = summary["halogen_bonds"]["count"]
            pi_count = summary["pi_interactions"]["count"]
            pi_pi_count = summary.get("pi_pi_interactions", {}).get("count", 0)
            carbonyl_count = summary.get("carbonyl_interactions", {}).get("count", 0)
            n_pi_count = summary.get("n_pi_interactions", {}).get("count", 0)
            chains_count = summary["cooperativity_chains"]["count"]

            print(
                f"Found {hb_count} hydrogen bonds, "
                f"{xb_count} halogen bonds, "
                f"{pi_count} π interactions, "
                f"{pi_pi_count} π-π stacking, "
                f"{carbonyl_count} carbonyl n→π*, "
                f"{n_pi_count} n→π* interactions, "
                f"{chains_count} cooperativity chains"
            )

        # Output results
        if args.output:
            try:
                output_format = detect_output_format(args.output)
                if output_format == "text":
                    print_progress(f"Writing results to {args.output}", verbose)
                    export_to_txt_single_file(analyzer, args.output)
                elif output_format == "json":
                    print_progress(f"Exporting to JSON: {args.output}", verbose)
                    export_to_json_single_file(analyzer, args.output, args.input)
            except ValueError as e:
                print_error(str(e))
                return 1

        if args.json:
            print_progress(f"Exporting to multiple JSON files: {args.json}", verbose)
            export_to_json_files(analyzer, args.json, args.input)

        if args.csv:
            print_progress(f"Exporting to multiple CSV files: {args.csv}", verbose)
            export_to_csv_files(analyzer, args.csv)

        # Print to stdout if no output files specified
        if not any([args.output, args.json, args.csv]) and not args.quiet:
            print("\n" + format_results_text(analyzer, args.input, args.summary_only))

        return 0

    except KeyboardInterrupt:
        print_error("Analysis interrupted by user")
        return 130
    except Exception as e:
        print_error(f"Analysis failed: {e}")
        if verbose:
            import traceback

            traceback.print_exc()
        return 1


def main() -> int:
    """Main CLI entry point for HBAT molecular interaction analysis.

    Parses command-line arguments and dispatches to appropriate functionality.
    Supports comprehensive analysis of molecular interactions including:

    - Hydrogen bonds (classical N-H···O, O-H···O)
    - Weak hydrogen bonds (C-H···O interactions)
    - Halogen bonds (C-X···A with default 150° angle cutoff)
    - π interactions with multiple subtypes:
      • Hydrogen-π: C-H···π, N-H···π, O-H···π, S-H···π
      • Halogen-π: C-Cl···π, C-Br···π, C-I···π
    - Cooperativity chains and interaction networks

    Includes built-in parameter presets and PDB structure fixing capabilities.

    :returns: Exit code (0 for success, non-zero for failure)
    :rtype: int
    """
    # Initialize HBAT environment first
    try:
        from ..core.app_config import initialize_hbat_environment

        initialize_hbat_environment(
            verbose=False
        )  # We'll handle verbosity based on args
    except ImportError:
        pass  # Continue without app config if import fails

    parser = create_parser()
    args = parser.parse_args()

    # Show HBAT environment info if verbose
    if hasattr(args, "verbose") and args.verbose:
        try:
            from ..core.app_config import get_hbat_config

            config = get_hbat_config()
            info = config.get_info()
            print(f"📁 HBAT data directory: {info['hbat_directory']}")
            if info["ccd_files_present"]:
                print(f"✅ CCD data available")
            else:
                print(f"⚠️  CCD data will be downloaded as needed")
        except ImportError:
            pass

    # Handle preset listing first
    if hasattr(args, "list_presets") and args.list_presets:
        list_available_presets()
        return 0

    # Validate that input file is provided for analysis
    if not args.input:
        print_error("Input PDB file is required for analysis")
        print_error(
            "Use --help for usage information or --list-presets to see available presets"
        )
        return 1

    # Handle conflicting options
    if args.verbose and args.quiet:
        print_error("Cannot use both --verbose and --quiet options")
        return 1

    return run_analysis(args)


if __name__ == "__main__":
    sys.exit(main())
