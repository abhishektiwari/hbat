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
from ..constants import AnalysisDefaults
from ..core.analysis import AnalysisParameters, HBondAnalyzer
from ..core.pdb_parser import PDBParser


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
  %(prog)s input.pdb                          # Basic analysis
  %(prog)s input.pdb -o results.txt          # Save results to file
  %(prog)s input.pdb --hb-distance 3.0       # Custom H-bond distance cutoff
  %(prog)s input.pdb --mode local             # Local interactions only
  %(prog)s input.pdb --json results.json     # Export to JSON format
  %(prog)s --list-presets                     # List available presets
  %(prog)s input.pdb --preset high_resolution # Use preset with custom overrides
  %(prog)s input.pdb --preset drug_design_strict --hb-distance 3.0
        """,
    )

    # Version
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    # Input file (optional when listing presets)
    parser.add_argument("input", nargs="?", help="Input PDB file")

    # Output options
    parser.add_argument("-o", "--output", help="Output text file for results")
    parser.add_argument("--json", help="Output JSON file for structured results")
    parser.add_argument("--csv", help="Output CSV file for tabular results")

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
        default=AnalysisDefaults.HB_DISTANCE_CUTOFF,
        help=f"Hydrogen bond H...A distance cutoff in Å (default: {AnalysisDefaults.HB_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--hb-angle",
        type=float,
        default=AnalysisDefaults.HB_ANGLE_CUTOFF,
        help=f"Hydrogen bond D-H...A angle cutoff in degrees (default: {AnalysisDefaults.HB_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--da-distance",
        type=float,
        default=AnalysisDefaults.HB_DA_DISTANCE,
        help=f"Donor-acceptor distance cutoff in Å (default: {AnalysisDefaults.HB_DA_DISTANCE})",
    )
    param_group.add_argument(
        "--xb-distance",
        type=float,
        default=AnalysisDefaults.XB_DISTANCE_CUTOFF,
        help=f"Halogen bond X...A distance cutoff in Å (default: {AnalysisDefaults.XB_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--xb-angle",
        type=float,
        default=AnalysisDefaults.XB_ANGLE_CUTOFF,
        help=f"Halogen bond C-X...A angle cutoff in degrees (default: {AnalysisDefaults.XB_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-distance",
        type=float,
        default=AnalysisDefaults.PI_DISTANCE_CUTOFF,
        help=f"π interaction H...π distance cutoff in Å (default: {AnalysisDefaults.PI_DISTANCE_CUTOFF})",
    )
    param_group.add_argument(
        "--pi-angle",
        type=float,
        default=AnalysisDefaults.PI_ANGLE_CUTOFF,
        help=f"π interaction D-H...π angle cutoff in degrees (default: {AnalysisDefaults.PI_ANGLE_CUTOFF})",
    )
    param_group.add_argument(
        "--covalent-factor",
        type=float,
        default=AnalysisDefaults.COVALENT_CUTOFF_FACTOR,
        help=f"Covalent bond detection factor (default: {AnalysisDefaults.COVALENT_CUTOFF_FACTOR})",
    )

    # Analysis mode
    param_group.add_argument(
        "--mode",
        choices=["complete", "local"],
        default=AnalysisDefaults.ANALYSIS_MODE,
        help="Analysis mode: complete (all interactions) or local (intra-residue only)",
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
        xb_params = params.get("halogen_bonds", {})
        pi_params = params.get("pi_interactions", {})
        general_params = params.get("general", {})

        return AnalysisParameters(
            hb_distance_cutoff=hb_params.get(
                "h_a_distance_cutoff", AnalysisDefaults.HB_DISTANCE_CUTOFF
            ),
            hb_angle_cutoff=hb_params.get(
                "dha_angle_cutoff", AnalysisDefaults.HB_ANGLE_CUTOFF
            ),
            hb_donor_acceptor_cutoff=hb_params.get(
                "d_a_distance_cutoff", AnalysisDefaults.HB_DA_DISTANCE
            ),
            xb_distance_cutoff=xb_params.get(
                "x_a_distance_cutoff", AnalysisDefaults.XB_DISTANCE_CUTOFF
            ),
            xb_angle_cutoff=xb_params.get(
                "cxa_angle_cutoff", AnalysisDefaults.XB_ANGLE_CUTOFF
            ),
            pi_distance_cutoff=pi_params.get(
                "h_pi_distance_cutoff", AnalysisDefaults.PI_DISTANCE_CUTOFF
            ),
            pi_angle_cutoff=pi_params.get(
                "dh_pi_angle_cutoff", AnalysisDefaults.PI_ANGLE_CUTOFF
            ),
            covalent_cutoff_factor=general_params.get(
                "covalent_cutoff_factor", AnalysisDefaults.COVALENT_CUTOFF_FACTOR
            ),
            analysis_mode=general_params.get(
                "analysis_mode", AnalysisDefaults.ANALYSIS_MODE
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
        if args.xb_distance != defaults.get("xb_distance"):
            params.xb_distance_cutoff = args.xb_distance
        if args.xb_angle != defaults.get("xb_angle"):
            params.xb_angle_cutoff = args.xb_angle
        if args.pi_distance != defaults.get("pi_distance"):
            params.pi_distance_cutoff = args.pi_distance
        if args.pi_angle != defaults.get("pi_angle"):
            params.pi_angle_cutoff = args.pi_angle
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
            xb_distance_cutoff=args.xb_distance,
            xb_angle_cutoff=args.xb_angle,
            pi_distance_cutoff=args.pi_distance,
            pi_angle_cutoff=args.pi_angle,
            covalent_cutoff_factor=args.covalent_factor,
            analysis_mode=args.mode,
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
    analyzer: HBondAnalyzer, input_file: str, summary_only: bool = False
) -> str:
    """Format analysis results as text.

    Creates a human-readable text report of analysis results,
    with options for summary or detailed output.

    :param analyzer: Analysis results to format
    :type analyzer: HBondAnalyzer
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
    stats = analyzer.get_statistics()
    lines.append("Summary:")
    lines.append(f"  Hydrogen bonds: {stats['hydrogen_bonds']}")
    lines.append(f"  Halogen bonds: {stats['halogen_bonds']}")
    lines.append(f"  π interactions: {stats['pi_interactions']}")
    lines.append(f"  Cooperativity chains: {stats.get('cooperativity_chains', 0)}")
    lines.append(f"  Total interactions: {stats['total_interactions']}")
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

    if analyzer.cooperativity_chains:
        lines.append("Cooperativity Chains:")
        lines.append("-" * 30)
        for i, chain in enumerate(analyzer.cooperativity_chains, 1):
            lines.append(f"{i:3d}. {chain}")
        lines.append("")

    return "\n".join(lines)


def export_to_json(analyzer: HBondAnalyzer, input_file: str, output_file: str) -> None:
    """Export results to JSON format.

    Exports complete analysis results to a structured JSON file
    with metadata, statistics, and detailed interaction data.

    :param analyzer: Analysis results to export
    :type analyzer: HBondAnalyzer
    :param input_file: Path to the input file analyzed
    :type input_file: str
    :param output_file: Path to the JSON output file
    :type output_file: str
    :returns: None
    :rtype: None
    """
    import math

    data: Dict[str, Any] = {
        "metadata": {
            "input_file": input_file,
            "analysis_time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "hbat_version": __version__,
        },
        "statistics": analyzer.get_statistics(),
        "hydrogen_bonds": [],
        "halogen_bonds": [],
        "pi_interactions": [],
        "cooperativity_chains": [],
    }

    # Convert hydrogen bonds
    for hb in analyzer.hydrogen_bonds:
        data["hydrogen_bonds"].append(
            {
                "donor_residue": hb.donor_residue,
                "donor_atom": hb.donor.name,
                "donor_coords": hb.donor.coords.to_list(),
                "hydrogen_atom": hb.hydrogen.name,
                "hydrogen_coords": hb.hydrogen.coords.to_list(),
                "acceptor_residue": hb.acceptor_residue,
                "acceptor_atom": hb.acceptor.name,
                "acceptor_coords": hb.acceptor.coords.to_list(),
                "distance": round(hb.distance, 3),
                "angle": round(math.degrees(hb.angle), 1),
                "donor_acceptor_distance": round(hb.donor_acceptor_distance, 3),
                "bond_type": hb.bond_type,
            }
        )

    # Convert halogen bonds
    for xb in analyzer.halogen_bonds:
        data["halogen_bonds"].append(
            {
                "halogen_residue": xb.halogen_residue,
                "halogen_atom": xb.halogen.name,
                "halogen_coords": xb.halogen.coords.to_list(),
                "acceptor_residue": xb.acceptor_residue,
                "acceptor_atom": xb.acceptor.name,
                "acceptor_coords": xb.acceptor.coords.to_list(),
                "distance": round(xb.distance, 3),
                "angle": round(math.degrees(xb.angle), 1),
                "bond_type": xb.bond_type,
            }
        )

    # Convert π interactions
    for pi in analyzer.pi_interactions:
        data["pi_interactions"].append(
            {
                "donor_residue": pi.donor_residue,
                "donor_atom": pi.donor.name,
                "donor_coords": pi.donor.coords.to_list(),
                "hydrogen_atom": pi.hydrogen.name,
                "hydrogen_coords": pi.hydrogen.coords.to_list(),
                "pi_residue": pi.pi_residue,
                "pi_center": pi.pi_center.to_list(),
                "distance": round(pi.distance, 3),
                "angle": round(math.degrees(pi.angle), 1),
            }
        )

    # Convert cooperativity chains
    for chain in analyzer.cooperativity_chains:
        chain_data: Dict[str, Any] = {
            "chain_length": chain.chain_length,
            "chain_type": chain.chain_type,
            "interactions": [],
        }

        for interaction in chain.interactions:
            interaction_data: Dict[str, Any] = {
                "interaction_type": interaction.interaction_type,
                "donor_residue": interaction.get_donor_residue(),
                "acceptor_residue": interaction.get_acceptor_residue(),
                "distance": round(interaction.distance, 3),
                "angle": round(math.degrees(interaction.angle), 1),
            }

            # Add type-specific fields
            if hasattr(interaction, "bond_type"):
                interaction_data["bond_type"] = interaction.bond_type

            donor_atom = interaction.get_donor_atom()
            if donor_atom:
                interaction_data["donor_atom"] = donor_atom.name

            acceptor_atom = interaction.get_acceptor_atom()
            if acceptor_atom:
                interaction_data["acceptor_atom"] = acceptor_atom.name

            chain_data["interactions"].append(interaction_data)

        data["cooperativity_chains"].append(chain_data)

    with open(output_file, "w") as f:
        json.dump(data, f, indent=2)


def export_to_csv(analyzer: HBondAnalyzer, output_file: str) -> None:
    """Export results to CSV format.

    Exports analysis results to a CSV file with separate sections
    for different interaction types.

    :param analyzer: Analysis results to export
    :type analyzer: HBondAnalyzer
    :param output_file: Path to the CSV output file
    :type output_file: str
    :returns: None
    :rtype: None
    """
    import csv
    import math

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)

        # Hydrogen bonds section
        if analyzer.hydrogen_bonds:
            writer.writerow(["# Hydrogen Bonds"])
            writer.writerow(
                [
                    "Donor_Residue",
                    "Donor_Atom",
                    "Hydrogen_Atom",
                    "Acceptor_Residue",
                    "Acceptor_Atom",
                    "Distance_A",
                    "Angle_deg",
                    "DA_Distance_A",
                    "Bond_Type",
                ]
            )

            for hb in analyzer.hydrogen_bonds:
                writer.writerow(
                    [
                        hb.donor_residue,
                        hb.donor.name,
                        hb.hydrogen.name,
                        hb.acceptor_residue,
                        hb.acceptor.name,
                        f"{hb.distance:.3f}",
                        f"{math.degrees(hb.angle):.1f}",
                        f"{hb.donor_acceptor_distance:.3f}",
                        hb.bond_type,
                    ]
                )
            writer.writerow([])  # Empty row

        # Halogen bonds section
        if analyzer.halogen_bonds:
            writer.writerow(["# Halogen Bonds"])
            writer.writerow(
                [
                    "Halogen_Residue",
                    "Halogen_Atom",
                    "Acceptor_Residue",
                    "Acceptor_Atom",
                    "Distance_A",
                    "Angle_deg",
                    "Bond_Type",
                ]
            )

            for xb in analyzer.halogen_bonds:
                writer.writerow(
                    [
                        xb.halogen_residue,
                        xb.halogen.name,
                        xb.acceptor_residue,
                        xb.acceptor.name,
                        f"{xb.distance:.3f}",
                        f"{math.degrees(xb.angle):.1f}",
                        xb.bond_type,
                    ]
                )
            writer.writerow([])  # Empty row

        # π interactions section
        if analyzer.pi_interactions:
            writer.writerow(["# Pi Interactions"])
            writer.writerow(
                [
                    "Donor_Residue",
                    "Donor_Atom",
                    "Hydrogen_Atom",
                    "Pi_Residue",
                    "Distance_A",
                    "Angle_deg",
                ]
            )

            for pi in analyzer.pi_interactions:
                writer.writerow(
                    [
                        pi.donor_residue,
                        pi.donor.name,
                        pi.hydrogen.name,
                        pi.pi_residue,
                        f"{pi.distance:.3f}",
                        f"{math.degrees(pi.angle):.1f}",
                    ]
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
        analyzer = HBondAnalyzer(parameters)

        # Run analysis
        start_time = time.time()
        success = analyzer.analyze_file(args.input)
        analysis_time = time.time() - start_time

        if not success:
            print_error("Analysis failed")
            return 1

        print_progress(f"Analysis completed in {analysis_time:.2f} seconds", verbose)

        # Get results
        stats = analyzer.get_statistics()

        if not args.quiet:
            print(
                f"Found {stats['hydrogen_bonds']} hydrogen bonds, "
                f"{stats['halogen_bonds']} halogen bonds, "
                f"{stats['pi_interactions']} π interactions, "
                f"{stats.get('cooperativity_chains', 0)} cooperativity chains"
            )

        # Output results
        if args.output:
            print_progress(f"Writing results to {args.output}", verbose)
            with open(args.output, "w") as f:
                f.write(format_results_text(analyzer, args.input, args.summary_only))

        if args.json:
            print_progress(f"Exporting to JSON: {args.json}", verbose)
            export_to_json(analyzer, args.input, args.json)

        if args.csv:
            print_progress(f"Exporting to CSV: {args.csv}", verbose)
            export_to_csv(analyzer, args.csv)

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
    """Main CLI entry point.

    Parses command-line arguments and dispatches to appropriate
    functionality (preset listing or analysis execution).

    :returns: Exit code (0 for success, non-zero for failure)
    :rtype: int
    """
    parser = create_parser()
    args = parser.parse_args()

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
