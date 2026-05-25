"""
Centralized export functions for HBAT analysis results.

This module provides functions to export molecular interaction analysis
results to CSV and JSON formats. These functions are used by both the
CLI and GUI interfaces.
"""

import csv
import io
import json
import math
from pathlib import Path
from typing import Optional, Dict, Any

from hbat import __version__
from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
from hbat.config import (
    INTERACTION_CONFIGS,
    extract_interaction_data,
    get_interaction_config,
)


def export_to_txt_single_file(
    analyzer: NPMolecularInteractionAnalyzer, output_file: str
) -> None:
    """Export all interactions to a single text file with human-readable format.

    Creates a text file containing a summary and detailed listing of all
    interaction types found in the analysis, including ligand interactions
    organized by ligand. Uses same organization as export_to_json_single_file(),
    just with text formatting.

    :param analyzer: Analyzer instance with interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param output_file: Path to the output text file
    :type output_file: str
    :returns: None
    :rtype: None
    """
    # Get organized interactions using shared helper
    organized = _organize_all_interactions_by_type(analyzer)

    with open(output_file, "w", encoding="utf-8") as f:
        # Write summary
        summary = analyzer.get_summary()
        f.write("Summary:\n")

        # Write summary counts from raw interactions
        for output_key, raw_data in organized["types_raw"].items():
            count = len(raw_data["interactions"])
            f.write(f"  {raw_data['label']}: {count}\n")

        if organized["ligands"]:
            total_ligand_interactions = sum(
                len(l["interactions_raw"]) + len(l["water_bridges_raw"])
                for l in organized["ligands"]
            )
            f.write(f"  Ligand interactions: {total_ligand_interactions}\n")
            f.write(f"  Unique ligands: {len(organized['ligands'])}\n")

        f.write(f"  Total interactions: {summary['total_interactions']}\n\n")

        # Write detailed results for each interaction type (using raw objects)
        for output_key, raw_data in organized["types_raw"].items():
            f.write(f"\n{raw_data['label']}:\n")
            f.write("-" * 30 + "\n")
            for interaction in raw_data["interactions"]:
                f.write(f"{interaction}\n")

        # Write ligand interactions organized by ligand
        if organized["ligands"]:
            f.write("\nLigand Interactions:\n")
            f.write("=" * 50 + "\n\n")
            for ligand_data in organized["ligands"]:
                ligand_res = ligand_data["ligand_residue"]
                ligand_info = ligand_data["ligand_info"]
                regular_interactions = ligand_data["interactions_raw"]
                water_bridges = ligand_data["water_bridges_raw"]

                f.write(f"{ligand_res} ({ligand_info['count']} interactions):\n")
                f.write("-" * 40 + "\n")

                # Write regular interactions
                if regular_interactions:
                    f.write("Regular Interactions:\n")
                    for interaction in regular_interactions:
                        f.write(f"  {interaction}\n")

                # Write water bridges
                if water_bridges:
                    if regular_interactions:
                        f.write("\n")
                    f.write("Water Bridges:\n")
                    for bridge in water_bridges:
                        f.write(f"  Start: {bridge.get_donor_residue()}\n")
                        f.write(f"  Water: {' → '.join(bridge.water_residues)}\n")
                        f.write(f"  End: {bridge.get_acceptor_residue()}\n")
                        f.write(f"  Hops: {bridge.bridge_length}\n")
                        f.write(
                            f"  Distance: {bridge.get_donor_acceptor_distance():.2f} Å\n\n"
                        )

                f.write("\n")


def export_to_csv_files(
    analyzer: NPMolecularInteractionAnalyzer, base_filename: str
) -> None:
    """Export all interaction types to separate CSV files.

    Creates one CSV file per interaction type with the naming pattern:
    {base_name}_interaction_type.csv

    Also creates one CSV file per ligand with ligand-specific interactions:
    {base_name}_ligand_{ligand_residue}.csv

    :param analyzer: Analyzer instance with interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param base_filename: Base filename (extension will be removed)
    :type base_filename: str
    :returns: None
    :rtype: None
    """
    base_path = Path(base_filename)
    base_name = base_path.stem
    directory = base_path.parent

    # Mapping of interaction types to filename patterns
    filename_patterns = {
        "hydrogen_bonds": f"{base_name}_h_bonds.csv",
        "halogen_bonds": f"{base_name}_x_bonds.csv",
        "pi_interactions": f"{base_name}_pi_interactions.csv",
        "pi_pi_interactions": f"{base_name}_pi_pi_interactions.csv",
        "carbonyl_interactions": f"{base_name}_carbonyl_interactions.csv",
        "n_pi_interactions": f"{base_name}_n_pi_interactions.csv",
        "cooperativity_chains": f"{base_name}_cooperativity_chains.csv",
        "water_bridges": f"{base_name}_water_bridges.csv",
    }

    # Export each interaction type using generic function
    for interaction_type, filename_pattern in filename_patterns.items():
        config = get_interaction_config(interaction_type)
        if config and hasattr(analyzer, config.analyzer_attr):
            interactions = getattr(analyzer, config.analyzer_attr)
            if interactions:
                output_file = directory / filename_pattern
                write_interaction_to_csv(analyzer, interaction_type, output_file)

    # Export ligand interactions - one CSV per ligand (regular interactions only)
    # and one CSV per ligand for water bridges
    if hasattr(analyzer, "ligand_interactions") and analyzer.ligand_interactions:
        for ligand_residue in analyzer.ligand_interactions.ligand_info.keys():
            # Regular interactions
            ligand_file = (
                directory / f"{base_name}_ligand_{ligand_residue.replace(':', '_')}.csv"
            )
            write_ligand_interactions_csv(analyzer, ligand_file, ligand_residue)

            # Water bridges
            wb_file = (
                directory
                / f"{base_name}_ligand_{ligand_residue.replace(':', '_')}_water_bridges.csv"
            )
            write_ligand_water_bridges_csv(analyzer, wb_file, ligand_residue)


def export_to_json_files(
    analyzer: NPMolecularInteractionAnalyzer,
    base_filename: str,
    input_file: Optional[str] = None,
) -> None:
    """Export all interaction types to separate JSON files.

    Creates one JSON file per interaction type with the naming pattern:
    {base_name}_interaction_type.json

    Also creates one JSON file per ligand with ligand-specific interactions:
    {base_name}_ligand_{ligand_residue}.json

    :param analyzer: Analyzer instance with interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param base_filename: Base filename (extension will be removed)
    :type base_filename: str
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    base_path = Path(base_filename)
    base_name = base_path.stem
    directory = base_path.parent

    # Mapping of interaction types to filename patterns
    filename_patterns = {
        "hydrogen_bonds": f"{base_name}_h_bonds.json",
        "halogen_bonds": f"{base_name}_x_bonds.json",
        "pi_interactions": f"{base_name}_pi_interactions.json",
        "pi_pi_interactions": f"{base_name}_pi_pi_interactions.json",
        "carbonyl_interactions": f"{base_name}_carbonyl_interactions.json",
        "n_pi_interactions": f"{base_name}_n_pi_interactions.json",
        "cooperativity_chains": f"{base_name}_cooperativity_chains.json",
        "water_bridges": f"{base_name}_water_bridges.json",
    }

    # Export each interaction type using generic function
    for interaction_type, filename_pattern in filename_patterns.items():
        config = get_interaction_config(interaction_type)
        if config and hasattr(analyzer, config.analyzer_attr):
            interactions = getattr(analyzer, config.analyzer_attr)
            if interactions:
                output_file = directory / filename_pattern
                write_interaction_to_json(analyzer, interaction_type, output_file, input_file)

    # Export ligand interactions - one JSON per ligand (regular interactions only)
    # and one JSON per ligand for water bridges
    if hasattr(analyzer, "ligand_interactions") and analyzer.ligand_interactions:
        for ligand_residue in analyzer.ligand_interactions.ligand_info.keys():
            # Regular interactions
            ligand_file = (
                directory
                / f"{base_name}_ligand_{ligand_residue.replace(':', '_')}.json"
            )
            write_ligand_interactions_json(
                analyzer, ligand_file, ligand_residue, input_file
            )

            # Water bridges
            wb_file = (
                directory
                / f"{base_name}_ligand_{ligand_residue.replace(':', '_')}_water_bridges.json"
            )
            write_ligand_water_bridges_json(
                analyzer, wb_file, ligand_residue, input_file
            )


def _organize_all_interactions_by_type(
    analyzer: NPMolecularInteractionAnalyzer,
) -> dict:
    """Organize all interactions by type and ligand.

    Helper shared by export_to_json_single_file() and export_to_txt_single_file()
    to avoid duplicating the organization logic.

    Returns both formatted (for JSON) and raw (for text) interaction lists.

    :param analyzer: Analyzer instance with interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :returns: Dict with organized interactions: {"types": {...}, "ligands": [...]}
    :rtype: dict
    """
    # Mapping of interaction types to output keys (handling key name differences)
    interaction_key_map = {
        "hydrogen_bonds": "hydrogen_bonds",
        "halogen_bonds": "halogen_bonds",
        "pi_interactions": "pi_interactions",
        "pi_pi_interactions": "pi_pi_stacking",
        "carbonyl_interactions": "carbonyl_interactions",
        "n_pi_interactions": "n_pi_interactions",
        "cooperativity_chains": "cooperativity_chains",
        "water_bridges": "water_bridges",
    }

    # Organize standard interaction types (both formatted and raw)
    organized = {
        "types": {},  # Formatted data for JSON
        "types_raw": {},  # Raw objects for text
    }

    for interaction_type, output_key in interaction_key_map.items():
        config = get_interaction_config(interaction_type)
        if config and hasattr(analyzer, config.analyzer_attr):
            interactions = getattr(analyzer, config.analyzer_attr, [])
            if interactions:
                # Store both formatted (for JSON) and raw (for text)
                organized["types"][output_key] = _extract_interaction_json_data(
                    analyzer, interaction_type
                )
                organized["types_raw"][output_key] = {
                    "label": config.label,
                    "interactions": interactions,
                }

    # Organize ligand interactions
    organized["ligands"] = []
    if hasattr(analyzer, "ligand_interactions") and analyzer.ligand_interactions:
        for ligand_res in sorted(analyzer.ligand_interactions.ligand_info.keys()):
            ligand_info = analyzer.ligand_interactions.ligand_info.get(ligand_res, {})

            ligand_data = {
                "ligand_residue": ligand_res,
                "ligand_info": ligand_info,
                # Formatted data (for JSON)
                "interactions_formatted": _build_ligand_interactions_data(
                    analyzer, ligand_res
                ),
                "water_bridges_formatted": _build_ligand_water_bridges_data(
                    analyzer, ligand_res
                ),
                # Raw objects (for text)
                "interactions_raw": _get_ligand_interactions_list(analyzer, ligand_res),
                "water_bridges_raw": _get_ligand_water_bridges_list(analyzer, ligand_res),
            }

            organized["ligands"].append(ligand_data)

    return organized


def export_to_json_single_file(
    analyzer: NPMolecularInteractionAnalyzer,
    output_file: str,
    input_file: Optional[str] = None,
) -> None:
    """Export all interaction types to a single comprehensive JSON file.

    Combines all per-type JSON exports (using same logic) into a single file with summary.

    :param analyzer: Analyzer instance with interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param output_file: Output JSON file path
    :type output_file: str
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    import time

    # Get organized interactions using shared helper
    organized = _organize_all_interactions_by_type(analyzer)

    # Build JSON structure using formatted data
    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "hbat_version": __version__,
        },
        "summary": analyzer.get_summary(),
    }

    # Add standard interaction types (formatted)
    for output_key, formatted_data in organized["types"].items():
        data[output_key] = formatted_data

    # Add ligand interactions (formatted)
    data["ligand_interactions"] = []
    for ligand_data in organized["ligands"]:
        data["ligand_interactions"].append(
            {
                "ligand_residue": ligand_data["ligand_residue"],
                "ligand_info": ligand_data["ligand_info"],
                "interactions": ligand_data["interactions_formatted"],
                "water_bridges": ligand_data["water_bridges_formatted"],
            }
        )

    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


# ============================================================================
# GENERIC WRITE FUNCTIONS (using centralized config)
# ============================================================================

def write_interaction_to_csv(
    analyzer: NPMolecularInteractionAnalyzer,
    interaction_type: str,
    filename: Path,
) -> None:
    """Generic CSV writer for any interaction type using centralized config.

    Dynamically generates CSV columns and data extraction from InteractionConfig.

    :param analyzer: Analyzer with interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param interaction_type: Interaction type ID (e.g., "hydrogen_bonds")
    :type interaction_type: str
    :param filename: Output CSV file path
    :type filename: Path
    :returns: None
    :rtype: None
    """
    config = get_interaction_config(interaction_type)
    if not config:
        return

    interactions = getattr(analyzer, config.analyzer_attr, [])
    if not interactions:
        return

    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)

        # Write headers from config
        headers = [col.csv_header for col in config.columns]
        writer.writerow(headers)

        # Write data rows using extract_interaction_data
        for interaction in interactions:
            data = extract_interaction_data(interaction, config)
            row = [data.get(col.name) for col in config.columns]
            writer.writerow(row)


def _extract_interaction_json_data(
    analyzer: NPMolecularInteractionAnalyzer,
    interaction_type: str,
) -> list:
    """Extract formatted JSON data for an interaction type without metadata.

    Helper function used by both per-file and single-file JSON exports.

    :param analyzer: Analyzer with interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param interaction_type: Interaction type ID (e.g., "hydrogen_bonds")
    :type interaction_type: str
    :returns: List of formatted interaction dictionaries
    :rtype: list
    """
    config = get_interaction_config(interaction_type)
    if not config:
        return []

    interactions = getattr(analyzer, config.analyzer_attr, [])
    if not interactions:
        return []

    result = []
    for interaction in interactions:
        row_data = extract_interaction_data(interaction, config)
        # Convert keys to use json_key names
        interaction_dict = {
            col.json_key: row_data.get(col.name)
            for col in config.columns
        }
        result.append(interaction_dict)

    return result


def write_interaction_to_json(
    analyzer: NPMolecularInteractionAnalyzer,
    interaction_type: str,
    filename: Path,
    input_file: Optional[str] = None,
) -> None:
    """Generic JSON writer for any interaction type using centralized config.

    Dynamically generates JSON structure from InteractionConfig.

    :param analyzer: Analyzer with interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param interaction_type: Interaction type ID (e.g., "hydrogen_bonds")
    :type interaction_type: str
    :param filename: Output JSON file path
    :type filename: Path
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    config = get_interaction_config(interaction_type)
    if not config:
        return

    interactions_data = _extract_interaction_json_data(analyzer, interaction_type)
    if not interactions_data:
        return

    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "interaction_type": config.label,
        },
        "interactions": interactions_data,
    }

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


# ============================================================================
# LIGAND-SPECIFIC WRITE FUNCTIONS
# ============================================================================

# Individual JSON write functions for ligand interactions


def write_ligand_interactions_json(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Path,
    ligand_residue: str,
    input_file: Optional[str] = None,
) -> None:
    """Write ligand interactions to JSON file for a specific ligand.

    :param analyzer: Analyzer with ligand interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output JSON file path
    :type filename: Path
    :param ligand_residue: Ligand residue identifier (e.g., "A:GTP:301")
    :type ligand_residue: str
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    ligand_info = analyzer.ligand_interactions.ligand_info.get(ligand_residue, {})

    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "ligand_residue": ligand_residue,
            "ligand_info": ligand_info,
        },
        "interactions": _build_ligand_interactions_data(analyzer, ligand_residue),
    }

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def _is_water_bridge(interaction) -> bool:
    """Check if an interaction is a water bridge.

    Water bridges have the water_residues attribute from the WaterBridge class.

    :param interaction: MolecularInteraction object
    :returns: True if interaction is a water bridge
    :rtype: bool
    """
    return hasattr(interaction, "water_residues")


def _format_water_bridge_row(water_bridge) -> dict:
    """Format a single water bridge for CSV/JSON output.

    Helper function to avoid duplication between export formats.

    :param water_bridge: WaterBridge object
    :returns: Dictionary with formatted water bridge data
    :rtype: dict
    """
    return {
        "donor_res": water_bridge.get_donor_residue(),
        "acceptor_res": water_bridge.get_acceptor_residue(),
        "bridge_length": water_bridge.bridge_length,
        "water_residues": "; ".join(water_bridge.water_residues),
        "distance": f"{water_bridge.get_donor_acceptor_distance():.2f}",
    }


def _format_ligand_interaction_row(interaction) -> dict:
    """Format a single ligand interaction for CSV/JSON output.

    Helper function to avoid duplication between export formats.

    :param interaction: MolecularInteraction object
    :returns: Dictionary with formatted interaction data
    :rtype: dict
    """
    donor_res = interaction.get_donor_residue()
    acceptor_res = interaction.get_acceptor_residue()
    donor_atom = interaction.get_donor()
    acceptor_atom = interaction.get_acceptor()

    # Determine interaction type label
    int_type = interaction.get_interaction_type()
    if "hydrogen" in int_type.lower() or "h-bond" in int_type.lower():
        type_label = "H-Bond"
    elif "halogen" in int_type.lower() or "x-bond" in int_type.lower():
        type_label = "Halogen Bond"
    elif "pi-pi" in int_type.lower() or "stacking" in int_type.lower():
        type_label = "π-π Stacking"
    elif "pi" in int_type.lower():
        type_label = "π-Interaction"
    elif "carbonyl" in int_type.lower():
        type_label = "Carbonyl"
    elif "n-pi" in int_type.lower():
        type_label = "n-π*"
    elif "water" in int_type.lower():
        type_label = "Water Bridge"
    else:
        type_label = int_type

    # Get distance/angle metric
    distance_str = "N/A"
    if hasattr(interaction, "distance"):
        distance_str = f"{interaction.distance:.2f}"
    elif hasattr(interaction, "_distance"):
        distance_str = f"{interaction._distance:.2f}"

    # Get angle metric if available
    angle_str = "N/A"
    if hasattr(interaction, "angle"):
        angle_str = f"{math.degrees(interaction.angle):.1f}"
    elif hasattr(interaction, "plane_angle"):
        angle_str = f"{interaction.plane_angle:.1f}"
    elif hasattr(interaction, "burgi_dunitz_angle"):
        angle_str = f"{interaction.burgi_dunitz_angle:.1f}"
    elif hasattr(interaction, "angle_to_plane"):
        angle_str = f"{interaction.angle_to_plane:.1f}"

    # Get properties
    properties = ""
    if hasattr(interaction, "donor_acceptor_properties"):
        properties = interaction.donor_acceptor_properties

    # Get atom names
    donor_atom_name = donor_atom.name if hasattr(donor_atom, "name") else "N/A"
    acceptor_atom_name = acceptor_atom.name if hasattr(acceptor_atom, "name") else "N/A"

    return {
        "type_label": type_label,
        "donor_res": donor_res,
        "donor_atom": donor_atom_name,
        "acceptor_res": acceptor_res,
        "acceptor_atom": acceptor_atom_name,
        "distance": distance_str,
        "angle": angle_str,
        "properties": properties,
    }


def _get_ligand_interactions_list(
    analyzer: NPMolecularInteractionAnalyzer,
    ligand_residue: str,
) -> list:
    """Get unformatted regular ligand interactions for a specific ligand.

    Helper used by all ligand export functions to avoid duplicating filter logic.

    :param analyzer: Analyzer with ligand interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param ligand_residue: Ligand residue identifier (e.g., "A:GTP:301")
    :type ligand_residue: str
    :returns: List of MolecularInteraction objects (regular interactions only)
    :rtype: list
    """
    all_interactions = analyzer.ligand_interactions.get_interactions_for_ligand(
        ligand_residue
    )
    return [i for i in all_interactions if not _is_water_bridge(i)]


def _get_ligand_water_bridges_list(
    analyzer: NPMolecularInteractionAnalyzer,
    ligand_residue: str,
) -> list:
    """Get unformatted water bridges for a specific ligand.

    Helper used by all ligand export functions to avoid duplicating filter logic.

    :param analyzer: Analyzer with ligand interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param ligand_residue: Ligand residue identifier (e.g., "A:GTP:301")
    :type ligand_residue: str
    :returns: List of WaterBridge objects
    :rtype: list
    """
    all_interactions = analyzer.ligand_interactions.get_interactions_for_ligand(
        ligand_residue
    )
    return [i for i in all_interactions if _is_water_bridge(i)]


def _build_ligand_interactions_data(
    analyzer: NPMolecularInteractionAnalyzer,
    ligand_residue: str,
) -> list:
    """Extract and format regular ligand interactions for a specific ligand.

    Helper shared by write_ligand_interactions_json() and export_to_json_single_file().

    :param analyzer: Analyzer with ligand interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param ligand_residue: Ligand residue identifier (e.g., "A:GTP:301")
    :type ligand_residue: str
    :returns: List of formatted interaction dictionaries
    :rtype: list
    """
    regular_interactions = _get_ligand_interactions_list(analyzer, ligand_residue)

    interactions_list = []
    for interaction in regular_interactions:
        try:
            row = _format_ligand_interaction_row(interaction)
            interaction_data = {
                "interaction_type": row["type_label"],
                "donor_residue": row["donor_res"],
                "donor_atom": row["donor_atom"],
                "acceptor_residue": row["acceptor_res"],
                "acceptor_atom": row["acceptor_atom"],
            }

            if row["distance"] != "N/A":
                interaction_data["distance_angstrom"] = float(row["distance"])
            if row["angle"] != "N/A":
                interaction_data["angle_or_metric_degrees"] = float(row["angle"])
            if row["properties"]:
                interaction_data["properties"] = row["properties"]

            interactions_list.append(interaction_data)
        except Exception:
            continue

    return interactions_list


def _build_ligand_water_bridges_data(
    analyzer: NPMolecularInteractionAnalyzer,
    ligand_residue: str,
) -> list:
    """Extract and format water bridges for a specific ligand.

    Helper shared by write_ligand_water_bridges_json() and export_to_json_single_file().

    :param analyzer: Analyzer with ligand interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param ligand_residue: Ligand residue identifier (e.g., "A:GTP:301")
    :type ligand_residue: str
    :returns: List of formatted water bridge dictionaries
    :rtype: list
    """
    water_bridges = _get_ligand_water_bridges_list(analyzer, ligand_residue)

    bridges_list = []
    for wb in water_bridges:
        try:
            water_bridge_data = {
                "start_residue": wb.get_donor_residue(),
                "water_molecules": wb.water_residues,
                "end_residue": wb.get_acceptor_residue(),
                "bridge_length": wb.bridge_length,
                "distance_angstrom": float(
                    f"{wb.get_donor_acceptor_distance():.2f}"
                ),
            }

            if hasattr(wb, "donor_acceptor_properties"):
                water_bridge_data["properties"] = wb.donor_acceptor_properties

            bridges_list.append(water_bridge_data)
        except Exception:
            continue

    return bridges_list


def write_ligand_interactions_csv(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Optional[Path] = None,
    ligand_residue: Optional[str] = None,
) -> Optional[str]:
    """Write ligand interactions to CSV file or return as string.

    If filename is provided, writes to file and returns None.
    If filename is None, returns CSV content as string.
    If ligand_residue is specified, filters to only that ligand.
    If ligand_residue is None, includes all ligands.

    :param analyzer: Analyzer with ligand interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output CSV file path, or None to return as string
    :type filename: Optional[Path]
    :param ligand_residue: Optional ligand residue identifier (e.g., "A:GTP:301") to filter by
    :type ligand_residue: Optional[str]
    :returns: CSV content as string (if filename is None), otherwise None
    :rtype: Optional[str]
    """
    # Use StringIO if returning as string, otherwise use file
    if filename is None:
        output = io.StringIO()
        use_string_io = True
    else:
        output = open(filename, "w", newline="", encoding="utf-8")
        use_string_io = False

    try:
        writer = csv.writer(output)

        # Write header
        writer.writerow(
            [
                "Interaction_Type",
                "Donor_Residue",
                "Donor_Atom",
                "Acceptor_Residue",
                "Acceptor_Atom",
                "Distance_Angstrom",
                "Angle_Or_Metric_Degrees",
                "Properties",
            ]
        )

        # Get interactions from analyzer's ligand_interactions container
        if ligand_residue:
            all_interactions = analyzer.ligand_interactions.get_interactions_for_ligand(
                ligand_residue
            )
        else:
            all_interactions = analyzer.ligand_interactions.interactions

        # Filter out water bridges (they go in separate files)
        interactions = [i for i in all_interactions if not _is_water_bridge(i)]

        # Write interactions
        for interaction in interactions:
            try:
                row = _format_ligand_interaction_row(interaction)
                writer.writerow(
                    [
                        row["type_label"],
                        row["donor_res"],
                        row["donor_atom"],
                        row["acceptor_res"],
                        row["acceptor_atom"],
                        row["distance"],
                        row["angle"],
                        row["properties"],
                    ]
                )
            except Exception:
                # Skip interactions that can't be processed
                continue

        if use_string_io:
            return output.getvalue()
        else:
            return None
    finally:
        if use_string_io:
            output.close()
        else:
            output.close()


def write_ligand_water_bridges_csv(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Optional[Path] = None,
    ligand_residue: Optional[str] = None,
) -> Optional[str]:
    """Write ligand water bridges to CSV file or return as string.

    If filename is provided, writes to file and returns None.
    If filename is None, returns CSV content as string.
    If ligand_residue is specified, filters to only that ligand.
    If ligand_residue is None, includes all ligands.

    :param analyzer: Analyzer with ligand interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output CSV file path, or None to return as string
    :type filename: Optional[Path]
    :param ligand_residue: Optional ligand residue identifier (e.g., "A:GTP:301") to filter by
    :type ligand_residue: Optional[str]
    :returns: CSV content as string (if filename is None), otherwise None
    :rtype: Optional[str]
    """
    # Use StringIO if returning as string, otherwise use file
    if filename is None:
        output = io.StringIO()
        use_string_io = True
    else:
        output = open(filename, "w", newline="", encoding="utf-8")
        use_string_io = False

    try:
        writer = csv.writer(output)

        # Write header for water bridges (consistent with write_water_bridges_csv)
        writer.writerow(
            [
                "Donor_Residue",
                "Water_Molecules",
                "Acceptor_Residue",
                "Hops",
                "Distance_Angstrom",
            ]
        )

        # Get interactions from analyzer's ligand_interactions container
        if ligand_residue:
            all_interactions = analyzer.ligand_interactions.get_interactions_for_ligand(
                ligand_residue
            )
        else:
            all_interactions = analyzer.ligand_interactions.interactions

        # Filter to only water bridges
        water_bridges = [i for i in all_interactions if _is_water_bridge(i)]

        # Write water bridges using shared formatting function
        for wb in water_bridges:
            try:
                row = _format_water_bridge_row(wb)
                writer.writerow(
                    [
                        row["donor_res"],
                        row["water_residues"],
                        row["acceptor_res"],
                        row["bridge_length"],
                        row["distance"],
                    ]
                )
            except Exception:
                # Skip water bridges that can't be processed
                continue

        if use_string_io:
            return output.getvalue()
        else:
            return None
    finally:
        if use_string_io:
            output.close()
        else:
            output.close()


def write_ligand_water_bridges_json(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Path,
    ligand_residue: str,
    input_file: Optional[str] = None,
) -> None:
    """Write ligand water bridges to JSON file for a specific ligand.

    :param analyzer: Analyzer with ligand interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output JSON file path
    :type filename: Path
    :param ligand_residue: Ligand residue identifier (e.g., "A:GTP:301")
    :type ligand_residue: str
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    ligand_info = analyzer.ligand_interactions.ligand_info.get(ligand_residue, {})

    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "ligand_residue": ligand_residue,
            "ligand_info": ligand_info,
        },
        "water_bridges": _build_ligand_water_bridges_data(analyzer, ligand_residue),
    }

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
