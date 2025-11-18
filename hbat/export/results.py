"""
Centralized export functions for HBAT analysis results.

This module provides functions to export molecular interaction analysis
results to CSV and JSON formats. These functions are used by both the
CLI and GUI interfaces.
"""

import csv
import json
import math
from pathlib import Path
from typing import Optional

from hbat import __version__
from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer


def export_to_txt_single_file(
    analyzer: NPMolecularInteractionAnalyzer, output_file: str
) -> None:
    """Export all interactions to a single text file with human-readable format.

    Creates a text file containing a summary and detailed listing of all
    interaction types found in the analysis.

    :param analyzer: Analyzer instance with interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param output_file: Path to the output text file
    :type output_file: str
    :returns: None
    :rtype: None
    """
    with open(output_file, "w", encoding="utf-8") as f:
        # Write summary
        summary = analyzer.get_summary()
        f.write("Summary:\n")
        f.write(f"  Hydrogen bonds: {summary['hydrogen_bonds']['count']}\n")
        f.write(f"  Halogen bonds: {summary['halogen_bonds']['count']}\n")
        f.write(f"  π interactions: {summary['pi_interactions']['count']}\n")
        f.write(
            f"  π-π stacking: {summary.get('pi_pi_stacking', {}).get('count', 0)}\n"
        )
        f.write(
            f"  Carbonyl interactions: {summary.get('carbonyl_interactions', {}).get('count', 0)}\n"
        )
        f.write(
            f"  n-π interactions: {summary.get('n_pi_interactions', {}).get('count', 0)}\n"
        )
        f.write(f"  Total interactions: {summary['total_interactions']}\n\n")

        # Write detailed results
        f.write("Hydrogen Bonds:\n")
        f.write("-" * 30 + "\n")
        for hb in analyzer.hydrogen_bonds:
            f.write(f"{hb}\n")

        f.write("\nHalogen Bonds:\n")
        f.write("-" * 30 + "\n")
        for xb in analyzer.halogen_bonds:
            f.write(f"{xb}\n")

        f.write("\nπ Interactions:\n")
        f.write("-" * 30 + "\n")
        for pi in analyzer.pi_interactions:
            f.write(f"{pi}\n")

        # Write π-π stacking interactions if available
        if hasattr(analyzer, "pi_pi_interactions") and analyzer.pi_pi_interactions:
            f.write("\nπ-π Stacking Interactions:\n")
            f.write("-" * 30 + "\n")
            for pi_pi in analyzer.pi_pi_interactions:
                f.write(f"{pi_pi}\n")

        # Write carbonyl interactions if available
        if (
            hasattr(analyzer, "carbonyl_interactions")
            and analyzer.carbonyl_interactions
        ):
            f.write("\nCarbonyl Interactions:\n")
            f.write("-" * 30 + "\n")
            for carbonyl in analyzer.carbonyl_interactions:
                f.write(f"{carbonyl}\n")

        # Write n-π interactions if available
        if hasattr(analyzer, "n_pi_interactions") and analyzer.n_pi_interactions:
            f.write("\nn-π Interactions:\n")
            f.write("-" * 30 + "\n")
            for n_pi in analyzer.n_pi_interactions:
                f.write(f"{n_pi}\n")

        # Write cooperativity chains if available
        if hasattr(analyzer, "cooperativity_chains") and analyzer.cooperativity_chains:
            f.write("\nCooperativity Chains:\n")
            f.write("-" * 30 + "\n")
            for chain in analyzer.cooperativity_chains:
                f.write(f"{chain}\n")


def export_to_csv_files(
    analyzer: NPMolecularInteractionAnalyzer, base_filename: str
) -> None:
    """Export all interaction types to separate CSV files.

    Creates one CSV file per interaction type with the naming pattern:
    {base_name}_interaction_type.csv

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

    # Export each interaction type
    if analyzer.hydrogen_bonds:
        hb_file = directory / f"{base_name}_h_bonds.csv"
        write_hydrogen_bonds_csv(analyzer, hb_file)

    if analyzer.halogen_bonds:
        xb_file = directory / f"{base_name}_x_bonds.csv"
        write_halogen_bonds_csv(analyzer, xb_file)

    if analyzer.pi_interactions:
        pi_file = directory / f"{base_name}_pi_interactions.csv"
        write_pi_interactions_csv(analyzer, pi_file)

    if hasattr(analyzer, "pi_pi_interactions") and analyzer.pi_pi_interactions:
        pi_pi_file = directory / f"{base_name}_pi_pi_interactions.csv"
        write_pi_pi_interactions_csv(analyzer, pi_pi_file)

    if hasattr(analyzer, "carbonyl_interactions") and analyzer.carbonyl_interactions:
        carbonyl_file = directory / f"{base_name}_carbonyl_interactions.csv"
        write_carbonyl_interactions_csv(analyzer, carbonyl_file)

    if hasattr(analyzer, "n_pi_interactions") and analyzer.n_pi_interactions:
        n_pi_file = directory / f"{base_name}_n_pi_interactions.csv"
        write_n_pi_interactions_csv(analyzer, n_pi_file)

    if hasattr(analyzer, "cooperativity_chains") and analyzer.cooperativity_chains:
        chains_file = directory / f"{base_name}_cooperativity_chains.csv"
        write_cooperativity_chains_csv(analyzer, chains_file)


def export_to_json_files(
    analyzer: NPMolecularInteractionAnalyzer,
    base_filename: str,
    input_file: Optional[str] = None,
) -> None:
    """Export all interaction types to separate JSON files.

    Creates one JSON file per interaction type with the naming pattern:
    {base_name}_interaction_type.json

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

    # Export each interaction type
    if analyzer.hydrogen_bonds:
        hb_file = directory / f"{base_name}_h_bonds.json"
        write_hydrogen_bonds_json(analyzer, hb_file, input_file)

    if analyzer.halogen_bonds:
        xb_file = directory / f"{base_name}_x_bonds.json"
        write_halogen_bonds_json(analyzer, xb_file, input_file)

    if analyzer.pi_interactions:
        pi_file = directory / f"{base_name}_pi_interactions.json"
        write_pi_interactions_json(analyzer, pi_file, input_file)

    if hasattr(analyzer, "pi_pi_interactions") and analyzer.pi_pi_interactions:
        pi_pi_file = directory / f"{base_name}_pi_pi_interactions.json"
        write_pi_pi_interactions_json(analyzer, pi_pi_file, input_file)

    if hasattr(analyzer, "carbonyl_interactions") and analyzer.carbonyl_interactions:
        carbonyl_file = directory / f"{base_name}_carbonyl_interactions.json"
        write_carbonyl_interactions_json(analyzer, carbonyl_file, input_file)

    if hasattr(analyzer, "n_pi_interactions") and analyzer.n_pi_interactions:
        n_pi_file = directory / f"{base_name}_n_pi_interactions.json"
        write_n_pi_interactions_json(analyzer, n_pi_file, input_file)

    if hasattr(analyzer, "cooperativity_chains") and analyzer.cooperativity_chains:
        chains_file = directory / f"{base_name}_cooperativity_chains.json"
        write_cooperativity_chains_json(analyzer, chains_file, input_file)


def export_to_json_single_file(
    analyzer: NPMolecularInteractionAnalyzer,
    output_file: str,
    input_file: Optional[str] = None,
) -> None:
    """Export all interaction types to a single JSON file.

    Creates a comprehensive JSON file with all interaction types.

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

    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "hbat_version": __version__,
        },
        "summary": analyzer.get_summary(),
        "hydrogen_bonds": [],
        "halogen_bonds": [],
        "pi_interactions": [],
        "pi_pi_stacking": [],
        "carbonyl_interactions": [],
        "n_pi_interactions": [],
        "cooperativity_chains": [],
    }

    # Hydrogen bonds
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
                "backbone_sidechain": hb.get_backbone_sidechain_interaction(),
                "donor_acceptor_properties": hb.donor_acceptor_properties,
            }
        )

    # Halogen bonds
    for xb in analyzer.halogen_bonds:
        data["halogen_bonds"].append(
            {
                "donor_residue": xb.donor_residue,
                "donor_atom": xb.donor.name,
                "halogen_atom": xb.halogen.name,
                "halogen_coords": xb.halogen.coords.to_list(),
                "acceptor_residue": xb.acceptor_residue,
                "acceptor_atom": xb.acceptor.name,
                "acceptor_coords": xb.acceptor.coords.to_list(),
                "distance": round(xb.distance, 3),
                "angle": round(math.degrees(xb.angle), 1),
                "bond_type": xb.bond_type,
                "backbone_sidechain": xb.get_backbone_sidechain_interaction(),
                "donor_acceptor_properties": xb.donor_acceptor_properties,
            }
        )

    # π interactions
    for pi in analyzer.pi_interactions:
        data["pi_interactions"].append(
            {
                "donor_residue": pi.donor_residue,
                "donor_atom": pi.donor.name,
                "hydrogen_atom": pi.hydrogen.name,
                "pi_residue": pi.pi_residue,
                "distance": round(pi.distance, 3),
                "angle": round(math.degrees(pi.angle), 1),
                "interaction_type": pi.get_interaction_type_display(),
                "backbone_sidechain": pi.get_backbone_sidechain_interaction(),
                "donor_acceptor_properties": pi.donor_acceptor_properties,
            }
        )

    # π-π stacking
    if hasattr(analyzer, "pi_pi_interactions"):
        for pi_pi in analyzer.pi_pi_interactions:
            data["pi_pi_stacking"].append(
                {
                    "ring1_residue": pi_pi.ring1_residue,
                    "ring1_type": pi_pi.ring1_type,
                    "ring2_residue": pi_pi.ring2_residue,
                    "ring2_type": pi_pi.ring2_type,
                    "distance": round(pi_pi.distance, 3),
                    "plane_angle": round(pi_pi.plane_angle, 1),
                    "offset": round(pi_pi.offset, 3),
                    "stacking_type": pi_pi.stacking_type,
                }
            )

    # Carbonyl interactions
    if hasattr(analyzer, "carbonyl_interactions"):
        for carbonyl in analyzer.carbonyl_interactions:
            data["carbonyl_interactions"].append(
                {
                    "donor_residue": carbonyl.donor_residue,
                    "donor_carbon": carbonyl.donor_carbon.name,
                    "donor_oxygen": carbonyl.donor_oxygen.name,
                    "acceptor_residue": carbonyl.acceptor_residue,
                    "acceptor_carbon": carbonyl.acceptor_carbon.name,
                    "acceptor_oxygen": carbonyl.acceptor_oxygen.name,
                    "distance": round(carbonyl.distance, 3),
                    "burgi_dunitz_angle": round(carbonyl.burgi_dunitz_angle, 1),
                    "interaction_type": carbonyl.interaction_classification,
                    "is_backbone": carbonyl.is_backbone,
                }
            )

    # n-π interactions
    if hasattr(analyzer, "n_pi_interactions"):
        for n_pi in analyzer.n_pi_interactions:
            pi_atom_name = n_pi.pi_atoms[0].name if n_pi.pi_atoms else "?"
            data["n_pi_interactions"].append(
                {
                    "donor_residue": n_pi.donor_residue,
                    "lone_pair_atom": n_pi.lone_pair_atom.name,
                    "lone_pair_element": n_pi.lone_pair_atom.element,
                    "acceptor_residue": n_pi.acceptor_residue,
                    "pi_atom": pi_atom_name,
                    "distance": round(n_pi.distance, 3),
                    "angle_to_plane": round(n_pi.angle_to_plane, 1),
                    "subtype": n_pi.subtype,
                }
            )

    # Cooperativity chains
    if hasattr(analyzer, "cooperativity_chains"):
        for i, chain in enumerate(analyzer.cooperativity_chains):
            chain_data = {
                "chain_id": i + 1,
                "chain_length": chain.chain_length,
                "chain_type": chain.chain_type,
                "interactions": [],
            }

            for interaction in chain.interactions:
                interaction_data = {
                    "donor_residue": interaction.get_donor_residue(),
                    "acceptor_residue": interaction.get_acceptor_residue(),
                    "interaction_type": interaction.get_interaction_type(),
                }

                donor_atom = interaction.get_donor_atom()
                if donor_atom:
                    interaction_data["donor_atom"] = donor_atom.name

                acceptor_atom = interaction.get_acceptor_atom()
                if acceptor_atom:
                    interaction_data["acceptor_atom"] = acceptor_atom.name

                chain_data["interactions"].append(interaction_data)

            data["cooperativity_chains"].append(chain_data)

    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


# Individual CSV write functions
def write_hydrogen_bonds_csv(
    analyzer: NPMolecularInteractionAnalyzer, filename: Path
) -> None:
    """Write hydrogen bonds to CSV file.

    :param analyzer: Analyzer with hydrogen bond results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output CSV file path
    :type filename: Path
    :returns: None
    :rtype: None
    """
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "Donor_Residue",
                "Donor_Atom",
                "Hydrogen_Atom",
                "Acceptor_Residue",
                "Acceptor_Atom",
                "Distance_Angstrom",
                "Angle_Degrees",
                "Donor_Acceptor_Distance_Angstrom",
                "Bond_Type",
                "B/S_Interaction",
                "D-A_Properties",
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
                    hb.get_backbone_sidechain_interaction(),
                    hb.donor_acceptor_properties,
                ]
            )


def write_halogen_bonds_csv(
    analyzer: NPMolecularInteractionAnalyzer, filename: Path
) -> None:
    """Write halogen bonds to CSV file.

    :param analyzer: Analyzer with halogen bond results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output CSV file path
    :type filename: Path
    :returns: None
    :rtype: None
    """
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "Halogen_Residue",
                "Halogen_Atom",
                "Acceptor_Residue",
                "Acceptor_Atom",
                "Distance_Angstrom",
                "Angle_Degrees",
                "Bond_Type",
                "B/S_Interaction",
                "D-A_Properties",
            ]
        )
        for xb in analyzer.halogen_bonds:
            writer.writerow(
                [
                    xb.donor_residue,
                    xb.halogen.name,
                    xb.acceptor_residue,
                    xb.acceptor.name,
                    f"{xb.distance:.3f}",
                    f"{math.degrees(xb.angle):.1f}",
                    xb.bond_type,
                    xb.get_backbone_sidechain_interaction(),
                    xb.donor_acceptor_properties,
                ]
            )


def write_pi_interactions_csv(
    analyzer: NPMolecularInteractionAnalyzer, filename: Path
) -> None:
    """Write π interactions to CSV file.

    :param analyzer: Analyzer with π interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output CSV file path
    :type filename: Path
    :returns: None
    :rtype: None
    """
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "Donor_Residue",
                "Donor_Atom",
                "Hydrogen_Atom",
                "Pi_Residue",
                "Distance_Angstrom",
                "Angle_Degrees",
                "Interaction_Type",
                "B/S_Interaction",
                "D-A_Properties",
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
                    pi.get_interaction_type_display(),
                    pi.get_backbone_sidechain_interaction(),
                    pi.donor_acceptor_properties,
                ]
            )


def write_pi_pi_interactions_csv(
    analyzer: NPMolecularInteractionAnalyzer, filename: Path
) -> None:
    """Write π-π stacking interactions to CSV file.

    :param analyzer: Analyzer with π-π stacking results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output CSV file path
    :type filename: Path
    :returns: None
    :rtype: None
    """
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "Ring1_Residue",
                "Ring1_Type",
                "Ring2_Residue",
                "Ring2_Type",
                "Distance_Angstrom",
                "Plane_Angle_Degrees",
                "Offset_Angstrom",
                "Stacking_Type",
            ]
        )
        for pi_pi in analyzer.pi_pi_interactions:
            writer.writerow(
                [
                    pi_pi.ring1_residue,
                    pi_pi.ring1_type,
                    pi_pi.ring2_residue,
                    pi_pi.ring2_type,
                    f"{pi_pi.distance:.3f}",
                    f"{pi_pi.plane_angle:.1f}",
                    f"{pi_pi.offset:.3f}",
                    pi_pi.stacking_type,
                ]
            )


def write_carbonyl_interactions_csv(
    analyzer: NPMolecularInteractionAnalyzer, filename: Path
) -> None:
    """Write carbonyl interactions to CSV file.

    :param analyzer: Analyzer with carbonyl interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output CSV file path
    :type filename: Path
    :returns: None
    :rtype: None
    """
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "Donor_Residue",
                "Donor_Carbonyl",
                "Acceptor_Residue",
                "Acceptor_Carbonyl",
                "Distance_Angstrom",
                "Burgi_Dunitz_Angle_Degrees",
                "Interaction_Type",
            ]
        )
        for carbonyl in analyzer.carbonyl_interactions:
            writer.writerow(
                [
                    carbonyl.donor_residue,
                    f"{carbonyl.donor_carbon.name}={carbonyl.donor_oxygen.name}",
                    carbonyl.acceptor_residue,
                    f"{carbonyl.acceptor_carbon.name}={carbonyl.acceptor_oxygen.name}",
                    f"{carbonyl.distance:.3f}",
                    f"{carbonyl.burgi_dunitz_angle:.1f}",
                    carbonyl.interaction_classification,
                ]
            )


def write_n_pi_interactions_csv(
    analyzer: NPMolecularInteractionAnalyzer, filename: Path
) -> None:
    """Write n-π interactions to CSV file.

    :param analyzer: Analyzer with n-π interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output CSV file path
    :type filename: Path
    :returns: None
    :rtype: None
    """
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "Donor_Residue",
                "Lone_Pair_Atom",
                "Acceptor_Residue",
                "Pi_System",
                "Distance_Angstrom",
                "Angle_To_Plane_Degrees",
                "Subtype",
            ]
        )
        for n_pi in analyzer.n_pi_interactions:
            pi_atom_name = n_pi.pi_atoms[0].name if n_pi.pi_atoms else "?"
            writer.writerow(
                [
                    n_pi.donor_residue,
                    n_pi.lone_pair_atom.name,
                    n_pi.acceptor_residue,
                    pi_atom_name,
                    f"{n_pi.distance:.3f}",
                    f"{n_pi.angle_to_plane:.1f}",
                    n_pi.subtype,
                ]
            )


def write_cooperativity_chains_csv(
    analyzer: NPMolecularInteractionAnalyzer, filename: Path
) -> None:
    """Write cooperativity chains to CSV file.

    :param analyzer: Analyzer with cooperativity chain results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output CSV file path
    :type filename: Path
    :returns: None
    :rtype: None
    """
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Chain_ID", "Chain_Length", "Chain_Type", "Interactions"])
        for i, chain in enumerate(analyzer.cooperativity_chains):
            interactions_str = " -> ".join(
                [
                    f"{interaction.get_donor_residue()}({interaction.get_donor_atom().name if interaction.get_donor_atom() else '?'})"
                    for interaction in chain.interactions
                ]
            )
            writer.writerow(
                [i + 1, chain.chain_length, chain.chain_type, interactions_str]
            )


# Individual JSON write functions
def write_hydrogen_bonds_json(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Path,
    input_file: Optional[str] = None,
) -> None:
    """Write hydrogen bonds to JSON file.

    :param analyzer: Analyzer with hydrogen bond results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output JSON file path
    :type filename: Path
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "interaction_type": "Hydrogen Bonds",
        },
        "interactions": [],
    }

    for hb in analyzer.hydrogen_bonds:
        data["interactions"].append(
            {
                "donor_residue": hb.donor_residue,
                "donor_atom": hb.donor.name,
                "hydrogen_atom": hb.hydrogen.name,
                "acceptor_residue": hb.acceptor_residue,
                "acceptor_atom": hb.acceptor.name,
                "distance_angstrom": round(hb.distance, 3),
                "angle_degrees": round(math.degrees(hb.angle), 1),
                "donor_acceptor_distance_angstrom": round(
                    hb.donor_acceptor_distance, 3
                ),
                "bond_type": hb.bond_type,
                "backbone_sidechain_interaction": hb.get_backbone_sidechain_interaction(),
                "donor_acceptor_properties": hb.donor_acceptor_properties,
            }
        )

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def write_halogen_bonds_json(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Path,
    input_file: Optional[str] = None,
) -> None:
    """Write halogen bonds to JSON file.

    :param analyzer: Analyzer with halogen bond results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output JSON file path
    :type filename: Path
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "interaction_type": "Halogen Bonds",
        },
        "interactions": [],
    }

    for xb in analyzer.halogen_bonds:
        data["interactions"].append(
            {
                "halogen_residue": xb.donor_residue,
                "halogen_atom": xb.halogen.name,
                "acceptor_residue": xb.acceptor_residue,
                "acceptor_atom": xb.acceptor.name,
                "distance_angstrom": round(xb.distance, 3),
                "angle_degrees": round(math.degrees(xb.angle), 1),
                "bond_type": xb.bond_type,
                "backbone_sidechain_interaction": xb.get_backbone_sidechain_interaction(),
                "donor_acceptor_properties": xb.donor_acceptor_properties,
            }
        )

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def write_pi_interactions_json(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Path,
    input_file: Optional[str] = None,
) -> None:
    """Write π interactions to JSON file.

    :param analyzer: Analyzer with π interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output JSON file path
    :type filename: Path
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "interaction_type": "π Interactions",
        },
        "interactions": [],
    }

    for pi in analyzer.pi_interactions:
        data["interactions"].append(
            {
                "donor_residue": pi.donor_residue,
                "donor_atom": pi.donor.name,
                "hydrogen_atom": pi.hydrogen.name,
                "pi_residue": pi.pi_residue,
                "distance_angstrom": round(pi.distance, 3),
                "angle_degrees": round(math.degrees(pi.angle), 1),
                "interaction_type": pi.get_interaction_type_display(),
                "backbone_sidechain_interaction": pi.get_backbone_sidechain_interaction(),
                "donor_acceptor_properties": pi.donor_acceptor_properties,
            }
        )

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def write_pi_pi_interactions_json(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Path,
    input_file: Optional[str] = None,
) -> None:
    """Write π-π stacking interactions to JSON file.

    :param analyzer: Analyzer with π-π stacking results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output JSON file path
    :type filename: Path
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "interaction_type": "π-π Stacking Interactions",
        },
        "interactions": [],
    }

    for pi_pi in analyzer.pi_pi_interactions:
        data["interactions"].append(
            {
                "ring1_residue": pi_pi.ring1_residue,
                "ring1_type": pi_pi.ring1_type,
                "ring2_residue": pi_pi.ring2_residue,
                "ring2_type": pi_pi.ring2_type,
                "distance_angstrom": round(pi_pi.distance, 3),
                "plane_angle_degrees": round(pi_pi.plane_angle, 1),
                "offset_angstrom": round(pi_pi.offset, 3),
                "stacking_type": pi_pi.stacking_type,
            }
        )

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def write_carbonyl_interactions_json(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Path,
    input_file: Optional[str] = None,
) -> None:
    """Write carbonyl interactions to JSON file.

    :param analyzer: Analyzer with carbonyl interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output JSON file path
    :type filename: Path
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "interaction_type": "Carbonyl Interactions",
        },
        "interactions": [],
    }

    for carbonyl in analyzer.carbonyl_interactions:
        data["interactions"].append(
            {
                "donor_residue": carbonyl.donor_residue,
                "donor_carbon": carbonyl.donor_carbon.name,
                "donor_oxygen": carbonyl.donor_oxygen.name,
                "acceptor_residue": carbonyl.acceptor_residue,
                "acceptor_carbon": carbonyl.acceptor_carbon.name,
                "acceptor_oxygen": carbonyl.acceptor_oxygen.name,
                "distance_angstrom": round(carbonyl.distance, 3),
                "burgi_dunitz_angle_degrees": round(carbonyl.burgi_dunitz_angle, 1),
                "interaction_type": carbonyl.interaction_classification,
                "is_backbone": carbonyl.is_backbone,
            }
        )

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def write_n_pi_interactions_json(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Path,
    input_file: Optional[str] = None,
) -> None:
    """Write n-π interactions to JSON file.

    :param analyzer: Analyzer with n-π interaction results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output JSON file path
    :type filename: Path
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "interaction_type": "n-π Interactions",
        },
        "interactions": [],
    }

    for n_pi in analyzer.n_pi_interactions:
        pi_atom_name = n_pi.pi_atoms[0].name if n_pi.pi_atoms else "?"

        data["interactions"].append(
            {
                "donor_residue": n_pi.donor_residue,
                "lone_pair_atom": n_pi.lone_pair_atom.name,
                "lone_pair_element": n_pi.lone_pair_atom.element,
                "acceptor_residue": n_pi.acceptor_residue,
                "pi_atom": pi_atom_name,
                "distance_angstrom": round(n_pi.distance, 3),
                "angle_to_plane_degrees": round(n_pi.angle_to_plane, 1),
                "subtype": n_pi.subtype,
            }
        )

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def write_cooperativity_chains_json(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Path,
    input_file: Optional[str] = None,
) -> None:
    """Write cooperativity chains to JSON file.

    :param analyzer: Analyzer with cooperativity chain results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output JSON file path
    :type filename: Path
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "interaction_type": "Cooperativity Chains",
        },
        "chains": [],
    }

    for i, chain in enumerate(analyzer.cooperativity_chains):
        chain_data = {
            "chain_id": i + 1,
            "chain_length": chain.chain_length,
            "chain_type": chain.chain_type,
            "interactions": [],
        }

        for interaction in chain.interactions:
            interaction_data = {
                "donor_residue": interaction.get_donor_residue(),
                "acceptor_residue": interaction.get_acceptor_residue(),
                "interaction_type": interaction.get_interaction_type(),
            }

            donor_atom = interaction.get_donor_atom()
            if donor_atom:
                interaction_data["donor_atom"] = donor_atom.name

            acceptor_atom = interaction.get_acceptor_atom()
            if acceptor_atom:
                interaction_data["acceptor_atom"] = acceptor_atom.name

            chain_data["interactions"].append(interaction_data)

        data["chains"].append(chain_data)

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
