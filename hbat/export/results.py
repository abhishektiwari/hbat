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
from typing import Optional

from hbat import __version__
from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer


def export_to_txt_single_file(
    analyzer: NPMolecularInteractionAnalyzer, output_file: str
) -> None:
    """Export all interactions to a single text file with human-readable format.

    Creates a text file containing a summary and detailed listing of all
    interaction types found in the analysis, including ligand interactions
    organized by ligand.

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
        f.write(
            f"  Water bridges: {summary.get('water_bridges', {}).get('count', 0)}\n"
        )
        if hasattr(analyzer, "ligand_interactions") and analyzer.ligand_interactions:
            f.write(f"  Ligand interactions: {len(analyzer.ligand_interactions.interactions)}\n")
            f.write(f"  Unique ligands: {len(analyzer.ligand_interactions.ligand_info)}\n")
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

        # Write water bridges if available
        if hasattr(analyzer, "water_bridges") and analyzer.water_bridges:
            f.write("\nWater Bridges:\n")
            f.write("-" * 30 + "\n")
            for bridge in analyzer.water_bridges:
                f.write(f"{bridge}\n")

        # Write ligand interactions organized by ligand if available
        if hasattr(analyzer, "ligand_interactions") and analyzer.ligand_interactions:
            f.write("\nLigand Interactions:\n")
            f.write("=" * 50 + "\n\n")
            for ligand_res in sorted(analyzer.ligand_interactions.ligand_info.keys()):
                ligand_info = analyzer.ligand_interactions.ligand_info[ligand_res]
                all_interactions = analyzer.ligand_interactions.get_interactions_for_ligand(ligand_res)

                # Separate regular interactions and water bridges
                regular_interactions = [i for i in all_interactions if not _is_water_bridge(i)]
                water_bridges = [i for i in all_interactions if _is_water_bridge(i)]

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
                        f.write(f"  Distance: {bridge.get_donor_acceptor_distance():.2f} Å\n\n")

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

    if hasattr(analyzer, "water_bridges") and analyzer.water_bridges:
        wb_file = directory / f"{base_name}_water_bridges.csv"
        write_water_bridges_csv(analyzer, wb_file)

    # Export ligand interactions - one CSV per ligand (regular interactions only)
    # and one CSV per ligand for water bridges
    if hasattr(analyzer, "ligand_interactions") and analyzer.ligand_interactions:
        for ligand_residue in analyzer.ligand_interactions.ligand_info.keys():
            # Regular interactions
            ligand_file = directory / f"{base_name}_ligand_{ligand_residue.replace(':', '_')}.csv"
            write_ligand_interactions_csv(analyzer, ligand_file, ligand_residue)

            # Water bridges
            wb_file = directory / f"{base_name}_ligand_{ligand_residue.replace(':', '_')}_water_bridges.csv"
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

    if hasattr(analyzer, "water_bridges") and analyzer.water_bridges:
        wb_file = directory / f"{base_name}_water_bridges.json"
        write_water_bridges_json(analyzer, wb_file, input_file)

    # Export ligand interactions - one JSON per ligand (regular interactions only)
    # and one JSON per ligand for water bridges
    if hasattr(analyzer, "ligand_interactions") and analyzer.ligand_interactions:
        for ligand_residue in analyzer.ligand_interactions.ligand_info.keys():
            # Regular interactions
            ligand_file = directory / f"{base_name}_ligand_{ligand_residue.replace(':', '_')}.json"
            write_ligand_interactions_json(analyzer, ligand_file, ligand_residue, input_file)

            # Water bridges
            wb_file = directory / f"{base_name}_ligand_{ligand_residue.replace(':', '_')}_water_bridges.json"
            write_ligand_water_bridges_json(analyzer, wb_file, ligand_residue, input_file)


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
        "ligand_interactions": [],
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

    # Ligand interactions organized by ligand
    if hasattr(analyzer, "ligand_interactions") and analyzer.ligand_interactions:
        for ligand_res in sorted(analyzer.ligand_interactions.ligand_info.keys()):
            ligand_info = analyzer.ligand_interactions.ligand_info[ligand_res]
            all_interactions = analyzer.ligand_interactions.get_interactions_for_ligand(ligand_res)

            # Separate regular interactions and water bridges
            regular_interactions = [i for i in all_interactions if not _is_water_bridge(i)]
            water_bridges = [i for i in all_interactions if _is_water_bridge(i)]

            ligand_data = {
                "ligand_residue": ligand_res,
                "ligand_info": ligand_info,
                "interactions": [],
                "water_bridges": [],
            }

            # Add regular interactions
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

                    ligand_data["interactions"].append(interaction_data)
                except Exception:
                    continue

            # Add water bridges
            for wb in water_bridges:
                try:
                    water_bridge_data = {
                        "start_residue": wb.get_donor_residue(),
                        "water_molecules": wb.water_residues,
                        "end_residue": wb.get_acceptor_residue(),
                        "bridge_length": wb.bridge_length,
                        "distance_angstrom": float(f"{wb.get_donor_acceptor_distance():.2f}"),
                    }

                    if hasattr(wb, 'donor_acceptor_properties'):
                        water_bridge_data["properties"] = wb.donor_acceptor_properties

                    ligand_data["water_bridges"].append(water_bridge_data)
                except Exception:
                    continue

            data["ligand_interactions"].append(ligand_data)

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


def write_water_bridges_csv(
    analyzer: NPMolecularInteractionAnalyzer, filename: Path
) -> None:
    """Write water bridges to CSV file.

    :param analyzer: Analyzer instance with water bridge results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output CSV file path
    :type filename: Path
    :returns: None
    :rtype: None
    """
    with open(filename, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "Donor_Residue",
                "Water_Molecules",
                "Acceptor_Residue",
                "Hops",
                "Distance_Angstrom",
            ]
        )

        for wb in analyzer.water_bridges:
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


def write_water_bridges_json(
    analyzer: NPMolecularInteractionAnalyzer,
    filename: Path,
    input_file: Optional[str] = None,
) -> None:
    """Write water bridges to JSON file.

    :param analyzer: Analyzer instance with water bridge results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param filename: Output JSON file path
    :type filename: Path
    :param input_file: Original input file path (for metadata)
    :type input_file: Optional[str]
    :returns: None
    :rtype: None
    """
    data = {
        "version": __version__,
        "input_file": input_file,
        "interaction_type": "water_bridges",
        "bridges": [],
    }

    for wb in analyzer.water_bridges:
        water_residues = [
            {
                "chain_id": parts[0],
                "res_name": parts[1],
                "res_seq": int(parts[2]),
            }
            for parts in [res.split(":") for res in wb.water_residues]
        ]
        bridge_data = {
            "donor_residue": wb.get_donor_residue(),
            "acceptor_residue": wb.get_acceptor_residue(),
            "bridge_length": wb.bridge_length,
            "water_residues": water_residues,
            "donor_acceptor_distance": f"{wb.get_donor_acceptor_distance():.2f}",
        }
        data["bridges"].append(bridge_data)

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


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
    all_interactions = analyzer.ligand_interactions.get_interactions_for_ligand(ligand_residue)

    # Filter out water bridges (they go in separate files)
    interactions = [i for i in all_interactions if not _is_water_bridge(i)]

    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "ligand_residue": ligand_residue,
            "ligand_info": ligand_info,
        },
        "interactions": [],
    }

    for interaction in interactions:
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

            data["interactions"].append(interaction_data)
        except Exception:
            continue

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def _is_water_bridge(interaction) -> bool:
    """Check if an interaction is a water bridge.

    Water bridges have the water_residues attribute from the WaterBridge class.

    :param interaction: MolecularInteraction object
    :returns: True if interaction is a water bridge
    :rtype: bool
    """
    return hasattr(interaction, 'water_residues')


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
    if hasattr(interaction, 'distance'):
        distance_str = f"{interaction.distance:.2f}"
    elif hasattr(interaction, '_distance'):
        distance_str = f"{interaction._distance:.2f}"

    # Get angle metric if available
    angle_str = "N/A"
    if hasattr(interaction, 'angle'):
        angle_str = f"{math.degrees(interaction.angle):.1f}"
    elif hasattr(interaction, 'plane_angle'):
        angle_str = f"{interaction.plane_angle:.1f}"
    elif hasattr(interaction, 'burgi_dunitz_angle'):
        angle_str = f"{interaction.burgi_dunitz_angle:.1f}"
    elif hasattr(interaction, 'angle_to_plane'):
        angle_str = f"{interaction.angle_to_plane:.1f}"

    # Get properties
    properties = ""
    if hasattr(interaction, 'donor_acceptor_properties'):
        properties = interaction.donor_acceptor_properties

    # Get atom names
    donor_atom_name = donor_atom.name if hasattr(donor_atom, 'name') else "N/A"
    acceptor_atom_name = acceptor_atom.name if hasattr(acceptor_atom, 'name') else "N/A"

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
        writer.writerow([
            "Interaction_Type",
            "Donor_Residue",
            "Donor_Atom",
            "Acceptor_Residue",
            "Acceptor_Atom",
            "Distance_Angstrom",
            "Angle_Or_Metric_Degrees",
            "Properties",
        ])

        # Get interactions from analyzer's ligand_interactions container
        if ligand_residue:
            all_interactions = analyzer.ligand_interactions.get_interactions_for_ligand(ligand_residue)
        else:
            all_interactions = analyzer.ligand_interactions.interactions

        # Filter out water bridges (they go in separate files)
        interactions = [i for i in all_interactions if not _is_water_bridge(i)]

        # Write interactions
        for interaction in interactions:
            try:
                row = _format_ligand_interaction_row(interaction)
                writer.writerow([
                    row["type_label"],
                    row["donor_res"],
                    row["donor_atom"],
                    row["acceptor_res"],
                    row["acceptor_atom"],
                    row["distance"],
                    row["angle"],
                    row["properties"],
                ])
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
        writer.writerow([
            "Donor_Residue",
            "Water_Molecules",
            "Acceptor_Residue",
            "Hops",
            "Distance_Angstrom",
        ])

        # Get interactions from analyzer's ligand_interactions container
        if ligand_residue:
            all_interactions = analyzer.ligand_interactions.get_interactions_for_ligand(ligand_residue)
        else:
            all_interactions = analyzer.ligand_interactions.interactions

        # Filter to only water bridges
        water_bridges = [i for i in all_interactions if _is_water_bridge(i)]

        # Write water bridges using shared formatting function
        for wb in water_bridges:
            try:
                row = _format_water_bridge_row(wb)
                writer.writerow([
                    row["donor_res"],
                    row["water_residues"],
                    row["acceptor_res"],
                    row["bridge_length"],
                    row["distance"],
                ])
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
    all_interactions = analyzer.ligand_interactions.get_interactions_for_ligand(ligand_residue)

    # Filter to only water bridges
    water_bridges = [i for i in all_interactions if _is_water_bridge(i)]

    data = {
        "metadata": {
            "input_file": input_file or "",
            "analysis_engine": "HBAT",
            "version": __version__,
            "ligand_residue": ligand_residue,
            "ligand_info": ligand_info,
        },
        "water_bridges": [],
    }

    for wb in water_bridges:
        try:
            water_bridge_data = {
                "start_residue": wb.get_donor_residue(),
                "water_molecules": wb.water_residues,
                "end_residue": wb.get_acceptor_residue(),
                "bridge_length": wb.bridge_length,
                "distance_angstrom": float(f"{wb.get_donor_acceptor_distance():.2f}"),
            }

            if hasattr(wb, 'donor_acceptor_properties'):
                water_bridge_data["properties"] = wb.donor_acceptor_properties

            data["water_bridges"].append(water_bridge_data)
        except Exception:
            continue

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


