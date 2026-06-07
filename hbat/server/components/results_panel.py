"""
Results display component with 3D visualization for HBAT web interface.

This module provides results display with integrated py3Dmol visualization
for molecular interactions.
"""

import math
from pathlib import Path
from typing import Optional, Callable, Dict, Any

import os
import tempfile
import zipfile
import random

from nicegui import ui

from ...core.analysis import NPMolecularInteractionAnalyzer
from ...config import (
    INTERACTION_CONFIGS,
    extract_interaction_data,
    get_interaction_config,
)
from ...visualization.minimal_pdb_extractor import (
    format_minimal_pdb,
    extract_water_bridge_pdb,
)
from ...visualization.chain_graph import render_chain_for_web
from ...visualization import export_interactions_to_pymol
from ...visualization.ligplot import LigplotGenerator
from ...visualization.pymol3d import (
    generate_carbonyl_interaction_viewer_js,
    generate_halogen_bond_viewer_js,
    generate_hydrogen_bond_viewer_js,
    generate_n_pi_interaction_viewer_js,
    generate_pi_interaction_viewer_js,
    generate_pi_pi_stacking_viewer_js,
    generate_water_bridge_viewer_js,
    generate_png_export_js,
    generate_ligand_interactions_viewer_js,
)


class WebResultsPanel:
    """Results panel component with 3D visualization."""

    def __init__(self, container, session_dir_callback=None):
        """Initialize the results panel.

        :param container: Container element for results display
        :param session_dir_callback: Callable that returns the current session directory
        """
        self.container = container
        self.analyzer: Optional[NPMolecularInteractionAnalyzer] = None
        self.current_file: Optional[str] = None
        self.session_dir_callback = session_dir_callback

        # Tab components
        self.tabs = None
        self.tab_panels = None

        # Tab panel placeholders
        self.summary_panel = None
        self.hydrogen_panel = None
        self.halogen_panel = None
        self.pi_panel = None
        self.pi_pi_panel = None
        self.carbonyl_panel = None
        self.n_pi_panel = None
        self.water_bridges_panel = None
        self.cooperativity_panel = None
        self.ligands_panel = None

    def _get_pdb_basename(self) -> str:
        """Get the base name of the current PDB file (without extension).

        :returns: Base name of PDB file
        :rtype: str
        """
        if self.current_file:
            return Path(self.current_file).stem
        return "structure"

    def _download_pymol_visualization(
        self, interaction_label: str, **interaction_kwargs
    ):
        """Generic helper to View in PyMOL visualization for any interaction type.

        :param interaction_label: Label for the interaction (e.g., "A:ALA:1_to_B:GLY:2")
        :param interaction_kwargs: Keyword arguments to pass to export_interactions_to_pymol
                                  (e.g., hydrogen_bonds=[...], water_bridges=[...])
        """
        try:
            # Determine which PDB file to use
            pdb_file_path = None
            if hasattr(
                self.analyzer, "_pdb_fixing_info"
            ) and self.analyzer._pdb_fixing_info.get("applied"):
                pdb_file_path = self.analyzer._pdb_fixing_info.get("fixed_file_path")
            else:
                pdb_file_path = self.analyzer._pdb_fixing_info.get("input_file_path")

            if not pdb_file_path or not os.path.exists(pdb_file_path):
                ui.notify(
                    "PDB file not found for visualization",
                    type="negative",
                    position="top-left",
                )
                return

            # Create temporary directory for PyMOL export
            with tempfile.TemporaryDirectory() as temp_dir:
                # Export interaction with provided kwargs
                success = export_interactions_to_pymol(
                    pdb_file_path=pdb_file_path,
                    parser=self.analyzer.parser,
                    output_dir=temp_dir,
                    **interaction_kwargs,
                )

                if not success:
                    ui.notify(
                        "Failed to generate PyMOL visualization",
                        type="negative",
                        position="top-left",
                    )
                    return

                # Create a zip file with visualization files
                pdb_base = self._get_pdb_basename()
                zip_filename = f"{pdb_base}_{interaction_label}_pymol.zip"

                basename = os.path.splitext(os.path.basename(pdb_file_path))[0]
                pdb_file = os.path.join(temp_dir, f"{basename}.pdb")
                script_file = os.path.join(temp_dir, f"{basename}.pml")
                readme_file = os.path.join(temp_dir, "README.txt")

                # Create zip in temp location then download
                zip_path = os.path.join(temp_dir, zip_filename)
                with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
                    if os.path.exists(pdb_file):
                        zipf.write(pdb_file, os.path.basename(pdb_file))
                    if os.path.exists(script_file):
                        zipf.write(script_file, os.path.basename(script_file))
                    if os.path.exists(readme_file):
                        zipf.write(readme_file, os.path.basename(readme_file))

                # Read and download the zip file
                with open(zip_path, "rb") as f:
                    ui.download(f.read(), filename=zip_filename)

                ui.notify(
                    "PyMOL visualization downloaded",
                    type="positive",
                    position="top-left",
                )

        except Exception as e:
            ui.notify(
                f"PyMOL export failed: {str(e)}",
                type="negative",
                position="top-left",
            )

    async def update_results(
        self, analyzer: NPMolecularInteractionAnalyzer, filename: str
    ):
        """Update results display with analysis results.

        :param analyzer: The analyzer with results
        :param filename: Name of the analyzed file
        """
        self.analyzer = analyzer
        self.current_file = filename

        # Clear container and recreate tabs
        self.container.clear()

        with self.container:
            # Create tabs and panels
            with ui.tabs().classes("w-full") as self.tabs:
                self._create_tabs()

            with ui.tab_panels(self.tabs, value="summary").classes(
                "w-full"
            ) as self.tab_panels:
                self._create_panels()

    def _get_interaction_configs(self):
        """Get configuration for all interaction types from centralized config.

        :returns: List of interaction configurations
        :rtype: list[dict]
        """
        # Mapping of config IDs to visualization handlers and display properties
        visualization_handlers = {
            "hydrogen_bonds": {
                "viz_method": "_show_hydrogen_bond_visualization",
                "icon": "link",
            },
            "water_bridges": {
                "viz_method": "_show_water_bridge_visualization",
                "icon": "opacity",
            },
            "halogen_bonds": {
                "viz_method": "_show_halogen_bond_visualization",
                "icon": "whatshot",
            },
            "pi_interactions": {
                "viz_method": "_show_pi_interaction_visualization",
                "icon": "fiber_manual_record",
            },
            "pi_pi_interactions": {
                "viz_method": "_show_pi_pi_stacking_visualization",
                "icon": "layers",
            },
            "carbonyl_interactions": {
                "viz_method": "_show_carbonyl_interaction_visualization",
                "icon": "architecture",
            },
            "n_pi_interactions": {
                "viz_method": "_show_n_pi_interaction_visualization",
                "icon": "arrow_forward",
            },
            "cooperativity_chains": {
                "viz_method": "_show_cooperativity_chain_visualization",
                "icon": "account_tree",
            },
        }

        configs = [
            {
                "id": "summary",
                "label": "Summary",
                "icon": "summarize",
                "attr": None,
                "always_show": True,
                "panel_attr": "summary_panel",
                "update_method": "_update_summary_panel",
                "config": None,
            },
            {
                "id": "ligands",
                "label": "Ligand Interactions",
                "icon": "donut_large",
                "attr": None,
                "panel_attr": "ligands_panel",
                "update_method": "_update_ligand_interactions_panel",
                "config": None,
            },
        ]

        # Add interaction type configs from centralized system
        for config_id, cfg in INTERACTION_CONFIGS.items():
            handler_info = visualization_handlers.get(config_id, {})
            configs.append(
                {
                    "id": config_id,
                    "label": cfg.label,
                    "icon": handler_info.get("icon", cfg.icon),
                    "attr": cfg.analyzer_attr,
                    "panel_attr": f"{config_id}_panel",
                    "update_method": "_update_interaction_panel",
                    "config": cfg,
                    "viz_method": handler_info.get("viz_method"),
                }
            )

        return configs

    def _should_show_interaction(self, config):
        """Check if interaction should be shown.

        :param config: Interaction configuration
        :returns: True if interaction should be shown
        :rtype: bool
        """
        if config.get("always_show"):
            return True

        # Special case for ligand interactions
        if config.get("id") == "ligands":
            return bool(self.analyzer.ligand_interactions)

        attr = config.get("attr")
        if not attr:
            return False

        if not hasattr(self.analyzer, attr):
            return False

        interactions = getattr(self.analyzer, attr)
        return interactions and len(interactions) > 0

    def _create_tabs(self):
        """Create tabs for all available interactions."""
        for config in self._get_interaction_configs():
            if self._should_show_interaction(config):
                attr = config.get("attr")
                label = config["label"]

                # Add count to label if not summary
                if attr and hasattr(self.analyzer, attr):
                    interactions = getattr(self.analyzer, attr)
                    if interactions:
                        label = f"{label} ({len(interactions)})"

                ui.tab(config["id"], label=label, icon=config["icon"])

    def _create_panels(self):
        """Create panels for all available interactions."""
        for config in self._get_interaction_configs():
            if self._should_show_interaction(config):
                with ui.tab_panel(config["id"]):
                    panel = ui.column().classes("w-full")
                    setattr(self, config["panel_attr"], panel)

                    # Call the update method (pass config for centralized interactions)
                    update_method = getattr(self, config["update_method"])
                    if config["update_method"] == "_update_interaction_panel":
                        update_method(config)
                    else:
                        update_method()

    def _update_summary_panel(self):
        """Update summary statistics panel."""
        with self.summary_panel:
            ui.label(f"Analysis Summary For: {self.current_file}").classes("text-h5")

            summary = self.analyzer.get_summary()

            # Create summary cards
            with ui.row().classes("w-full gap-2"):
                with ui.card().classes("flex-1"):
                    ui.label("Hydrogen Bonds").classes("text-h6")
                    ui.label(str(summary["hydrogen_bonds"]["count"])).classes(
                        "text-h4 text-primary"
                    )

                with ui.card().classes("flex-1"):
                    ui.label("Halogen Bonds").classes("text-h6")
                    ui.label(str(summary["halogen_bonds"]["count"])).classes(
                        "text-h4 text-orange"
                    )

                with ui.card().classes("flex-1"):
                    ui.label("π Inter").classes("text-h6")
                    ui.label(str(summary["pi_interactions"]["count"])).classes(
                        "text-h4 text-purple"
                    )

            with ui.row().classes("w-full gap-2"):
                with ui.card().classes("flex-1"):
                    ui.label("π-π Stacking").classes("text-h6")
                    ui.label(
                        str(summary.get("pi_pi_interactions", {}).get("count", 0))
                    ).classes("text-h4 text-indigo")

                with ui.card().classes("flex-1"):
                    ui.label("Carbonyl n→π*").classes("text-h6")
                    ui.label(
                        str(summary.get("carbonyl_interactions", {}).get("count", 0))
                    ).classes("text-h4 text-teal")

                with ui.card().classes("flex-1"):
                    ui.label("n→π* Inter").classes("text-h6")
                    ui.label(
                        str(summary.get("n_pi_interactions", {}).get("count", 0))
                    ).classes("text-h4 text-green")

                # Ligand Interactions count
                with ui.card().classes("flex-1"):
                    ui.label("Ligand Inter").classes("text-h6")
                    ui.label(
                        str(summary.get("ligand_interactions", {}).get("count", 0))
                    ).classes("text-h4 text-amber")

            with ui.row().classes("w-full gap-2"):
                with ui.card().classes("flex-1"):
                    ui.label("Potential Chains").classes("text-h6")
                    ui.label(str(summary["cooperativity_chains"]["count"])).classes(
                        "text-h4 text-blue"
                    )

                pdb_fix = summary["pdb_fixing"]
                original = pdb_fix.get("original_atoms", 0)
                fixed = pdb_fix.get("fixed_atoms", 0)
                added_h = pdb_fix.get("added_hydrogens", 0)
                with ui.card().classes("flex-1"):
                    ui.label("Atoms Fixed").classes("text-h6")
                    ui.label(str(f"{fixed - original}")).classes("text-h4 text-blue")

                with ui.card().classes("flex-1"):
                    ui.label("Hydrogen Added").classes("text-h6")
                    ui.label(str(f"{added_h}")).classes("text-h4 text-blue")

                with ui.card().classes("flex-1"):
                    ui.label("Water Bridges").classes("text-h6")
                    ui.label(str(summary["water_bridges"]["count"])).classes(
                        "text-h4 text-cyan"
                    )

            # Bond Detection Information
            if "bond_detection" in summary:
                bond_det = summary["bond_detection"]
                breakdown = bond_det.get("breakdown", {})

                with ui.row().classes("w-full gap-2 q-mt-md"):
                    # Total Bonds Detected
                    with ui.card().classes("flex-1"):
                        ui.label("Bonds Detected").classes("text-h6")
                        ui.label(str(bond_det.get("total_bonds", 0))).classes(
                            "text-h4 text-primary"
                        )

                    # CONECT Records
                    with ui.card().classes("flex-1"):
                        ui.label("CONECT").classes("text-h6")
                        if "conect_records" in breakdown:
                            conect = breakdown["conect_records"]
                            ui.label(f"{conect['percentage']}%").classes(
                                "text-h4 text-orange"
                            )
                        else:
                            ui.label("0.0%").classes("text-h4 text-orange")

                    # Residue Lookup
                    with ui.card().classes("flex-1"):
                        ui.label("Residue Lookup").classes("text-h6")
                        if "residue_lookup" in breakdown:
                            residue = breakdown["residue_lookup"]
                            ui.label(f"{residue['percentage']}%").classes(
                                "text-h4 text-purple"
                            )
                        else:
                            ui.label("0.0%").classes("text-h4 text-purple")

                    # Distance Based
                    with ui.card().classes("flex-1"):
                        ui.label("Distance Based").classes("text-h6")
                        if "distance_based" in breakdown:
                            distance = breakdown["distance_based"]
                            ui.label(f"{distance['percentage']}%").classes(
                                "text-h4 text-teal"
                            )
                        else:
                            ui.label("0.0%").classes("text-h4 text-teal")

    def _update_interaction_panel(self, config: Dict[str, Any]):
        """Generic panel updater for all interaction types using centralized config.

        Uses the InteractionConfig from INTERACTION_CONFIGS to dynamically
        generate tables with columns and data extraction.

        :param config: Configuration dict with 'config', 'panel_attr', 'id', 'viz_method'
        """
        interaction_config = config.get("config")
        if not interaction_config:
            return

        panel = getattr(self, config["panel_attr"])
        interaction_type = config["id"]
        analyzer_attr = interaction_config.analyzer_attr
        interactions = getattr(self.analyzer, analyzer_attr, [])

        if not interactions:
            return

        with panel:
            # Title with count
            ui.label(f"{interaction_config.label} ({len(interactions)})").classes(
                "text-h5"
            )

            # Build columns from config
            columns = [
                {
                    "name": "visualize",
                    "label": "3D",
                    "field": "visualize",
                    "align": "center",
                }
            ]

            # Add columns from config (skip any that have accessor, they'll be computed)
            for col in interaction_config.columns:
                columns.append(
                    {
                        "name": col.name,
                        "label": col.label,
                        "field": col.name,
                        "align": "right" if col.precision is not None else "left",
                    }
                )

            # Build rows using extract_interaction_data
            rows = []
            for idx, interaction in enumerate(interactions):
                row = {"id": idx, "visualize": ""}
                data = extract_interaction_data(interaction, interaction_config, idx)
                row.update(data)
                rows.append(row)

            # Create filter input
            filter_input = (
                ui.input(placeholder="Filter interactions (residue, atom, type, etc.)")
                .props("clearable outlined dense")
                .classes("w-full mb-2")
            )

            # Create table
            table = ui.table(columns=columns, rows=rows, row_key="id").classes("w-full")
            filter_input.bind_value(table, "filter")

            # Add visualization button column (only if viz method exists)
            viz_method = config.get("viz_method")
            if viz_method and hasattr(self, viz_method):
                table.add_slot(
                    "body-cell-visualize",
                    """
                    <q-td :props="props">
                        <q-btn size="sm" color="primary" round dense icon="visibility"
                               @click="$parent.$emit('visualize', props.row)" />
                    </q-td>
                """,
                )

                # Create visualization handler
                def show_visualization(e, viz_func=getattr(self, viz_method)):
                    idx = e.args["id"]
                    interaction = interactions[idx]
                    viz_func(interaction)

                table.on("visualize", show_visualization)

    def _show_hydrogen_bond_visualization(self, hb):
        """Show 3D visualization of a hydrogen bond in a dialog.

        :param hb: Hydrogen bond interaction
        """
        # Create unique viewer ID
        import random

        viewer_id = f"hb_viewer_{random.randint(1000, 9999)}"

        # Extract minimal PDB early so it can be used for download
        minimal_pdb = format_minimal_pdb(self.analyzer.parser, [hb])

        def export_png():
            """Export viewer as PNG."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{hb.get_donor_residue()}_to_{hb.get_acceptor_residue()}_hbond.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        def download_pdb():
            """Download extracted PDB file."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{hb.get_donor_residue()}_to_{hb.get_acceptor_residue()}_hbond.pdb".replace(
                ":", "_"
            )
            ui.download(minimal_pdb.encode(), filename=filename)

        def download_pymol():
            """Download hydrogen bond as PyMOL visualization package."""
            hbond_label = (
                f"{hb.get_donor_residue()}_to_{hb.get_acceptor_residue()}".replace(
                    ":", "_"
                )
            )
            self._download_pymol_visualization(hbond_label, hydrogen_bonds=[hb])

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
                    "outline color=secondary"
                )
                ui.button("Download PDB", icon="download", on_click=download_pdb).props(
                    "outline color=secondary"
                )
                ui.button("View in PyMOL", icon="movie", on_click=download_pymol).props(
                    "outline color=secondary"
                )
                ui.button("Close", on_click=dialog.close).props("color=primary")

            with ui.card_section().classes("q-pa-md"):
                ui.label(
                    f"{hb.get_donor_residue()} → {hb.get_acceptor_residue()}"
                ).classes("text-subtitle1 q-mb-md")

                # Create div for viewer - container must have actual width
                ui.html(
                    f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>',
                    sanitize=False,
                )

        dialog.open()

        # Initialize viewer after dialog is opened and rendered
        ui.timer(
            0.1,
            lambda: self._initialize_hydrogen_bond_viewer(viewer_id, hb, minimal_pdb),
            once=True,
        )

    def _initialize_hydrogen_bond_viewer(self, viewer_id: str, hb, minimal_pdb: str):
        """Initialize py3Dmol viewer for hydrogen bond visualization.

        :param viewer_id: Unique ID for the viewer div
        :param hb: Hydrogen bond interaction
        :param minimal_pdb: Pre-extracted minimal PDB content
        """
        javascript = generate_hydrogen_bond_viewer_js(hb, minimal_pdb, viewer_id)
        ui.run_javascript(javascript)

    def _show_halogen_bond_visualization(self, xb):
        """Show 3D visualization of a halogen bond in a dialog."""
        # Create unique viewer ID
        import random

        viewer_id = f"xb_viewer_{random.randint(1000, 9999)}"

        # Extract minimal PDB early so it can be used for download
        minimal_pdb = format_minimal_pdb(self.analyzer.parser, [xb])

        def export_png():
            """Export viewer as PNG."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{xb.get_donor_residue()}_to_{xb.get_acceptor_residue()}_xbond.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        def download_pdb():
            """Download extracted PDB file."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{xb.get_donor_residue()}_to_{xb.get_acceptor_residue()}_xbond.pdb".replace(
                ":", "_"
            )
            ui.download(minimal_pdb.encode(), filename=filename)

        def download_pymol():
            """Download halogen bond as PyMOL visualization package."""
            xbond_label = (
                f"{xb.get_donor_residue()}_to_{xb.get_acceptor_residue()}".replace(
                    ":", "_"
                )
            )
            self._download_pymol_visualization(xbond_label, halogen_bonds=[xb])

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
                    "outline color=secondary"
                )
                ui.button("Download PDB", icon="download", on_click=download_pdb).props(
                    "outline color=secondary"
                )
                ui.button("View in PyMOL", icon="movie", on_click=download_pymol).props(
                    "outline color=secondary"
                )
                ui.button("Close", on_click=dialog.close).props("color=primary")

            with ui.card_section().classes("q-pa-md"):
                ui.label(
                    f"{xb.get_donor_residue()} → {xb.get_acceptor_residue()}"
                ).classes("text-subtitle1 q-mb-md")

                # Create div for viewer - container must have actual width
                ui.html(
                    f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>',
                    sanitize=False,
                )

        dialog.open()

        # Initialize viewer after dialog is opened and rendered
        ui.timer(
            0.1,
            lambda: self._initialize_halogen_bond_viewer(viewer_id, xb, minimal_pdb),
            once=True,
        )

    def _initialize_halogen_bond_viewer(self, viewer_id: str, xb, minimal_pdb: str):
        """Initialize py3Dmol viewer for halogen bond visualization.

        :param viewer_id: Unique ID for the viewer div
        :param xb: Halogen bond interaction
        :param minimal_pdb: Pre-extracted minimal PDB content
        """
        javascript = generate_halogen_bond_viewer_js(xb, minimal_pdb, viewer_id)
        ui.run_javascript(javascript)

    def _show_pi_interaction_visualization(self, pi):
        """Show 3D visualization of a π interaction in a dialog."""
        import random

        viewer_id = f"pi_viewer_{random.randint(1000, 9999)}"

        # Extract minimal PDB early so it can be used for download
        minimal_pdb = format_minimal_pdb(self.analyzer.parser, [pi])

        def export_png():
            """Export viewer as PNG."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{pi.get_donor_residue()}_to_{pi.get_acceptor_residue()}_pi.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        def download_pdb():
            """Download extracted PDB file."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{pi.get_donor_residue()}_to_{pi.get_acceptor_residue()}_pi.pdb".replace(
                ":", "_"
            )
            ui.download(minimal_pdb.encode(), filename=filename)

        def download_pymol():
            """Download pi interaction as PyMOL visualization package."""
            pi_label = (
                f"{pi.get_donor_residue()}_to_{pi.get_acceptor_residue()}".replace(
                    ":", "_"
                )
            )
            self._download_pymol_visualization(pi_label, pi_interactions=[pi])

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
                    "outline color=secondary"
                )
                ui.button("Download PDB", icon="download", on_click=download_pdb).props(
                    "outline color=secondary"
                )
                ui.button("View in PyMOL", icon="movie", on_click=download_pymol).props(
                    "outline color=secondary"
                )
                ui.button("Close", on_click=dialog.close).props("color=primary")

            with ui.card_section().classes("q-pa-md"):
                ui.label(
                    f"{pi.get_donor_residue()} → {pi.get_acceptor_residue()}"
                ).classes("text-subtitle1 q-mb-md")

                ui.html(
                    f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>',
                    sanitize=False,
                )

        dialog.open()
        ui.timer(
            0.1,
            lambda: self._initialize_pi_interaction_viewer(viewer_id, pi, minimal_pdb),
            once=True,
        )

    def _initialize_pi_interaction_viewer(self, viewer_id: str, pi, minimal_pdb: str):
        """Initialize py3Dmol viewer for π interaction visualization.

        :param viewer_id: Unique ID for the viewer div
        :param pi: Pi interaction
        :param minimal_pdb: Pre-extracted minimal PDB content
        """
        javascript = generate_pi_interaction_viewer_js(pi, minimal_pdb, viewer_id)
        ui.run_javascript(javascript)

    def _show_water_bridge_visualization(self, wb):
        """Show 3D visualization of a water bridge in a dialog.

        :param wb: Water bridge interaction
        """
        import random

        viewer_id = f"wb_viewer_{random.randint(1000, 9999)}"

        # Extract minimal PDB with only water bridge residues (donor, water, acceptor)
        minimal_pdb = extract_water_bridge_pdb(self.analyzer.parser, wb)

        def export_png():
            """Export viewer as PNG."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{wb.get_donor_residue()}_to_{wb.get_acceptor_residue()}_wbridge.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        def download_pdb():
            """Download extracted PDB file."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{wb.get_donor_residue()}_to_{wb.get_acceptor_residue()}_wbridge.pdb".replace(
                ":", "_"
            )
            ui.download(minimal_pdb.encode(), filename=filename)

        def download_pymol():
            """Download water bridge as PyMOL visualization package."""
            wbridge_label = (
                f"{wb.get_donor_residue()}_to_{wb.get_acceptor_residue()}".replace(
                    ":", "_"
                )
            )
            self._download_pymol_visualization(wbridge_label, water_bridges=[wb])

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 1000px; max-width: 95vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
                    "outline color=secondary"
                )
                ui.button("Download PDB", icon="download", on_click=download_pdb).props(
                    "outline color=secondary"
                )
                ui.button("View in PyMOL", icon="movie", on_click=download_pymol).props(
                    "outline color=secondary"
                )
                ui.button("Close", on_click=dialog.close).props("color=primary")

            with ui.card_section().classes("q-pa-md"):
                water_residues_display = ", ".join(wb.water_residues)
                ui.label(
                    f"{wb.get_donor_residue()} → "
                    f"({water_residues_display}) → "
                    f"{wb.get_acceptor_residue()}"
                ).classes("text-subtitle1 q-mb-md")

                ui.label(
                    f"Bridge Length: {wb.bridge_length} hop(s) | "
                    f"Distance: {wb.get_donor_acceptor_distance():.2f}Å"
                ).classes("text-caption text-grey q-mb-md")

                # Create div for viewer - container must have actual width
                ui.html(
                    f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>',
                    sanitize=False,
                )

        dialog.open()

        # Initialize viewer after dialog is opened and rendered
        ui.timer(
            0.1,
            lambda: self._initialize_water_bridge_viewer(viewer_id, wb, minimal_pdb),
            once=True,
        )

    def _initialize_water_bridge_viewer(self, viewer_id: str, wb, minimal_pdb: str):
        """Initialize py3Dmol viewer for water bridge visualization.

        :param viewer_id: Unique ID for the viewer div
        :param wb: Water bridge interaction
        :param minimal_pdb: Extracted PDB content as string
        """
        javascript = generate_water_bridge_viewer_js(wb, minimal_pdb, viewer_id)
        ui.run_javascript(javascript)

    def _show_pi_pi_stacking_visualization(self, pi_pi):
        """Show 3D visualization of π-π stacking in a dialog."""
        import random

        viewer_id = f"pipi_viewer_{random.randint(1000, 9999)}"

        # Extract minimal PDB early so it can be used for download
        minimal_pdb = format_minimal_pdb(self.analyzer.parser, [pi_pi])

        def export_png():
            """Export viewer as PNG."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{pi_pi.ring1_residue}_to_{pi_pi.ring2_residue}_pipi.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        def download_pdb():
            """Download extracted PDB file."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{pi_pi.ring1_residue}_to_{pi_pi.ring2_residue}_pipi.pdb".replace(
                ":", "_"
            )
            ui.download(minimal_pdb.encode(), filename=filename)

        def download_pymol():
            """Download pi-pi stacking as PyMOL visualization package."""
            pipi_label = f"{pi_pi.ring1_residue}_to_{pi_pi.ring2_residue}".replace(
                ":", "_"
            )
            self._download_pymol_visualization(pipi_label, pi_pi_stacking=[pi_pi])

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
                    "outline color=secondary"
                )
                ui.button("Download PDB", icon="download", on_click=download_pdb).props(
                    "outline color=secondary"
                )
                ui.button("View in PyMOL", icon="movie", on_click=download_pymol).props(
                    "outline color=secondary"
                )
                ui.button("Close", on_click=dialog.close).props("color=primary")

            with ui.card_section().classes("q-pa-md"):
                ui.label(f"{pi_pi.ring1_residue} ⇄ {pi_pi.ring2_residue}").classes(
                    "text-subtitle1 q-mb-md"
                )

                ui.html(
                    f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>',
                    sanitize=False,
                )

        dialog.open()
        ui.timer(
            0.1,
            lambda: self._initialize_pi_pi_stacking_viewer(
                viewer_id, pi_pi, minimal_pdb
            ),
            once=True,
        )

    def _initialize_pi_pi_stacking_viewer(
        self, viewer_id: str, pi_pi, minimal_pdb: str
    ):
        """Initialize py3Dmol viewer for π-π stacking visualization.

        :param viewer_id: Unique ID for the viewer div
        :param pi_pi: Pi-pi stacking interaction
        :param minimal_pdb: Pre-extracted minimal PDB content
        """
        javascript = generate_pi_pi_stacking_viewer_js(pi_pi, minimal_pdb, viewer_id)
        ui.run_javascript(javascript)

    def _show_carbonyl_interaction_visualization(self, carbonyl):
        """Show 3D visualization of carbonyl n→π* interaction in a dialog."""
        import random

        viewer_id = f"carbonyl_viewer_{random.randint(1000, 9999)}"

        # Extract minimal PDB early so it can be used for download
        minimal_pdb = format_minimal_pdb(self.analyzer.parser, [carbonyl])

        def export_png():
            """Export viewer as PNG."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{carbonyl.get_donor_residue()}_to_{carbonyl.get_acceptor_residue()}_carbonyl.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        def download_pdb():
            """Download extracted PDB file."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{carbonyl.get_donor_residue()}_to_{carbonyl.get_acceptor_residue()}_carbonyl.pdb".replace(
                ":", "_"
            )
            ui.download(minimal_pdb.encode(), filename=filename)

        def download_pymol():
            """Download carbonyl interaction as PyMOL visualization package."""
            carbonyl_label = f"{carbonyl.get_donor_residue()}_to_{carbonyl.get_acceptor_residue()}".replace(
                ":", "_"
            )
            self._download_pymol_visualization(
                carbonyl_label, carbonyl_interactions=[carbonyl]
            )

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
                    "outline color=secondary"
                )
                ui.button("Download PDB", icon="download", on_click=download_pdb).props(
                    "outline color=secondary"
                )
                ui.button("View in PyMOL", icon="movie", on_click=download_pymol).props(
                    "outline color=secondary"
                )
                ui.button("Close", on_click=dialog.close).props("color=primary")

            with ui.card_section().classes("q-pa-md"):
                ui.label(
                    f"{carbonyl.get_donor_residue()} → {carbonyl.get_acceptor_residue()}"
                ).classes("text-subtitle1 q-mb-md")

                ui.html(
                    f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>',
                    sanitize=False,
                )

        dialog.open()
        ui.timer(
            0.1,
            lambda: self._initialize_carbonyl_interaction_viewer(
                viewer_id, carbonyl, minimal_pdb
            ),
            once=True,
        )

    def _initialize_carbonyl_interaction_viewer(
        self, viewer_id: str, carbonyl, minimal_pdb: str
    ):
        """Initialize py3Dmol viewer for carbonyl n→π* visualization.

        :param viewer_id: Unique ID for the viewer div
        :param carbonyl: Carbonyl interaction
        :param minimal_pdb: Pre-extracted minimal PDB content
        """
        javascript = generate_carbonyl_interaction_viewer_js(
            carbonyl, minimal_pdb, viewer_id
        )
        ui.run_javascript(javascript)

    def _show_n_pi_interaction_visualization(self, n_pi):
        """Show 3D visualization of n→π* interaction in a dialog."""
        import random

        viewer_id = f"npi_viewer_{random.randint(1000, 9999)}"

        # Extract minimal PDB early so it can be used for download
        minimal_pdb = format_minimal_pdb(self.analyzer.parser, [n_pi])

        def export_png():
            """Export viewer as PNG."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{n_pi.get_donor_residue()}_to_{n_pi.get_acceptor_residue()}_npi.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        def download_pdb():
            """Download extracted PDB file."""
            pdb_base = self._get_pdb_basename()
            filename = f"{pdb_base}_{n_pi.get_donor_residue()}_to_{n_pi.get_acceptor_residue()}_npi.pdb".replace(
                ":", "_"
            )
            ui.download(minimal_pdb.encode(), filename=filename)

        def download_pymol():
            """Download N-pi interaction as PyMOL visualization package."""
            npi_label = (
                f"{n_pi.get_donor_residue()}_to_{n_pi.get_acceptor_residue()}".replace(
                    ":", "_"
                )
            )
            self._download_pymol_visualization(npi_label, n_pi_interactions=[n_pi])

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
                    "outline color=secondary"
                )
                ui.button("Download PDB", icon="download", on_click=download_pdb).props(
                    "outline color=secondary"
                )
                ui.button("View in PyMOL", icon="movie", on_click=download_pymol).props(
                    "outline color=secondary"
                )
                ui.button("Close", on_click=dialog.close).props("color=primary")

            with ui.card_section().classes("q-pa-md"):
                ui.label(
                    f"{n_pi.get_donor_residue()} → {n_pi.get_acceptor_residue()}"
                ).classes("text-subtitle1 q-mb-md")

                ui.html(
                    f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>',
                    sanitize=False,
                )

        dialog.open()
        ui.timer(
            0.1,
            lambda: self._initialize_n_pi_interaction_viewer(
                viewer_id, n_pi, minimal_pdb
            ),
            once=True,
        )

    def _initialize_n_pi_interaction_viewer(
        self, viewer_id: str, n_pi, minimal_pdb: str
    ):
        """Initialize py3Dmol viewer for n→π* interaction visualization.

        :param viewer_id: Unique ID for the viewer div
        :param n_pi: N-pi interaction
        :param minimal_pdb: Pre-extracted minimal PDB content
        """
        javascript = generate_n_pi_interaction_viewer_js(n_pi, minimal_pdb, viewer_id)
        ui.run_javascript(javascript)

    def _show_cooperativity_chain_visualization(self, chain):
        """Show 2D graph visualization of cooperativity chain in a dialog."""
        try:
            # Get session directory from callback or fallback to SESSIONS_BASE_DIR
            if self.session_dir_callback:
                session_dir = self.session_dir_callback()
            else:
                from ...server.app import SESSIONS_BASE_DIR

                session_dir = SESSIONS_BASE_DIR

            # Generate both SVG (for display) and PNG (for export) using reusable function
            # Replace spaces and special chars in chain type for clean filenames
            pdb_base = self._get_pdb_basename()
            safe_chain_type = chain.chain_type.replace(" ", "_").replace("-", "_")
            filename_prefix = f"{pdb_base}_chain_{safe_chain_type}_{chain.chain_length}"
            svg_content, png_path = render_chain_for_web(
                chain,
                output_dir=session_dir,
                filename_prefix=filename_prefix,
                engine="dot",
                rankdir="LR",
            )

            if svg_content:

                def export_png():
                    """Export chain graph as PNG."""
                    if png_path and png_path.exists():
                        ui.download(str(png_path))
                        ui.notify(
                            "Downloaded chain graph as PNG",
                            type="positive",
                            position="top-left",
                        )
                    else:
                        ui.notify(
                            "PNG file not found", type="warning", position="top-left"
                        )

                with (
                    ui.dialog().props("persistent") as dialog,
                    ui.card().style("width: 95vw; max-width: 1400px;"),
                ):
                    with ui.card_actions().classes("justify-end"):
                        ui.button(
                            "Export PNG", icon="download", on_click=export_png
                        ).props("outline color=secondary")
                        ui.button("Close", on_click=dialog.close).props("color=primary")

                    with ui.card_section().classes("q-pa-md"):
                        ui.label(
                            f"Chain: {chain.chain_type} (Length: {chain.chain_length})"
                        ).classes("text-subtitle1 q-mb-md")

                        # Display the SVG
                        ui.html(svg_content, sanitize=False).style(
                            "width: 100%; overflow: auto;"
                        )

                dialog.open()
            else:
                ui.notify(
                    "Failed to generate chain graph",
                    type="negative",
                    position="top-left",
                )

        except ImportError:
            ui.notify(
                "Graphviz is not installed. Install with: pip install graphviz",
                type="warning",
                position="top-left",
            )
        except Exception as e:
            ui.notify(
                f"Error generating chain graph: {str(e)}",
                type="negative",
                position="top-left",
            )

    def _update_ligand_interactions_panel(self):
        """Update ligand interactions panel with separate sections for regular and water bridge interactions."""
        import random

        # Use pre-computed ligand interactions from the analyzer
        ligand_container = self.analyzer.ligand_interactions

        if not ligand_container:
            with self.ligands_panel:
                ui.label("No ligand interactions found").classes(
                    "text-subtitle1 q-pa-lg"
                )
            return

        # Get ligand interactions and info from the container
        ligand_interactions = ligand_container.interactions
        ligand_info = ligand_container.ligand_info

        # Separate regular interactions from water bridges
        regular_interactions = [
            i for i in ligand_interactions if not hasattr(i, "water_residues")
        ]
        water_bridge_interactions = [
            i for i in ligand_interactions if hasattr(i, "water_residues")
        ]

        with self.ligands_panel:
            # Title
            ui.label(f"Ligand Interactions ({len(ligand_interactions)})").classes(
                "text-h5"
            )

            # Dropdown selector - show each ligand with interaction count
            # NiceGUI expects {value: label} format
            ligand_options = {}
            first_ligand = None
            # Sort by count (descending)
            for ligand_res, info in sorted(
                ligand_info.items(), key=lambda x: -x[1]["count"]
            ):
                label = f"{ligand_res} ({info['count']})"
                ligand_options[ligand_res] = label
                if first_ligand is None:
                    first_ligand = ligand_res

            # Store interactions in dicts for retrieval (can't serialize to JSON)
            # Separate stores to avoid ID collisions between regular and water bridge interactions
            interaction_store = {}
            water_bridge_store = {}

            # Build table data for regular interactions
            def build_regular_interactions_data(selected_value):
                """Build table data for regular interactions of the selected ligand."""
                rows = []
                interaction_index = 0

                for interaction in regular_interactions:
                    donor_res = interaction.get_donor_residue()
                    acceptor_res = interaction.get_acceptor_residue()

                    # Show only interactions involving the selected ligand
                    if donor_res != selected_value and acceptor_res != selected_value:
                        continue

                    # Determine interaction type
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

                    # Get distance metric
                    if hasattr(interaction, "distance"):
                        distance = f"{interaction.distance:.2f} Å"
                    elif hasattr(interaction, "_distance"):
                        distance = f"{interaction._distance:.2f} Å"
                    else:
                        distance = "N/A"

                    # Get properties
                    properties = ""
                    if hasattr(interaction, "donor_acceptor_properties"):
                        properties = interaction.donor_acceptor_properties

                    # Store the interaction object in the dict (not serialized)
                    interaction_store[interaction_index] = interaction

                    rows.append(
                        {
                            "id": interaction_index,
                            "type": type_label,
                            "donor_res": donor_res,
                            "donor_atom": interaction.get_donor().name
                            if hasattr(interaction.get_donor(), "name")
                            else "N/A",
                            "acceptor_res": acceptor_res,
                            "acceptor_atom": interaction.get_acceptor().name
                            if hasattr(interaction.get_acceptor(), "name")
                            else "N/A",
                            "distance": distance,
                            "properties": properties,
                            "int_type": int_type,
                        }
                    )
                    interaction_index += 1

                return rows

            # Build table data for water bridges
            def build_water_bridges_data(selected_value):
                """Build table data for water bridges of the selected ligand."""
                rows = []
                interaction_index = 0

                for interaction in water_bridge_interactions:
                    donor_res = interaction.get_donor_residue()
                    acceptor_res = interaction.get_acceptor_residue()

                    # Show only interactions involving the selected ligand
                    if donor_res != selected_value and acceptor_res != selected_value:
                        continue

                    # Water bridges have unique structure (use semicolon separator like main Water Bridges tab)
                    water_residues_str = "; ".join(interaction.water_residues)

                    # Store the interaction object in water bridge store (separate from regular interactions)
                    water_bridge_store[interaction_index] = interaction

                    rows.append(
                        {
                            "id": interaction_index,
                            "donor_res": donor_res,
                            "acceptor_res": acceptor_res,
                            "water_residues": water_residues_str,
                            "bridge_length": interaction.bridge_length,
                            "distance": f"{interaction.get_donor_acceptor_distance():.2f} Å",
                            "int_type": interaction.get_interaction_type(),
                        }
                    )
                    interaction_index += 1

                return rows

            # Table columns for regular interactions
            columns = [
                {
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
                {"name": "type", "label": "Type", "field": "type", "align": "left"},
                {
                    "name": "donor_res",
                    "label": "Donor",
                    "field": "donor_res",
                    "align": "left",
                },
                {
                    "name": "donor_atom",
                    "label": "Donor Atom",
                    "field": "donor_atom",
                    "align": "left",
                },
                {
                    "name": "acceptor_res",
                    "label": "Acceptor",
                    "field": "acceptor_res",
                    "align": "left",
                },
                {
                    "name": "acceptor_atom",
                    "label": "Acceptor Atom",
                    "field": "acceptor_atom",
                    "align": "left",
                },
                {
                    "name": "distance",
                    "label": "Distance",
                    "field": "distance",
                    "align": "left",
                },
                {
                    "name": "properties",
                    "label": "Properties",
                    "field": "properties",
                    "align": "left",
                },
            ]

            # Table columns for water bridges (consistent with main Water Bridges tab)
            water_bridge_columns = [
                {
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
                {
                    "name": "donor_res",
                    "label": "Start Residue",
                    "field": "donor_res",
                    "align": "left",
                },
                {
                    "name": "acceptor_res",
                    "label": "End Residue",
                    "field": "acceptor_res",
                    "align": "left",
                },
                {
                    "name": "bridge_length",
                    "label": "Hops",
                    "field": "bridge_length",
                    "align": "center",
                },
                {
                    "name": "water_residues",
                    "label": "Water Residues",
                    "field": "water_residues",
                    "align": "left",
                },
                {
                    "name": "distance",
                    "label": "Distance (Å)",
                    "field": "distance",
                    "align": "right",
                },
            ]
            columns = [
                {
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
                {"name": "type", "label": "Type", "field": "type", "align": "left"},
                {
                    "name": "donor_res",
                    "label": "Donor",
                    "field": "donor_res",
                    "align": "left",
                },
                {
                    "name": "donor_atom",
                    "label": "Donor Atom",
                    "field": "donor_atom",
                    "align": "left",
                },
                {
                    "name": "acceptor_res",
                    "label": "Acceptor",
                    "field": "acceptor_res",
                    "align": "left",
                },
                {
                    "name": "acceptor_atom",
                    "label": "Acceptor Atom",
                    "field": "acceptor_atom",
                    "align": "left",
                },
                {
                    "name": "distance",
                    "label": "Distance",
                    "field": "distance",
                    "align": "left",
                },
                {
                    "name": "properties",
                    "label": "Properties",
                    "field": "properties",
                    "align": "left",
                },
            ]

            # Create dropdown selector (above table)
            selected_ligand = (
                ui.select(ligand_options, label="Select Ligand", value=first_ligand)
                .classes("w-full mb-4")
                .props("outlined dense")
            )

            # Filter input for searching (above table)
            filter_input = (
                ui.input(placeholder="Filter interactions (residue, atom, type, etc.)")
                .props("clearable outlined dense")
                .classes("w-full mb-4")
            )

            # Ligplot visualization (2D structure with highlighted atoms)
            ligplot_container = ui.html().classes("w-full mb-4")

            def update_ligplot(ligand_res):
                """Update ligplot when ligand selection changes."""
                try:
                    ligand_name = ligand_res.split(":")[1]
                    gen = LigplotGenerator(
                        ligand_name, self.analyzer, residue_id=ligand_res
                    )
                    html = gen.generate_svg(width=400, height=300)
                    ligplot_container.content = html
                except Exception as e:
                    ligplot_container.content = (
                        f"<p>Could not generate ligplot: {str(e)}</p>"
                    )

            # Initial ligplot
            if first_ligand:
                update_ligplot(first_ligand)

            # Update ligplot when selection changes
            def on_ligand_change(selected_value):
                update_ligplot(selected_value.value)

            selected_ligand.on_value_change(on_ligand_change)

            # Action buttons for regular interactions
            def show_all_interactions_viewer(selected_ligand_res):
                """Show 3D viewer dialog with all regular interactions for selected ligand."""
                # Get structured data for this ligand
                lig_info = ligand_info.get(selected_ligand_res)
                if not lig_info:
                    ui.notify("Error: Ligand information not found", type="negative")
                    return

                # Get all regular interactions for this ligand (exclude water bridges)
                selected_interactions = []
                for interaction in regular_interactions:
                    donor_res = interaction.get_donor_residue()
                    acceptor_res = interaction.get_acceptor_residue()
                    if (
                        donor_res == selected_ligand_res
                        or acceptor_res == selected_ligand_res
                    ):
                        selected_interactions.append(interaction)

                if not selected_interactions:
                    ui.notify("No interactions found for this ligand", type="warning")
                    return

                # Extract interaction data (donor/acceptor atom pairs, with π-center for π-interactions)
                interactions_data = []
                for interaction in selected_interactions:
                    try:
                        donor_atom = interaction.get_donor()
                        acceptor_atom = interaction.get_acceptor()
                        donor_res = interaction.get_donor_residue()
                        acceptor_res = interaction.get_acceptor_residue()
                        int_type = interaction.get_interaction_type().lower()

                        interaction_dict = {
                            "donor_atom": donor_atom.name
                            if hasattr(donor_atom, "name")
                            else "",
                            "donor_res": donor_res,
                            "acceptor_atom": acceptor_atom.name
                            if hasattr(acceptor_atom, "name")
                            else "",
                            "acceptor_res": acceptor_res,
                            "type": int_type,
                        }

                        # For π-interactions, include the ring center coordinates
                        if ("π" in int_type or "pi" in int_type) and hasattr(
                            interaction, "pi_center"
                        ):
                            interaction_dict["pi_center"] = {
                                "x": interaction.pi_center.x,
                                "y": interaction.pi_center.y,
                                "z": interaction.pi_center.z,
                            }

                        interactions_data.append(interaction_dict)
                    except Exception as e:
                        # Skip interactions that can't be extracted
                        continue

                # Generate minimal PDB with all interactions
                minimal_pdb = format_minimal_pdb(
                    self.analyzer.parser, selected_interactions
                )
                viewer_id = f"lig_all_viewer_{random.randint(1000, 9999)}"

                with (
                    ui.dialog().props("persistent") as dialog,
                    ui.card().style("width: 900px; max-width: 90vw;"),
                ):
                    with ui.card_actions().classes("justify-end"):
                        ui.button(
                            "Export PNG",
                            icon="download",
                            on_click=lambda: ui.run_javascript(
                                generate_png_export_js(
                                    viewer_id,
                                    f"{self._get_pdb_basename()}_{selected_ligand_res.replace(':', '_')}_all.png",
                                )
                            ),
                        ).props("outline color=secondary")
                        ui.button(
                            "Download PDB",
                            icon="download",
                            on_click=lambda: ui.download(
                                minimal_pdb.encode(),
                                filename=f"{self._get_pdb_basename()}_{selected_ligand_res.replace(':', '_')}_all.pdb",
                            ),
                        ).props("outline color=secondary")
                        ui.button(
                            "View in PyMOL",
                            icon="movie",
                            on_click=lambda: self._download_pymol_visualization(
                                selected_ligand_res.replace(":", "_") + "_all",
                                hydrogen_bonds=[
                                    i
                                    for i in selected_interactions
                                    if "h-bond" in i.get_interaction_type().lower()
                                ],
                                halogen_bonds=[
                                    i
                                    for i in selected_interactions
                                    if "x-bond" in i.get_interaction_type().lower()
                                ],
                                pi_interactions=[
                                    i
                                    for i in selected_interactions
                                    if (
                                        "π" in i.get_interaction_type().lower()
                                        or "pi" in i.get_interaction_type().lower()
                                    )
                                    and "pi-pi" not in i.get_interaction_type().lower()
                                    and "n-pi" not in i.get_interaction_type().lower()
                                ],
                                pi_pi_stacking=[
                                    i
                                    for i in selected_interactions
                                    if "pi-pi" in i.get_interaction_type().lower()
                                ],
                                carbonyl_interactions=[
                                    i
                                    for i in selected_interactions
                                    if "carbonyl" in i.get_interaction_type().lower()
                                ],
                                n_pi_interactions=[
                                    i
                                    for i in selected_interactions
                                    if "n-pi" in i.get_interaction_type().lower()
                                ],
                                water_bridges=[
                                    i
                                    for i in selected_interactions
                                    if "water_bridge"
                                    in i.get_interaction_type().lower()
                                ],
                            ),
                        ).props("outline color=secondary")
                        ui.button("Close", on_click=dialog.close).props("color=primary")

                    with ui.card_section().classes("q-pa-md"):
                        ui.label(f"All Interactions - {selected_ligand_res}").classes(
                            "text-subtitle1 q-mb-md"
                        )

                        # Create div for viewer
                        ui.html(
                            f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>',
                            sanitize=False,
                        )

                dialog.open()

                # Initialize viewer after dialog is opened and rendered
                # Pass interaction data, PDB content, viewer ID, and ligand residue name (following standard convention)
                ui.timer(
                    0.1,
                    lambda: ui.run_javascript(
                        generate_ligand_interactions_viewer_js(
                            interactions_data, minimal_pdb, viewer_id, lig_info["name"]
                        )
                    ),
                    once=True,
                )

            def download_all_pymol(selected_ligand_res):
                """Download all interactions as PyMOL visualization package."""
                # Get all interactions for this ligand
                selected_interactions = []
                for interaction in ligand_interactions:
                    donor_res = interaction.get_donor_residue()
                    acceptor_res = interaction.get_acceptor_residue()
                    if (
                        donor_res == selected_ligand_res
                        or acceptor_res == selected_ligand_res
                    ):
                        selected_interactions.append(interaction)

                if not selected_interactions:
                    ui.notify("No interactions found for this ligand", type="warning")
                    return

                # Organize interactions by type
                self._download_pymol_visualization(
                    selected_ligand_res.replace(":", "_") + "_all",
                    hydrogen_bonds=[
                        i
                        for i in selected_interactions
                        if "h-bond" in i.get_interaction_type().lower()
                    ],
                    halogen_bonds=[
                        i
                        for i in selected_interactions
                        if "x-bond" in i.get_interaction_type().lower()
                    ],
                    pi_interactions=[
                        i
                        for i in selected_interactions
                        if (
                            "π" in i.get_interaction_type().lower()
                            or "pi" in i.get_interaction_type().lower()
                        )
                        and "pi-pi" not in i.get_interaction_type().lower()
                        and "n-pi" not in i.get_interaction_type().lower()
                    ],
                    pi_pi_stacking=[
                        i
                        for i in selected_interactions
                        if "pi-pi" in i.get_interaction_type().lower()
                    ],
                    carbonyl_interactions=[
                        i
                        for i in selected_interactions
                        if "carbonyl" in i.get_interaction_type().lower()
                    ],
                    n_pi_interactions=[
                        i
                        for i in selected_interactions
                        if "n-pi" in i.get_interaction_type().lower()
                    ],
                    water_bridges=[
                        i
                        for i in selected_interactions
                        if "water_bridge" in i.get_interaction_type().lower()
                    ],
                )

            def download_csv(selected_ligand_res):
                """Download interactions as ZIP containing CSV files for selected ligand."""
                import io
                import zipfile
                from ...export.results import (
                    write_ligand_interactions_csv,
                    write_ligand_water_bridges_csv,
                )

                # Create zip file in memory
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                    # Add regular interactions CSV
                    csv_content = write_ligand_interactions_csv(
                        self.analyzer, filename=None, ligand_residue=selected_ligand_res
                    )
                    if (
                        csv_content and csv_content.strip()
                    ):  # Only add if there are interactions
                        zip_file.writestr(
                            f"{self._get_pdb_basename()}_{selected_ligand_res.replace(':', '_')}_interactions.csv",
                            csv_content,
                        )

                    # Add water bridges CSV
                    wb_content = write_ligand_water_bridges_csv(
                        self.analyzer, filename=None, ligand_residue=selected_ligand_res
                    )
                    if (
                        wb_content and wb_content.strip()
                    ):  # Only add if there are water bridges
                        zip_file.writestr(
                            f"{self._get_pdb_basename()}_{selected_ligand_res.replace(':', '_')}_water_bridges.csv",
                            wb_content,
                        )

                ui.download(
                    zip_buffer.getvalue(),
                    filename=f"{self._get_pdb_basename()}_{selected_ligand_res.replace(':', '_')}_interactions.zip",
                )

            # Action buttons container
            with ui.row().classes("w-full gap-2 q-mb-4"):
                ui.button(
                    "3D View All Interactions",
                    icon="visibility",
                    on_click=lambda: show_all_interactions_viewer(
                        selected_ligand.value
                    ),
                ).props("color=primary")
                ui.button(
                    "View All Interactions in PyMOL",
                    icon="movie",
                    on_click=lambda: download_all_pymol(selected_ligand.value),
                ).props("color=secondary")
                ui.button(
                    "Download Interactions as CSV",
                    icon="download",
                    on_click=lambda: download_csv(selected_ligand.value),
                ).props("color=info")

            # Handle ligand selection changes
            def on_ligand_selected(e):
                if e.value:
                    # Update tables with filtered interactions
                    regular_rows = build_regular_interactions_data(e.value)
                    table.update_rows(regular_rows, clear_selection=True)

                    if water_bridges_enabled:
                        water_bridge_rows = build_water_bridges_data(e.value)
                        water_bridge_table.update_rows(
                            water_bridge_rows, clear_selection=True
                        )

            selected_ligand.on_value_change(on_ligand_selected)

            # Create table for regular interactions
            table = ui.table(
                columns=columns,
                rows=build_regular_interactions_data(first_ligand)
                if first_ligand
                else [],
                row_key="id",
            )
            table.classes("w-full")

            # Flag to track if water bridges section is shown
            water_bridges_enabled = bool(water_bridge_interactions)

            # Bind filter input to table
            filter_input.bind_value(table, "filter")

            # Add visualize button slot with color based on interaction type
            table.add_slot(
                "body-cell-visualize",
                """
                <q-td :props="props">
                    <q-btn
                        flat dense round
                        icon="visibility"
                        :color="
                            props.row.int_type.toLowerCase().includes('hydrogen') ? 'primary' :
                            props.row.int_type.toLowerCase().includes('halogen') ? 'orange' :
                            props.row.int_type.toLowerCase().includes('pi-pi') || props.row.int_type.toLowerCase().includes('stacking') ? 'purple' :
                            props.row.int_type.toLowerCase().includes('pi') ? 'green' :
                            props.row.int_type.toLowerCase().includes('carbonyl') ? 'red' :
                            props.row.int_type.toLowerCase().includes('n-pi') ? 'teal' :
                            props.row.int_type.toLowerCase().includes('water') ? 'cyan' :
                            'blue'
                        "
                        size="sm"
                        @click="$parent.$emit('visualize', props.row)"
                    />
                </q-td>
            """,
            )

            # Create water bridges section if there are any
            if water_bridges_enabled:
                ui.label("Water Bridges").classes("text-h6 q-mt-lg")
                water_bridge_table = ui.table(
                    columns=water_bridge_columns,
                    rows=build_water_bridges_data(first_ligand) if first_ligand else [],
                    row_key="id",
                )
                water_bridge_table.classes("w-full")

                # Add visualize button slot for water bridges
                water_bridge_table.add_slot(
                    "body-cell-visualize",
                    """
                    <q-td :props="props">
                        <q-btn
                            flat dense round
                            icon="visibility"
                            color="cyan"
                            size="sm"
                            @click="$parent.$emit('visualize', props.row)"
                        />
                    </q-td>
                """,
                )

                # Handle water bridge visualization
                def show_water_bridge_visualization(e):
                    """Show 3D visualization for water bridge from ligand interactions."""
                    row_data = e.args if hasattr(e, "args") else {}
                    row_id = row_data.get("id")
                    if row_id is None:
                        ui.notify("Error: Water bridge ID not found", type="negative")
                        return

                    interaction = water_bridge_store.get(row_id)
                    if not interaction:
                        ui.notify("Error: Water bridge not found", type="negative")
                        return

                    # Use extract_water_bridge_pdb for water bridges to include water molecules
                    minimal_pdb = extract_water_bridge_pdb(
                        self.analyzer.parser, interaction
                    )
                    viewer_id = f"lig_wbridge_viewer_{random.randint(1000, 9999)}"

                    def export_png():
                        pdb_base = self._get_pdb_basename()
                        donor = row_data.get("donor_res", "donor").replace(":", "_")
                        acceptor = row_data.get("acceptor_res", "acceptor").replace(
                            ":", "_"
                        )
                        filename = f"{pdb_base}_{donor}_to_{acceptor}_wbridge.png"
                        ui.run_javascript(generate_png_export_js(viewer_id, filename))

                    def download_pdb():
                        pdb_base = self._get_pdb_basename()
                        donor = row_data.get("donor_res", "donor").replace(":", "_")
                        acceptor = row_data.get("acceptor_res", "acceptor").replace(
                            ":", "_"
                        )
                        filename = f"{pdb_base}_{donor}_to_{acceptor}_wbridge.pdb"
                        ui.download(minimal_pdb.encode(), filename=filename)

                    def download_pymol():
                        donor = row_data.get("donor_res", "donor").replace(":", "_")
                        acceptor = row_data.get("acceptor_res", "acceptor").replace(
                            ":", "_"
                        )
                        label = f"{donor}_to_{acceptor}_wbridge"
                        self._download_pymol_visualization(
                            label, water_bridges=[interaction]
                        )

                    with (
                        ui.dialog().props("persistent") as dialog,
                        ui.card().style("width: 900px; max-width: 90vw;"),
                    ):
                        with ui.card_actions().classes("justify-end"):
                            ui.button(
                                "Export PNG", icon="download", on_click=export_png
                            ).props("outline color=secondary")
                            ui.button(
                                "Download PDB", icon="download", on_click=download_pdb
                            ).props("outline color=secondary")
                            ui.button(
                                "View in PyMOL", icon="movie", on_click=download_pymol
                            ).props("outline color=secondary")
                            ui.button("Close", on_click=dialog.close).props(
                                "color=primary"
                            )

                        with ui.card_section().classes("q-pa-md"):
                            viewer_id_div = f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>'
                            ui.html(viewer_id_div, sanitize=False)

                    dialog.open()
                    ui.timer(
                        0.1,
                        lambda: ui.run_javascript(
                            generate_water_bridge_viewer_js(
                                interaction, minimal_pdb, viewer_id
                            )
                        ),
                        once=True,
                    )

                water_bridge_table.on("visualize", show_water_bridge_visualization)

            # Handle visualization requests
            def show_interaction_visualization(row_data):
                # row_data can be either a dict (from event args) or an ID
                row_id = row_data.get("id") if isinstance(row_data, dict) else row_data

                # Retrieve interaction from store using row ID
                interaction = interaction_store.get(row_id)
                if not interaction:
                    ui.notify("Error: Interaction not found", type="negative")
                    return
                int_type = (
                    row_data.get("int_type", "").lower()
                    if isinstance(row_data, dict)
                    else interaction.get_interaction_type().lower()
                )

                # Determine which viewer to use based on actual type strings
                viewer_gen = None
                # Check for hydrogen bonds
                if (
                    "hydrogen" in int_type
                    or "h-bond" in int_type
                    or "h bond" in int_type
                ):
                    viewer_id = f"lig_hb_viewer_{random.randint(1000, 9999)}"
                    viewer_gen = generate_hydrogen_bond_viewer_js
                # Check for halogen bonds
                elif (
                    "halogen" in int_type
                    or "x-bond" in int_type
                    or "x bond" in int_type
                ):
                    viewer_id = f"lig_xb_viewer_{random.randint(1000, 9999)}"
                    viewer_gen = generate_halogen_bond_viewer_js
                # Check for π-π stacking (with actual π symbol or spelled out)
                elif (
                    ("π" in int_type and "π" in int_type and "stack" in int_type)
                    or "pi-pi" in int_type
                    or "pi pi" in int_type
                    or "stacking" in int_type
                ):
                    viewer_id = f"lig_pipi_viewer_{random.randint(1000, 9999)}"
                    viewer_gen = generate_pi_pi_stacking_viewer_js
                # Check for π-interactions (with actual π symbol or spelled out)
                elif (
                    ("π" in int_type and "inter" in int_type)
                    or "pi inter" in int_type
                    or "pi-inter" in int_type
                ):
                    viewer_id = f"lig_pi_viewer_{random.randint(1000, 9999)}"
                    viewer_gen = generate_pi_interaction_viewer_js
                # Check for carbonyl interactions
                elif "carbonyl" in int_type or "n→π" in int_type or "npi" in int_type:
                    viewer_id = f"lig_carbonyl_viewer_{random.randint(1000, 9999)}"
                    viewer_gen = generate_carbonyl_interaction_viewer_js
                # Check for n-π* interactions
                elif (
                    "n-pi" in int_type
                    or "n pi" in int_type
                    or ("n" in int_type and "π" in int_type)
                ):
                    viewer_id = f"lig_npi_viewer_{random.randint(1000, 9999)}"
                    viewer_gen = generate_n_pi_interaction_viewer_js
                # Check for water bridges
                elif "water" in int_type or "bridge" in int_type:
                    viewer_id = f"lig_wbridge_viewer_{random.randint(1000, 9999)}"
                    viewer_gen = generate_water_bridge_viewer_js

                if not viewer_gen:
                    ui.notify(
                        f"Error: Unknown interaction type '{int_type}'", type="negative"
                    )
                    return

                # Generate minimal PDB for this interaction
                # Water bridges need special handling to include water molecules
                if hasattr(interaction, "water_residues"):
                    minimal_pdb = extract_water_bridge_pdb(
                        self.analyzer.parser, interaction
                    )
                else:
                    minimal_pdb = format_minimal_pdb(
                        self.analyzer.parser, [interaction]
                    )

                # Create dialog with proper structure matching other interaction visualizations
                def export_png():
                    """Export viewer as PNG."""
                    pdb_base = self._get_pdb_basename()
                    donor = row_data.get("donor_res", "donor").replace(":", "_")
                    acceptor = row_data.get("acceptor_res", "acceptor").replace(
                        ":", "_"
                    )
                    filename = f"{pdb_base}_{donor}_to_{acceptor}_ligand.png"
                    ui.run_javascript(generate_png_export_js(viewer_id, filename))

                def download_pdb():
                    """Download extracted PDB file."""
                    pdb_base = self._get_pdb_basename()
                    donor = row_data.get("donor_res", "donor").replace(":", "_")
                    acceptor = row_data.get("acceptor_res", "acceptor").replace(
                        ":", "_"
                    )
                    filename = f"{pdb_base}_{donor}_to_{acceptor}_ligand.pdb"
                    ui.download(minimal_pdb.encode(), filename=filename)

                def download_pymol():
                    """Download as PyMOL visualization package."""
                    donor = row_data.get("donor_res", "donor").replace(":", "_")
                    acceptor = row_data.get("acceptor_res", "acceptor").replace(
                        ":", "_"
                    )
                    label = f"{donor}_to_{acceptor}"
                    if "hydrogen" in int_type:
                        self._download_pymol_visualization(
                            label, hydrogen_bonds=[interaction]
                        )
                    elif "halogen" in int_type:
                        self._download_pymol_visualization(
                            label, halogen_bonds=[interaction]
                        )
                    elif "pi-pi" in int_type:
                        self._download_pymol_visualization(
                            label, pi_pi_stacking=[interaction]
                        )
                    elif "pi" in int_type:
                        self._download_pymol_visualization(
                            label, pi_interactions=[interaction]
                        )
                    elif "carbonyl" in int_type:
                        self._download_pymol_visualization(
                            label, carbonyl_interactions=[interaction]
                        )
                    elif "n-pi" in int_type:
                        self._download_pymol_visualization(
                            label, n_pi_interactions=[interaction]
                        )
                    elif "water" in int_type:
                        self._download_pymol_visualization(
                            label, water_bridges=[interaction]
                        )

                with (
                    ui.dialog().props("persistent") as dialog,
                    ui.card().style("width: 900px; max-width: 90vw;"),
                ):
                    with ui.card_actions().classes("justify-end"):
                        ui.button(
                            "Export PNG", icon="download", on_click=export_png
                        ).props("outline color=secondary")
                        ui.button(
                            "Download PDB", icon="download", on_click=download_pdb
                        ).props("outline color=secondary")
                        ui.button(
                            "View in PyMOL", icon="movie", on_click=download_pymol
                        ).props("outline color=secondary")
                        ui.button("Close", on_click=dialog.close).props("color=primary")

                    with ui.card_section().classes("q-pa-md"):
                        ui.label(
                            f"{row_data.get('donor_res', '')} → {row_data.get('acceptor_res', '')}"
                        ).classes("text-subtitle1 q-mb-md")

                        # Create div for viewer - container must have actual width
                        ui.html(
                            f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>',
                            sanitize=False,
                        )

                dialog.open()

                # Initialize viewer after dialog is opened and rendered
                ui.timer(
                    0.1,
                    lambda: ui.run_javascript(
                        viewer_gen(interaction, minimal_pdb, viewer_id)
                    ),
                    once=True,
                )

            # Handle table row visualization
            def handle_visualize(e):
                # e.args should contain the row data
                row_data = e.args if hasattr(e, "args") else {}
                if isinstance(row_data, dict) and "id" in row_data:
                    show_interaction_visualization(row_data)
                else:
                    ui.notify("Error: Could not retrieve row data", type="negative")

            table.on("visualize", handle_visualize)

    def reset(self):
        """Reset the results panel and clear all analysis results."""
        self.analyzer = None
        self.current_file = None
        self.container.clear()

        with self.container:
            ui.label("No analysis results yet").classes("text-h6 text-grey q-pa-lg")
