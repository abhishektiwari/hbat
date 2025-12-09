"""
Results display component with 3D visualization for HBAT web interface.

This module provides results display with integrated py3Dmol visualization
for molecular interactions.
"""

import math
from typing import Optional

from nicegui import ui

from ...core.analysis import NPMolecularInteractionAnalyzer
from ...visualization.chain_graph import render_chain_for_web
from ...visualization.pymol3d import (
    generate_carbonyl_interaction_viewer_js,
    generate_halogen_bond_viewer_js,
    generate_hydrogen_bond_viewer_js,
    generate_n_pi_interaction_viewer_js,
    generate_pi_interaction_viewer_js,
    generate_pi_pi_stacking_viewer_js,
    generate_png_export_js,
)


class WebResultsPanel:
    """Results panel component with 3D visualization."""

    def __init__(self, container):
        """Initialize the results panel.

        :param container: Container element for results display
        """
        self.container = container
        self.analyzer: Optional[NPMolecularInteractionAnalyzer] = None
        self.pdb_content: Optional[str] = None
        self.current_file: Optional[str] = None

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
        self.cooperativity_panel = None

    async def update_results(
        self, analyzer: NPMolecularInteractionAnalyzer, pdb_content: str, filename: str
    ):
        """Update results display with analysis results.

        :param analyzer: The analyzer with results
        :param pdb_content: PDB file content for visualization
        :param filename: Name of the analyzed file
        """
        self.analyzer = analyzer
        self.pdb_content = pdb_content
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
        """Get configuration for all interaction types.

        :returns: List of interaction configurations
        :rtype: list[dict]
        """
        return [
            {
                "id": "summary",
                "label": "Summary",
                "icon": "summarize",
                "attr": None,
                "always_show": True,
                "panel_attr": "summary_panel",
                "update_method": "_update_summary_panel",
            },
            {
                "id": "hydrogen",
                "label": "Hydrogen Bonds",
                "icon": "link",
                "attr": "hydrogen_bonds",
                "panel_attr": "hydrogen_panel",
                "update_method": "_update_hydrogen_bonds_panel",
            },
            {
                "id": "halogen",
                "label": "Halogen Bonds",
                "icon": "whatshot",
                "attr": "halogen_bonds",
                "panel_attr": "halogen_panel",
                "update_method": "_update_halogen_bonds_panel",
            },
            {
                "id": "pi",
                "label": "π Interactions",
                "icon": "fiber_manual_record",
                "attr": "pi_interactions",
                "panel_attr": "pi_panel",
                "update_method": "_update_pi_interactions_panel",
            },
            {
                "id": "pi_pi",
                "label": "π-π Stacking",
                "icon": "layers",
                "attr": "pi_pi_interactions",
                "panel_attr": "pi_pi_panel",
                "update_method": "_update_pi_pi_stacking_panel",
            },
            {
                "id": "carbonyl",
                "label": "Carbonyl n→π*",
                "icon": "architecture",
                "attr": "carbonyl_interactions",
                "panel_attr": "carbonyl_panel",
                "update_method": "_update_carbonyl_interactions_panel",
            },
            {
                "id": "n_pi",
                "label": "n→π*",
                "icon": "arrow_forward",
                "attr": "n_pi_interactions",
                "panel_attr": "n_pi_panel",
                "update_method": "_update_n_pi_interactions_panel",
            },
            {
                "id": "cooperativity",
                "label": "Cooperativity",
                "icon": "account_tree",
                "attr": "cooperativity_chains",
                "panel_attr": "cooperativity_panel",
                "update_method": "_update_cooperativity_chains_panel",
            },
        ]

    def _should_show_interaction(self, config):
        """Check if interaction should be shown.

        :param config: Interaction configuration
        :returns: True if interaction should be shown
        :rtype: bool
        """
        if config.get("always_show"):
            return True

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

                    # Call the update method
                    update_method = getattr(self, config["update_method"])
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
                    ui.label("π Interactions").classes("text-h6")
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
                    ui.label("n→π* Interactions").classes("text-h6")
                    ui.label(
                        str(summary.get("n_pi_interactions", {}).get("count", 0))
                    ).classes("text-h4 text-green")

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

    def _update_hydrogen_bonds_panel(self):
        """Update hydrogen bonds panel with table and visualization."""
        with self.hydrogen_panel:
            ui.label(f"Hydrogen Bonds ({len(self.analyzer.hydrogen_bonds)})").classes(
                "text-h5"
            )

            # Create table
            columns = [
                {
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
                {
                    "name": "donor_res",
                    "label": "Donor Residue",
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
                    "name": "hydrogen",
                    "label": "Hydrogen Atom",
                    "field": "hydrogen",
                    "align": "left",
                },
                {
                    "name": "acceptor_res",
                    "label": "Acceptor Residue",
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
                    "label": "H...A (Å)",
                    "field": "distance",
                    "align": "right",
                },
                {
                    "name": "angle",
                    "label": "Angle (°)",
                    "field": "angle",
                    "align": "right",
                },
                {
                    "name": "da_distance",
                    "label": "D...A (Å)",
                    "field": "da_distance",
                    "align": "right",
                },
                {"name": "type", "label": "Type", "field": "type", "align": "left"},
                {
                    "name": "da_props",
                    "label": "D-A Props",
                    "field": "da_props",
                    "align": "left",
                },
                {
                    "name": "bs_int",
                    "label": "B/S",
                    "field": "bs_int",
                    "align": "center",
                },
            ]

            rows = []
            for idx, hb in enumerate(self.analyzer.hydrogen_bonds):
                bs_int = (
                    "B"
                    if (hb.get_donor_residue() != hb.get_acceptor_residue())
                    else "S"
                )
                rows.append(
                    {
                        "id": idx,
                        "donor_res": hb.get_donor_residue(),
                        "donor_atom": hb.donor.name,
                        "hydrogen": hb.hydrogen.name,
                        "acceptor_res": hb.get_acceptor_residue(),
                        "acceptor_atom": hb.acceptor.name,
                        "distance": f"{hb.distance:.2f}",
                        "angle": f"{math.degrees(hb.angle):.1f}",
                        "da_distance": f"{hb.donor_acceptor_distance:.2f}",
                        "type": hb.bond_type,
                        "da_props": hb.donor_acceptor_properties,
                        "bs_int": bs_int,
                    }
                )

            # Add filter input
            filter_input = (
                ui.input(placeholder="Filter interactions (residue, atom, type, etc.)")
                .props("clearable outlined dense")
                .classes("w-full mb-2")
            )

            table = ui.table(columns=columns, rows=rows, row_key="id").classes("w-full")
            filter_input.bind_value(table, "filter")

            table.add_slot(
                "body-cell-visualize",
                """
                <q-td :props="props">
                    <q-btn size="sm" color="primary" round dense icon="visibility"
                           @click="$parent.$emit('visualize', props.row)" />
                </q-td>
            """,
            )

            # Handle visualization button clicks
            def show_hydrogen_bond(e):
                idx = e.args["id"]
                hb = self.analyzer.hydrogen_bonds[idx]
                self._show_hydrogen_bond_visualization(hb)

            table.on("visualize", show_hydrogen_bond)

    def _show_hydrogen_bond_visualization(self, hb):
        """Show 3D visualization of a hydrogen bond in a dialog.

        :param hb: Hydrogen bond interaction
        """
        # Create unique viewer ID
        import random

        viewer_id = f"hb_viewer_{random.randint(1000, 9999)}"

        def export_png():
            """Export viewer as PNG."""
            filename = f"{hb.get_donor_residue()}_to_{hb.get_acceptor_residue()}_hbond.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
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
            0.1, lambda: self._initialize_hydrogen_bond_viewer(viewer_id, hb), once=True
        )

    def _initialize_hydrogen_bond_viewer(self, viewer_id: str, hb):
        """Initialize py3Dmol viewer for hydrogen bond visualization.

        :param viewer_id: Unique ID for the viewer div
        :param hb: Hydrogen bond interaction
        """
        javascript = generate_hydrogen_bond_viewer_js(hb, self.pdb_content, viewer_id)
        ui.run_javascript(javascript)

    def _update_halogen_bonds_panel(self):
        """Update halogen bonds panel."""
        with self.halogen_panel:
            ui.label(f"Halogen Bonds ({len(self.analyzer.halogen_bonds)})").classes(
                "text-h5"
            )

            # Similar to hydrogen bonds but for halogen bonds
            columns = [
                {
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
                {
                    "name": "halogen_res",
                    "label": "Halogen Residue",
                    "field": "halogen_res",
                    "align": "left",
                },
                {
                    "name": "donor_atom",
                    "label": "Donor Atom",
                    "field": "donor_atom",
                    "align": "left",
                },
                {
                    "name": "halogen_atom",
                    "label": "Halogen Atom",
                    "field": "halogen_atom",
                    "align": "left",
                },
                {
                    "name": "acceptor_res",
                    "label": "Acceptor Residue",
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
                    "label": "X...A (Å)",
                    "field": "distance",
                    "align": "right",
                },
                {
                    "name": "angle",
                    "label": "Angle (°)",
                    "field": "angle",
                    "align": "right",
                },
                {"name": "type", "label": "Type", "field": "type", "align": "left"},
                {
                    "name": "bs_interaction",
                    "label": "B/S Interaction",
                    "field": "bs_interaction",
                    "align": "center",
                },
                {
                    "name": "da_properties",
                    "label": "D-A Properties",
                    "field": "da_properties",
                    "align": "left",
                },
            ]

            rows = []
            for idx, xb in enumerate(self.analyzer.halogen_bonds):
                bs_int = (
                    "B"
                    if (xb.get_donor_residue() != xb.get_acceptor_residue())
                    else "S"
                )
                rows.append(
                    {
                        "id": idx,
                        "halogen_res": xb.get_donor_residue(),
                        "donor_atom": (
                            xb.donor.name if hasattr(xb, "donor") and xb.donor else ""
                        ),
                        "halogen_atom": xb.halogen.name,
                        "acceptor_res": xb.get_acceptor_residue(),
                        "acceptor_atom": xb.acceptor.name,
                        "distance": f"{xb.distance:.2f}",
                        "angle": f"{math.degrees(xb.angle):.1f}",
                        "type": xb.bond_type,
                        "bs_interaction": bs_int,
                        "da_properties": xb.donor_acceptor_properties,
                    }
                )

            # Add filter input
            filter_input = (
                ui.input(placeholder="Filter interactions (residue, atom, type, etc.)")
                .props("clearable outlined dense")
                .classes("w-full mb-2")
            )

            table = ui.table(columns=columns, rows=rows, row_key="id").classes("w-full")
            filter_input.bind_value(table, "filter")

            table.add_slot(
                "body-cell-visualize",
                """
                <q-td :props="props">
                    <q-btn size="sm" color="orange" round dense icon="visibility"
                           @click="$parent.$emit('visualize', props.row)" />
                </q-td>
            """,
            )

            def show_halogen_bond(e):
                idx = e.args["id"]
                xb = self.analyzer.halogen_bonds[idx]
                self._show_halogen_bond_visualization(xb)

            table.on("visualize", show_halogen_bond)

    def _show_halogen_bond_visualization(self, xb):
        """Show 3D visualization of a halogen bond in a dialog."""
        # Create unique viewer ID
        import random

        viewer_id = f"xb_viewer_{random.randint(1000, 9999)}"

        def export_png():
            """Export viewer as PNG."""
            filename = f"{xb.get_donor_residue()}_to_{xb.get_acceptor_residue()}_xbond.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
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
            0.1, lambda: self._initialize_halogen_bond_viewer(viewer_id, xb), once=True
        )

    def _initialize_halogen_bond_viewer(self, viewer_id: str, xb):
        """Initialize py3Dmol viewer for halogen bond visualization.

        :param viewer_id: Unique ID for the viewer div
        :param xb: Halogen bond interaction
        """
        javascript = generate_halogen_bond_viewer_js(xb, self.pdb_content, viewer_id)
        ui.run_javascript(javascript)

    def _show_pi_interaction_visualization(self, pi):
        """Show 3D visualization of a π interaction in a dialog."""
        import random

        viewer_id = f"pi_viewer_{random.randint(1000, 9999)}"

        def export_png():
            """Export viewer as PNG."""
            filename = f"{pi.get_donor_residue()}_to_{pi.get_acceptor_residue()}_pi.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
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
            lambda: self._initialize_pi_interaction_viewer(viewer_id, pi),
            once=True,
        )

    def _initialize_pi_interaction_viewer(self, viewer_id: str, pi):
        """Initialize py3Dmol viewer for π interaction visualization.

        :param viewer_id: Unique ID for the viewer div
        :param pi: Pi interaction
        """
        javascript = generate_pi_interaction_viewer_js(pi, self.pdb_content, viewer_id)
        ui.run_javascript(javascript)

    def _update_pi_interactions_panel(self):
        """Update π interactions panel."""
        with self.pi_panel:
            ui.label(f"π Interactions ({len(self.analyzer.pi_interactions)})").classes(
                "text-h5"
            )

            # Create table
            columns = [
                {
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
                {
                    "name": "donor_res",
                    "label": "Donor Residue",
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
                    "name": "pi_res",
                    "label": "π Residue",
                    "field": "pi_res",
                    "align": "left",
                },
                {
                    "name": "distance",
                    "label": "H...π (Å)",
                    "field": "distance",
                    "align": "right",
                },
                {
                    "name": "angle",
                    "label": "Angle (°)",
                    "field": "angle",
                    "align": "right",
                },
                {"name": "type", "label": "Type", "field": "type", "align": "left"},
                {
                    "name": "da_props",
                    "label": "D-A Props",
                    "field": "da_props",
                    "align": "left",
                },
                {
                    "name": "bs_int",
                    "label": "B/S",
                    "field": "bs_int",
                    "align": "center",
                },
            ]

            rows = []
            for idx, pi in enumerate(self.analyzer.pi_interactions):
                # Determine subtype from interaction atom
                x_atom = pi.hydrogen.element
                donor_atom = pi.donor.element
                if x_atom == "H":
                    if donor_atom == "C":
                        subtype = "C-H...π"
                    elif donor_atom == "N":
                        subtype = "N-H...π"
                    elif donor_atom == "O":
                        subtype = "O-H...π"
                    elif donor_atom == "S":
                        subtype = "S-H...π"
                    else:
                        subtype = "H...π"
                elif x_atom == "CL":
                    subtype = "C-Cl...π"
                elif x_atom == "BR":
                    subtype = "C-Br...π"
                elif x_atom == "I":
                    subtype = "C-I...π"
                else:
                    subtype = f"{donor_atom}-{x_atom}...π"

                bs_int = (
                    "B"
                    if (pi.get_donor_residue() != pi.get_acceptor_residue())
                    else "S"
                )

                rows.append(
                    {
                        "id": idx,
                        "donor_res": pi.get_donor_residue(),
                        "donor_atom": pi.donor.name,
                        "pi_res": pi.get_acceptor_residue(),
                        "distance": f"{pi.distance:.2f}",
                        "angle": f"{math.degrees(pi.angle):.1f}",
                        "type": subtype,
                        "da_props": pi.donor_acceptor_properties,
                        "bs_int": bs_int,
                    }
                )

            # Add filter input
            filter_input = (
                ui.input(placeholder="Filter interactions (residue, atom, type, etc.)")
                .props("clearable outlined dense")
                .classes("w-full mb-2")
            )

            # Create table with custom styling
            table = (
                ui.table(columns=columns, rows=rows, row_key="id")
                .classes("w-full")
                .props("dense")
            )
            filter_input.bind_value(table, "filter")

            # Add custom slot for visualize column to show button
            table.add_slot(
                "body-cell-visualize",
                """
                <q-td :props="props">
                    <q-btn size="sm" color="green" round dense icon="visibility"
                           @click="$parent.$emit('visualize', props.row)" />
                </q-td>
            """,
            )

            def show_pi_interaction(e):
                pi_idx = e.args["id"]
                pi = self.analyzer.pi_interactions[pi_idx]
                self._show_pi_interaction_visualization(pi)

            table.on("visualize", show_pi_interaction)

    def _update_pi_pi_stacking_panel(self):
        """Update π-π stacking panel."""
        with self.pi_pi_panel:
            ui.label(f"π-π Stacking ({len(self.analyzer.pi_pi_interactions)})").classes(
                "text-h5"
            )

            # Create table
            columns = [
                {
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
                {
                    "name": "ring1_res",
                    "label": "Ring 1 Residue",
                    "field": "ring1_res",
                    "align": "left",
                },
                {
                    "name": "ring1_atoms",
                    "label": "Ring 1 Atoms",
                    "field": "ring1_atoms",
                    "align": "left",
                },
                {
                    "name": "ring2_res",
                    "label": "Ring 2 Residue",
                    "field": "ring2_res",
                    "align": "left",
                },
                {
                    "name": "ring2_atoms",
                    "label": "Ring 2 Atoms",
                    "field": "ring2_atoms",
                    "align": "left",
                },
                {
                    "name": "distance",
                    "label": "Distance (Å)",
                    "field": "distance",
                    "align": "right",
                },
                {
                    "name": "plane_angle",
                    "label": "Plane Angle (°)",
                    "field": "plane_angle",
                    "align": "right",
                },
                {
                    "name": "offset",
                    "label": "Offset (Å)",
                    "field": "offset",
                    "align": "right",
                },
                {
                    "name": "stacking_type",
                    "label": "Stacking Type",
                    "field": "stacking_type",
                    "align": "left",
                },
                {
                    "name": "bs_int",
                    "label": "B/S",
                    "field": "bs_int",
                    "align": "center",
                },
            ]

            rows = []
            for idx, pi_pi in enumerate(self.analyzer.pi_pi_interactions):
                ring1_atoms = ",".join([atom.name for atom in pi_pi.ring1_atoms[:3]])
                if len(pi_pi.ring1_atoms) > 3:
                    ring1_atoms += "..."
                ring2_atoms = ",".join([atom.name for atom in pi_pi.ring2_atoms[:3]])
                if len(pi_pi.ring2_atoms) > 3:
                    ring2_atoms += "..."

                bs_int = "B" if pi_pi.is_between_residues else "S"

                rows.append(
                    {
                        "id": idx,
                        "ring1_res": pi_pi.ring1_residue,
                        "ring1_atoms": ring1_atoms,
                        "ring2_res": pi_pi.ring2_residue,
                        "ring2_atoms": ring2_atoms,
                        "distance": f"{pi_pi._distance:.2f}",
                        "plane_angle": f"{pi_pi.plane_angle:.1f}",
                        "offset": f"{pi_pi.offset:.2f}",
                        "stacking_type": pi_pi.stacking_type.capitalize(),
                        "bs_int": bs_int,
                    }
                )

            # Add filter input
            filter_input = (
                ui.input(placeholder="Filter interactions (residue, atom, type, etc.)")
                .props("clearable outlined dense")
                .classes("w-full mb-2")
            )

            table = (
                ui.table(columns=columns, rows=rows, row_key="id")
                .classes("w-full")
                .props("dense")
            )
            filter_input.bind_value(table, "filter")

            # Add custom slot for visualize column to show button
            table.add_slot(
                "body-cell-visualize",
                """
                <q-td :props="props">
                    <q-btn size="sm" color="purple" round dense icon="visibility"
                           @click="$parent.$emit('visualize', props.row)" />
                </q-td>
            """,
            )

            def show_pi_pi_stacking(e):
                idx = e.args["id"]
                pi_pi = self.analyzer.pi_pi_interactions[idx]
                self._show_pi_pi_stacking_visualization(pi_pi)

            table.on("visualize", show_pi_pi_stacking)

    def _update_carbonyl_interactions_panel(self):
        """Update carbonyl interactions panel."""
        with self.carbonyl_panel:
            ui.label(
                f"Carbonyl n→π* ({len(self.analyzer.carbonyl_interactions)})"
            ).classes("text-h5")

            # Create table
            columns = [
                {
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
                {
                    "name": "acceptor_res",
                    "label": "Acceptor Residue",
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
                    "name": "carbonyl_res",
                    "label": "Carbonyl Residue",
                    "field": "carbonyl_res",
                    "align": "left",
                },
                {
                    "name": "carbonyl_atoms",
                    "label": "Carbonyl C=O",
                    "field": "carbonyl_atoms",
                    "align": "left",
                },
                {
                    "name": "distance",
                    "label": "O···C Distance (Å)",
                    "field": "distance",
                    "align": "right",
                },
                {
                    "name": "angle",
                    "label": "Bürgi-Dunitz Angle (°)",
                    "field": "angle",
                    "align": "right",
                },
                {
                    "name": "carbonyl_type",
                    "label": "Carbonyl Type",
                    "field": "carbonyl_type",
                    "align": "left",
                },
                {
                    "name": "bs_int",
                    "label": "B/S",
                    "field": "bs_int",
                    "align": "center",
                },
            ]

            rows = []
            for idx, carbonyl in enumerate(self.analyzer.carbonyl_interactions):
                interaction_type = "Backbone" if carbonyl.is_backbone else "Sidechain"
                bs_int = "B" if carbonyl.is_between_residues else "S"
                carbonyl_atoms = (
                    f"{carbonyl.donor_carbon.name}={carbonyl.donor_oxygen.name}"
                )

                rows.append(
                    {
                        "id": idx,
                        "acceptor_res": carbonyl.get_acceptor_residue(),
                        "acceptor_atom": carbonyl.acceptor_carbon.name,
                        "carbonyl_res": carbonyl.get_donor_residue(),
                        "carbonyl_atoms": carbonyl_atoms,
                        "distance": f"{carbonyl.distance:.2f}",
                        "angle": f"{carbonyl.burgi_dunitz_angle:.1f}",
                        "carbonyl_type": interaction_type,
                        "bs_int": bs_int,
                    }
                )

            # Add filter input
            filter_input = (
                ui.input(placeholder="Filter interactions (residue, atom, type, etc.)")
                .props("clearable outlined dense")
                .classes("w-full mb-2")
            )

            table = (
                ui.table(columns=columns, rows=rows, row_key="id")
                .classes("w-full")
                .props("dense")
            )
            filter_input.bind_value(table, "filter")

            # Add custom slot for visualize column to show button
            table.add_slot(
                "body-cell-visualize",
                """
                <q-td :props="props">
                    <q-btn size="sm" color="red" round dense icon="visibility"
                           @click="$parent.$emit('visualize', props.row)" />
                </q-td>
            """,
            )

            def show_carbonyl(e):
                idx = e.args["id"]
                carbonyl = self.analyzer.carbonyl_interactions[idx]
                self._show_carbonyl_interaction_visualization(carbonyl)

            table.on("visualize", show_carbonyl)

    def _update_n_pi_interactions_panel(self):
        """Update n→π* interactions panel."""
        with self.n_pi_panel:
            ui.label(
                f"n→π* Interactions ({len(self.analyzer.n_pi_interactions)})"
            ).classes("text-h5")

            # Create table
            columns = [
                {
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
                {
                    "name": "donor_res",
                    "label": "Donor Residue",
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
                    "name": "pi_res",
                    "label": "π Residue",
                    "field": "pi_res",
                    "align": "left",
                },
                {
                    "name": "distance",
                    "label": "Distance (Å)",
                    "field": "distance",
                    "align": "right",
                },
                {
                    "name": "angle",
                    "label": "Angle (°)",
                    "field": "angle",
                    "align": "right",
                },
                {
                    "name": "donor_element",
                    "label": "Donor Element",
                    "field": "donor_element",
                    "align": "left",
                },
                {
                    "name": "bs_int",
                    "label": "B/S",
                    "field": "bs_int",
                    "align": "center",
                },
            ]

            rows = []
            for idx, n_pi in enumerate(self.analyzer.n_pi_interactions):
                bs_int = "B" if n_pi.is_between_residues else "S"
                rows.append(
                    {
                        "id": idx,
                        "donor_res": n_pi.get_donor_residue(),
                        "donor_atom": n_pi.lone_pair_atom.name,
                        "pi_res": n_pi.get_acceptor_residue(),
                        "distance": f"{n_pi.distance:.2f}",
                        "angle": f"{n_pi.angle_to_plane:.1f}",
                        "donor_element": n_pi.donor_element,
                        "bs_int": bs_int,
                    }
                )

            # Add filter input
            filter_input = (
                ui.input(placeholder="Filter interactions (residue, atom, type, etc.)")
                .props("clearable outlined dense")
                .classes("w-full mb-2")
            )

            table = (
                ui.table(columns=columns, rows=rows, row_key="id")
                .classes("w-full")
                .props("dense")
            )
            filter_input.bind_value(table, "filter")

            # Add custom slot for visualize column to show button
            table.add_slot(
                "body-cell-visualize",
                """
                <q-td :props="props">
                    <q-btn size="sm" color="teal" round dense icon="visibility"
                           @click="$parent.$emit('visualize', props.row)" />
                </q-td>
            """,
            )

            def show_n_pi(e):
                idx = e.args["id"]
                n_pi = self.analyzer.n_pi_interactions[idx]
                self._show_n_pi_interaction_visualization(n_pi)

            table.on("visualize", show_n_pi)

    def _update_cooperativity_chains_panel(self):
        """Update cooperativity chains panel."""
        with self.cooperativity_panel:
            ui.label(
                f"Cooperativity Chains ({len(self.analyzer.cooperativity_chains)})"
            ).classes("text-h5")

            # Display each chain as a card showing the chain structure
            for idx, chain in enumerate(self.analyzer.cooperativity_chains):
                with ui.card().classes("w-full q-mb-md"):
                    # Chain header
                    with ui.card_section().classes("bg-grey-2 w-full"):
                        with ui.row().classes("w-full items-center justify-between"):
                            with ui.column():
                                ui.label(f"Chain {idx + 1}").classes("text-h6")
                                ui.label(
                                    f"Length: {chain.chain_length} interactions"
                                ).classes("text-caption")
                            ui.button(
                                "View Graph",
                                icon="account_tree",
                                on_click=lambda c=chain: self._show_cooperativity_chain_visualization(
                                    c
                                ),
                            ).props("outline color=primary")

                    # Chain details
                    with ui.card_section():
                        ui.label("Interaction Chain:").classes("text-subtitle2 q-mb-sm")

                        # Display each interaction in the chain
                        for i, interaction in enumerate(chain.interactions):
                            with ui.row().classes("items-center q-mb-xs"):
                                # Interaction number
                                ui.label(f"{i + 1}.").classes("text-bold")

                                # Donor
                                ui.label(interaction.get_donor_residue()).classes(
                                    "text-primary"
                                )

                                # Arrow with interaction type
                                int_type = interaction.get_interaction_type()
                                if (
                                    "hydrogen" in int_type.lower()
                                    or "h-bond" in int_type.lower()
                                ):
                                    arrow_icon = "→"
                                    int_label = "H-Bond"
                                    color = "blue"
                                elif (
                                    "halogen" in int_type.lower()
                                    or "x-bond" in int_type.lower()
                                ):
                                    arrow_icon = "⇢"
                                    int_label = "X-Bond"
                                    color = "orange"
                                elif "pi" in int_type.lower():
                                    arrow_icon = "⤏"
                                    int_label = "π-Int"
                                    color = "green"
                                else:
                                    arrow_icon = "→"
                                    int_label = int_type
                                    color = "grey"

                                ui.label(arrow_icon).classes(f"text-{color}")
                                ui.label(f"[{int_label}]").classes(
                                    f"text-caption text-{color}"
                                )
                                ui.label(arrow_icon).classes(f"text-{color}")

                                # Acceptor
                                ui.label(interaction.get_acceptor_residue()).classes(
                                    "text-primary"
                                )

                                # Distance
                                distance = interaction.get_donor_acceptor_distance()
                                ui.label(f"({distance:.2f} Å)").classes(
                                    "text-caption text-grey"
                                )

                        # Chain summary
                        ui.separator().classes("q-my-sm")
                        with ui.row().classes(
                            "w-full justify-between text-caption text-grey"
                        ):
                            ui.label(f"Start: {chain.get_donor_residue()}")
                            ui.label(f"End: {chain.get_acceptor_residue()}")
                            ui.label(
                                f"Total span: {chain.get_donor_acceptor_distance():.2f} Å"
                            )

    def _show_pi_pi_stacking_visualization(self, pi_pi):
        """Show 3D visualization of π-π stacking in a dialog."""
        import random

        viewer_id = f"pipi_viewer_{random.randint(1000, 9999)}"

        def export_png():
            """Export viewer as PNG."""
            filename = (
                f"{pi_pi.ring1_residue}_to_{pi_pi.ring2_residue}_pipi.png".replace(
                    ":", "_"
                )
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
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
            lambda: self._initialize_pi_pi_stacking_viewer(viewer_id, pi_pi),
            once=True,
        )

    def _initialize_pi_pi_stacking_viewer(self, viewer_id: str, pi_pi):
        """Initialize py3Dmol viewer for π-π stacking visualization.

        :param viewer_id: Unique ID for the viewer div
        :param pi_pi: Pi-pi stacking interaction
        """
        javascript = generate_pi_pi_stacking_viewer_js(
            pi_pi, self.pdb_content, viewer_id
        )
        ui.run_javascript(javascript)

    def _show_carbonyl_interaction_visualization(self, carbonyl):
        """Show 3D visualization of carbonyl n→π* interaction in a dialog."""
        import random

        viewer_id = f"carbonyl_viewer_{random.randint(1000, 9999)}"

        def export_png():
            """Export viewer as PNG."""
            filename = f"{carbonyl.get_donor_residue()}_to_{carbonyl.get_acceptor_residue()}_carbonyl.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
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
            lambda: self._initialize_carbonyl_interaction_viewer(viewer_id, carbonyl),
            once=True,
        )

    def _initialize_carbonyl_interaction_viewer(self, viewer_id: str, carbonyl):
        """Initialize py3Dmol viewer for carbonyl n→π* visualization.

        :param viewer_id: Unique ID for the viewer div
        :param carbonyl: Carbonyl interaction
        """
        javascript = generate_carbonyl_interaction_viewer_js(
            carbonyl, self.pdb_content, viewer_id
        )
        ui.run_javascript(javascript)

    def _show_n_pi_interaction_visualization(self, n_pi):
        """Show 3D visualization of n→π* interaction in a dialog."""
        import random

        viewer_id = f"npi_viewer_{random.randint(1000, 9999)}"

        def export_png():
            """Export viewer as PNG."""
            filename = f"{n_pi.get_donor_residue()}_to_{n_pi.get_acceptor_residue()}_npi.png".replace(
                ":", "_"
            )
            ui.run_javascript(generate_png_export_js(viewer_id, filename))

        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 900px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Export PNG", icon="download", on_click=export_png).props(
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
            lambda: self._initialize_n_pi_interaction_viewer(viewer_id, n_pi),
            once=True,
        )

    def _initialize_n_pi_interaction_viewer(self, viewer_id: str, n_pi):
        """Initialize py3Dmol viewer for n→π* interaction visualization.

        :param viewer_id: Unique ID for the viewer div
        :param n_pi: N-pi interaction
        """
        javascript = generate_n_pi_interaction_viewer_js(
            n_pi, self.pdb_content, viewer_id
        )
        ui.run_javascript(javascript)

    def _show_cooperativity_chain_visualization(self, chain):
        """Show 2D graph visualization of cooperativity chain in a dialog."""
        try:
            from ...server.app import UPLOADS_DIR

            # Generate both SVG (for display) and PNG (for export) using reusable function
            # Replace spaces and special chars in chain type for clean filenames
            safe_chain_type = chain.chain_type.replace(" ", "_").replace("-", "_")
            filename_prefix = f"chain_{safe_chain_type}_{chain.chain_length}"
            svg_content, png_path = render_chain_for_web(
                chain,
                output_dir=UPLOADS_DIR,
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

                        # Fix SVG dimensions to match viewBox for consistency across Graphviz versions
                        import re

                        # Extract viewBox dimensions
                        viewbox_match = re.search(
                            r'viewBox="([^"]+)"', svg_content
                        )
                        if viewbox_match:
                            viewbox = viewbox_match.group(1).split()
                            if len(viewbox) == 4:
                                vb_width = float(viewbox[2])
                                vb_height = float(viewbox[3])

                                # Set width and height to match viewBox dimensions
                                fixed_svg = re.sub(
                                    r'width="[^"]*"',
                                    f'width="{vb_width}pt"',
                                    svg_content
                                )
                                fixed_svg = re.sub(
                                    r'height="[^"]*"',
                                    f'height="{vb_height}pt"',
                                    fixed_svg
                                )
                                svg_content = fixed_svg

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

    def reset(self):
        """Reset the results panel and clear all analysis results."""
        self.analyzer = None
        self.pdb_content = None
        self.current_file = None
        self.container.clear()

        with self.container:
            ui.label("No analysis results yet").classes("text-h6 text-grey q-pa-lg")
