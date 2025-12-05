"""
Results display component with 3D visualization for HBAT web interface.

This module provides results display with integrated py3Dmol visualization
for molecular interactions.
"""

import math
import tempfile
from pathlib import Path
from typing import Optional

from nicegui import ui

from ...core.analysis import NPMolecularInteractionAnalyzer
from ...visualization.chain_graph import render_chain_graphviz


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
            # Create tabs dynamically based on available interactions
            with ui.tabs().classes("w-full") as self.tabs:
                # Summary tab - always shown
                ui.tab("summary", label="Summary", icon="summarize")

                # Conditionally add interaction type tabs
                if self.analyzer.hydrogen_bonds:
                    ui.tab("hydrogen", label=f"Hydrogen Bonds ({len(self.analyzer.hydrogen_bonds)})", icon="link")

                if self.analyzer.halogen_bonds:
                    ui.tab("halogen", label=f"Halogen Bonds ({len(self.analyzer.halogen_bonds)})", icon="whatshot")

                if self.analyzer.pi_interactions:
                    ui.tab("pi", label=f"π Interactions ({len(self.analyzer.pi_interactions)})", icon="fiber_manual_record")

                if hasattr(self.analyzer, "pi_pi_interactions") and self.analyzer.pi_pi_interactions:
                    ui.tab("pi_pi", label=f"π-π Stacking ({len(self.analyzer.pi_pi_interactions)})", icon="layers")

                if hasattr(self.analyzer, "carbonyl_interactions") and self.analyzer.carbonyl_interactions:
                    ui.tab("carbonyl", label=f"Carbonyl n→π* ({len(self.analyzer.carbonyl_interactions)})", icon="architecture")

                if hasattr(self.analyzer, "n_pi_interactions") and self.analyzer.n_pi_interactions:
                    ui.tab("n_pi", label=f"n→π* ({len(self.analyzer.n_pi_interactions)})", icon="arrow_forward")

                if hasattr(self.analyzer, "cooperativity_chains") and self.analyzer.cooperativity_chains:
                    ui.tab("cooperativity", label=f"Cooperativity ({len(self.analyzer.cooperativity_chains)})", icon="account_tree")

            # Create tab panels
            with ui.tab_panels(self.tabs, value="summary").classes("w-full") as self.tab_panels:
                # Summary panel - always shown
                with ui.tab_panel("summary"):
                    self.summary_panel = ui.column().classes("w-full")
                    self._update_summary_panel()

                # Conditionally add interaction panels
                if self.analyzer.hydrogen_bonds:
                    with ui.tab_panel("hydrogen"):
                        self.hydrogen_panel = ui.column().classes("w-full")
                        self._update_hydrogen_bonds_panel()

                if self.analyzer.halogen_bonds:
                    with ui.tab_panel("halogen"):
                        self.halogen_panel = ui.column().classes("w-full")
                        self._update_halogen_bonds_panel()

                if self.analyzer.pi_interactions:
                    with ui.tab_panel("pi"):
                        self.pi_panel = ui.column().classes("w-full")
                        self._update_pi_interactions_panel()

                if hasattr(self.analyzer, "pi_pi_interactions") and self.analyzer.pi_pi_interactions:
                    with ui.tab_panel("pi_pi"):
                        self.pi_pi_panel = ui.column().classes("w-full")
                        self._update_pi_pi_stacking_panel()

                if hasattr(self.analyzer, "carbonyl_interactions") and self.analyzer.carbonyl_interactions:
                    with ui.tab_panel("carbonyl"):
                        self.carbonyl_panel = ui.column().classes("w-full")
                        self._update_carbonyl_interactions_panel()

                if hasattr(self.analyzer, "n_pi_interactions") and self.analyzer.n_pi_interactions:
                    with ui.tab_panel("n_pi"):
                        self.n_pi_panel = ui.column().classes("w-full")
                        self._update_n_pi_interactions_panel()

                if hasattr(self.analyzer, "cooperativity_chains") and self.analyzer.cooperativity_chains:
                    with ui.tab_panel("cooperativity"):
                        self.cooperativity_panel = ui.column().classes("w-full")
                        self._update_cooperativity_chains_panel()

    def _update_summary_panel(self):
        """Update summary statistics panel."""
        with self.summary_panel:
            ui.label("Analysis Summary").classes("text-h5")
            ui.label(f"File: {self.current_file}").classes("text-subtitle1")

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

            with ui.card().classes("w-full"):
                ui.label("Cooperativity Chains").classes("text-h6")
                ui.label(str(summary["cooperativity_chains"]["count"])).classes(
                    "text-h4 text-blue"
                )

    def _update_hydrogen_bonds_panel(self):
        """Update hydrogen bonds panel with table and visualization."""
        with self.hydrogen_panel:
            ui.label(f"Hydrogen Bonds ({len(self.analyzer.hydrogen_bonds)})").classes(
                "text-h5"
            )

            # Create table
            columns = [
                {"name": "donor", "label": "Donor", "field": "donor", "align": "left"},
                {
                    "name": "hydrogen",
                    "label": "Hydrogen",
                    "field": "hydrogen",
                    "align": "left",
                },
                {
                    "name": "acceptor",
                    "label": "Acceptor",
                    "field": "acceptor",
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
                {"name": "type", "label": "Type", "field": "type", "align": "left"},
                {
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
            ]

            rows = []
            for idx, hb in enumerate(self.analyzer.hydrogen_bonds):
                rows.append(
                    {
                        "id": idx,
                        "donor": f"{hb.get_donor_residue()}",
                        "hydrogen": hb.hydrogen.name,
                        "acceptor": f"{hb.get_acceptor_residue()}",
                        "distance": f"{hb.distance:.2f}",
                        "angle": f"{hb.angle:.1f}",
                        "type": hb.bond_type,
                    }
                )

            table = ui.table(columns=columns, rows=rows, row_key="id").classes(
                "w-full"
            )
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
        with ui.dialog().props("persistent") as dialog, ui.card().style("width: 900px; max-width: 90vw;"):
            with ui.card_actions().classes("justify-end"):
                ui.button("Close", on_click=dialog.close).props("color=primary")

            with ui.card_section().classes("q-pa-md"):
                ui.label(
                    f"{hb.donor.res_name}{hb.donor.res_seq} → {hb.acceptor.res_name}{hb.acceptor.res_seq}"
                ).classes("text-subtitle1 q-mb-md")

                # Create unique viewer ID
                import random
                viewer_id = f"hb_viewer_{random.randint(1000, 9999)}"

                # Create div for viewer - container must have actual width
                ui.html(f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>', sanitize=False)


        dialog.open()

        # Initialize viewer after dialog is opened and rendered
        ui.timer(0.1, lambda: self._initialize_hydrogen_bond_viewer(viewer_id, hb), once=True)

    def _initialize_hydrogen_bond_viewer(self, viewer_id: str, hb):
        """Initialize py3Dmol viewer for hydrogen bond visualization.

        :param viewer_id: Unique ID for the viewer div
        :param hb: Hydrogen bond interaction
        """
        donor_chain = hb.donor.chain_id
        donor_resi = hb.donor.res_seq
        acceptor_chain = hb.acceptor.chain_id
        acceptor_resi = hb.acceptor.res_seq

        # Escape PDB content for JavaScript - escape backslashes first, then quotes
        pdb_escaped = self.pdb_content.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')

        javascript = f"""
        (function() {{
            // Wait for both 3Dmol to load and div to be ready
            function init3Dmol() {{
                if (typeof $3Dmol === 'undefined') {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                let element = document.getElementById("{viewer_id}");
                if (!element) {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                try {{
                    let viewer = $3Dmol.createViewer("{viewer_id}", {{backgroundColor: 'white'}});
                    let pdbData = `{pdb_escaped}`;

                    viewer.addModel(pdbData, "pdb");
                    viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                    // Show donor residue
                    viewer.addStyle({{chain: '{donor_chain}', resi: {donor_resi}}},
                                   {{stick: {{colorscheme: 'cyanCarbon', radius: 0.25}}}});

                    // Show acceptor residue
                    viewer.addStyle({{chain: '{acceptor_chain}', resi: {acceptor_resi}}},
                                   {{stick: {{colorscheme: 'orangeCarbon', radius: 0.25}}}});

                    // Add dashed line for hydrogen bond
                    viewer.addCylinder({{
                        start: {{x: {hb.hydrogen.coords.x}, y: {hb.hydrogen.coords.y}, z: {hb.hydrogen.coords.z}}},
                        end: {{x: {hb.acceptor.coords.x}, y: {hb.acceptor.coords.y}, z: {hb.acceptor.coords.z}}},
                        radius: 0.15,
                        color: 'yellow',
                        dashed: true
                    }});

                    // Add labels
                    viewer.addLabel("{hb.donor.res_name}{hb.donor.res_seq}",
                                   {{position: {{x: {hb.donor.coords.x}, y: {hb.donor.coords.y}, z: {hb.donor.coords.z}}},
                                    backgroundColor: 'cyan', fontColor: 'black', fontSize: 14}});

                    viewer.addLabel("{hb.acceptor.res_name}{hb.acceptor.res_seq}",
                                   {{position: {{x: {hb.acceptor.coords.x}, y: {hb.acceptor.coords.y}, z: {hb.acceptor.coords.z}}},
                                    backgroundColor: 'orange', fontColor: 'white', fontSize: 14}});

                    viewer.zoomTo({{chain: ['{donor_chain}', '{acceptor_chain}'], resi: [{donor_resi}, {acceptor_resi}]}});
                    viewer.render();
                    viewer.zoom(1.2);

                    // Force multiple resizes to ensure canvas is sized correctly
                    setTimeout(function() {{
                        viewer.resize();
                        viewer.render();
                    }}, 100);
                    setTimeout(function() {{
                        viewer.resize();
                        viewer.render();
                    }}, 300);
                    setTimeout(function() {{
                        viewer.resize();
                        viewer.render();
                    }}, 500);
                }} catch (error) {{
                    console.error("Error creating viewer:", error);
                }}
            }}

            // Start with a delay to ensure div is in DOM and sized
            setTimeout(init3Dmol, 400);
        }})();
        """
        ui.run_javascript(javascript)

    def _update_halogen_bonds_panel(self):
        """Update halogen bonds panel."""
        with self.halogen_panel:
            ui.label(f"Halogen Bonds ({len(self.analyzer.halogen_bonds)})").classes(
                "text-h5"
            )

            # Similar to hydrogen bonds but for halogen bonds
            columns = [
                {"name": "donor", "label": "Donor", "field": "donor", "align": "left"},
                {
                    "name": "halogen",
                    "label": "Halogen",
                    "field": "halogen",
                    "align": "left",
                },
                {
                    "name": "acceptor",
                    "label": "Acceptor",
                    "field": "acceptor",
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
                    "name": "visualize",
                    "label": "3D View",
                    "field": "visualize",
                    "align": "center",
                },
            ]

            rows = []
            for idx, xb in enumerate(self.analyzer.halogen_bonds):
                rows.append(
                    {
                        "id": idx,
                        "donor": f"{xb.get_donor_residue()}",
                        "halogen": xb.halogen.name,
                        "acceptor": f"{xb.get_acceptor_residue()}",
                        "distance": f"{xb.distance:.2f}",
                        "angle": f"{xb.angle:.1f}",
                    }
                )

            table = ui.table(columns=columns, rows=rows, row_key="id").classes(
                "w-full"
            )
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
        with ui.dialog().props("persistent") as dialog, ui.card().style("width: 900px; max-width: 90vw;"):
            with ui.card_section().classes("bg-orange text-white"):
                with ui.row().classes("w-full items-center"):
                    ui.label("3D Visualization - Halogen Bond").classes("text-h6")
                    ui.space()
                    ui.button(icon="close", on_click=dialog.close).props("flat round dense color=white")

            with ui.card_section().classes("q-pa-md"):
                ui.label(
                    f"{xb.donor_atom.res_name}{xb.donor_atom.res_seq} → {xb.acceptor.res_name}{xb.acceptor.res_seq}"
                ).classes("text-subtitle1 q-mb-md")

                # Create unique viewer ID
                import random
                viewer_id = f"xb_viewer_{random.randint(1000, 9999)}"

                # Create div for viewer - container must have actual width
                ui.html(f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>', sanitize=False)

            with ui.card_actions().classes("justify-end"):
                ui.button("Close", on_click=dialog.close).props("color=orange")

        dialog.open()

        # Initialize viewer after dialog is opened and rendered
        ui.timer(0.1, lambda: self._initialize_halogen_bond_viewer(viewer_id, xb), once=True)

    def _initialize_halogen_bond_viewer(self, viewer_id: str, xb):
        """Initialize py3Dmol viewer for halogen bond visualization.

        :param viewer_id: Unique ID for the viewer div
        :param xb: Halogen bond interaction
        """
        donor_chain = xb.donor_atom.chain_id
        donor_resi = xb.donor_atom.res_seq
        acceptor_chain = xb.acceptor.chain_id
        acceptor_resi = xb.acceptor.res_seq

        # Escape PDB content for JavaScript - escape backslashes first, then quotes
        pdb_escaped = self.pdb_content.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')

        javascript = f"""
        (function() {{
            // Wait for both 3Dmol to load and div to be ready
            function init3Dmol() {{
                if (typeof $3Dmol === 'undefined') {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                let element = document.getElementById("{viewer_id}");
                if (!element) {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                try {{
                    let viewer = $3Dmol.createViewer("{viewer_id}", {{backgroundColor: 'white'}});
                    let pdbData = `{pdb_escaped}`;

                    viewer.addModel(pdbData, "pdb");
                    viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                    viewer.addStyle({{chain: '{donor_chain}', resi: {donor_resi}}},
                                   {{stick: {{colorscheme: 'purpleCarbon', radius: 0.25}}}});

                    viewer.addStyle({{chain: '{acceptor_chain}', resi: {acceptor_resi}}},
                                   {{stick: {{colorscheme: 'orangeCarbon', radius: 0.25}}}});

                    viewer.addCylinder({{
                        start: {{x: {xb.halogen.coords.x}, y: {xb.halogen.coords.y}, z: {xb.halogen.coords.z}}},
                        end: {{x: {xb.acceptor.coords.x}, y: {xb.acceptor.coords.y}, z: {xb.acceptor.coords.z}}},
                        radius: 0.15,
                        color: 'orange',
                        dashed: true
                    }});

                    viewer.zoomTo({{chain: ['{donor_chain}', '{acceptor_chain}'], resi: [{donor_resi}, {acceptor_resi}]}});
                    viewer.render();
                    viewer.zoom(1.2);

                    // Force multiple resizes to ensure canvas is sized correctly
                    setTimeout(function() {{
                        viewer.resize();
                        viewer.render();
                    }}, 100);
                    setTimeout(function() {{
                        viewer.resize();
                        viewer.render();
                    }}, 300);
                    setTimeout(function() {{
                        viewer.resize();
                        viewer.render();
                    }}, 500);
                }} catch (error) {{
                    console.error("Error creating viewer:", error);
                }}
            }}

            // Start with a delay to ensure div is in DOM and sized
            setTimeout(init3Dmol, 400);
        }})();
        """
        ui.run_javascript(javascript)

    def _show_pi_interaction_visualization(self, pi):
        """Show 3D visualization of a π interaction in a dialog."""
        import random

        with ui.dialog().props("persistent") as dialog, ui.card().style("width: 900px; max-width: 90vw;"):
            with ui.card_section().classes("bg-primary text-white"):
                with ui.row().classes("w-full items-center"):
                    ui.label("3D Visualization - π Interaction").classes("text-h6")
                    ui.space()
                    ui.button(icon="close", on_click=dialog.close).props("flat round dense color=white")

            with ui.card_section().classes("q-pa-md"):
                ui.label(f"{pi.get_donor_residue()} → {pi.get_acceptor_residue()}").classes("text-subtitle1 q-mb-md")

                viewer_id = f"pi_viewer_{random.randint(1000, 9999)}"
                ui.html(f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>', sanitize=False)

            with ui.card_actions().classes("justify-end"):
                ui.button("Close", on_click=dialog.close).props("color=primary")

        dialog.open()
        ui.timer(0.1, lambda: self._initialize_pi_interaction_viewer(viewer_id, pi), once=True)

    def _initialize_pi_interaction_viewer(self, viewer_id: str, pi):
        """Initialize py3Dmol viewer for π interaction visualization."""
        javascript = f"""
        (function() {{
            function init3Dmol() {{
                if (typeof $3Dmol === 'undefined') {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                let element = document.getElementById("{viewer_id}");
                if (!element) {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                try {{
                    let viewer = $3Dmol.createViewer("{viewer_id}", {{backgroundColor: 'white'}});

                    // Add the entire structure
                    let pdb_data = `{self.pdb_content}`;
                    viewer.addModel(pdb_data, "pdb");

                    // Style: show only relevant residues
                    viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                    // Highlight donor residue (C-H donor) in cyan
                    viewer.addStyle({{resi: {pi.donor_res_seq}, chain: '{pi.donor_chain_id}'}},
                                   {{stick: {{colorscheme: 'cyanCarbon', radius: 0.25}}}});

                    // Highlight π acceptor residue (aromatic ring) in magenta
                    viewer.addStyle({{resi: {pi.acceptor_res_seq}, chain: '{pi.acceptor_chain_id}'}},
                                   {{stick: {{colorscheme: 'magentaCarbon', radius: 0.3}}}});

                    // Add sphere at π center to show the interaction point
                    viewer.addSphere({{
                        center: {{x: {pi.pi_center.x}, y: {pi.pi_center.y}, z: {pi.pi_center.z}}},
                        radius: 0.3,
                        color: 'magenta',
                        alpha: 0.7
                    }});

                    // Highlight the hydrogen/X-atom explicitly as a sphere
                    viewer.addStyle({{serial: {pi.hydrogen.serial}}},
                                   {{sphere: {{color: 'yellow', radius: 0.5}}}});

                    // Add residue labels
                    viewer.addLabel('{pi.get_donor_residue()}',
                                  {{position: {{x: {pi.donor.coords.x}, y: {pi.donor.coords.y}, z: {pi.donor.coords.z}}},
                                   backgroundColor: 'cyan', fontColor: 'black', fontSize: 10}});
                    viewer.addLabel('{pi.get_acceptor_residue()}',
                                  {{position: {{x: {pi.pi_center.x}, y: {pi.pi_center.y}, z: {pi.pi_center.z}}},
                                   backgroundColor: 'magenta', fontColor: 'white', fontSize: 10}});

                    // Add distance label
                    viewer.addLabel('{pi.distance:.2f} Å',
                                  {{position: {{x: ({pi.hydrogen.coords.x} + {pi.pi_center.x})/2,
                                              y: ({pi.hydrogen.coords.y} + {pi.pi_center.y})/2,
                                              z: ({pi.hydrogen.coords.z} + {pi.pi_center.z})/2}},
                                   backgroundColor: 'black', fontColor: 'white', fontSize: 12}});

                    // Add dashed line from X-atom to π center
                    viewer.addCylinder({{
                        start: {{x: {pi.hydrogen.coords.x}, y: {pi.hydrogen.coords.y}, z: {pi.hydrogen.coords.z}}},
                        end: {{x: {pi.pi_center.x}, y: {pi.pi_center.y}, z: {pi.pi_center.z}}},
                        radius: 0.1,
                        color: 'yellow',
                        dashed: true
                    }});

                    // Center on interaction
                    viewer.zoomTo({{chain: ['{pi.donor_chain_id}', '{pi.acceptor_chain_id}'], resi: [{pi.donor_res_seq}, {pi.acceptor_res_seq}]}});
                    viewer.render();

                    // Multiple resize calls for proper sizing
                    setTimeout(function() {{
                        viewer.resize();
                        viewer.render();
                    }}, 100);
                    setTimeout(function() {{
                        viewer.resize();
                        viewer.render();
                    }}, 300);
                    setTimeout(function() {{
                        viewer.resize();
                        viewer.render();
                    }}, 500);
                }} catch (error) {{
                    console.error("Error creating viewer:", error);
                }}
            }}

            // Start with a delay to ensure div is in DOM and sized
            setTimeout(init3Dmol, 400);
        }})();
        """
        ui.run_javascript(javascript)

    def _update_pi_interactions_panel(self):
        """Update π interactions panel."""
        with self.pi_panel:
            ui.label(f"π Interactions ({len(self.analyzer.pi_interactions)})").classes(
                "text-h5"
            )

            # Create table
            columns = [
                {"name": "donor", "label": "Donor", "field": "donor", "align": "left"},
                {"name": "interaction", "label": "X-atom", "field": "interaction", "align": "left"},
                {"name": "pi_residue", "label": "π Residue", "field": "pi_residue", "align": "left"},
                {"name": "distance", "label": "Distance (Å)", "field": "distance", "align": "right"},
                {"name": "angle", "label": "Angle (°)", "field": "angle", "align": "right"},
                {"name": "type", "label": "Subtype", "field": "type", "align": "left"},
                {"name": "visualize", "label": "3D View", "field": "visualize", "align": "center"},
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

                rows.append({
                    "id": idx,
                    "donor": f"{pi.get_donor_residue()}:{pi.donor.name}",
                    "interaction": pi.hydrogen.name,
                    "pi_residue": pi.get_acceptor_residue(),
                    "distance": f"{pi.distance:.2f}",
                    "angle": f"{math.degrees(pi.angle):.1f}",
                    "type": subtype,
                })

            # Create table with custom styling
            table = ui.table(columns=columns, rows=rows, row_key="id").classes("w-full").props("dense")

            # Add custom slot for visualize column to show button
            table.add_slot('body-cell-visualize', '''
                <q-td :props="props">
                    <q-btn size="sm" color="green" round dense icon="visibility"
                           @click="$parent.$emit('visualize', props.row)" />
                </q-td>
            ''')

            def show_pi_interaction(e):
                pi_idx = e.args["id"]
                pi = self.analyzer.pi_interactions[pi_idx]
                self._show_pi_interaction_visualization(pi)

            table.on("visualize", show_pi_interaction)

    def _update_pi_pi_stacking_panel(self):
        """Update π-π stacking panel."""
        with self.pi_pi_panel:
            ui.label(
                f"π-π Stacking ({len(self.analyzer.pi_pi_interactions)})"
            ).classes("text-h5")

            # Create table
            columns = [
                {"name": "ring1", "label": "Ring 1", "field": "ring1", "align": "left"},
                {"name": "ring2", "label": "Ring 2", "field": "ring2", "align": "left"},
                {"name": "distance", "label": "Distance (Å)", "field": "distance", "align": "right"},
                {"name": "angle", "label": "Plane Angle (°)", "field": "angle", "align": "right"},
                {"name": "offset", "label": "Offset (Å)", "field": "offset", "align": "right"},
                {"name": "type", "label": "Type", "field": "type", "align": "left"},
                {"name": "visualize", "label": "3D View", "field": "visualize", "align": "center"},
            ]

            rows = []
            for idx, pi_pi in enumerate(self.analyzer.pi_pi_interactions):
                rows.append({
                    "id": idx,
                    "ring1": f"{pi_pi.ring1_residue}",
                    "ring2": f"{pi_pi.ring2_residue}",
                    "distance": f"{pi_pi._distance:.2f}",
                    "angle": f"{pi_pi.plane_angle:.1f}",
                    "offset": f"{pi_pi.offset:.2f}",
                    "type": pi_pi.stacking_type.capitalize(),
                })

            table = ui.table(columns=columns, rows=rows, row_key="id").classes("w-full").props("dense")

            # Add custom slot for visualize column to show button
            table.add_slot('body-cell-visualize', '''
                <q-td :props="props">
                    <q-btn size="sm" color="purple" round dense icon="visibility"
                           @click="$parent.$emit('visualize', props.row)" />
                </q-td>
            ''')

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
                {"name": "donor", "label": "Donor C=O", "field": "donor", "align": "left"},
                {"name": "acceptor", "label": "Acceptor C=O", "field": "acceptor", "align": "left"},
                {"name": "distance", "label": "O···C Distance (Å)", "field": "distance", "align": "right"},
                {"name": "angle", "label": "Bürgi-Dunitz Angle (°)", "field": "angle", "align": "right"},
                {"name": "type", "label": "Type", "field": "type", "align": "left"},
                {"name": "visualize", "label": "3D View", "field": "visualize", "align": "center"},
            ]

            rows = []
            for idx, carbonyl in enumerate(self.analyzer.carbonyl_interactions):
                interaction_type = "Backbone" if carbonyl.is_backbone else "Sidechain"
                rows.append({
                    "id": idx,
                    "donor": f"{carbonyl.get_donor_residue()}",
                    "acceptor": f"{carbonyl.get_acceptor_residue()}",
                    "distance": f"{carbonyl.distance:.2f}",
                    "angle": f"{carbonyl.burgi_dunitz_angle:.1f}",
                    "type": interaction_type,
                })

            table = ui.table(columns=columns, rows=rows, row_key="id").classes("w-full").props("dense")

            # Add custom slot for visualize column to show button
            table.add_slot('body-cell-visualize', '''
                <q-td :props="props">
                    <q-btn size="sm" color="red" round dense icon="visibility"
                           @click="$parent.$emit('visualize', props.row)" />
                </q-td>
            ''')

            def show_carbonyl(e):
                idx = e.args["id"]
                carbonyl = self.analyzer.carbonyl_interactions[idx]
                self._show_carbonyl_interaction_visualization(carbonyl)

            table.on("visualize", show_carbonyl)

    def _update_n_pi_interactions_panel(self):
        """Update n→π* interactions panel."""
        with self.n_pi_panel:
            ui.label(f"n→π* Interactions ({len(self.analyzer.n_pi_interactions)})").classes(
                "text-h5"
            )

            # Create table
            columns = [
                {"name": "donor", "label": "Lone Pair Atom", "field": "donor", "align": "left"},
                {"name": "pi_residue", "label": "π System", "field": "pi_residue", "align": "left"},
                {"name": "distance", "label": "Distance (Å)", "field": "distance", "align": "right"},
                {"name": "angle", "label": "Angle to Plane (°)", "field": "angle", "align": "right"},
                {"name": "subtype", "label": "Subtype", "field": "subtype", "align": "left"},
                {"name": "visualize", "label": "3D View", "field": "visualize", "align": "center"},
            ]

            rows = []
            for idx, n_pi in enumerate(self.analyzer.n_pi_interactions):
                rows.append({
                    "id": idx,
                    "donor": f"{n_pi.get_donor_residue()}:{n_pi.lone_pair_atom.name}",
                    "pi_residue": f"{n_pi.get_acceptor_residue()}",
                    "distance": f"{n_pi.distance:.2f}",
                    "angle": f"{n_pi.angle_to_plane:.1f}",
                    "subtype": n_pi.subtype,
                })

            table = ui.table(columns=columns, rows=rows, row_key="id").classes("w-full").props("dense")

            # Add custom slot for visualize column to show button
            table.add_slot('body-cell-visualize', '''
                <q-td :props="props">
                    <q-btn size="sm" color="teal" round dense icon="visibility"
                           @click="$parent.$emit('visualize', props.row)" />
                </q-td>
            ''')

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
                    with ui.card_section().classes("bg-grey-2"):
                        with ui.row().classes("w-full items-center justify-between"):
                            with ui.column():
                                ui.label(f"Chain {idx + 1}").classes("text-h6")
                                ui.label(f"Length: {chain.chain_length} interactions").classes("text-caption")
                            ui.button(
                                "View Graph",
                                icon="account_tree",
                                on_click=lambda c=chain: self._show_cooperativity_chain_visualization(c)
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
                                ui.label(interaction.get_donor_residue()).classes("text-primary")

                                # Arrow with interaction type
                                int_type = interaction.get_interaction_type()
                                if "hydrogen" in int_type.lower() or "h-bond" in int_type.lower():
                                    arrow_icon = "→"
                                    int_label = "H-Bond"
                                    color = "blue"
                                elif "halogen" in int_type.lower() or "x-bond" in int_type.lower():
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
                                ui.label(f"[{int_label}]").classes(f"text-caption text-{color}")
                                ui.label(arrow_icon).classes(f"text-{color}")

                                # Acceptor
                                ui.label(interaction.get_acceptor_residue()).classes("text-primary")

                                # Distance
                                distance = interaction.get_donor_acceptor_distance()
                                ui.label(f"({distance:.2f} Å)").classes("text-caption text-grey")

                        # Chain summary
                        ui.separator().classes("q-my-sm")
                        with ui.row().classes("w-full justify-between text-caption text-grey"):
                            ui.label(f"Start: {chain.get_donor_residue()}")
                            ui.label(f"End: {chain.get_acceptor_residue()}")
                            ui.label(f"Total span: {chain.get_donor_acceptor_distance():.2f} Å")

    def _show_pi_pi_stacking_visualization(self, pi_pi):
        """Show 3D visualization of π-π stacking in a dialog."""
        import random

        with ui.dialog().props("persistent") as dialog, ui.card().style("width: 900px; max-width: 90vw;"):
            with ui.card_section().classes("bg-primary text-white"):
                with ui.row().classes("w-full items-center"):
                    ui.label("3D Visualization - π-π Stacking").classes("text-h6")
                    ui.space()
                    ui.button(icon="close", on_click=dialog.close).props("flat round dense color=white")

            with ui.card_section().classes("q-pa-md"):
                ui.label(f"{pi_pi.ring1_residue} ⇄ {pi_pi.ring2_residue}").classes("text-subtitle1 q-mb-md")

                viewer_id = f"pipi_viewer_{random.randint(1000, 9999)}"
                ui.html(f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>', sanitize=False)

            with ui.card_actions().classes("justify-end"):
                ui.button("Close", on_click=dialog.close).props("color=primary")

        dialog.open()
        ui.timer(0.1, lambda: self._initialize_pi_pi_stacking_viewer(viewer_id, pi_pi), once=True)

    def _initialize_pi_pi_stacking_viewer(self, viewer_id: str, pi_pi):
        """Initialize py3Dmol viewer for π-π stacking visualization."""
        # Get residue info from ring atoms
        ring1_res = pi_pi.ring1_atoms[0].res_seq if pi_pi.ring1_atoms else ""
        ring1_chain = pi_pi.ring1_atoms[0].chain_id if pi_pi.ring1_atoms else ""
        ring2_res = pi_pi.ring2_atoms[0].res_seq if pi_pi.ring2_atoms else ""
        ring2_chain = pi_pi.ring2_atoms[0].chain_id if pi_pi.ring2_atoms else ""

        javascript = f"""
        (function() {{
            function init3Dmol() {{
                if (typeof $3Dmol === 'undefined') {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                let element = document.getElementById("{viewer_id}");
                if (!element) {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                try {{
                    let viewer = $3Dmol.createViewer("{viewer_id}", {{backgroundColor: 'white'}});

                    let pdb_data = `{self.pdb_content}`;
                    viewer.addModel(pdb_data, "pdb");

                    viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                    // Highlight both aromatic rings
                    viewer.setStyle({{resi: '{ring1_res}', chain: '{ring1_chain}'}}, {{stick: {{colorscheme: 'cyanCarbon'}}}});
                    viewer.setStyle({{resi: '{ring2_res}', chain: '{ring2_chain}'}}, {{stick: {{colorscheme: 'magentaCarbon'}}}});

                    // Add distance label
                    viewer.addLabel('{pi_pi._distance:.2f} Å ({pi_pi.stacking_type})',
                                  {{position: {{x: {pi_pi.midpoint.x}, y: {pi_pi.midpoint.y}, z: {pi_pi.midpoint.z}}},
                                   backgroundColor: 'black', fontColor: 'white', fontSize: 12}});

                    // Add line between ring centers
                    viewer.addCylinder({{
                        start: {{x: {pi_pi.ring1_center.x}, y: {pi_pi.ring1_center.y}, z: {pi_pi.ring1_center.z}}},
                        end: {{x: {pi_pi.ring2_center.x}, y: {pi_pi.ring2_center.y}, z: {pi_pi.ring2_center.z}}},
                        radius: 0.1,
                        color: 'purple',
                        dashed: true
                    }});

                    viewer.zoomTo({{chain: ['{ring1_chain}', '{ring2_chain}'], resi: ['{ring1_res}', '{ring2_res}']}});
                    viewer.render();

                    setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 100);
                    setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 300);
                    setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 500);
                }} catch (error) {{
                    console.error("Error creating viewer:", error);
                }}
            }}

            setTimeout(init3Dmol, 400);
        }})();
        """
        ui.run_javascript(javascript)

    def _show_carbonyl_interaction_visualization(self, carbonyl):
        """Show 3D visualization of carbonyl n→π* interaction in a dialog."""
        import random

        with ui.dialog().props("persistent") as dialog, ui.card().style("width: 900px; max-width: 90vw;"):
            with ui.card_section().classes("bg-primary text-white"):
                with ui.row().classes("w-full items-center"):
                    ui.label("3D Visualization - Carbonyl n→π*").classes("text-h6")
                    ui.space()
                    ui.button(icon="close", on_click=dialog.close).props("flat round dense color=white")

            with ui.card_section().classes("q-pa-md"):
                ui.label(f"{carbonyl.get_donor_residue()} → {carbonyl.get_acceptor_residue()}").classes("text-subtitle1 q-mb-md")

                viewer_id = f"carbonyl_viewer_{random.randint(1000, 9999)}"
                ui.html(f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>', sanitize=False)

            with ui.card_actions().classes("justify-end"):
                ui.button("Close", on_click=dialog.close).props("color=primary")

        dialog.open()
        ui.timer(0.1, lambda: self._initialize_carbonyl_interaction_viewer(viewer_id, carbonyl), once=True)

    def _initialize_carbonyl_interaction_viewer(self, viewer_id: str, carbonyl):
        """Initialize py3Dmol viewer for carbonyl n→π* visualization."""
        javascript = f"""
        (function() {{
            function init3Dmol() {{
                if (typeof $3Dmol === 'undefined') {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                let element = document.getElementById("{viewer_id}");
                if (!element) {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                try {{
                    let viewer = $3Dmol.createViewer("{viewer_id}", {{backgroundColor: 'white'}});

                    let pdb_data = `{self.pdb_content}`;
                    viewer.addModel(pdb_data, "pdb");

                    viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                    // Highlight donor and acceptor residues
                    viewer.setStyle({{resi: '{carbonyl.donor_oxygen.res_seq}', chain: '{carbonyl.donor_oxygen.chain_id}'}},
                                  {{stick: {{colorscheme: 'redCarbon'}}}});
                    viewer.setStyle({{resi: '{carbonyl.acceptor_carbon.res_seq}', chain: '{carbonyl.acceptor_carbon.chain_id}'}},
                                  {{stick: {{colorscheme: 'blueCarbon'}}}});

                    // Add distance label
                    viewer.addLabel('{carbonyl.distance:.2f} Å',
                                  {{position: {{x: {carbonyl.donor_oxygen.coords.x}, y: {carbonyl.donor_oxygen.coords.y}, z: {carbonyl.donor_oxygen.coords.z}}},
                                   backgroundColor: 'black', fontColor: 'white', fontSize: 12}});

                    // Add line from donor O to acceptor C
                    viewer.addCylinder({{
                        start: {{x: {carbonyl.donor_oxygen.coords.x}, y: {carbonyl.donor_oxygen.coords.y}, z: {carbonyl.donor_oxygen.coords.z}}},
                        end: {{x: {carbonyl.acceptor_carbon.coords.x}, y: {carbonyl.acceptor_carbon.coords.y}, z: {carbonyl.acceptor_carbon.coords.z}}},
                        radius: 0.1,
                        color: 'orange',
                        dashed: true
                    }});

                    viewer.zoomTo({{chain: ['{carbonyl.donor_oxygen.chain_id}', '{carbonyl.acceptor_carbon.chain_id}'], resi: ['{carbonyl.donor_oxygen.res_seq}', '{carbonyl.acceptor_carbon.res_seq}']}});
                    viewer.render();

                    setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 100);
                    setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 300);
                    setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 500);
                }} catch (error) {{
                    console.error("Error creating viewer:", error);
                }}
            }}

            setTimeout(init3Dmol, 400);
        }})();
        """
        ui.run_javascript(javascript)

    def _show_n_pi_interaction_visualization(self, n_pi):
        """Show 3D visualization of n→π* interaction in a dialog."""
        import random

        with ui.dialog().props("persistent") as dialog, ui.card().style("width: 900px; max-width: 90vw;"):
            with ui.card_section().classes("bg-primary text-white"):
                with ui.row().classes("w-full items-center"):
                    ui.label("3D Visualization - n→π* Interaction").classes("text-h6")
                    ui.space()
                    ui.button(icon="close", on_click=dialog.close).props("flat round dense color=white")

            with ui.card_section().classes("q-pa-md"):
                ui.label(f"{n_pi.get_donor_residue()} → {n_pi.get_acceptor_residue()}").classes("text-subtitle1 q-mb-md")

                viewer_id = f"npi_viewer_{random.randint(1000, 9999)}"
                ui.html(f'<div id="{viewer_id}" style="width: 100%; height: 600px; min-width: 800px; position: relative;"></div>', sanitize=False)

            with ui.card_actions().classes("justify-end"):
                ui.button("Close", on_click=dialog.close).props("color=primary")

        dialog.open()
        ui.timer(0.1, lambda: self._initialize_n_pi_interaction_viewer(viewer_id, n_pi), once=True)

    def _initialize_n_pi_interaction_viewer(self, viewer_id: str, n_pi):
        """Initialize py3Dmol viewer for n→π* interaction visualization."""
        # Get the π system residue sequence number from the first pi atom
        pi_res_seq = n_pi.pi_atoms[0].res_seq if n_pi.pi_atoms else ""
        pi_chain_id = n_pi.pi_atoms[0].chain_id if n_pi.pi_atoms else ""

        javascript = f"""
        (function() {{
            function init3Dmol() {{
                if (typeof $3Dmol === 'undefined') {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                let element = document.getElementById("{viewer_id}");
                if (!element) {{
                    setTimeout(init3Dmol, 100);
                    return;
                }}

                try {{
                    let viewer = $3Dmol.createViewer("{viewer_id}", {{backgroundColor: 'white'}});

                    let pdb_data = `{self.pdb_content}`;
                    viewer.addModel(pdb_data, "pdb");

                    viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                    // Highlight lone pair atom residue
                    viewer.setStyle({{resi: '{n_pi.lone_pair_atom.res_seq}', chain: '{n_pi.lone_pair_atom.chain_id}'}},
                                  {{stick: {{colorscheme: 'orangeCarbon'}}}});

                    // Highlight π system residue
                    viewer.setStyle({{resi: '{pi_res_seq}', chain: '{pi_chain_id}'}},
                                  {{stick: {{colorscheme: 'tealCarbon'}}}});

                    // Add distance label
                    viewer.addLabel('{n_pi.distance:.2f} Å',
                                  {{position: {{x: {n_pi.lone_pair_atom.coords.x}, y: {n_pi.lone_pair_atom.coords.y}, z: {n_pi.lone_pair_atom.coords.z}}},
                                   backgroundColor: 'black', fontColor: 'white', fontSize: 12}});

                    // Add line from lone pair atom to π center
                    viewer.addCylinder({{
                        start: {{x: {n_pi.lone_pair_atom.coords.x}, y: {n_pi.lone_pair_atom.coords.y}, z: {n_pi.lone_pair_atom.coords.z}}},
                        end: {{x: {n_pi.pi_center.x}, y: {n_pi.pi_center.y}, z: {n_pi.pi_center.z}}},
                        radius: 0.1,
                        color: 'green',
                        dashed: true
                    }});

                    viewer.zoomTo({{chain: ['{n_pi.lone_pair_atom.chain_id}', '{pi_chain_id}'], resi: ['{n_pi.lone_pair_atom.res_seq}', '{pi_res_seq}']}});
                    viewer.render();

                    setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 100);
                    setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 300);
                    setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 500);
                }} catch (error) {{
                    console.error("Error creating viewer:", error);
                }}
            }}

            setTimeout(init3Dmol, 400);
        }})();
        """
        ui.run_javascript(javascript)

    def _show_cooperativity_chain_visualization(self, chain):
        """Show 2D graph visualization of cooperativity chain in a dialog."""
        try:
            # Generate SVG graph using chain_graph module
            with tempfile.TemporaryDirectory() as tmpdir:
                output_path = Path(tmpdir) / "chain_graph"

                # Generate the graph and get the SVG content
                dot = render_chain_graphviz(
                    chain,
                    engine="dot",
                    rankdir="LR",
                    filename=str(output_path),
                    format="svg",
                    view=False
                )

                # Read the SVG file
                svg_file = Path(f"{output_path}.svg")
                if svg_file.exists():
                    svg_content = svg_file.read_text()
                else:
                    svg_content = None

            if svg_content:
                with ui.dialog().props("persistent") as dialog, ui.card().style("width: 900px; max-width: 90vw;"):
                    with ui.card_section().classes("bg-primary text-white"):
                        with ui.row().classes("w-full items-center"):
                            ui.label("Cooperativity Chain Graph").classes("text-h6")
                            ui.space()
                            ui.button(icon="close", on_click=dialog.close).props("flat round dense color=white")

                    with ui.card_section().classes("q-pa-md"):
                        ui.label(f"Chain: {chain.chain_type} (Length: {chain.chain_length})").classes("text-subtitle1 q-mb-md")

                        # Display the SVG
                        ui.html(svg_content, sanitize=False).style("width: 100%; overflow: auto;")

                    with ui.card_actions().classes("justify-end"):
                        ui.button("Close", on_click=dialog.close).props("color=primary")

                dialog.open()
            else:
                ui.notify("Failed to generate chain graph", type="negative")

        except ImportError:
            ui.notify("Graphviz is not installed. Install with: pip install graphviz", type="warning")
        except Exception as e:
            ui.notify(f"Error generating chain graph: {str(e)}", type="negative")
