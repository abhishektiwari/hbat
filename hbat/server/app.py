"""
Main web application for HBAT using NiceGUI.

This module provides a web-based interface for molecular interaction analysis
with integrated 3D visualization using py3Dmol.
"""

import asyncio
import os
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Optional

from nicegui import app, ui

from ..core.analysis import NPMolecularInteractionAnalyzer
from ..export import (
    export_to_csv_files,
    export_to_json_single_file,
    export_to_txt_single_file,
)
from .components.parameter_panel import ParameterPanel
from .components.results_panel import WebResultsPanel
from .components.upload_panel import UploadPanel
from .session import SessionManager

# Base uploads directory (can be mounted as Docker volume)
UPLOADS_DIR = Path("uploads")
UPLOADS_DIR.mkdir(exist_ok=True)

# Sessions directory inside uploads for easy Docker volume management
SESSIONS_BASE_DIR = UPLOADS_DIR / "sessions"
SESSIONS_BASE_DIR.mkdir(exist_ok=True)

# Global session manager
session_manager = SessionManager(SESSIONS_BASE_DIR, session_timeout_hours=24)


class HBATWebApp:
    """Main HBAT web application class."""

    def __init__(self):
        """Initialize the HBAT web application."""
        self.analyzer: Optional[NPMolecularInteractionAnalyzer] = None
        self.current_file: Optional[str] = None
        self.current_file_path: Optional[Path] = None
        self.analysis_running = False
        self.pdb_content: Optional[str] = None

        # Session management
        self.session_id: Optional[str] = None
        self.session_dir: Optional[Path] = None

        # UI components
        self.upload_panel: Optional[UploadPanel] = None
        self.parameter_panel: Optional[ParameterPanel] = None
        self.results_panel: Optional[WebResultsPanel] = None

        # Status tracking
        self.status_label: Optional[ui.label] = None
        self.analyze_button: Optional[ui.button] = None
        self.download_structure_button: Optional[ui.button] = None

        # Stepper
        self.stepper: Optional[ui.stepper] = None

        # Navigation items
        self.nav_configure: Optional[ui.item] = None
        self.nav_results: Optional[ui.item] = None
        self.nav_export: Optional[ui.item] = None

    def _update_status_label(self, message: str = "Ready"):
        """Update status label with current file information.

        :param message: Status message to display
        """
        if self.status_label:
            if self.current_file and message == "Ready":
                self.status_label.text = f"Ready - File: {self.current_file}"
            else:
                self.status_label.text = message

    def _show_citation_dialog(self):
        """Show citation dialog."""
        with (
            ui.dialog().props("persistent") as dialog,
            ui.card().style("width: 700px; max-width: 90vw;"),
        ):
            with ui.card_actions().classes("justify-end"):
                ui.button("Close", on_click=dialog.close).props("color=primary")

            with ui.card_section():
                ui.label("Citation").classes("text-h6 q-mb-md")
                ui.label("If you use HBAT 2 in your research, please cite:").classes(
                    "text-body2 q-mb-md text-bold"
                )

                # HBAT 2 citation
                ui.label(
                    "Tiwari, A. (2026). HBAT 2: A Python Package to analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures. arXiv. https://doi.org/10.48550/arXiv.2602.17712"
                ).classes("text-body2 q-mb-sm")

                with ui.expansion("BibTeX", icon="code").classes("w-full q-mb-md"):
                    ui.markdown(
                        """
```bibtex
@article{tiwari_2026_hbat_arxiv,
  author       = {Tiwari, Abhishek},
  title        = {HBAT 2: A Python Package to analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures},
  year         = 2026,
  publisher    = {arXiv},
  doi          = {10.48550/arXiv.2602.17712},
  url          = {https://arxiv.org/abs/2602.17712}, 
}
```
                    """
                    )

                # HBAT 2 citation
                ui.label(
                    "Tiwari, A. (2026). HBAT 2: A Python Package to analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures. ChemRxiv. https://chemrxiv.org/doi/abs/10.26434/chemrxiv.15000141/v1"
                ).classes("text-body2 q-mb-sm")

                with ui.expansion("BibTeX", icon="code").classes("w-full q-mb-md"):
                    ui.markdown(
                        """
```bibtex
@article{tiwari_2026_hbat_chemrxiv,
  author = {Abhishek Tiwari },
  title = {HBAT 2: A Python Package to Analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures},
  publisher = {ChemRxiv},
  year = {2026},
  doi = {10.26434/chemrxiv.15000141/v1},
  URL = {https://chemrxiv.org/doi/abs/10.26434/chemrxiv.15000141/v1},
  eprint = {https://chemrxiv.org/doi/pdf/10.26434/chemrxiv.15000141/v1},
}
```
                    """
                    )

                ui.label("If you use HBAT 1.0 or 1.1 in your research, please cite:").classes(
                    "text-body2 q-mb-md text-bold"
                )
                # Original HBAT citation
                ui.label(
                    "Tiwari, A., & Panigrahi, S. K. (2007). HBAT: A Complete Package for Analysing Strong and Weak Hydrogen Bonds in Macromolecular Crystal Structures. In Silico Biology, 7(6). https://doi.org/10.3233/ISI-2007-00337"
                ).classes("text-body2 q-mb-sm")

                with ui.expansion("BibTeX", icon="code").classes("w-full"):
                    ui.markdown(
                        """
```bibtex
@article{tiwari2007hbat,
   author  = {Tiwari, Abhishek and Panigrahi, Sunil Kumar},
   doi     = {10.3233/ISI-2007-00337},
   journal = {In Silico Biology},
   month   = dec,
   number  = {6},
   title   = {{HBAT: A Complete Package for Analysing Strong and Weak Hydrogen Bonds in Macromolecular Crystal Structures}},
   volume  = {7},
   year    = {2007}
}
```
                    """
                    )

        dialog.open()

    def create_ui(self):
        """Create the main user interface."""
        # Header with menu button
        with ui.header().classes("items-center"):
            ui.button(icon="menu", on_click=lambda: left_drawer.toggle()).props(
                "flat color=white"
            )
            ui.image("/static/hbat.svg").classes("w-10 h-10 q-ml-md")
            ui.label("HBAT 2 - Web Server").classes("text-h5 q-ml-sm")
            ui.space()
            ui.button(
                "Cite", icon="format_quote", on_click=self._show_citation_dialog
            ).props("flat color=yellow")
            with ui.link(
                target="http://hbat.abhishek-tiwari.com/web", new_tab=True
            ).classes("text-white no-underline"):
                ui.button("Docs", icon="description").props("flat color=white")
            with ui.link(
                target="https://github.com/abhishektiwari/hbat", new_tab=True
            ).classes("text-white no-underline"):
                ui.button("GitHub", icon="code").props("flat color=white")

        # Left drawer (set fixed=True to prevent JavaScript state queries during heavy operations)
        with ui.left_drawer(fixed=True).props("bordered") as left_drawer:
            with ui.scroll_area().classes("w-full h-full"):
                ui.label("Navigation").classes("text-h6 q-pa-md")
                ui.separator()

                with ui.list().props("dense"):
                    # Step 1: Always enabled
                    with (
                        ui.item()
                        .props("clickable")
                        .on("click", lambda: self.stepper.set_value("upload"))
                    ):
                        with ui.item_section().props("avatar"):
                            ui.icon("upload_file", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 1: Upload PDB")
                            ui.item_label("Upload or download PDB file").props(
                                "caption"
                            )

                    # Step 2: Enabled when file is uploaded
                    self.nav_configure = (
                        ui.item()
                        .props("clickable")
                        .on("click", lambda: self._navigate_to_step("configure"))
                    )
                    with self.nav_configure:
                        with ui.item_section().props("avatar"):
                            ui.icon("settings", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 2: Configure")
                            ui.item_label("Set analysis parameters").props("caption")
                    self.nav_configure.bind_enabled_from(
                        self, "current_file", lambda x: x is not None
                    )

                    # Step 3: Enabled when analysis is complete
                    self.nav_results = (
                        ui.item()
                        .props("clickable")
                        .on("click", lambda: self._navigate_to_step("results"))
                    )
                    with self.nav_results:
                        with ui.item_section().props("avatar"):
                            ui.icon("analytics", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 3: Results")
                            ui.item_label("View analysis results").props("caption")
                    self.nav_results.bind_enabled_from(
                        self, "analyzer", lambda x: x is not None
                    )

                    # Step 4: Enabled when analysis is complete
                    self.nav_export = (
                        ui.item()
                        .props("clickable")
                        .on("click", lambda: self._navigate_to_step("export"))
                    )
                    with self.nav_export:
                        with ui.item_section().props("avatar"):
                            ui.icon("download", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 4: Export")
                            ui.item_label("Download results").props("caption")
                    self.nav_export.bind_enabled_from(
                        self, "analyzer", lambda x: x is not None
                    )

                ui.separator()

                ui.label("Quick Actions").classes("text-subtitle2 q-pa-md")
                with ui.list().props("dense"):
                    with ui.item().props("clickable").on("click", self._start_over):
                        with ui.item_section().props("avatar"):
                            ui.icon("refresh", color="warning")
                        with ui.item_section():
                            ui.item_label("Start Over")
                            ui.item_label("Reset all analysis and settings").props(
                                "caption"
                            )

                    with (
                        ui.item()
                        .props("clickable")
                        .on("click", self._show_citation_dialog)
                    ):
                        with ui.item_section().props("avatar"):
                            ui.icon("format_quote", color="grey")
                        with ui.item_section():
                            ui.item_label("Cite HBAT")

        # Footer
        with ui.footer(fixed=False).classes("items-center"):
            with ui.row().classes("w-full items-center justify-between q-px-md"):
                ui.label(f"© {datetime.now().year} Abhishek Tiwari ").classes(
                    "text-caption"
                )

        # Right drawer for parameter editing
        with ui.right_drawer(fixed=False).props("bordered width=400") as param_drawer:
            param_drawer.hide()
            with ui.scroll_area().classes("w-full h-full"):
                param_drawer_content = ui.column().classes("w-full q-pa-md")

        # Main content area with stepper
        with ui.column().classes("w-full q-pa-md"):
            with ui.stepper().classes("w-full") as self.stepper:
                # Step 1: Upload PDB File
                with ui.step("upload", title="Upload PDB File", icon="upload_file"):
                    self.upload_panel = UploadPanel(on_file_upload=self._on_file_upload)
                    self.upload_panel.create_ui()

                    with ui.stepper_navigation():
                        ui.button(
                            "Next",
                            icon="arrow_forward",
                            on_click=self._handle_upload_next,
                        ).props("primary")

                # Step 2: Configure Parameters and Run Analysis
                with ui.step(
                    "configure", title="Configure Parameters & Run", icon="settings"
                ):
                    self.parameter_panel = ParameterPanel(
                        param_drawer, param_drawer_content
                    )
                    self.parameter_panel.create_ui()

                    ui.separator().classes("q-my-md")

                    with ui.column().classes("items-center w-full"):
                        self.analyze_button = ui.button(
                            "Analyze",
                            on_click=self._run_analysis,
                            icon="play_arrow",
                        ).props("color=primary size=lg")
                        self.analyze_button.bind_enabled_from(
                            self, "current_file", lambda x: x is not None
                        )
                        self.status_label = ui.label("Ready").classes(
                            "text-caption q-mt-sm text-grey"
                        )

                    with ui.stepper_navigation():
                        ui.button(
                            "Back",
                            icon="arrow_back",
                            on_click=lambda: self.stepper.previous(),
                        ).props("primary")
                        ui.button(
                            "Next",
                            icon="arrow_forward",
                            on_click=lambda: self.stepper.next(),
                        ).props("primary").bind_enabled_from(
                            self, "analyzer", lambda x: x is not None
                        )

                # Step 3: View Results
                with ui.step("results", title="View Results", icon="analytics"):
                    # Container for dynamically created results
                    results_container = ui.column().classes("w-full")
                    self.results_panel = WebResultsPanel(
                        results_container,
                        session_dir_callback=lambda: self.session_dir
                    )

                    with ui.stepper_navigation():
                        ui.button(
                            "Back",
                            icon="arrow_back",
                            on_click=lambda: self.stepper.previous(),
                        ).props("primary")
                        ui.button(
                            "Next",
                            icon="arrow_forward",
                            on_click=lambda: self.stepper.next(),
                        ).props("primary")

                # Step 4: Export Results
                with ui.step("export", title="Export Results", icon="download"):
                    with ui.card().classes("w-full q-pa-md mt-5"):
                        ui.label("Choose export format:").classes("text-h6 q-mb-md")

                        with ui.row().classes("gap-4"):
                            ui.button(
                                "Export JSON",
                                icon="code",
                                on_click=self._export_json,
                            ).props("color=primary")
                            ui.button(
                                "Export CSV",
                                icon="table_chart",
                                on_click=self._export_csv,
                            ).props("color=primary")
                            ui.button(
                                "Export TXT",
                                icon="description",
                                on_click=self._export_txt,
                            ).props("color=primary")
                        # Download Structure section
                        ui.separator().classes("q-my-md")
                        with ui.row().classes("gap-md"):
                            self.download_structure_button = ui.button(
                                "Download Original Structure",
                                icon="download",
                                on_click=self._export_fixed_structure,
                            ).props("color=primary")

                    with ui.stepper_navigation():
                        ui.button(
                            "Back",
                            icon="arrow_back",
                            on_click=lambda: self.stepper.previous(),
                        ).props("primary")
                        ui.button(
                            "Start Over", icon="refresh", on_click=self._start_over
                        ).props("flat color=primary")

            # Info card below stepper (responsive: stacks on mobile, horizontal on desktop)
            with ui.card().classes("w-full q-mt-md q-pa-md"):
                with ui.row().classes(
                    "w-full gap-4 flex-wrap items-center md:flex-nowrap md:items-start"
                ):
                    # Image section - 100% width on mobile, fixed size on desktop
                    with ui.column().classes("w-full md:w-auto flex-shrink-0"):
                        ui.image(
                            "https://static.abhishek-tiwari.com/hbat/B_3WH_501_to_B_MET_146_xbond.png"
                        ).classes("w-full md:w-48 md:h-48 rounded")

                    # Text section - full width on mobile, flexible on desktop
                    with ui.column().classes(
                        "w-full md:flex-1 items-center md:items-start"
                    ):
                        ui.label(
                            "HBAT 2: Analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures"
                        ).classes("text-h6 q-mb-sm text-center md:text-left")
                        ui.label(
                            "HBAT 2 is open-source software licensed under MIT. HBAT 2 provides multiple interfaces including web server, desktop GUI, command-line (CLI), and ready-to-use Jupyter/Colab notebooks. The web server version is freely accessible to all users, including for commercial use."
                        ).classes("text-body2 text-grey-7 text-center md:text-left")
                        with ui.row().classes(
                            "gap-2 q-mt-sm flex-wrap justify-center md:justify-start"
                        ):
                            ui.button(
                                "Cite",
                                icon="format_quote",
                                on_click=self._show_citation_dialog,
                            ).props("flat color=primary")
                            with ui.link(
                                target="http://hbat.abhishek-tiwari.com", new_tab=True
                            ).classes("no-underline"):
                                ui.button("Docs", icon="open_in_new").props(
                                    "flat color=primary"
                                )
                            with ui.link(
                                target="https://github.com/abhishektiwari/hbat",
                                new_tab=True,
                            ).classes("no-underline"):
                                ui.button("GitHub", icon="code").props(
                                    "flat color=primary"
                                )

            # Badges card below info card
            with ui.card().classes("w-full q-mt-md q-pa-md"):
                ui.html(
                    """
                    <div style="display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin-bottom: 8px;">
                        <a href="https://github.com/abhishektiwari/hbat/releases" target="_blank">
                            <img src="https://img.shields.io/github/v/release/abhishektiwari/hbat" alt="GitHub Release">
                        </a>
                        <a href="https://github.com/abhishektiwari/hbat/actions/workflows/test.yml" target="_blank">
                            <img src="https://img.shields.io/github/actions/workflow/status/abhishektiwari/hbat/test.yml?label=tests" alt="Tests">
                        </a>
                        <img src="https://img.shields.io/pypi/v/hbat" alt="PyPI Version">
                        <img src="https://img.shields.io/pypi/wheel/hbat" alt="PyPI Wheel">
                        <img src="https://img.shields.io/pypi/pyversions/hbat" alt="Python Versions">
                        <img src="https://img.shields.io/pypi/status/hbat" alt="PyPI Status">
                        <img src="https://img.shields.io/conda/v/hbat/hbat" alt="Conda Version">
                        <img src="https://img.shields.io/github/license/abhishektiwari/hbat" alt="License">
                        <img src="https://img.shields.io/github/last-commit/abhishektiwari/hbat/main" alt="Last Commit">
                        <a href="https://github.com/abhishektiwari/hbat/releases" target="_blank">
                            <img src="https://img.shields.io/github/downloads/abhishektiwari/hbat/total?label=GitHub%20Downloads" alt="GitHub Downloads">
                        </a>
                        <a href="https://sourceforge.net/projects/hbat/files/" target="_blank">
                            <img src="https://img.shields.io/sourceforge/dt/hbat?label=SourceForge%20Downloads" alt="SourceForge Downloads">
                        </a>
                        <a href="https://pypi.org/project/hbat/" target="_blank">
                            <img src="https://img.shields.io/pepy/dt/hbat?label=PyPI%20Downloads" alt="PyPI Downloads">
                        </a>
                        <a href="https://codecov.io/gh/abhishektiwari/hbat" target="_blank">
                            <img src="https://codecov.io/gh/abhishektiwari/hbat/graph/badge.svg?token=QSKYLB3M1V" alt="Code Coverage">
                        </a>
                        <a href="https://www.codefactor.io/repository/github/abhishektiwari/hbat/overview/main" target="_blank">
                            <img src="https://www.codefactor.io/repository/github/abhishektiwari/hbat/badge/main" alt="CodeFactor">
                        </a>
                        <a href="https://socket.dev/pypi/package/hbat" target="_blank">
                            <img src="https://socket.dev/api/badge/pypi/package/hbat/2.2.11?artifact_id=py3-none-any-whl" alt="Socket Security">
                        </a>
                        <a href="https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Mb7eYKYAAAAJ&citation_for_view=Mb7eYKYAAAAJ:u-x6o8ySG0sC" target="_blank">
                            <img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.3233%2FISI-2007-00337" alt="Citations HBAT 1.0 and 1.1">
                        </a>

                        <a href="https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Mb7eYKYAAAAJ&citation_for_view=Mb7eYKYAAAAJ:3bvyWxjaHKcC" target="_blank">
                            <img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.26434%2Fchemrxiv.15000141%2Fv1&link=https%3A%2F%2Fscholar.google.com%2Fcitations%3Fview_op%3Dview_citation%26hl%3Den%26user%3DMb7eYKYAAAAJ%26citation_for_view%3DMb7eYKYAAAAJ%3A3bvyWxjaHKcC" alt="Citations HBAT 2">
                        </a>
                        
                    </div>
                    <div style="display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin-bottom: 8px;">
                        <span class="__dimensions_badge_embed__" data-doi="10.3233/isi-2007-00337" data-legend="always" data-style="small_circle"></span>
                    </div>
                """,
                    sanitize=False,
                )

    async def _handle_upload_next(self):
        """Handle Next button in upload step - download PDB if ID provided."""
        if self.upload_panel:
            # Check if file already uploaded
            if self.upload_panel.file_uploaded:
                self.stepper.next()
                return

            # Try to download if PDB ID is provided
            success = await self.upload_panel.download_if_pdb_id_provided()
            if success:
                self.stepper.next()
            else:
                ui.notify(
                    "Please upload a PDB file or enter a PDB ID to download",
                    type="warning",
                    position="top-left",
                )
        else:
            self.stepper.next()

    def _on_file_upload(self, filename: str, content: bytes):
        """Handle file upload and save to session directory."""
        self.current_file = filename
        self.pdb_content = content.decode("utf-8")

        # Create a new session if one doesn't exist
        if not self.session_id:
            self.session_id = session_manager.create_session()
            self.session_dir = session_manager.get_session_dir(self.session_id)

        # Save file to session directory
        self.current_file_path = session_manager.get_session_file_path(
            self.session_id, filename
        )
        with open(self.current_file_path, "w") as f:
            f.write(self.pdb_content)

        # ui.notify(f"File loaded: {filename}", type="positive", position="top-left")

        # Update status label with filename
        self._update_status_label("Ready")

        # Allow user to proceed to next step
        self.stepper.next()

    async def _run_analysis(self):
        """Run the molecular interaction analysis."""
        if self.analysis_running:
            ui.notify("Analysis already running", type="warning", position="top-left")
            return

        if not self.current_file_path or not self.current_file_path.exists():
            ui.notify(
                "Please upload a PDB file first", type="warning", position="top-left"
            )
            return

        self.analysis_running = True
        self.analyze_button.disable()
        self.status_label.text = "Running analysis..."
        self.status_label.classes(replace="text-caption q-mt-sm text-warning text-bold")

        # Create notification for progress feedback
        notification = ui.notification(timeout=None)
        notification.message = "Analyzing molecular interactions..."
        notification.spinner = True

        try:
            # Get parameters from UI
            params = self.parameter_panel.get_parameters()

            # Create analyzer and run analysis
            self.analyzer = NPMolecularInteractionAnalyzer(parameters=params)

            # Run analysis in executor to avoid blocking
            loop = asyncio.get_event_loop()
            success = await loop.run_in_executor(
                None, self.analyzer.analyze_file, str(self.current_file_path)
            )

            if success:
                self.status_label.text = (
                    f"Analysis completed! File: {self.current_file}"
                )
                self.status_label.classes(replace="text-caption q-mt-sm text-positive")

                notification.message = "Analysis completed!"
                notification.spinner = False
                await asyncio.sleep(1)
                notification.dismiss()

                # Update download button label based on whether fixing was applied
                if (
                    hasattr(self.analyzer, "_pdb_fixing_info")
                    and self.analyzer._pdb_fixing_info.get("applied")
                ):
                    button_label = "Download Fixed Structure"
                else:
                    button_label = "Download Original Structure"

                if self.download_structure_button:
                    self.download_structure_button.text = button_label

                # Update results display
                await self.results_panel.update_results(
                    self.analyzer, self.current_file
                )

                # Auto-advance to results step
                self.stepper.next()
            else:
                self.status_label.text = "Analysis failed"
                self.status_label.classes(replace="text-caption q-mt-sm text-negative")

                notification.message = "Analysis failed. Please check the PDB file."
                notification.spinner = False

        except Exception as e:
            self.status_label.text = f"Error: {str(e)}"
            self.status_label.classes(replace="text-caption q-mt-sm text-negative")

            notification.message = f"Error: {str(e)}"
            notification.spinner = False

        finally:
            self.analysis_running = False
            self.analyze_button.enable()

    def _export_json(self):
        """Export results as JSON."""
        if not self.analyzer:
            ui.notify(
                "No analysis results to export", type="warning", position="top-left"
            )
            return

        base_name = Path(self.current_file).stem
        output_file = session_manager.get_session_file_path(
            self.session_id, f"{base_name}_results.json"
        )
        export_to_json_single_file(
            self.analyzer, str(output_file), input_file=self.current_file
        )
        ui.download(str(output_file))
        ui.notify(
            f"Exported to {output_file.name}", type="positive", position="top-left"
        )

    def _export_csv(self):
        """Export results as CSV (creates a zip file with all CSV files and fixed PDB)."""
        if not self.analyzer:
            ui.notify(
                "No analysis results to export", type="warning", position="top-left"
            )
            return

        # Remove extension from current_file to get base name
        base_name = Path(self.current_file).stem
        base_filename = self.session_dir / base_name
        export_to_csv_files(self.analyzer, str(base_filename))

        # Get all generated CSV files
        csv_files = list(self.session_dir.glob(f"{base_name}_*.csv"))

        # Check for fixed PDB file
        fixed_pdb_path = None
        if hasattr(
            self.analyzer, "_pdb_fixing_info"
        ) and self.analyzer._pdb_fixing_info.get("applied"):
            fixed_file_path = self.analyzer._pdb_fixing_info.get("fixed_file_path")
            if fixed_file_path and os.path.exists(fixed_file_path):
                # Copy fixed PDB to session directory
                method = self.analyzer._pdb_fixing_info.get("method", "unknown")
                fixed_pdb_path = self.session_dir / f"{base_name}_fixed_{method}.pdb"
                with open(fixed_file_path, "r") as src:
                    with open(fixed_pdb_path, "w") as dst:
                        dst.write(src.read())

        if csv_files:
            # Create a zip file containing all CSV files and fixed PDB (if available)
            zip_path = self.session_dir / f"{base_name}_csv_export.zip"
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
                for csv_file in csv_files:
                    zipf.write(csv_file, csv_file.name)

                # Add fixed PDB file if available
                if fixed_pdb_path and fixed_pdb_path.exists():
                    zipf.write(fixed_pdb_path, fixed_pdb_path.name)

            # Download the zip file
            ui.download(str(zip_path))

            # Build notification message
            files_count = len(csv_files)
            if fixed_pdb_path and fixed_pdb_path.exists():
                files_count += 1
                ui.notify(
                    f"Exported {len(csv_files)} CSV file(s) + fixed PDB as ZIP",
                    type="positive",
                    position="top-left",
                )
            else:
                ui.notify(
                    f"Exported {len(csv_files)} CSV file(s) as ZIP",
                    type="positive",
                    position="top-left",
                )

            # Clean up individual CSV files (optional - keep them for now)
            # for csv_file in csv_files:
            #     csv_file.unlink()
        else:
            ui.notify(
                "No CSV files were generated", type="warning", position="top-left"
            )

    def _export_txt(self):
        """Export results as TXT."""
        if not self.analyzer:
            ui.notify(
                "No analysis results to export", type="warning", position="top-left"
            )
            return

        base_name = Path(self.current_file).stem
        output_file = session_manager.get_session_file_path(
            self.session_id, f"{base_name}_results.txt"
        )
        export_to_txt_single_file(self.analyzer, str(output_file))
        ui.download(str(output_file))
        ui.notify(
            f"Exported to {output_file.name}", type="positive", position="top-left"
        )

    def _export_fixed_structure(self):
        """Download the fixed structure (or original if no fixing was applied)."""
        if not self.analyzer:
            ui.notify(
                "No analysis results to export", type="warning", position="top-left"
            )
            return

        # Check if PDB fixing was applied
        if hasattr(
            self.analyzer, "_pdb_fixing_info"
        ) and self.analyzer._pdb_fixing_info.get("applied"):
            # Download fixed structure
            fixed_file_path = self.analyzer._pdb_fixing_info.get("fixed_file_path")
            file_type = "fixed structure"
        else:
            # No fixing applied, download original structure
            if not hasattr(
                self.analyzer, "_pdb_fixing_info"
            ) or not self.analyzer._pdb_original_info.get("input_file_path"):
                ui.notify(
                    "Original structure file not found",
                    type="negative",
                    position="top-left",
                )
                return
            fixed_file_path = self.analyzer._pdb_original_info.get("input_file_path")
            file_type = "original structure"

        try:
            if not os.path.exists(fixed_file_path):
                ui.notify(
                    f"{file_type.capitalize()} file not found",
                    type="negative",
                    position="top-left",
                )
                return

            # Download the file directly
            ui.download(fixed_file_path)
            ui.notify(
                f"Downloaded {file_type}",
                type="positive",
                position="top-left",
            )

        except Exception as e:
            ui.notify(
                f"Download failed: {str(e)}", type="negative", position="top-left"
            )


    def _navigate_to_step(self, step: str):
        """Navigate to a specific step with validation.

        :param step: Step name to navigate to
        :type step: str
        """
        # Validate navigation based on step requirements
        if step == "configure" and self.current_file is None:
            ui.notify(
                "Please upload a PDB file first", type="warning", position="top-left"
            )
            return

        if step in ("results", "export") and self.analyzer is None:
            ui.notify("Please run analysis first", type="warning", position="top-left")
            return

        # Navigate to the requested step
        if self.stepper:
            self.stepper.set_value(step)

    def _start_over(self):
        """Reset the application to initial state."""
        # Clean up old session if it exists
        if self.session_id:
            session_manager.delete_session(self.session_id)

        # Reset analyzer and file state
        self.analyzer = None
        self.current_file = None
        self.current_file_path = None
        self.analysis_running = False
        self.pdb_content = None
        self.session_id = None
        self.session_dir = None

        # Reset all panels
        if self.upload_panel:
            self.upload_panel.reset()
        if self.parameter_panel:
            self.parameter_panel.reset()
        if self.results_panel:
            self.results_panel.reset()

        # Reset status label
        self._update_status_label("Ready")
        if self.status_label:
            self.status_label.classes(replace="text-caption text-grey")

        # Note: Analyze button and navigation items are automatically
        # disabled via bindings when current_file and analyzer are set to None

        # Go back to upload step
        if self.stepper:
            self.stepper.set_value("upload")

        ui.notify(
            "Application reset to initial state", type="positive", position="top-left"
        )


def create_app():
    """Create and configure the HBAT web application.

    :returns: None
    """
    # Cleanup expired sessions on startup
    cleaned = session_manager.cleanup_expired_sessions()
    if cleaned > 0:
        print(f"Cleaned up {cleaned} expired session(s) on startup")

    # Schedule periodic cleanup (every 6 hours)
    async def periodic_cleanup():
        """Periodically clean up expired sessions."""
        while True:
            await asyncio.sleep(6 * 3600)  # 6 hours
            cleaned = session_manager.cleanup_expired_sessions()
            if cleaned > 0:
                print(f"Periodic cleanup: removed {cleaned} expired session(s)")

    # Start background cleanup task
    app.on_startup(lambda: asyncio.create_task(periodic_cleanup()))

    # Configure static files BEFORE page routes
    static_dir = Path(__file__).parent / "static"
    if static_dir.exists():
        app.add_static_files("/static", str(static_dir))

    @ui.page("/")
    def index():
        """Main page route."""
        # Configure Quasar color theme
        ui.colors(
            primary="#20c997",
            secondary="#6c757d",
            accent="#9C27B0",
            dark="#1d1d1d",
            positive="#198754",
            negative="#C10015",
            info="#31CCEC",
            warning="#ffc107",
        )

        # Load 3Dmol library for this page
        ui.add_head_html(
            '<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>'
        )

        # Add meta tags for SEO and social sharing
        ui.add_head_html(
            """
            <meta name="description" content="HBAT 2 (Hydrogen Bond Analysis Tool 2) for analyzing molecular interactions including hydrogen bonds, halogen bonds, π interactions, and cooperativity chains in protein structures">
            <meta name="keywords" content="hydrogen bond analysis, hydrogen bond calculator, protein-ligand interactions, molecular interactions, protein structure, halogen bonds, pi interactions, cooperativity chains, structural biology, bioinformatics, PDB analysis">
            <meta name="author" content="Abhishek Tiwari">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">

            <!-- Open Graph / Facebook -->
            <meta property="og:type" content="website">
            <meta property="og:url" content="https://hbat.abhishek-tiwari.com/">
            <meta property="og:title" content="HBAT - Hydrogen Bond Analysis Tool">
            <meta property="og:description" content="Comprehensive tool for analyzing molecular interactions in protein structures including hydrogen bonds, halogen bonds, π interactions, and cooperativity chains">
            <meta property="og:image" content="/static/hbat.svg">

            <!-- Twitter -->
            <meta property="twitter:card" content="summary_large_image">
            <meta property="twitter:url" content="https://hbat.abhishek-tiwari.com/">
            <meta property="twitter:title" content="HBAT 2 (Hydrogen Bond Analysis Tool 2)">
            <meta property="twitter:description" content="Comprehensive tool for analyzing molecular interactions in protein structures including hydrogen bonds, halogen bonds, π interactions, and cooperativity chains">
            <meta property="twitter:image" content="/static/hbat.svg">

            <!-- Custom Styling -->
            <style>
                /* Custom styling for cards */
                .q-card {
                    border-radius: 8px;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                }

                /* Smooth transitions */
                .q-btn, .q-stepper__tab {
                    transition: all 0.3s ease;
                }

                /* Enhance button appearance */
                .q-btn {
                    text-transform: none;
                }

                /* Notifications at top */
                .q-notifications__list--top {
                    top: 70px;
                }

                /* Parameter drawer input fields - ensure full width and labels don't wrap */
                .q-drawer--right .q-field {
                    width: 50%;
                }

                .q-drawer--right .q-field__label {
                    white-space: normal;
                    word-wrap: break-word;
                    max-width: 100%;
                    font-weight: bold;
                }
            </style>
            <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
        """
        )

        hbat_app = HBATWebApp()
        hbat_app.create_ui()

    # Check if running in production/Docker environment
    is_production = os.getenv("HBAT_ENV") == "production"
    # Bind to 0.0.0.0 in production/Docker (intentional for containerized deployments)
    host = os.getenv("HBAT_HOST", "0.0.0.0" if is_production else "127.0.0.1")  # nosec B104
    port = int(os.getenv("HBAT_PORT", "8080"))

    # Auto-reload disabled by default to preserve analysis state during active use.
    # Only enable if explicitly requested via HBAT_RELOAD=true environment variable.
    # This prevents loss of analysis results when source files change during development.
    reload = os.getenv("HBAT_RELOAD", "false").lower() == "true"

    # Set favicon path
    favicon_path = static_dir / "favicon.ico" if static_dir.exists() else None

    ui.run(
        title="HBAT 2 - Web Server",
        favicon=str(favicon_path) if favicon_path else "🧬",
        dark=False,
        reload=reload,  # Disabled by default to preserve state; enable with HBAT_RELOAD=true
        show=not is_production,  # Don't open browser in production/Docker
        host=host,
        port=port,
        reconnect_timeout=3600.0,  # 60 minutes - increased for long-running analysis
    )


if __name__ in {"__main__", "__mp_main__"}:
    create_app()
