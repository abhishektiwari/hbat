"""
Main web application for HBAT using NiceGUI.

This module provides a web-based interface for molecular interaction analysis
with integrated 3D visualization using py3Dmol.
"""

import asyncio
from datetime import datetime
import os
import tempfile
import zipfile
from pathlib import Path
from typing import Optional

from nicegui import app, context, ui

from ..constants import APP_NAME
from ..core.analysis import AnalysisParameters, NPMolecularInteractionAnalyzer
from ..export import export_to_csv_files, export_to_json_single_file, export_to_txt_single_file
from .components.parameter_panel import ParameterPanel
from .components.results_panel import WebResultsPanel
from .components.upload_panel import UploadPanel

# Uploads directory
UPLOADS_DIR = Path("uploads")
UPLOADS_DIR.mkdir(exist_ok=True)


class HBATWebApp:
    """Main HBAT web application class."""

    def __init__(self):
        """Initialize the HBAT web application."""
        self.analyzer: Optional[NPMolecularInteractionAnalyzer] = None
        self.current_file: Optional[str] = None
        self.current_file_path: Optional[Path] = None
        self.analysis_running = False
        self.pdb_content: Optional[str] = None

        # UI components
        self.upload_panel: Optional[UploadPanel] = None
        self.parameter_panel: Optional[ParameterPanel] = None
        self.results_panel: Optional[WebResultsPanel] = None

        # Status tracking
        self.status_label: Optional[ui.label] = None
        self.analyze_button: Optional[ui.button] = None

        # Stepper
        self.stepper: Optional[ui.stepper] = None

        # Navigation items
        self.nav_configure: Optional[ui.item] = None
        self.nav_results: Optional[ui.item] = None
        self.nav_export: Optional[ui.item] = None

    def _show_citation_dialog(self):
        """Show citation dialog."""
        with ui.dialog().props("persistent") as dialog, ui.card().style("width: 700px; max-width: 90vw;"):
            with ui.card_actions().classes("justify-end"):
                ui.button("Close", on_click=dialog.close).props("color=primary")

            with ui.card_section():
                ui.label("Citation").classes("text-h6 q-mb-md")
                ui.label("If you use HBAT in your research, please cite:").classes("text-body2 q-mb-md")

                # HBAT 2 citation
                ui.label("Tiwari, A. (2025). HBAT 2: A Python Package to analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures. Zenodo. https://doi.org/10.5281/zenodo.17645377").classes("text-body2 q-mb-sm")

                with ui.expansion("BibTeX", icon="code").classes("w-full q-mb-md"):
                    ui.markdown("""
```bibtex
@misc{tiwari_2025_17645321,
   author    = {Tiwari, Abhishek},
   title     = {HBAT 2: A Python Package to analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures},
   month     = nov,
   year      = 2025,
   publisher = {Zenodo},
   doi       = {10.5281/zenodo.17645321},
   url       = {https://doi.org/10.5281/zenodo.17645321},
}
```
                    """)

                # Original HBAT citation
                ui.label("Tiwari, A., & Panigrahi, S. K. (2007). HBAT: A Complete Package for Analysing Strong and Weak Hydrogen Bonds in Macromolecular Crystal Structures. In Silico Biology, 7(6). https://doi.org/10.3233/ISI-2007-00337").classes("text-body2 q-mb-sm")

                with ui.expansion("BibTeX", icon="code").classes("w-full"):
                    ui.markdown("""
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
                    """)

        dialog.open()

    def create_ui(self):
        """Create the main user interface."""
        # Header with menu button
        with ui.header().classes("items-center"):
            ui.button(icon="menu", on_click=lambda: left_drawer.toggle()).props("flat color=white")
            ui.image("/static/hbat.svg").classes("w-10 h-10 q-ml-md")
            ui.label(f"HBAT 2 - Web Server").classes("text-h5 q-ml-sm")
            ui.space()
            ui.button("Cite", icon="format_quote", on_click=self._show_citation_dialog).props("flat color=yellow")
            with ui.link(target="http://hbat.abhishek-tiwari.com", new_tab=True).classes("text-white no-underline"):
                ui.button("Docs", icon="description").props("flat color=white")
            with ui.link(target="https://github.com/abhishektiwari/hbat", new_tab=True).classes("text-white no-underline"):
                ui.button("Code", icon="code").props("flat color=white")

        # Left drawer (set fixed=True to prevent JavaScript state queries during heavy operations)
        with ui.left_drawer(fixed=True).props("bordered") as left_drawer:
            with ui.scroll_area().classes("w-full h-full"):
                ui.label("Navigation").classes("text-h6 q-pa-md")
                ui.separator()

                with ui.list().props("dense"):
                    # Step 1: Always enabled
                    with ui.item().props("clickable").on("click", lambda: self.stepper.set_value("upload")):
                        with ui.item_section().props("avatar"):
                            ui.icon("upload_file", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 1: Upload PDB")
                            ui.item_label("Upload or download PDB file").props("caption")

                    # Step 2: Enabled when file is uploaded
                    self.nav_configure = ui.item().props("clickable").on("click", lambda: self._navigate_to_step("configure"))
                    with self.nav_configure:
                        with ui.item_section().props("avatar"):
                            ui.icon("settings", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 2: Configure")
                            ui.item_label("Set analysis parameters").props("caption")
                    self.nav_configure.bind_enabled_from(self, "current_file", lambda x: x is not None)

                    # Step 3: Enabled when analysis is complete
                    self.nav_results = ui.item().props("clickable").on("click", lambda: self._navigate_to_step("results"))
                    with self.nav_results:
                        with ui.item_section().props("avatar"):
                            ui.icon("analytics", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 3: Results")
                            ui.item_label("View analysis results").props("caption")
                    self.nav_results.bind_enabled_from(self, "analyzer", lambda x: x is not None)

                    # Step 4: Enabled when analysis is complete
                    self.nav_export = ui.item().props("clickable").on("click", lambda: self._navigate_to_step("export"))
                    with self.nav_export:
                        with ui.item_section().props("avatar"):
                            ui.icon("download", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 4: Export")
                            ui.item_label("Download results").props("caption")
                    self.nav_export.bind_enabled_from(self, "analyzer", lambda x: x is not None)

                ui.separator()

                ui.label("Quick Actions").classes("text-subtitle2 q-pa-md")
                with ui.list().props("dense"):
                    with ui.item().props("clickable").on("click", self._start_over):
                        with ui.item_section().props("avatar"):
                            ui.icon("refresh", color="warning")
                        with ui.item_section():
                            ui.item_label("Start Over")
                            ui.item_label("Reset all analysis and settings").props("caption")

                    with ui.item().props("clickable").on("click", self._show_citation_dialog):
                        with ui.item_section().props("avatar"):
                            ui.icon("format_quote", color="grey")
                        with ui.item_section():
                            ui.item_label("Cite HBAT")

        # Footer
        with ui.footer(fixed=False).classes("items-center"):
            with ui.row().classes("w-full items-center justify-between q-px-md"):
                ui.label(f"© {datetime.now().year} Abhishek Tiwari ").classes("text-caption")

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
                        ui.button("Next", icon="arrow_forward", on_click=self._handle_upload_next).props(
                            "primary"
                        )

                # Step 2: Configure Parameters and Run Analysis
                with ui.step(
                    "configure", title="Configure Parameters & Run", icon="settings"
                ):
                    self.parameter_panel = ParameterPanel(param_drawer, param_drawer_content)
                    self.parameter_panel.create_ui()

                    ui.separator().classes("q-my-md")

                    with ui.column().classes("items-center w-full"):
                        self.analyze_button = ui.button(
                            "Analyze",
                            on_click=self._run_analysis,
                            icon="play_arrow",
                        ).props("color=primary size=lg")
                        self.analyze_button.bind_enabled_from(self, "current_file", lambda x: x is not None)
                        self.status_label = ui.label("Ready").classes("text-caption q-mt-sm text-grey")

                    with ui.stepper_navigation():
                        ui.button("Back", icon="arrow_back", on_click=lambda: self.stepper.previous()).props(
                            "primary"
                        )
                        ui.button("Next", icon="arrow_forward", on_click=lambda: self.stepper.next()).props(
                            "primary"
                        ).bind_enabled_from(self, "analyzer", lambda x: x is not None)

                # Step 3: View Results
                with ui.step("results", title="View Results", icon="analytics"):

                    # Container for dynamically created results
                    results_container = ui.column().classes("w-full")
                    self.results_panel = WebResultsPanel(results_container)

                    with ui.stepper_navigation():
                        ui.button("Back", icon="arrow_back", on_click=lambda: self.stepper.previous()).props(
                            "primary"
                        )
                        ui.button("Next", icon="arrow_forward", on_click=lambda: self.stepper.next()).props(
                            "primary"
                        )

                # Step 4: Export Results
                with ui.step("export", title="Export Results", icon="download"):
                    ui.label("Export analysis results").classes("text-subtitle1 q-mb-md")

                    with ui.card().classes("w-full q-pa-md"):
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

                    with ui.stepper_navigation():
                        ui.button("Back", icon="arrow_back", on_click=lambda: self.stepper.previous()).props(
                            "primary"
                        )
                        ui.button(
                            "Start Over", icon="refresh", on_click=self._start_over
                        ).props("flat color=primary")

            # Info card below stepper
            with ui.card().classes("w-full q-mt-md"):
                with ui.row().classes("w-full items-center gap-4"):
                    # Image section
                    with ui.column().classes("flex-shrink-0"):
                        ui.image("https://static.abhishek-tiwari.com/hbat/B_3WH_501_to_B_MET_146_xbond.png").classes("w-48 h-48 rounded")

                    # Text section
                    with ui.column().classes("flex-1"):
                        ui.label("HBAT 2: Analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures").classes("text-h6 q-mb-sm")
                        ui.label("HBAT 2 is open-source software licensed under MIT. HBAT 2 provides multiple interfaces including web server, desktop GUI, command-line (CLI), and ready-to-use Jupyter/Colab notebooks. The web server version is freely accessible to all users, including for commercial use.").classes("text-body2 text-grey-7")
                        with ui.row().classes("gap-2 q-mt-sm"):
                            with ui.link(target="http://hbat.abhishek-tiwari.com", new_tab=True).classes("no-underline"):
                                ui.button("HBAT Documentation", icon="open_in_new").props("flat color=primary")
                            with ui.link(target="https://github.com/abhishektiwari/hbat", new_tab=True).classes("no-underline"):
                                ui.button("HBAT 2 Code", icon="code").props("flat color=primary")
                            ui.button("Cite HBAT 2", icon="format_quote", on_click=self._show_citation_dialog).props("flat color=primary")

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
                    position="top-left"
                )
        else:
            self.stepper.next()

    def _on_file_upload(self, filename: str, content: bytes):
        """Handle file upload and save to uploads directory."""
        self.current_file = filename
        self.pdb_content = content.decode("utf-8")

        # Save file to uploads directory
        self.current_file_path = UPLOADS_DIR / filename
        with open(self.current_file_path, "w") as f:
            f.write(self.pdb_content)

        ui.notify(f"File loaded: {filename}", type="positive", position="top-left")
        # Allow user to proceed to next step
        self.stepper.next()

    async def _run_analysis(self):
        """Run the molecular interaction analysis."""
        if self.analysis_running:
            ui.notify("Analysis already running", type="warning", position="top-left")
            return

        if not self.current_file_path or not self.current_file_path.exists():
            ui.notify("Please upload a PDB file first", type="warning", position="top-left")
            return

        self.analysis_running = True
        self.analyze_button.disable()
        self.status_label.text = "Running analysis..."
        self.status_label.classes(replace="text-caption q-mt-sm text-warning text-bold")

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
                self.status_label.text = f"Analysis completed! File: uploads/{self.current_file}"
                self.status_label.classes(replace="text-caption q-mt-sm text-positive")
                ui.notify("Analysis completed!", type="positive", position="top-left")

                # Update results display
                await self.results_panel.update_results(
                    self.analyzer, self.pdb_content, self.current_file
                )

                # Auto-advance to results step
                self.stepper.next()
            else:
                self.status_label.text = "Analysis failed"
                self.status_label.classes(replace="text-caption q-mt-sm text-negative")
                ui.notify("Analysis failed. Please check the PDB file.", type="negative", position="top-left")

        except Exception as e:
            self.status_label.text = f"Error: {str(e)}"
            self.status_label.classes(replace="text-caption q-mt-sm text-negative")
            ui.notify(f"Error: {str(e)}", type="negative", position="top-left")

        finally:
            self.analysis_running = False
            self.analyze_button.enable()

    def _export_json(self):
        """Export results as JSON."""
        if not self.analyzer:
            ui.notify("No analysis results to export", type="warning", position="top-left")
            return

        base_name = Path(self.current_file).stem
        output_file = UPLOADS_DIR / f"{base_name}_results.json"
        export_to_json_single_file(
            self.analyzer, str(output_file), input_file=self.current_file
        )
        ui.download(str(output_file))
        ui.notify(f"Exported to {output_file.name}", type="positive", position="top-left")

    def _export_csv(self):
        """Export results as CSV (creates a zip file with all CSV files)."""
        if not self.analyzer:
            ui.notify("No analysis results to export", type="warning", position="top-left")
            return

        # Remove extension from current_file to get base name
        base_name = Path(self.current_file).stem
        base_filename = UPLOADS_DIR / base_name
        export_to_csv_files(self.analyzer, str(base_filename))

        # Get all generated CSV files
        csv_files = list(UPLOADS_DIR.glob(f"{base_name}_*.csv"))
        if csv_files:
            # Create a zip file containing all CSV files
            zip_path = UPLOADS_DIR / f"{base_name}_csv_export.zip"
            with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
                for csv_file in csv_files:
                    zipf.write(csv_file, csv_file.name)

            # Download the zip file
            ui.download(str(zip_path))
            ui.notify(f"Exported {len(csv_files)} CSV file(s) as ZIP", type="positive", position="top-left")

            # Clean up individual CSV files (optional - keep them for now)
            # for csv_file in csv_files:
            #     csv_file.unlink()
        else:
            ui.notify("No CSV files were generated", type="warning", position="top-left")

    def _export_txt(self):
        """Export results as TXT."""
        if not self.analyzer:
            ui.notify("No analysis results to export", type="warning", position="top-left")
            return

        base_name = Path(self.current_file).stem
        output_file = UPLOADS_DIR / f"{base_name}_results.txt"
        export_to_txt_single_file(self.analyzer, str(output_file))
        ui.download(str(output_file))
        ui.notify(f"Exported to {output_file.name}", type="positive", position="top-left")

    def _navigate_to_step(self, step: str):
        """Navigate to a specific step with validation.

        :param step: Step name to navigate to
        :type step: str
        """
        # Validate navigation based on step requirements
        if step == "configure" and self.current_file is None:
            ui.notify("Please upload a PDB file first", type="warning", position="top-left")
            return

        if step in ("results", "export") and self.analyzer is None:
            ui.notify("Please run analysis first", type="warning", position="top-left")
            return

        # Navigate to the requested step
        if self.stepper:
            self.stepper.set_value(step)

    def _start_over(self):
        """Reset the application to initial state."""
        # Reset analyzer and file state
        self.analyzer = None
        self.current_file = None
        self.current_file_path = None
        self.analysis_running = False
        self.pdb_content = None

        # Reset all panels
        if self.upload_panel:
            self.upload_panel.reset()
        if self.parameter_panel:
            self.parameter_panel.reset()
        if self.results_panel:
            self.results_panel.reset()

        # Reset status label
        if self.status_label:
            self.status_label.text = "Ready"
            self.status_label.classes(replace="text-caption text-grey")

        # Note: Analyze button and navigation items are automatically
        # disabled via bindings when current_file and analyzer are set to None

        # Go back to upload step
        if self.stepper:
            self.stepper.set_value("upload")

        ui.notify("Application reset to initial state", type="positive", position="top-left")


def create_app():
    """Create and configure the HBAT web application.

    :returns: None
    """
    # Configure static files BEFORE page routes
    static_dir = Path(__file__).parent / "static"
    if static_dir.exists():
        app.add_static_files("/static", str(static_dir))

    @ui.page("/")
    def index():
        """Main page route."""
        # Configure Quasar color theme
        ui.colors(
            primary='#20c997',
            secondary='#6c757d',
            accent='#9C27B0',
            dark='#1d1d1d',
            positive='#198754',
            negative='#C10015',
            info='#31CCEC',
            warning='#ffc107'
        )

        # Load 3Dmol library for this page
        ui.add_head_html('<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>')

        # Add meta tags for SEO and social sharing
        ui.add_head_html('''
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
            </style>
        ''')

        hbat_app = HBATWebApp()
        hbat_app.create_ui()

    # Check if running in production/Docker environment
    is_production = os.getenv("HBAT_ENV") == "production"
    # Bind to 0.0.0.0 in production/Docker (intentional for containerized deployments)
    host = os.getenv("HBAT_HOST", "0.0.0.0" if is_production else "127.0.0.1")  # nosec B104
    port = int(os.getenv("HBAT_PORT", "8080"))

    # Set favicon path
    favicon_path = static_dir / "favicon.ico" if static_dir.exists() else None

    ui.run(
        title=f"HBAT 2 - Web Server",
        favicon=str(favicon_path) if favicon_path else "🧬",
        dark=False,
        reload=False,
        show=not is_production,  # Don't open browser in production/Docker
        host=host,
        port=port,
        reconnect_timeout=300.0,  # 5 minutes - increased for long-running analysis
    )


if __name__ == "__main__":
    create_app()
