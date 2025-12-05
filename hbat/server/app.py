"""
Main web application for HBAT using NiceGUI.

This module provides a web-based interface for molecular interaction analysis
with integrated 3D visualization using py3Dmol.
"""

import asyncio
import os
import tempfile
from pathlib import Path
from typing import Optional

from nicegui import app, ui

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

    def create_ui(self):
        """Create the main user interface."""
        # Header with menu button
        with ui.header().classes("items-center"):
            ui.button(icon="menu", on_click=lambda: left_drawer.toggle()).props("flat color=white")
            ui.label(f"{APP_NAME}").classes("text-h5 q-ml-md")
            ui.space()
            ui.link("Documentation", "https://hbat.abhishek-tiwari.com/", new_tab=True).classes(
                "text-white"
            )

        # Left drawer
        with ui.left_drawer(fixed=False).props("bordered") as left_drawer:
            with ui.scroll_area().classes("w-full h-full"):
                ui.label("Navigation").classes("text-h6 q-pa-md")
                ui.separator()

                with ui.list().props("dense"):
                    with ui.item().props("clickable").on("click", lambda: self.stepper.set_value("upload")):
                        with ui.item_section().props("avatar"):
                            ui.icon("upload_file", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 1: Upload PDB")
                            ui.item_label("Upload or download PDB file").props("caption")

                    with ui.item().props("clickable").on("click", lambda: self.stepper.set_value("configure")):
                        with ui.item_section().props("avatar"):
                            ui.icon("settings", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 2: Configure")
                            ui.item_label("Set analysis parameters").props("caption")

                    with ui.item().props("clickable").on("click", lambda: self.stepper.set_value("results")):
                        with ui.item_section().props("avatar"):
                            ui.icon("analytics", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 3: Results")
                            ui.item_label("View analysis results").props("caption")

                    with ui.item().props("clickable").on("click", lambda: self.stepper.set_value("export")):
                        with ui.item_section().props("avatar"):
                            ui.icon("download", color="primary")
                        with ui.item_section():
                            ui.item_label("Step 4: Export")
                            ui.item_label("Download results").props("caption")

                ui.separator()

                ui.label("Quick Actions").classes("text-subtitle2 q-pa-md")
                with ui.list().props("dense"):
                    with ui.item().props("clickable").on("click", lambda: ui.open("https://www.rcsb.org/", new_tab=True)):
                        with ui.item_section().props("avatar"):
                            ui.icon("public", color="grey")
                        with ui.item_section():
                            ui.item_label("RCSB PDB")

                    with ui.item().props("clickable").on("click", lambda: ui.open("https://hbat.abhishek-tiwari.com/", new_tab=True)):
                        with ui.item_section().props("avatar"):
                            ui.icon("help", color="grey")
                        with ui.item_section():
                            ui.item_label("Help & Docs")

        # Footer
        with ui.footer(fixed=False).classes("items-center"):
            with ui.row().classes("w-full items-center justify-between q-px-md"):
                ui.label(f"© 2025 {APP_NAME} ").classes("text-caption")
                ui.link("GitHub", "https://github.com/abhishektiwari/hbat", new_tab=True).classes("text-white")

        # Right drawer for parameter editing
        with ui.right_drawer(fixed=False).props("bordered width=400") as param_drawer:
            param_drawer.hide()
            with ui.scroll_area().classes("w-full h-full"):
                param_drawer_content = ui.column().classes("w-full q-pa-md")

        # Main content area with stepper
        with ui.column().classes("w-full q-pa-md"):
            with ui.stepper().props("vertical").classes("w-full") as self.stepper:
                # Step 1: Upload PDB File
                with ui.step("upload", title="Upload PDB File", icon="upload_file"):
                    ui.label("Upload a PDB file or enter a PDB ID to download").classes(
                        "text-subtitle1 q-mb-md"
                    )
                    self.upload_panel = UploadPanel(on_file_upload=self._on_file_upload)
                    self.upload_panel.create_ui()

                    with ui.stepper_navigation():
                        ui.button("Next", on_click=lambda: self.stepper.next()).props(
                            "flat"
                        )

                # Step 2: Configure Parameters and Run Analysis
                with ui.step(
                    "configure", title="Configure Parameters & Run", icon="settings"
                ):
                    ui.label("Configure analysis parameters").classes(
                        "text-subtitle1 q-mb-md"
                    )
                    self.parameter_panel = ParameterPanel(param_drawer, param_drawer_content)
                    self.parameter_panel.create_ui()

                    ui.separator().classes("q-my-md")

                    ui.label("Run Analysis").classes("text-h6 q-mt-md")
                    self.analyze_button = ui.button(
                        "Analyze",
                        on_click=self._run_analysis,
                        icon="play_arrow",
                    ).props("color=primary size=lg")
                    self.status_label = ui.label("").classes("text-caption q-mt-sm")

                    with ui.stepper_navigation():
                        ui.button("Back", on_click=lambda: self.stepper.previous()).props(
                            "flat"
                        )
                        ui.button("Next", on_click=lambda: self.stepper.next()).props(
                            "flat"
                        ).bind_enabled_from(self, "analyzer", lambda x: x is not None)

                # Step 3: View Results
                with ui.step("results", title="View Results", icon="analytics"):
                    ui.label("Analysis Results").classes("text-subtitle1 q-mb-md")

                    # Container for dynamically created results
                    results_container = ui.column().classes("w-full")
                    self.results_panel = WebResultsPanel(results_container)

                    with ui.stepper_navigation():
                        ui.button("Back", on_click=lambda: self.stepper.previous()).props(
                            "flat"
                        )
                        ui.button("Next", on_click=lambda: self.stepper.next()).props(
                            "flat"
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
                        ui.button("Back", on_click=lambda: self.stepper.previous()).props(
                            "flat"
                        )
                        ui.button(
                            "Start Over", on_click=lambda: self.stepper.set_value("upload")
                        ).props("flat color=primary")

    def _on_file_upload(self, filename: str, content: bytes):
        """Handle file upload and save to uploads directory."""
        self.current_file = filename
        self.pdb_content = content.decode("utf-8")

        # Save file to uploads directory
        self.current_file_path = UPLOADS_DIR / filename
        with open(self.current_file_path, "w") as f:
            f.write(self.pdb_content)

        ui.notify(f"File loaded: {filename}", type="positive")
        # Allow user to proceed to next step
        self.stepper.next()

    async def _run_analysis(self):
        """Run the molecular interaction analysis."""
        if self.analysis_running:
            ui.notify("Analysis already running", type="warning")
            return

        if not self.current_file_path or not self.current_file_path.exists():
            ui.notify("Please upload a PDB file first", type="warning")
            return

        self.analysis_running = True
        self.analyze_button.disable()
        self.status_label.text = "Running analysis..."

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
                ui.notify("Analysis completed!", type="positive")

                # Update results display
                await self.results_panel.update_results(
                    self.analyzer, self.pdb_content, self.current_file
                )

                # Auto-advance to results step
                self.stepper.next()
            else:
                self.status_label.text = "Analysis failed"
                ui.notify("Analysis failed. Please check the PDB file.", type="negative")

        except Exception as e:
            self.status_label.text = f"Error: {str(e)}"
            ui.notify(f"Error: {str(e)}", type="negative")

        finally:
            self.analysis_running = False
            self.analyze_button.enable()

    def _export_json(self):
        """Export results as JSON."""
        if not self.analyzer:
            ui.notify("No analysis results to export", type="warning")
            return

        base_name = Path(self.current_file).stem
        output_file = UPLOADS_DIR / f"{base_name}_results.json"
        export_to_json_single_file(
            self.analyzer, str(output_file), input_file=self.current_file
        )
        ui.download(str(output_file))
        ui.notify(f"Exported to {output_file.name}", type="positive")

    def _export_csv(self):
        """Export results as CSV."""
        if not self.analyzer:
            ui.notify("No analysis results to export", type="warning")
            return

        # Remove extension from current_file to get base name
        base_name = Path(self.current_file).stem
        base_filename = UPLOADS_DIR / base_name
        export_to_csv_files(self.analyzer, str(base_filename))

        # Download all generated CSV files
        csv_files = list(UPLOADS_DIR.glob(f"{base_name}_*.csv"))
        if csv_files:
            for csv_file in csv_files:
                ui.download(str(csv_file))
            ui.notify(f"Exported {len(csv_files)} CSV file(s)", type="positive")
        else:
            ui.notify("No CSV files were generated", type="warning")

    def _export_txt(self):
        """Export results as TXT."""
        if not self.analyzer:
            ui.notify("No analysis results to export", type="warning")
            return

        base_name = Path(self.current_file).stem
        output_file = UPLOADS_DIR / f"{base_name}_results.txt"
        export_to_txt_single_file(self.analyzer, str(output_file))
        ui.download(str(output_file))
        ui.notify(f"Exported to {output_file.name}", type="positive")


def create_app():
    """Create and configure the HBAT web application.

    :returns: None
    """

    @ui.page("/")
    def index():
        """Main page route."""
        # Load 3Dmol library for this page
        ui.add_head_html('<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>')

        hbat_app = HBATWebApp()
        hbat_app.create_ui()

    # Configure static files if directory exists
    static_dir = Path(__file__).parent / "static"
    if static_dir.exists():
        app.add_static_files("/static", str(static_dir))

    ui.run(
        title=f"{APP_NAME} Web",
        favicon="🧬",
        dark=False,
        reload=False,
        show=True,
    )


if __name__ == "__main__":
    create_app()
