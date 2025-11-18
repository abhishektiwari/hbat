"""
Main GUI window for HBAT application.

This module provides the main tkinter interface for the HBAT application,
allowing users to load PDB files, configure analysis parameters, and view results.
"""

import asyncio
import csv
import json
import os
import tkinter as tk
import webbrowser
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Any, Dict, Optional

import tk_async_execute as tae

from ..constants import APP_NAME, APP_VERSION, GUIDefaults
from ..constants.parameters import ParametersDefault
from ..core.analysis import (
    AnalysisParameters,
    NPMolecularInteractionAnalyzer,
)
from ..export.results import (
    export_to_csv_files,
    export_to_json_files,
    export_to_json_single_file,
    export_to_txt_single_file,
)
from .geometry_cutoffs_dialog import GeometryCutoffsDialog
from .pdb_fixing_dialog import PDBFixingDialog
from .results_panel import ResultsPanel


class MainWindow:
    """Main application window for HBAT.

    This class provides the primary GUI interface for HBAT, including
    file loading, parameter configuration, analysis execution, and
    results visualization.

    :param None: This class takes no parameters during initialization
    """

    def __init__(self) -> None:
        """Initialize the main window.

        Sets up the complete GUI interface including menus, toolbar,
        main content area, and status bar.

        :returns: None
        :rtype: None
        """
        # Initialize HBAT environment first
        try:
            from ..core.app_config import get_hbat_config, initialize_hbat_environment

            initialize_hbat_environment(verbose=True)
            self.hbat_config = get_hbat_config()
        except ImportError:
            self.hbat_config = None

        self.root = tk.Tk()
        self.root.title(f"{APP_NAME} v{APP_VERSION}")
        self.root.geometry(f"{GUIDefaults.WINDOW_WIDTH}x{GUIDefaults.WINDOW_HEIGHT}")
        self.root.minsize(GUIDefaults.MIN_WINDOW_WIDTH, GUIDefaults.MIN_WINDOW_HEIGHT)

        # Analysis components
        self.analyzer: Optional[NPMolecularInteractionAnalyzer] = None
        self.current_file: Optional[str] = None
        self.analysis_running = False

        # Create UI components
        self._create_menu()
        self._create_toolbar()
        self._create_main_content()
        self._create_status_bar()

        # Set up event handlers
        self.root.protocol("WM_DELETE_WINDOW", self._on_closing)

        # Initialize async executor
        tae.start()

    def _create_menu(self) -> None:
        """Create the menu bar.

        Sets up the application menu with File, Analysis, Tools, and Help menus,
        including keyboard shortcuts and event bindings.

        :returns: None
        :rtype: None
        """
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        self.menubar = menubar  # Store reference for state updates

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(
            label="Open PDB File...", accelerator="Ctrl+O", command=self._open_file
        )
        file_menu.add_separator()
        file_menu.add_command(
            label="Save Results...", accelerator="Ctrl+S", command=self._save_results
        )
        file_menu.add_command(
            label="Save Fixed PDB...",
            accelerator="Ctrl+Shift+S",
            command=self._save_fixed_pdb,
        )
        file_menu.add_separator()
        file_menu.add_command(
            label="Export All...", accelerator="Ctrl+E", command=self._export_all
        )
        file_menu.add_separator()
        file_menu.add_command(
            label="Exit", accelerator="Ctrl+Q", command=self._on_closing
        )

        # Analysis menu
        analysis_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)
        self.analysis_menu = analysis_menu  # Store reference
        analysis_menu.add_command(
            label="Run Analysis",
            accelerator="F5",
            command=self._run_analysis,
            state=tk.DISABLED,
        )
        self.run_analysis_index = 0  # Index of "Run Analysis" menu item
        analysis_menu.add_command(label="Clear Results", command=self._clear_results)

        # Settings menu
        settings_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Settings", menu=settings_menu)
        settings_menu.add_command(
            label="Geometry Cutoffs...", command=self._open_parameters_window
        )
        settings_menu.add_command(
            label="PDB Fixing...", command=self._open_pdb_fixing_window
        )
        settings_menu.add_separator()
        settings_menu.add_command(
            label="Manage Presets...", command=self._open_preset_manager
        )
        settings_menu.add_separator()
        settings_menu.add_command(
            label="Reset Parameters", command=self._reset_parameters
        )

        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="About", command=self._show_about)
        help_menu.add_command(label="User Guide", command=self._show_help)

        # Bind keyboard shortcuts
        self.root.bind("<Control-o>", lambda e: self._open_file())
        self.root.bind("<Control-s>", lambda e: self._save_results())
        self.root.bind("<Control-e>", lambda e: self._export_all())
        self.root.bind("<Control-q>", lambda e: self._on_closing())
        self.root.bind("<F5>", lambda e: self._run_analysis())
        self.root.bind("<Control-Shift-S>", lambda e: self._save_fixed_pdb())

    def _create_toolbar(self) -> None:
        """Create the toolbar.

        Creates a toolbar with progress bar and performance indicator.
        The toolbar is hidden by default and only shown during operations.

        :returns: None
        :rtype: None
        """
        self.toolbar = ttk.Frame(self.root)
        # Don't pack the toolbar initially - it will be shown when needed

        # Performance indicator
        self.performance_label = ttk.Label(
            self.toolbar,
            text="⚡",
            foreground="green",
        )

        # Progress bar
        self.progress_var = tk.DoubleVar(master=self.root)
        self.progress_bar = ttk.Progressbar(
            self.toolbar, variable=self.progress_var, mode="indeterminate"
        )

    def _create_main_content(self) -> None:
        """Create the main content area.

        Sets up the main interface with a vertical paned window containing PDB file
        content on top (30%) and results display area below (70%).

        :returns: None
        :rtype: None
        """
        # Create main paned window with vertical orientation
        main_paned = ttk.PanedWindow(self.root, orient=tk.VERTICAL)
        main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Top panel - PDB File content (30% of height)
        top_frame = ttk.Frame(main_paned)
        main_paned.add(top_frame, weight=3)  # 30% weight ratio

        # Create notebook for PDB file tabs
        self.left_notebook = ttk.Notebook(top_frame)
        self.left_notebook.pack(fill=tk.BOTH, expand=True)

        # File content tab
        file_frame = ttk.Frame(self.left_notebook)
        self.left_notebook.add(file_frame, text="PDB File")

        # Create text widget with both vertical and horizontal scrollbars
        text_frame = ttk.Frame(file_frame)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.file_text = tk.Text(text_frame, wrap=tk.NONE, font=("Courier", 12))
        file_v_scrollbar = ttk.Scrollbar(
            text_frame, orient=tk.VERTICAL, command=self.file_text.yview
        )
        file_h_scrollbar = ttk.Scrollbar(
            text_frame, orient=tk.HORIZONTAL, command=self.file_text.xview
        )
        self.file_text.configure(
            yscrollcommand=file_v_scrollbar.set, xscrollcommand=file_h_scrollbar.set
        )

        self.file_text.grid(row=0, column=0, sticky="nsew")
        file_v_scrollbar.grid(row=0, column=1, sticky="ns")
        file_h_scrollbar.grid(row=1, column=0, sticky="ew")

        text_frame.grid_rowconfigure(0, weight=1)
        text_frame.grid_columnconfigure(0, weight=1)

        # Fixed PDB content tab
        fixed_file_frame = ttk.Frame(self.left_notebook)
        self.left_notebook.add(fixed_file_frame, text="Fixed PDB")

        # Create text widget with both vertical and horizontal scrollbars
        fixed_text_frame = ttk.Frame(fixed_file_frame)
        fixed_text_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.fixed_file_text = tk.Text(
            fixed_text_frame, wrap=tk.NONE, font=("Courier", 12)
        )
        fixed_v_scrollbar = ttk.Scrollbar(
            fixed_text_frame, orient=tk.VERTICAL, command=self.fixed_file_text.yview
        )
        fixed_h_scrollbar = ttk.Scrollbar(
            fixed_text_frame, orient=tk.HORIZONTAL, command=self.fixed_file_text.xview
        )
        self.fixed_file_text.configure(
            yscrollcommand=fixed_v_scrollbar.set, xscrollcommand=fixed_h_scrollbar.set
        )

        self.fixed_file_text.grid(row=0, column=0, sticky="nsew")
        fixed_v_scrollbar.grid(row=0, column=1, sticky="ns")
        fixed_h_scrollbar.grid(row=1, column=0, sticky="ew")

        fixed_text_frame.grid_rowconfigure(0, weight=1)
        fixed_text_frame.grid_columnconfigure(0, weight=1)

        # Add context menu for Fixed PDB tab
        self._create_fixed_pdb_context_menu()

        # Store parameters for session persistence
        self.session_parameters = None

        # Bottom panel - Results (70% of height)
        bottom_frame = ttk.Frame(main_paned)
        main_paned.add(bottom_frame, weight=7)  # 70% weight ratio

        self.results_panel = ResultsPanel(bottom_frame)

        # Set initial pane positions (30% for top, 70% for bottom)
        # Position the sash at 30% of the total height
        main_paned.update_idletasks()  # Ensure the window is rendered
        total_height = self.root.winfo_height() - 60  # Account for padding and menu
        if total_height > 0:
            sash_position = int(total_height * 0.3)
            main_paned.sashpos(0, sash_position)

    def _create_status_bar(self) -> None:
        """Create the status bar.

        Creates a status bar at the bottom of the window to display
        application state and progress information.

        :returns: None
        :rtype: None
        """
        self.status_var = tk.StringVar(master=self.root)
        self.status_var.set("Ready")

        status_frame = ttk.Frame(self.root)
        status_frame.pack(fill=tk.X, side=tk.BOTTOM)

        ttk.Label(
            status_frame, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W
        ).pack(fill=tk.X, padx=2, pady=1)

    def _open_file(self) -> None:
        """Open a PDB file.

        Displays a file dialog to select a PDB file, loads its content,
        and enables analysis functionality.

        :returns: None
        :rtype: None
        """
        filename = filedialog.askopenfilename(
            title="Open PDB File",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
        )

        if filename:
            try:
                self.current_file = filename
                self._load_file_content(filename)
                # Enable "Run Analysis" menu item
                self.analysis_menu.entryconfig(self.run_analysis_index, state=tk.NORMAL)
                self.status_var.set(f"Loaded: {os.path.basename(filename)}")
                self._clear_results()
                self._clear_fixed_pdb_content()

            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")
                self.status_var.set("Error loading file")

    def _load_file_content(self, filename: str) -> None:
        """Load and display file content.

        Reads the PDB file content and displays it in the text widget
        with syntax highlighting for PDB record types.

        :param filename: Path to the PDB file to load
        :type filename: str
        :returns: None
        :rtype: None
        :raises Exception: If file cannot be read
        """
        try:
            # Show progress for large files
            self._show_loading_progress("Loading file...")

            # Load file in chunks to prevent GUI freezing
            self._load_file_in_chunks(filename)

            # Add to recent files if config is available
            if self.hbat_config:
                self.hbat_config.add_recent_file(filename)

        except Exception as e:
            raise Exception(f"Cannot read file: {e}")
        finally:
            self._hide_loading_progress()

    def _load_file_in_chunks(self, filename: str) -> None:
        """Load file content in chunks to prevent GUI freezing.

        :param filename: Path to the PDB file to load
        :type filename: str
        :returns: None
        :rtype: None
        """
        chunk_size = 50000  # Process in 50KB chunks
        self.file_text.delete(1.0, tk.END)

        try:
            with open(filename, "r") as file:
                while True:
                    chunk = file.read(chunk_size)
                    if not chunk:
                        break

                    # Insert chunk and update display
                    self.file_text.insert(tk.END, chunk)
                    self.root.update_idletasks()  # Process pending GUI events

            # Apply syntax highlighting after loading
            self.root.after_idle(self._highlight_pdb_records)

        except Exception as e:
            raise Exception(f"Cannot read file: {e}")

    def _show_loading_progress(self, message: str) -> None:
        """Show loading progress indicator.

        :param message: Status message to display
        :type message: str
        :returns: None
        :rtype: None
        """
        self.status_var.set(message)
        if not self.toolbar.winfo_ismapped():
            self.toolbar.pack(fill=tk.X, padx=5, pady=2)
        self.progress_bar.pack(fill=tk.BOTH, padx=5, expand=True)
        self.progress_bar.config(mode="indeterminate")
        self.progress_bar.start(GUIDefaults.PROGRESS_BAR_INTERVAL)
        self.root.update_idletasks()

    def _hide_loading_progress(self) -> None:
        """Hide loading progress indicator.

        :returns: None
        :rtype: None
        """
        self.progress_bar.stop()
        self.progress_bar.pack_forget()
        if self.toolbar.winfo_ismapped():
            self.toolbar.pack_forget()

    def _highlight_pdb_records(self) -> None:
        """Highlight important PDB record types.

        Applies color coding to different PDB record types (ATOM, HETATM,
        HEADER, etc.) for better readability.

        :returns: None
        :rtype: None
        """
        # Configure text tags
        self.file_text.tag_configure("atom", foreground="blue")
        self.file_text.tag_configure("hetatm", foreground="red")
        self.file_text.tag_configure(
            "header", foreground="green", font=("Courier", 12, "bold")
        )

        content = self.file_text.get(1.0, tk.END)
        lines = content.split("\n")

        for i, line in enumerate(lines):
            line_start = f"{i+1}.0"
            line_end = f"{i+1}.end"

            if line.startswith("ATOM"):
                self.file_text.tag_add("atom", line_start, line_end)
            elif line.startswith("HETATM"):
                self.file_text.tag_add("hetatm", line_start, line_end)
            elif line.startswith(("HEADER", "TITLE", "COMPND")):
                self.file_text.tag_add("header", line_start, line_end)

    def _run_analysis(self) -> None:
        """Run the molecular interaction analysis.

        Initiates analysis using async/await pattern to keep GUI responsive.

        :returns: None
        :rtype: None
        """
        if not self.current_file:
            messagebox.showwarning("Warning", "Please open a PDB file first.")
            return

        if self.analysis_running:
            messagebox.showinfo("Info", "Analysis is already running.")
            return

        # Get parameters from session storage or use defaults
        if self.session_parameters:
            params = self.session_parameters
        else:
            # Use default parameters if none have been set
            params = AnalysisParameters()
            self.session_parameters = params

        # Apply PDB fixing parameters if available
        # Note: Default AnalysisParameters() should have PDB fixing enabled by default
        if (
            hasattr(self, "session_pdb_fixing_params")
            and self.session_pdb_fixing_params
        ):
            # Update the params object with session PDB fixing settings
            params.fix_pdb_enabled = self.session_pdb_fixing_params["enabled"]
            params.fix_pdb_method = self.session_pdb_fixing_params["method"]
            params.fix_pdb_add_hydrogens = self.session_pdb_fixing_params[
                "add_hydrogens"
            ]
            params.fix_pdb_add_heavy_atoms = self.session_pdb_fixing_params[
                "add_heavy_atoms"
            ]
            params.fix_pdb_replace_nonstandard = self.session_pdb_fixing_params[
                "replace_nonstandard"
            ]
            params.fix_pdb_remove_heterogens = self.session_pdb_fixing_params[
                "remove_heterogens"
            ]
            params.fix_pdb_keep_water = self.session_pdb_fixing_params["keep_water"]

        # Start async analysis without popup window
        tae.async_execute(
            self._perform_analysis_async(params), visible=False, show_exceptions=False
        )

    async def _perform_analysis_async(self, params: AnalysisParameters) -> None:
        """Perform the analysis asynchronously to keep GUI responsive.

        :param params: Analysis parameters to use
        :type params: AnalysisParameters
        :returns: None
        :rtype: None
        """
        try:
            # Set up UI for analysis
            self.analysis_running = True
            self.analysis_menu.entryconfig(self.run_analysis_index, state=tk.DISABLED)

            # Show toolbar and progress bar
            if not self.toolbar.winfo_ismapped():
                self.toolbar.pack(fill=tk.X, padx=5, pady=2)
            self.performance_label.pack(side=tk.LEFT, padx=5)
            self.progress_bar.pack(fill=tk.BOTH, padx=5, expand=True)
            self.progress_bar.config(mode="indeterminate")
            self.progress_bar.start(GUIDefaults.PROGRESS_BAR_INTERVAL)
            self.status_var.set("Running analysis...")

            # Create analyzer
            self.analyzer = NPMolecularInteractionAnalyzer(params)

            # Set up progress callback for direct GUI updates
            def progress_callback(message: str) -> None:
                # Update progress directly without creating new async task
                self.root.after(0, lambda: self.status_var.set(message))

            self.analyzer.progress_callback = progress_callback

            # Run analysis in executor to avoid blocking
            success = await asyncio.get_event_loop().run_in_executor(
                None, self.analyzer.analyze_file, self.current_file
            )

            if success:
                await self._analysis_complete_async()
            else:
                await self._analysis_error_async("Analysis failed")

        except Exception as e:
            await self._analysis_error_async(str(e))

    async def _update_progress_async(self, message: str) -> None:
        """Update progress message asynchronously.

        :param message: Progress message to display
        :type message: str
        :returns: None
        :rtype: None
        """
        self.status_var.set(message)

    async def _analysis_complete_async(self) -> None:
        """Handle successful analysis completion asynchronously.

        :returns: None
        :rtype: None
        """
        self.analysis_running = False
        self.progress_bar.stop()
        self.progress_bar.config(mode="determinate")
        self.progress_var.set(0)
        # Hide progress bar, performance label, and toolbar
        self.progress_bar.pack_forget()
        self.performance_label.pack_forget()
        self.toolbar.pack_forget()
        self.analysis_menu.entryconfig(self.run_analysis_index, state=tk.NORMAL)

        # Update results panel
        self.results_panel.update_results(self.analyzer)

        # Update Fixed PDB tab if PDB fixing was applied
        self._update_fixed_pdb_content()

        # Update status
        summary = self.analyzer.get_summary()
        self.status_var.set(
            f"Analysis complete - H-bonds: {summary['hydrogen_bonds']['count']}, "
            f"X-bonds: {summary['halogen_bonds']['count']}, π-interactions: {summary['pi_interactions']['count']}"
        )

        messagebox.showinfo("Success", "Analysis completed successfully!")

    async def _analysis_error_async(self, error_msg: str) -> None:
        """Handle analysis error asynchronously.

        :param error_msg: Error message to display
        :type error_msg: str
        :returns: None
        :rtype: None
        """
        self.analysis_running = False
        self.progress_bar.stop()
        self.progress_bar.config(mode="determinate")
        self.progress_var.set(0)
        # Hide progress bar, performance label, and toolbar
        self.progress_bar.pack_forget()
        self.performance_label.pack_forget()
        self.toolbar.pack_forget()
        self.analysis_menu.entryconfig(self.run_analysis_index, state=tk.NORMAL)
        self.status_var.set("Analysis failed")
        messagebox.showerror("Analysis Error", f"Analysis failed:\n{error_msg}")

    def _clear_results(self) -> None:
        """Clear analysis results.

        Clears all analysis results from the interface and resets
        the analyzer state.

        :returns: None
        :rtype: None
        """
        self.results_panel.clear_results()
        self.analyzer = None
        self._clear_fixed_pdb_content()
        self.status_var.set("Results cleared")

    def _save_results(self) -> None:
        """Save analysis results to file.

        Displays a file dialog to save analysis results in CSV or JSON format.
        Requires completed analysis to function.

        :returns: None
        :rtype: None
        """
        if not self.analyzer:
            messagebox.showwarning("Warning", "No results to save. Run analysis first.")
            return

        filename = filedialog.asksaveasfilename(
            title="Save Results",
            defaultextension=".csv",
            filetypes=[
                ("CSV files", "*.csv"),
                ("JSON files", "*.json"),
            ],
        )

        if filename:
            try:
                file_path = Path(filename)
                extension = file_path.suffix.lower()

                self._export_results_to_file(filename)

                # Show appropriate success message based on format
                if extension in [".csv", ".json"]:
                    base_name = file_path.stem
                    directory = file_path.parent
                    messagebox.showinfo(
                        "Success",
                        f"Results saved as separate {extension.upper()} files:\n"
                        f"- {base_name}_h_bonds{extension}\n"
                        f"- {base_name}_x_bonds{extension}\n"
                        f"- {base_name}_pi_interactions{extension}\n"
                        f"- {base_name}_pi_pi_interactions{extension}\n"
                        f"- {base_name}_carbonyl_interactions{extension}\n"
                        f"- {base_name}_n_pi_interactions{extension}\n"
                        f"- {base_name}_cooperativity_chains{extension}\n\n"
                        f"Location: {directory}",
                    )
                    self.status_var.set(
                        f"Results saved as multiple {extension.upper()} files"
                    )
                else:
                    messagebox.showinfo("Success", f"Results saved to {filename}")
                    self.status_var.set(
                        f"Results saved to {os.path.basename(filename)}"
                    )
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save results:\n{str(e)}")

    def _export_results_to_file(self, filename: str) -> None:
        """Export results to a file.

        Supports multiple formats:
        - .txt: Single text file with all results
        - .csv: Separate CSV files for each interaction type
        - .json: Separate JSON files for each interaction type

        :param filename: Path to the output file
        :type filename: str
        :returns: None
        :rtype: None
        """
        file_path = Path(filename)
        extension = file_path.suffix.lower()

        if extension == ".txt":
            self._export_text_file(filename)
        elif extension == ".json":
            export_to_json_files(self.analyzer, str(file_path), self.current_file)
        elif extension == ".csv":
            export_to_csv_files(self.analyzer, str(file_path))
        else:
            # Default to CSV format for unknown extensions
            export_to_csv_files(self.analyzer, str(file_path))

    def _export_text_file(self, filename: str) -> None:
        """Export results to a single text file.

        :param filename: Path to the output file
        :type filename: str
        """
        export_to_txt_single_file(self.analyzer, filename)

    def _export_all(self) -> None:
        """Export all results to a single file.

        Exports all analysis results to a single file in either text (.txt)
        or JSON (.json) format.

        :returns: None
        :rtype: None
        """
        if not self.analyzer:
            messagebox.showwarning(
                "Warning", "No results to export. Run analysis first."
            )
            return

        filename = filedialog.asksaveasfilename(
            title="Export All Results",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("JSON files", "*.json")],
        )

        if filename:
            try:
                file_path = Path(filename)
                extension = file_path.suffix.lower()

                if extension == ".json":
                    export_to_json_single_file(
                        self.analyzer, filename, self.current_file
                    )
                else:
                    # Default to text format
                    export_to_txt_single_file(self.analyzer, filename)

                messagebox.showinfo("Success", f"Results exported to {filename}")
                self.status_var.set(f"Results exported to {os.path.basename(filename)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to export results:\n{str(e)}")

    def _reset_parameters(self) -> None:
        """Reset analysis parameters to defaults.

        Restores all analysis parameters to their default values
        as defined in the application constants.

        :returns: None
        :rtype: None
        """
        # Reset session parameters to defaults
        self.session_parameters = AnalysisParameters()
        self.status_var.set("Parameters reset to defaults")

    def _show_about(self) -> None:
        """Show about dialog.

        Displays application information including version, authors,
        and institutional affiliation.

        :returns: None
        :rtype: None
        """
        about_text = f"""
{APP_NAME} v{APP_VERSION}

A high-performance tool for analyzing molecular interactions in protein structures.

Features:
• Hydrogen bond detection
• Halogen bond analysis  
• π-interaction identification
• High-performance analysis engine
• Comprehensive visualization and export options

Author: Abhishek Tiwari
        """
        messagebox.showinfo("About HBAT", about_text.strip())

    def _open_parameters_window(self) -> None:
        """Open geometry cutoffs configuration dialog.

        Creates a modal dialog for configuring geometry cutoff parameters,
        preserving any existing parameter values for the session.

        :returns: None
        :rtype: None
        """
        # Get current session parameters or defaults
        current_params = self.session_parameters or AnalysisParameters()

        # Create dialog
        dialog = GeometryCutoffsDialog(self.root, current_params)

        # Get results
        result = dialog.get_result()

        if result:
            # Store the updated parameters in session
            self.session_parameters = result
            self.status_var.set("Geometry cutoffs updated")

    def _open_pdb_fixing_window(self) -> None:
        """Open PDB fixing settings dialog.

        Creates a modal dialog for configuring PDB fixing parameters.

        :returns: None
        :rtype: None
        """
        # Create dialog
        dialog = PDBFixingDialog(self.root)

        # Set current parameters if available
        if (
            hasattr(self, "session_pdb_fixing_params")
            and self.session_pdb_fixing_params
        ):
            dialog.set_parameters(self.session_pdb_fixing_params)

        # Get results
        result = dialog.get_parameters()

        if result:
            # Store the PDB fixing parameters
            self.session_pdb_fixing_params = result

            # Update the status bar if appropriate
            if hasattr(self, "status_bar"):
                fixing_status = "enabled" if result["enabled"] else "disabled"
                self.status_bar.config(
                    text=f"PDB fixing {fixing_status} ({result['method']})"
                )
                # Clear status after 3 seconds
                self.root.after(3000, lambda: self.status_bar.config(text="Ready"))

    def _show_help(self) -> None:
        """Show help dialog.

        Opens the HBAT documentation website in the default web browser.

        :returns: None
        :rtype: None
        """
        webbrowser.open("https://hbat.abhishek-tiwari.com")

    def _update_fixed_pdb_content(self) -> None:
        """Update the Fixed PDB tab with the processed structure content.

        Shows the PDB structure after any fixing has been applied,
        or indicates that no fixing was done.

        :returns: None
        :rtype: None
        """
        if not self.analyzer:
            self._clear_fixed_pdb_content()
            return

        summary = self.analyzer.get_summary()
        pdb_info = summary.get("pdb_fixing", {})

        self.fixed_file_text.delete(1.0, tk.END)

        if pdb_info.get("applied", False):
            # PDB fixing was applied - show the fixed structure from saved file
            try:
                fixed_file_path = pdb_info.get("fixed_file_path")
                if fixed_file_path and os.path.exists(fixed_file_path):
                    # Read content from the saved fixed PDB file
                    with open(fixed_file_path, "r") as f:
                        fixed_content = f.read()
                    self.fixed_file_text.insert(1.0, fixed_content)
                    self._highlight_pdb_records_in_widget(self.fixed_file_text)
                else:
                    # Fallback to generating content from atoms
                    fixed_content = self._generate_pdb_content_from_atoms()
                    self.fixed_file_text.insert(1.0, fixed_content)
                    self._highlight_pdb_records_in_widget(self.fixed_file_text)

                # Show tab title indicating changes
                self.left_notebook.tab(1, text="Fixed PDB ✓")

            except Exception as e:
                self.fixed_file_text.insert(
                    tk.END, f"Error loading fixed PDB content: {e}\n"
                )
                self.fixed_file_text.insert(
                    tk.END, "Original structure was used for analysis."
                )
        else:
            # No PDB fixing applied
            if "error" in pdb_info:
                self.fixed_file_text.insert(tk.END, "PDB Fixing Status: Failed\n")
                self.fixed_file_text.insert(tk.END, f"Error: {pdb_info['error']}\n\n")
                self.fixed_file_text.insert(
                    tk.END, "The original structure was used for analysis.\n"
                )
                self.fixed_file_text.insert(
                    tk.END, "Consider enabling PDB fixing in the analysis parameters."
                )
            else:
                self.fixed_file_text.insert(
                    tk.END, "PDB Fixing Status: Not Applied\n\n"
                )
                self.fixed_file_text.insert(
                    tk.END,
                    "The original structure was used for analysis without modifications.\n\n",
                )
                self.fixed_file_text.insert(tk.END, "To apply PDB fixing:\n")
                self.fixed_file_text.insert(
                    tk.END, "1. Open Settings → Geometry Cutoffs\n"
                )
                self.fixed_file_text.insert(tk.END, "2. Enable 'Fix PDB' option\n")
                self.fixed_file_text.insert(
                    tk.END, "3. Select fixing method (OpenBabel or PDBFixer)\n"
                )
                self.fixed_file_text.insert(tk.END, "4. Re-run the analysis")

            # Show tab title indicating no changes
            self.left_notebook.tab(1, text="Fixed PDB")

    def _clear_fixed_pdb_content(self) -> None:
        """Clear the Fixed PDB tab content.

        :returns: None
        :rtype: None
        """
        self.fixed_file_text.delete(1.0, tk.END)
        self.fixed_file_text.insert(tk.END, "No analysis results available.\n\n")
        self.fixed_file_text.insert(
            tk.END,
            "Please load a PDB file and run analysis to see the processed structure.",
        )
        self.left_notebook.tab(1, text="Fixed PDB")

    def _generate_pdb_content_from_atoms(self) -> str:
        """Generate PDB format content from the analyzer's atoms.

        Creates PDB format text from the processed atom list, which may
        include atoms added during PDB fixing.

        :returns: PDB format content
        :rtype: str
        :raises Exception: If atoms cannot be converted to PDB format
        """
        if not self.analyzer or not self.analyzer.parser.atoms:
            return "No atoms available in processed structure."

        lines = []

        # Add header information
        lines.append("REMARK   1 PROCESSED BY HBAT")
        lines.append(
            "REMARK   1 THIS STRUCTURE MAY INCLUDE MODIFICATIONS FROM PDB FIXING"
        )
        lines.append("REMARK   1")

        # Add PDB fixing information to header
        summary = self.analyzer.get_summary()
        pdb_info = summary.get("pdb_fixing", {})

        if pdb_info.get("applied", False):
            lines.append(f"REMARK   2 PDB FIXING METHOD: {pdb_info['method'].upper()}")
            lines.append(f"REMARK   2 ORIGINAL ATOMS: {pdb_info['original_atoms']}")
            lines.append(f"REMARK   2 FIXED ATOMS: {pdb_info['fixed_atoms']}")
            if pdb_info.get("added_hydrogens", 0) > 0:
                lines.append(
                    f"REMARK   2 ADDED HYDROGENS: {pdb_info['added_hydrogens']}"
                )
            lines.append(f"REMARK   2 REDETECTED BONDS: {pdb_info['redetected_bonds']}")

        lines.append("REMARK   2")

        # Convert atoms to PDB format
        for atom in self.analyzer.parser.atoms:
            line = self._atom_to_pdb_line(atom)
            lines.append(line)

        lines.append("END")

        return "\n".join(lines)

    def _atom_to_pdb_line(self, atom) -> str:
        """Convert an Atom object to PDB format line.

        :param atom: Atom object to convert
        :type atom: Atom
        :returns: PDB format line
        :rtype: str
        """
        # PDB format: ATOM/HETATM with specific column positions
        return (
            f"{atom.record_type:<6}"  # Record type (ATOM/HETATM)
            f"{atom.serial:>5} "  # Atom serial number
            f"{atom.name:<4}"  # Atom name
            f"{atom.alt_loc:>1}"  # Alternate location
            f"{atom.res_name:>3} "  # Residue name
            f"{atom.chain_id:>1}"  # Chain ID
            f"{atom.res_seq:>4}"  # Residue sequence number
            f"{atom.i_code:>1}   "  # Insertion code
            f"{atom.coords.x:>8.3f}"  # X coordinate
            f"{atom.coords.y:>8.3f}"  # Y coordinate
            f"{atom.coords.z:>8.3f}"  # Z coordinate
            f"{atom.occupancy:>6.2f}"  # Occupancy
            f"{atom.temp_factor:>6.2f}"  # Temperature factor
            f"          "  # Blank spaces
            f"{atom.element:>2}"  # Element symbol
            f"{atom.charge:>2}"  # Charge
        )

    def _highlight_pdb_records_in_widget(self, text_widget) -> None:
        """Apply PDB record highlighting to a specific text widget.

        :param text_widget: Text widget to apply highlighting to
        :type text_widget: tk.Text
        :returns: None
        :rtype: None
        """
        # Configure text tags
        text_widget.tag_configure("atom", foreground="blue")
        text_widget.tag_configure("hetatm", foreground="red")
        text_widget.tag_configure(
            "remark", foreground="green", font=("Courier", 12, "bold")
        )

        content = text_widget.get(1.0, tk.END)
        lines = content.split("\n")

        for i, line in enumerate(lines):
            line_start = f"{i+1}.0"
            line_end = f"{i+1}.end"

            if line.startswith("ATOM"):
                text_widget.tag_add("atom", line_start, line_end)
            elif line.startswith("HETATM"):
                text_widget.tag_add("hetatm", line_start, line_end)
            elif line.startswith("REMARK"):
                text_widget.tag_add("remark", line_start, line_end)

    def _create_fixed_pdb_context_menu(self) -> None:
        """Create context menu for the Fixed PDB text widget.

        :returns: None
        :rtype: None
        """
        # Create context menu
        self.fixed_pdb_context_menu = tk.Menu(self.root, tearoff=0)
        self.fixed_pdb_context_menu.add_command(
            label="Save Fixed PDB...", command=self._save_fixed_pdb
        )
        self.fixed_pdb_context_menu.add_separator()
        self.fixed_pdb_context_menu.add_command(
            label="Select All",
            command=lambda: self.fixed_file_text.tag_add(tk.SEL, "1.0", tk.END),
        )
        self.fixed_pdb_context_menu.add_command(
            label="Copy",
            command=lambda: self.root.clipboard_clear()
            or self.root.clipboard_append(self.fixed_file_text.selection_get()),
        )

        # Bind right-click to show context menu
        self.fixed_file_text.bind("<Button-3>", self._show_fixed_pdb_context_menu)

    def _show_fixed_pdb_context_menu(self, event) -> None:
        """Show context menu for Fixed PDB tab.

        :param event: Mouse event
        :type event: tkinter.Event
        :returns: None
        :rtype: None
        """
        try:
            self.fixed_pdb_context_menu.tk_popup(event.x_root, event.y_root)
        finally:
            self.fixed_pdb_context_menu.grab_release()

    def _save_fixed_pdb(self) -> None:
        """Save the Fixed PDB content to a file.

        :returns: None
        :rtype: None
        """
        if not self.analyzer:
            messagebox.showwarning("Warning", "No analysis results available to save.")
            return

        # Check if PDB fixing was actually applied
        summary = self.analyzer.get_summary()
        pdb_info = summary.get("pdb_fixing", {})

        if not pdb_info.get("applied", False):
            result = messagebox.askyesno(
                "No PDB Fixing Applied",
                "PDB fixing was not applied to this structure. "
                "The saved file will be identical to the original PDB.\n\n"
                "Do you want to continue?",
            )
            if not result:
                return

        # Get save filename
        filename = filedialog.asksaveasfilename(
            title="Save Fixed PDB",
            defaultextension=".pdb",
            filetypes=[
                ("PDB files", "*.pdb"),
                ("Text files", "*.txt"),
                ("All files", "*.*"),
            ],
            initialname=f"{os.path.splitext(os.path.basename(self.current_file or 'structure'))[0]}_fixed.pdb",
        )

        if filename:
            try:
                # Try to copy from the saved fixed file first
                fixed_file_path = pdb_info.get("fixed_file_path")
                if (
                    pdb_info.get("applied", False)
                    and fixed_file_path
                    and os.path.exists(fixed_file_path)
                ):
                    # Copy the existing fixed file
                    import shutil

                    shutil.copy2(fixed_file_path, filename)
                else:
                    # Fallback to saving from text widget content
                    content = self.fixed_file_text.get(1.0, tk.END)
                    with open(filename, "w") as f:
                        f.write(content)

                messagebox.showinfo("Success", f"Fixed PDB saved to {filename}")
                self.status_var.set(f"Fixed PDB saved to {os.path.basename(filename)}")

            except Exception as e:
                messagebox.showerror("Error", f"Failed to save fixed PDB:\n{str(e)}")

    def _open_preset_manager(self) -> None:
        """Open the preset manager dialog.

        Creates a standalone preset manager that can load/save presets
        and apply them to the current session parameters.

        :returns: None
        :rtype: None
        """
        from .preset_manager_dialog import PresetManagerDialog

        # Get current session parameters or defaults
        current_params = self.session_parameters or AnalysisParameters()

        # Open preset manager
        dialog = PresetManagerDialog(self.root, current_params)
        result = dialog.get_result()

        if result:
            # Apply the loaded preset to session parameters
            self._apply_preset_to_session(result)
            self.status_var.set("Preset loaded and applied to session")

    def _apply_preset_to_session(self, preset_data: Dict[str, Any]) -> None:
        """Apply preset data to session parameters.

        :param preset_data: Preset data to apply
        :type preset_data: Dict[str, Any]
        :returns: None
        :rtype: None
        """
        if "parameters" not in preset_data:
            messagebox.showerror(
                "Error", "Invalid preset format: missing 'parameters' section"
            )
            return

        params_data = preset_data["parameters"]

        # Create new AnalysisParameters with preset values
        kwargs = {}

        # Apply hydrogen bond parameters
        if "hydrogen_bonds" in params_data:
            hb = params_data["hydrogen_bonds"]
            kwargs.update(
                {
                    "hb_distance_cutoff": hb.get(
                        "h_a_distance_cutoff", ParametersDefault.HB_DISTANCE_CUTOFF
                    ),
                    "hb_angle_cutoff": hb.get(
                        "dha_angle_cutoff", ParametersDefault.HB_ANGLE_CUTOFF
                    ),
                    "hb_donor_acceptor_cutoff": hb.get(
                        "d_a_distance_cutoff", ParametersDefault.HB_DA_DISTANCE
                    ),
                }
            )

        # Apply halogen bond parameters
        if "halogen_bonds" in params_data:
            xb = params_data["halogen_bonds"]
            kwargs.update(
                {
                    "xb_distance_cutoff": xb.get(
                        "x_a_distance_cutoff", ParametersDefault.XB_DISTANCE_CUTOFF
                    ),
                    "xb_angle_cutoff": xb.get(
                        "dxa_angle_cutoff", ParametersDefault.XB_ANGLE_CUTOFF
                    ),
                }
            )

        # Apply π interaction parameters
        if "pi_interactions" in params_data:
            pi = params_data["pi_interactions"]
            kwargs.update(
                {
                    "pi_distance_cutoff": pi.get(
                        "h_pi_distance_cutoff", ParametersDefault.PI_DISTANCE_CUTOFF
                    ),
                    "pi_angle_cutoff": pi.get(
                        "dh_pi_angle_cutoff", ParametersDefault.PI_ANGLE_CUTOFF
                    ),
                }
            )

        # Apply general parameters
        if "general" in params_data:
            gen = params_data["general"]
            kwargs.update(
                {
                    "covalent_cutoff_factor": gen.get(
                        "covalent_cutoff_factor",
                        ParametersDefault.COVALENT_CUTOFF_FACTOR,
                    ),
                    "analysis_mode": gen.get(
                        "analysis_mode", ParametersDefault.ANALYSIS_MODE
                    ),
                }
            )

        # Create new parameters object and store in session
        self.session_parameters = AnalysisParameters(**kwargs)

    def _on_closing(self) -> None:
        """Handle window closing event.

        Handles application shutdown, checking for running analysis
        and prompting user confirmation if needed.

        :returns: None
        :rtype: None
        """
        if self.analysis_running:
            result = messagebox.askyesno(
                "Confirm Exit", "Analysis is running. Are you sure you want to exit?"
            )
            if not result:
                return

        # Stop async executor
        tae.stop()
        self.root.destroy()

    def run(self) -> None:
        """Start the GUI application.

        Enters the main GUI event loop to begin accepting user interactions.
        This method blocks until the application is closed.

        :returns: None
        :rtype: None
        """
        self.root.mainloop()
