"""
Main GUI window for HBAT application.

This module provides the main tkinter interface for the HBAT application,
allowing users to load PDB files, configure analysis parameters, and view results.
"""

import os
import threading
import tkinter as tk
import webbrowser
from tkinter import filedialog, messagebox, scrolledtext, ttk
from typing import Optional

from ..constants import APP_NAME, APP_VERSION, GUIDefaults
from ..core.analysis import (
    AnalysisParameters,
    NPMolecularInteractionAnalyzer,
)
from .parameter_panel import ParameterPanel
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
        self.analysis_thread: Optional[threading.Thread] = None

        # Create UI components
        self._create_menu()
        self._create_toolbar()
        self._create_main_content()
        self._create_status_bar()

        # Set up event handlers
        self.root.protocol("WM_DELETE_WINDOW", self._on_closing)

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
            label="Edit Parameters...", command=self._open_parameters_window
        )
        settings_menu.add_command(
            label="Reset Parameters", command=self._reset_parameters
        )
        settings_menu.add_separator()

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

        Sets up the main interface with a paned window containing file content,
        parameter panels, and results display areas.

        :returns: None
        :rtype: None
        """
        # Create main paned window
        main_paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left panel - File content and parameters
        left_frame = ttk.Frame(main_paned)
        main_paned.add(left_frame, weight=1)

        # Create notebook for left panel
        left_notebook = ttk.Notebook(left_frame)
        left_notebook.pack(fill=tk.BOTH, expand=True)

        # File content tab
        file_frame = ttk.Frame(left_notebook)
        left_notebook.add(file_frame, text="PDB File")

        self.file_text = scrolledtext.ScrolledText(
            file_frame, wrap=tk.NONE, font=("Courier", 12)
        )
        self.file_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Fixed PDB content tab
        fixed_file_frame = ttk.Frame(left_notebook)
        left_notebook.add(fixed_file_frame, text="Fixed PDB")

        self.fixed_file_text = scrolledtext.ScrolledText(
            fixed_file_frame, wrap=tk.NONE, font=("Courier", 12)
        )
        self.fixed_file_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Add context menu for Fixed PDB tab
        self._create_fixed_pdb_context_menu()

        # Store reference to the notebook for later updates
        self.left_notebook = left_notebook

        # Store parameters separately (no longer in tab)
        self.parameter_panel = None
        self.parameters_window = None
        self.session_parameters = None  # Store parameters for session persistence

        # Right panel - Results
        right_frame = ttk.Frame(main_paned)
        main_paned.add(right_frame, weight=2)

        self.results_panel = ResultsPanel(right_frame)

        # Set initial pane positions
        main_paned.sashpos(0, GUIDefaults.LEFT_PANEL_WIDTH)

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
            with open(filename, "r") as file:
                content = file.read()

            self.file_text.delete(1.0, tk.END)
            self.file_text.insert(1.0, content)

            # Highlight ATOM and HETATM lines
            self._highlight_pdb_records()

            # Add to recent files if config is available
            if self.hbat_config:
                self.hbat_config.add_recent_file(filename)

        except Exception as e:
            raise Exception(f"Cannot read file: {e}")

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

        Initiates analysis in a separate thread using current parameters
        and loaded PDB file. Updates UI to show progress.

        :returns: None
        :rtype: None
        """
        if not self.current_file:
            messagebox.showwarning("Warning", "Please open a PDB file first.")
            return

        if self.analysis_thread and self.analysis_thread.is_alive():
            messagebox.showinfo("Info", "Analysis is already running.")
            return

        # Get parameters from the parameter panel or session storage
        if self.parameter_panel:
            params = self.parameter_panel.get_parameters()
        elif self.session_parameters:
            params = self.session_parameters
        else:
            # Use default parameters if none have been set
            from ..core.analysis import AnalysisParameters

            params = AnalysisParameters()
            self.session_parameters = params

        # Start analysis in a separate thread
        self.analysis_thread = threading.Thread(
            target=self._perform_analysis, args=(params,)
        )
        self.analysis_thread.daemon = True
        self.analysis_thread.start()

        # Update UI
        self.analysis_menu.entryconfig(self.run_analysis_index, state=tk.DISABLED)
        # Show toolbar and progress bar
        if not self.toolbar.winfo_ismapped():
            self.toolbar.pack(fill=tk.X, padx=5, pady=2)
        self.performance_label.pack(side=tk.LEFT, padx=5)
        self.progress_bar.pack(fill=tk.BOTH, padx=5, expand=True)
        self.progress_bar.config(mode="indeterminate")
        self.progress_bar.start(GUIDefaults.PROGRESS_BAR_INTERVAL)
        self.status_var.set("Running analysis...")

    def _perform_analysis(self, params: AnalysisParameters) -> None:
        """Perform the analysis in a separate thread.

        Executes the molecular interaction analysis using the provided
        parameters and updates the UI upon completion or error.

        :param params: Analysis parameters to use
        :type params: AnalysisParameters
        :returns: None
        :rtype: None
        """
        try:
            # Create analyzer
            self.analyzer = NPMolecularInteractionAnalyzer(params)

            # Update status
            self.root.after(0, lambda: self.status_var.set("Running analysis..."))

            # Run analysis
            success = self.analyzer.analyze_file(self.current_file)

            if success:
                # Update results on main thread
                self.root.after(0, self._analysis_complete)
            else:
                self.root.after(0, self._analysis_error, "Analysis failed")

        except Exception as e:
            self.root.after(0, self._analysis_error, str(e))

    def _analysis_complete(self) -> None:
        """Handle successful analysis completion.

        Updates the UI after successful analysis, stops progress indication,
        displays results, and shows completion notification.

        :returns: None
        :rtype: None
        """
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

        messagebox.showinfo("Success", f"Analysis completed successfully!")

    def _analysis_error(self, error_msg: str) -> None:
        """Handle analysis error.

        Updates the UI after analysis failure, stops progress indication,
        and displays error message to the user.

        :param error_msg: Error message to display
        :type error_msg: str
        :returns: None
        :rtype: None
        """
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

        Displays a file dialog to save analysis results in text format.
        Requires completed analysis to function.

        :returns: None
        :rtype: None
        """
        if not self.analyzer:
            messagebox.showwarning("Warning", "No results to save. Run analysis first.")
            return

        filename = filedialog.asksaveasfilename(
            title="Save Results",
            defaultextension=".txt",
            filetypes=[
                ("Text files", "*.txt"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ],
        )

        if filename:
            try:
                self._export_results_to_file(filename)
                messagebox.showinfo("Success", f"Results saved to {filename}")
                self.status_var.set(f"Results saved to {os.path.basename(filename)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save results:\n{str(e)}")

    def _export_results_to_file(self, filename: str) -> None:
        """Export results to a file.

        Writes complete analysis results to the specified file in
        human-readable text format.

        :param filename: Path to the output file
        :type filename: str
        :returns: None
        :rtype: None
        """
        with open(filename, "w") as f:
            f.write("HBAT Analysis Results\n")
            f.write("=" * 50 + "\n\n")

            if self.current_file:
                f.write(f"Input file: {self.current_file}\n")

            f.write(f"Analysis engine: HBAT\n")
            f.write(f"HBAT version: {APP_VERSION}\n\n")

            # Write summary
            summary = self.analyzer.get_summary()
            f.write("Summary:\n")
            f.write(f"  Hydrogen bonds: {summary['hydrogen_bonds']['count']}\n")
            f.write(f"  Halogen bonds: {summary['halogen_bonds']['count']}\n")
            f.write(f"  π interactions: {summary['pi_interactions']['count']}\n")
            f.write(f"  Total interactions: {summary['total_interactions']}\n\n")

            # Write detailed results
            f.write("Hydrogen Bonds:\n")
            f.write("-" * 30 + "\n")
            for hb in self.analyzer.hydrogen_bonds:
                f.write(f"{hb}\n")

            f.write("\nHalogen Bonds:\n")
            f.write("-" * 30 + "\n")
            for xb in self.analyzer.halogen_bonds:
                f.write(f"{xb}\n")

            f.write("\nπ Interactions:\n")
            f.write("-" * 30 + "\n")
            for pi in self.analyzer.pi_interactions:
                f.write(f"{pi}\n")

    def _export_all(self) -> None:
        """Export all results in multiple formats.

        Exports analysis results to a directory in multiple file formats
        for comprehensive data preservation.

        :returns: None
        :rtype: None
        """
        if not self.analyzer:
            messagebox.showwarning(
                "Warning", "No results to export. Run analysis first."
            )
            return

        directory = filedialog.askdirectory(title="Select Export Directory")
        if directory:
            try:
                base_name = (
                    os.path.splitext(os.path.basename(self.current_file))[0]
                    if self.current_file
                    else "hbat_results"
                )

                # Export text summary
                self._export_results_to_file(
                    os.path.join(directory, f"{base_name}_summary.txt")
                )

                # Export Fixed PDB if available
                summary = self.analyzer.get_summary()
                pdb_info = summary.get("pdb_fixing", {})
                if pdb_info.get("applied", False):
                    fixed_pdb_path = os.path.join(directory, f"{base_name}_fixed.pdb")
                    try:
                        fixed_content = self._generate_pdb_content_from_atoms()
                        with open(fixed_pdb_path, "w") as f:
                            f.write(fixed_content)
                    except Exception as e:
                        print(f"Warning: Could not export fixed PDB: {e}")

                messagebox.showinfo("Success", f"Results exported to {directory}")
                self.status_var.set("All results exported")
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
        from ..core.analysis import AnalysisParameters

        self.session_parameters = AnalysisParameters()

        # If parameters window is open, update it too
        if self.parameter_panel:
            self.parameter_panel.reset_to_defaults()

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
        """Open parameters configuration in a popup window.

        Creates a popup window containing the parameter panel, preserving
        any existing parameter values for the session.

        :returns: None
        :rtype: None
        """
        if self.parameters_window and self.parameters_window.winfo_exists():
            # Bring existing window to front
            self.parameters_window.lift()
            self.parameters_window.focus_force()
            return

        # Create new parameters window
        self.parameters_window = tk.Toplevel(self.root)
        self.parameters_window.title("Analysis Parameters")
        self.parameters_window.geometry("800x800")
        self.parameters_window.resizable(True, True)

        # Create parameter panel in popup window
        self.parameter_panel = ParameterPanel(self.parameters_window)
        self.parameter_panel.frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Restore previous session parameters if they exist
        if self.session_parameters:
            self.parameter_panel.set_parameters(self.session_parameters)

        # Add Apply and Cancel buttons
        button_frame = ttk.Frame(self.parameters_window)
        button_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

        ttk.Button(button_frame, text="Apply", command=self._apply_parameters).pack(
            side=tk.RIGHT, padx=(5, 0)
        )

        ttk.Button(button_frame, text="Cancel", command=self._cancel_parameters).pack(
            side=tk.RIGHT
        )

        # Center the window
        self.parameters_window.transient(self.root)
        self.parameters_window.grab_set()

        # Handle window closing
        self.parameters_window.protocol("WM_DELETE_WINDOW", self._on_parameters_close)

    def _apply_parameters(self) -> None:
        """Apply parameter changes and close the window.

        Saves parameters to session storage for persistence.

        :returns: None
        :rtype: None
        """
        if self.parameter_panel:
            # Save current parameters to session storage
            self.session_parameters = self.parameter_panel.get_parameters()
            self.status_var.set("Parameters updated")
        self._close_parameters_window()

    def _cancel_parameters(self) -> None:
        """Cancel parameter changes and revert to session values.

        Restores the last saved session parameters if they exist.

        :returns: None
        :rtype: None
        """
        if self.parameter_panel and self.session_parameters:
            # Revert to last saved session parameters
            self.parameter_panel.set_parameters(self.session_parameters)
        self._close_parameters_window()

    def _close_parameters_window(self) -> None:
        """Close the parameters window.

        :returns: None
        :rtype: None
        """
        if self.parameters_window:
            self.parameters_window.destroy()
            self.parameters_window = None
            self.parameter_panel = None

    def _on_parameters_close(self) -> None:
        """Handle parameters window closing via window manager.

        Treats window closing as a cancel operation.

        :returns: None
        :rtype: None
        """
        self._cancel_parameters()

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
                    tk.END, "1. Open Settings → Edit Parameters\n"
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

    def _on_closing(self) -> None:
        """Handle window closing event.

        Handles application shutdown, checking for running analysis
        and prompting user confirmation if needed.

        :returns: None
        :rtype: None
        """
        if self.analysis_thread and self.analysis_thread.is_alive():
            result = messagebox.askyesno(
                "Confirm Exit", "Analysis is running. Are you sure you want to exit?"
            )
            if not result:
                return

        self.root.destroy()

    def run(self) -> None:
        """Start the GUI application.

        Enters the main GUI event loop to begin accepting user interactions.
        This method blocks until the application is closed.

        :returns: None
        :rtype: None
        """
        self.root.mainloop()
