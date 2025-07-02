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
from ..core.analysis import AnalysisParameters, HBondAnalyzer
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
        self.root = tk.Tk()
        self.root.title(f"{APP_NAME} v{APP_VERSION}")
        self.root.geometry(f"{GUIDefaults.WINDOW_WIDTH}x{GUIDefaults.WINDOW_HEIGHT}")
        self.root.minsize(GUIDefaults.MIN_WINDOW_WIDTH, GUIDefaults.MIN_WINDOW_HEIGHT)

        # Analysis components
        self.analyzer: Optional[HBondAnalyzer] = None
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

    def _create_toolbar(self) -> None:
        """Create the toolbar.

        Creates a toolbar with only a progress bar for showing operation status.
        The progress bar is hidden by default and only shown during operations.

        :returns: None
        :rtype: None
        """
        self.toolbar = ttk.Frame(self.root)
        # Don't pack the toolbar initially - it will be shown when needed

        # Progress bar
        self.progress_var = tk.DoubleVar()
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
        self.status_var = tk.StringVar()
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
            self.analyzer = HBondAnalyzer(params)

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
        # Hide progress bar and toolbar
        self.progress_bar.pack_forget()
        self.toolbar.pack_forget()
        self.analysis_menu.entryconfig(self.run_analysis_index, state=tk.NORMAL)

        # Update results panel
        self.results_panel.update_results(self.analyzer)

        # Update status
        stats = self.analyzer.get_statistics()
        self.status_var.set(
            f"Analysis complete - H-bonds: {stats['hydrogen_bonds']}, "
            f"X-bonds: {stats['halogen_bonds']}, π-interactions: {stats['pi_interactions']}"
        )

        messagebox.showinfo("Success", "Analysis completed successfully!")

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
        # Hide progress bar and toolbar
        self.progress_bar.pack_forget()
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
                f.write(f"Input file: {self.current_file}\n\n")

            # Write summary
            stats = self.analyzer.get_statistics()
            f.write("Summary:\n")
            f.write(f"  Hydrogen bonds: {stats['hydrogen_bonds']}\n")
            f.write(f"  Halogen bonds: {stats['halogen_bonds']}\n")
            f.write(f"  π interactions: {stats['pi_interactions']}\n")
            f.write(f"  Total interactions: {stats['total_interactions']}\n\n")

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
        self.parameters_window.geometry("600x700")
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
