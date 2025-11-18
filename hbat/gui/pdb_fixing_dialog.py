"""
PDB Fixing configuration dialog for HBAT GUI.

This module provides a dedicated dialog for configuring PDB fixing parameters,
separate from the main geometry parameter settings.
"""

import tkinter as tk
from tkinter import ttk
from typing import Any, Dict, Optional

from ..constants.parameters import ParametersDefault


class PDBFixingDialog:
    """Dialog for configuring PDB fixing parameters."""

    def __init__(self, parent: tk.Tk):
        """Initialize PDB fixing dialog.

        :param parent: Parent window
        :type parent: tk.Tk
        """
        self.parent = parent
        self.result = None

        # Create dialog window
        self.dialog = tk.Toplevel(parent)
        self.dialog.title("PDB Fixing Settings")
        self.dialog.geometry("800x600")
        self.dialog.resizable(False, False)

        # Make dialog modal
        self.dialog.transient(parent)
        self.dialog.grab_set()

        # Initialize variables
        self._init_variables()

        # Create widgets
        self._create_widgets()

        # Center the dialog
        self.dialog.update_idletasks()
        x = (self.dialog.winfo_screenwidth() // 2) - (self.dialog.winfo_width() // 2)
        y = (self.dialog.winfo_screenheight() // 2) - (self.dialog.winfo_height() // 2)
        self.dialog.geometry(f"+{x}+{y}")

        # Handle window closing
        self.dialog.protocol("WM_DELETE_WINDOW", self._cancel)

    def _init_variables(self):
        """Initialize tkinter variables with default values."""
        self.fix_pdb_enabled = tk.BooleanVar(value=ParametersDefault.FIX_PDB_ENABLED)
        self.fix_pdb_method = tk.StringVar(value=ParametersDefault.FIX_PDB_METHOD)
        self.fix_pdb_add_hydrogens = tk.BooleanVar(
            value=ParametersDefault.FIX_PDB_ADD_HYDROGENS
        )
        self.fix_pdb_add_heavy_atoms = tk.BooleanVar(
            value=ParametersDefault.FIX_PDB_ADD_HEAVY_ATOMS
        )
        self.fix_pdb_replace_nonstandard = tk.BooleanVar(
            value=ParametersDefault.FIX_PDB_REPLACE_NONSTANDARD
        )
        self.fix_pdb_remove_heterogens = tk.BooleanVar(
            value=ParametersDefault.FIX_PDB_REMOVE_HETEROGENS
        )
        self.fix_pdb_keep_water = tk.BooleanVar(
            value=ParametersDefault.FIX_PDB_KEEP_WATER
        )

    def _create_widgets(self):
        """Create dialog widgets."""
        # Main frame with padding
        main_frame = ttk.Frame(self.dialog, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Title
        title_label = ttk.Label(
            main_frame,
            text="PDB Structure Fixing Settings",
            font=("TkDefaultFont", 12, "bold"),
        )
        title_label.pack(pady=(0, 20))

        # Enable PDB fixing checkbox
        self.enable_check = ttk.Checkbutton(
            main_frame,
            text="Enable PDB structure fixing",
            variable=self.fix_pdb_enabled,
            command=self._on_enable_changed,
        )
        self.enable_check.pack(anchor=tk.W, pady=(0, 15))

        # Method selection frame
        method_frame = ttk.LabelFrame(main_frame, text="Fixing Method", padding="10")
        method_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Radiobutton(
            method_frame,
            text="PDBFixer (recommended for proteins)",
            variable=self.fix_pdb_method,
            value="pdbfixer",
            command=self._on_method_changed,
        ).pack(anchor=tk.W)

        ttk.Radiobutton(
            method_frame,
            text="OpenBabel",
            variable=self.fix_pdb_method,
            value="openbabel",
            command=self._on_method_changed,
        ).pack(anchor=tk.W)

        # Operations frame
        operations_frame = ttk.LabelFrame(main_frame, text="Operations", padding="10")
        operations_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # Add hydrogens (both methods)
        self.add_hydrogens_check = ttk.Checkbutton(
            operations_frame,
            text="Add missing hydrogen atoms",
            variable=self.fix_pdb_add_hydrogens,
        )
        self.add_hydrogens_check.pack(anchor=tk.W, pady=2)

        # PDBFixer-only operations
        self.add_heavy_atoms_check = ttk.Checkbutton(
            operations_frame,
            text="Add missing heavy atoms (PDBFixer only)",
            variable=self.fix_pdb_add_heavy_atoms,
        )
        self.add_heavy_atoms_check.pack(anchor=tk.W, pady=2)

        self.replace_nonstandard_check = ttk.Checkbutton(
            operations_frame,
            text="Replace nonstandard residues (PDBFixer only)",
            variable=self.fix_pdb_replace_nonstandard,
        )
        self.replace_nonstandard_check.pack(anchor=tk.W, pady=2)

        self.remove_heterogens_check = ttk.Checkbutton(
            operations_frame,
            text="Remove heterogens (PDBFixer only)",
            variable=self.fix_pdb_remove_heterogens,
            command=self._on_remove_heterogens_changed,
        )
        self.remove_heterogens_check.pack(anchor=tk.W, pady=2)

        # Keep water sub-option
        self.keep_water_check = ttk.Checkbutton(
            operations_frame,
            text="    Keep water molecules",
            variable=self.fix_pdb_keep_water,
        )
        self.keep_water_check.pack(anchor=tk.W, pady=2)

        # Store PDBFixer-only widgets
        self.pdbfixer_only_widgets = [
            self.add_heavy_atoms_check,
            self.replace_nonstandard_check,
            self.remove_heterogens_check,
            self.keep_water_check,
        ]

        # Button frame
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=(20, 0))

        ttk.Button(button_frame, text="OK", command=self._ok).pack(
            side=tk.RIGHT, padx=(5, 0)
        )

        ttk.Button(button_frame, text="Cancel", command=self._cancel).pack(
            side=tk.RIGHT
        )

        ttk.Button(
            button_frame, text="Reset to Defaults", command=self._reset_defaults
        ).pack(side=tk.LEFT, padx=(0, 5))

        ttk.Button(
            button_frame, text="Manage Presets...", command=self._open_preset_manager
        ).pack(side=tk.LEFT, padx=5)

        # Initialize widget states
        self._update_widget_states()

    def _on_enable_changed(self):
        """Handle enable checkbox change."""
        self._update_widget_states()

    def _on_method_changed(self):
        """Handle method selection change."""
        self._update_widget_states()

    def _on_remove_heterogens_changed(self):
        """Handle remove heterogens checkbox change."""
        self._update_widget_states()

    def _update_widget_states(self):
        """Update widget enable/disable states based on selections."""
        enabled = self.fix_pdb_enabled.get()
        method = self.fix_pdb_method.get()

        # Enable/disable based on main checkbox
        for widget in [self.add_hydrogens_check] + self.pdbfixer_only_widgets:
            widget.configure(state="normal" if enabled else "disabled")

        # Further disable PDBFixer-only options if OpenBabel selected
        if enabled and method == "openbabel":
            for widget in self.pdbfixer_only_widgets:
                widget.configure(state="disabled")

        # Enable keep water only if remove heterogens is checked
        if enabled and method == "pdbfixer" and self.fix_pdb_remove_heterogens.get():
            self.keep_water_check.configure(state="normal")
        else:
            self.keep_water_check.configure(state="disabled")

    def _reset_defaults(self):
        """Reset all values to defaults."""
        self.fix_pdb_enabled.set(ParametersDefault.FIX_PDB_ENABLED)
        self.fix_pdb_method.set(ParametersDefault.FIX_PDB_METHOD)
        self.fix_pdb_add_hydrogens.set(ParametersDefault.FIX_PDB_ADD_HYDROGENS)
        self.fix_pdb_add_heavy_atoms.set(ParametersDefault.FIX_PDB_ADD_HEAVY_ATOMS)
        self.fix_pdb_replace_nonstandard.set(
            ParametersDefault.FIX_PDB_REPLACE_NONSTANDARD
        )
        self.fix_pdb_remove_heterogens.set(ParametersDefault.FIX_PDB_REMOVE_HETEROGENS)
        self.fix_pdb_keep_water.set(ParametersDefault.FIX_PDB_KEEP_WATER)
        self._update_widget_states()

    def _ok(self):
        """Handle OK button - save settings and close."""
        self.result = {
            "enabled": self.fix_pdb_enabled.get(),
            "method": self.fix_pdb_method.get(),
            "add_hydrogens": self.fix_pdb_add_hydrogens.get(),
            "add_heavy_atoms": self.fix_pdb_add_heavy_atoms.get(),
            "replace_nonstandard": self.fix_pdb_replace_nonstandard.get(),
            "remove_heterogens": self.fix_pdb_remove_heterogens.get(),
            "keep_water": self.fix_pdb_keep_water.get(),
        }
        self.dialog.destroy()

    def _cancel(self):
        """Handle Cancel button - close without saving."""
        self.result = None
        self.dialog.destroy()

    def get_parameters(self) -> Optional[Dict[str, Any]]:
        """Get the configured parameters.

        :returns: Dictionary of PDB fixing parameters or None if cancelled
        :rtype: Optional[Dict[str, Any]]
        """
        self.dialog.wait_window()
        return self.result

    def set_parameters(self, params: Dict[str, Any]):
        """Set parameters from a dictionary.

        :param params: Dictionary of PDB fixing parameters
        :type params: Dict[str, Any]
        """
        if "enabled" in params:
            self.fix_pdb_enabled.set(params["enabled"])
        if "method" in params:
            self.fix_pdb_method.set(params["method"])
        if "add_hydrogens" in params:
            self.fix_pdb_add_hydrogens.set(params["add_hydrogens"])
        if "add_heavy_atoms" in params:
            self.fix_pdb_add_heavy_atoms.set(params["add_heavy_atoms"])
        if "replace_nonstandard" in params:
            self.fix_pdb_replace_nonstandard.set(params["replace_nonstandard"])
        if "remove_heterogens" in params:
            self.fix_pdb_remove_heterogens.set(params["remove_heterogens"])
        if "keep_water" in params:
            self.fix_pdb_keep_water.set(params["keep_water"])
        self._update_widget_states()

    def _open_preset_manager(self):
        """Open the preset manager dialog."""
        from tkinter import messagebox

        from ..core.analysis import AnalysisParameters
        from .preset_manager_dialog import PresetManagerDialog

        # Create a dummy AnalysisParameters object with current PDB fixing settings
        # The preset manager expects AnalysisParameters, so we'll use default geometry params
        # and the user can focus on PDB fixing presets
        current_params = AnalysisParameters()

        # Open preset manager
        dialog = PresetManagerDialog(self.dialog, current_params)
        result = dialog.get_result()

        if result and "parameters" in result and "pdb_fixing" in result["parameters"]:
            # Apply only the PDB fixing part of the loaded preset
            pdb_params = result["parameters"]["pdb_fixing"]
            self.set_parameters(pdb_params)
            messagebox.showinfo("Success", "PDB fixing preset loaded successfully")
        elif result:
            messagebox.showwarning(
                "Warning", "Selected preset does not contain PDB fixing parameters"
            )
