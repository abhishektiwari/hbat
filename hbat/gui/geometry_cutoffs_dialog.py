"""
Geometry cutoffs configuration dialog for HBAT GUI.

This module provides a dialog for configuring molecular interaction 
analysis parameters (distances, angles, etc.) without PDB fixing options.
"""

import tkinter as tk
from tkinter import messagebox, ttk
from typing import Any, Dict, Optional

from ..constants.parameters import (
    AnalysisModes,
    AnalysisParameters,
    ParameterRanges,
    ParametersDefault,
)


class GeometryCutoffsDialog:
    """Dialog for configuring geometry cutoffs parameters."""

    def __init__(self, parent: tk.Tk, current_params: Optional[AnalysisParameters] = None):
        """Initialize geometry cutoffs dialog.

        :param parent: Parent window
        :type parent: tk.Tk
        :param current_params: Current analysis parameters
        :type current_params: Optional[AnalysisParameters]
        """
        self.parent = parent
        self.current_params = current_params or AnalysisParameters()
        self.result = None
        
        # Create dialog window
        self.dialog = tk.Toplevel(parent)
        self.dialog.title("Geometry Cutoffs")
        self.dialog.geometry("800x600")
        self.dialog.resizable(True, True)
        
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
        
        # Set initial values
        self.set_parameters(self.current_params)

    def _init_variables(self):
        """Initialize tkinter variables with default values."""
        self.analysis_mode = tk.StringVar(value=ParametersDefault.ANALYSIS_MODE)
        self.covalent_factor = tk.DoubleVar(value=ParametersDefault.COVALENT_CUTOFF_FACTOR)
        self.hb_distance = tk.DoubleVar(value=ParametersDefault.HB_DISTANCE_CUTOFF)
        self.hb_angle = tk.DoubleVar(value=ParametersDefault.HB_ANGLE_CUTOFF)
        self.da_distance = tk.DoubleVar(value=ParametersDefault.HB_DA_DISTANCE)
        self.xb_distance = tk.DoubleVar(value=ParametersDefault.XB_DISTANCE_CUTOFF)
        self.xb_angle = tk.DoubleVar(value=ParametersDefault.XB_ANGLE_CUTOFF)
        self.pi_distance = tk.DoubleVar(value=ParametersDefault.PI_DISTANCE_CUTOFF)
        self.pi_angle = tk.DoubleVar(value=ParametersDefault.PI_ANGLE_CUTOFF)

    def _create_widgets(self) -> None:
        """Create and layout all parameter widgets.

        :returns: None
        :rtype: None
        """
        # Main frame with padding
        main_frame = ttk.Frame(self.dialog, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Create container for scrollable content
        content_frame = ttk.Frame(main_frame)
        content_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # Create scrollable area
        canvas = tk.Canvas(content_frame)
        scrollbar = ttk.Scrollbar(content_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        # Grid layout for canvas and scrollbar
        canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
        
        content_frame.grid_rowconfigure(0, weight=1)
        content_frame.grid_columnconfigure(0, weight=1)

        # Create parameter groups in order
        self._create_general_parameters(scrollable_frame)
        self._create_hydrogen_bond_parameters(scrollable_frame)
        self._create_halogen_bond_parameters(scrollable_frame)
        self._create_pi_interaction_parameters(scrollable_frame)

        # Buttons at bottom - separate from scrollable content
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, side=tk.BOTTOM)
        
        ttk.Button(
            button_frame, 
            text="OK", 
            command=self._ok
        ).pack(side=tk.RIGHT, padx=(5, 0))
        
        ttk.Button(
            button_frame, 
            text="Cancel", 
            command=self._cancel
        ).pack(side=tk.RIGHT)

        ttk.Button(
            button_frame, text="Reset to Defaults", command=self._set_defaults
        ).pack(side=tk.LEFT, padx=(0, 5))

        ttk.Button(button_frame, text="Manage Presets...", command=self._open_preset_manager).pack(
            side=tk.LEFT, padx=5
        )

    def _create_general_parameters(self, parent):
        """Create general analysis parameters."""
        group = ttk.LabelFrame(parent, text="General Parameters", padding=10)
        group.pack(fill=tk.X, padx=10, pady=5)

        # Analysis mode
        ttk.Label(group, text="Analysis Mode:").grid(
            row=0, column=0, sticky=tk.W, pady=2
        )
        self.analysis_mode = tk.StringVar(value=ParametersDefault.ANALYSIS_MODE)
        mode_frame = ttk.Frame(group)
        mode_frame.grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        ttk.Radiobutton(
            mode_frame,
            text="Complete PDB Analysis",
            variable=self.analysis_mode,
            value="complete",
        ).pack(anchor=tk.W)
        ttk.Radiobutton(
            mode_frame,
            text="Local Interactions Only",
            variable=self.analysis_mode,
            value="local",
        ).pack(anchor=tk.W)

        # Covalent bond cutoff factor
        ttk.Label(group, text="Covalent Bond Factor:").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        self.covalent_factor = tk.DoubleVar(
            value=ParametersDefault.COVALENT_CUTOFF_FACTOR
        )
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_COVALENT_FACTOR,
            to=ParameterRanges.MAX_COVALENT_FACTOR,
            variable=self.covalent_factor,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        # Value display
        factor_label = ttk.Label(group, text="")
        factor_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_factor_label(*args):
            factor_label.config(text=f"{self.covalent_factor.get():.2f}")

        self.covalent_factor.trace("w", update_factor_label)
        update_factor_label()

    def _create_hydrogen_bond_parameters(self, parent):
        """Create hydrogen bond parameter controls."""
        group = ttk.LabelFrame(parent, text="Hydrogen Bond Parameters", padding=10)
        group.pack(fill=tk.X, padx=10, pady=5)

        # H...A distance
        ttk.Label(group, text="H...A Distance (Å):").grid(
            row=0, column=0, sticky=tk.W, pady=2
        )
        self.hb_distance = tk.DoubleVar(value=ParametersDefault.HB_DISTANCE_CUTOFF)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.hb_distance,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        hb_dist_label = ttk.Label(group, text="")
        hb_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_hb_dist(*args):
            hb_dist_label.config(text=f"{self.hb_distance.get():.1f}")

        self.hb_distance.trace("w", update_hb_dist)
        update_hb_dist()

        # D-H...A angle
        ttk.Label(group, text="D-H...A Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        self.hb_angle = tk.DoubleVar(value=ParametersDefault.HB_ANGLE_CUTOFF)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.hb_angle,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        hb_angle_label = ttk.Label(group, text="")
        hb_angle_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_hb_angle(*args):
            hb_angle_label.config(text=f"{self.hb_angle.get():.0f}")

        self.hb_angle.trace("w", update_hb_angle)
        update_hb_angle()

        # D...A distance
        ttk.Label(group, text="D...A Distance (Å):").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )
        self.da_distance = tk.DoubleVar(value=ParametersDefault.HB_DA_DISTANCE)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.da_distance,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=2, column=1, sticky=tk.W, padx=10, pady=2)

        da_dist_label = ttk.Label(group, text="")
        da_dist_label.grid(row=2, column=2, sticky=tk.W, padx=5, pady=2)

        def update_da_dist(*args):
            da_dist_label.config(text=f"{self.da_distance.get():.1f}")

        self.da_distance.trace("w", update_da_dist)
        update_da_dist()

    def _create_halogen_bond_parameters(self, parent):
        """Create halogen bond parameter controls."""
        group = ttk.LabelFrame(parent, text="Halogen Bond Parameters", padding=10)
        group.pack(fill=tk.X, padx=10, pady=5)

        # X...A distance
        ttk.Label(group, text="X...A Distance (Å):").grid(
            row=0, column=0, sticky=tk.W, pady=2
        )
        self.xb_distance = tk.DoubleVar(value=ParametersDefault.XB_DISTANCE_CUTOFF)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.xb_distance,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        xb_dist_label = ttk.Label(group, text="")
        xb_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_xb_dist(*args):
            xb_dist_label.config(text=f"{self.xb_distance.get():.1f}")

        self.xb_distance.trace("w", update_xb_dist)
        update_xb_dist()

        # C-X...A angle
        ttk.Label(group, text="C-X...A Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        self.xb_angle = tk.DoubleVar(value=ParametersDefault.XB_ANGLE_CUTOFF)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.xb_angle,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        xb_angle_label = ttk.Label(group, text="")
        xb_angle_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_xb_angle(*args):
            xb_angle_label.config(text=f"{self.xb_angle.get():.0f}")

        self.xb_angle.trace("w", update_xb_angle)
        update_xb_angle()

    def _create_pi_interaction_parameters(self, parent):
        """Create π interaction parameter controls."""
        group = ttk.LabelFrame(parent, text="π Interaction Parameters", padding=10)
        group.pack(fill=tk.X, padx=10, pady=5)

        # H...π distance
        ttk.Label(group, text="H...π Distance (Å):").grid(
            row=0, column=0, sticky=tk.W, pady=2
        )
        self.pi_distance = tk.DoubleVar(value=ParametersDefault.PI_DISTANCE_CUTOFF)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.pi_distance,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        pi_dist_label = ttk.Label(group, text="")
        pi_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_dist(*args):
            pi_dist_label.config(text=f"{self.pi_distance.get():.1f}")

        self.pi_distance.trace("w", update_pi_dist)
        update_pi_dist()

        # D-H...π angle
        ttk.Label(group, text="D-H...π Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        self.pi_angle = tk.DoubleVar(value=ParametersDefault.PI_ANGLE_CUTOFF)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.pi_angle,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        pi_angle_label = ttk.Label(group, text="")
        pi_angle_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_angle(*args):
            pi_angle_label.config(text=f"{self.pi_angle.get():.0f}")

        self.pi_angle.trace("w", update_pi_angle)
        update_pi_angle()

    def get_parameters(self) -> AnalysisParameters:
        """Get current parameter values.

        :returns: Current analysis parameters
        :rtype: AnalysisParameters
        """
        return AnalysisParameters(
            hb_distance_cutoff=self.hb_distance.get(),
            hb_angle_cutoff=self.hb_angle.get(),
            hb_donor_acceptor_cutoff=self.da_distance.get(),
            xb_distance_cutoff=self.xb_distance.get(),
            xb_angle_cutoff=self.xb_angle.get(),
            pi_distance_cutoff=self.pi_distance.get(),
            pi_angle_cutoff=self.pi_angle.get(),
            covalent_cutoff_factor=self.covalent_factor.get(),
            analysis_mode=self.analysis_mode.get(),
        )

    def set_parameters(self, params: AnalysisParameters) -> None:
        """Set parameter values from AnalysisParameters object.

        :param params: Analysis parameters to set
        :type params: AnalysisParameters
        """
        self.hb_distance.set(params.hb_distance_cutoff)
        self.hb_angle.set(params.hb_angle_cutoff)
        self.da_distance.set(params.hb_donor_acceptor_cutoff)
        self.xb_distance.set(params.xb_distance_cutoff)
        self.xb_angle.set(params.xb_angle_cutoff)
        self.pi_distance.set(params.pi_distance_cutoff)
        self.pi_angle.set(params.pi_angle_cutoff)
        self.covalent_factor.set(params.covalent_cutoff_factor)
        self.analysis_mode.set(params.analysis_mode)

    def _set_defaults(self):
        """Reset all parameters to default values."""
        default_params = AnalysisParameters()
        self.set_parameters(default_params)

    def reset_to_defaults(self) -> None:
        """Public method to reset parameters to defaults."""
        self._set_defaults()

    def _open_preset_manager(self):
        """Open the preset manager dialog."""
        from .preset_manager_dialog import PresetManagerDialog
        
        # Get current parameters
        current_params = self.get_parameters()
        
        # Open preset manager
        dialog = PresetManagerDialog(self.dialog, current_params)
        result = dialog.get_result()
        
        if result:
            # Apply the loaded preset
            self._apply_preset_data(result)
            messagebox.showinfo("Success", "Preset loaded successfully")
    
    def _apply_preset_data(self, data: Dict[str, Any]) -> None:
        """Apply preset data to parameters."""
        if "parameters" not in data:
            raise ValueError("Invalid preset format: missing 'parameters' section")
            
        params = data["parameters"]
        
        # Apply hydrogen bond parameters
        if "hydrogen_bonds" in params:
            hb = params["hydrogen_bonds"]
            self.hb_distance.set(hb.get("h_a_distance_cutoff", ParametersDefault.HB_DISTANCE_CUTOFF))
            self.hb_angle.set(hb.get("dha_angle_cutoff", ParametersDefault.HB_ANGLE_CUTOFF))
            self.da_distance.set(hb.get("d_a_distance_cutoff", ParametersDefault.HB_DA_DISTANCE))
            
        # Apply halogen bond parameters
        if "halogen_bonds" in params:
            xb = params["halogen_bonds"]
            self.xb_distance.set(xb.get("x_a_distance_cutoff", ParametersDefault.XB_DISTANCE_CUTOFF))
            self.xb_angle.set(xb.get("cxa_angle_cutoff", ParametersDefault.XB_ANGLE_CUTOFF))
            
        # Apply π interaction parameters
        if "pi_interactions" in params:
            pi = params["pi_interactions"]
            self.pi_distance.set(pi.get("h_pi_distance_cutoff", ParametersDefault.PI_DISTANCE_CUTOFF))
            self.pi_angle.set(pi.get("dh_pi_angle_cutoff", ParametersDefault.PI_ANGLE_CUTOFF))
            
        # Apply general parameters
        if "general" in params:
            gen = params["general"]
            self.covalent_factor.set(gen.get("covalent_cutoff_factor", ParametersDefault.COVALENT_CUTOFF_FACTOR))
            self.analysis_mode.set(gen.get("analysis_mode", ParametersDefault.ANALYSIS_MODE))

    def _ok(self):
        """Handle OK button - save settings and close."""
        self.result = self.get_parameters()
        self.dialog.destroy()
        
    def _cancel(self):
        """Handle Cancel button - close without saving.""" 
        self.result = None
        self.dialog.destroy()
        
    def get_result(self) -> Optional[AnalysisParameters]:
        """Get the configured parameters.
        
        :returns: Analysis parameters or None if cancelled
        :rtype: Optional[AnalysisParameters]
        """
        self.dialog.wait_window()
        return self.result