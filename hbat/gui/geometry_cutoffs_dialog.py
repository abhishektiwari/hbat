"""
Geometry cutoffs configuration dialog for HBAT GUI.

This module provides a dialog for configuring molecular interaction
analysis parameters (distances, angles, etc.) without PDB fixing options.
"""

import tkinter as tk
from tkinter import messagebox, ttk
from typing import Any, Dict, Optional


class ToolTip:
    """Simple tooltip widget for providing help text on hover."""

    def __init__(self, widget, text: str, delay: int = 500):
        self.widget = widget
        self.text = text
        self.delay = delay
        self.tooltip = None
        self.id = None
        self.widget.bind("<Enter>", self.enter, "+")
        self.widget.bind("<Leave>", self.leave, "+")
        self.widget.bind("<Motion>", self.motion, "+")

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hide()

    def motion(self, event=None):
        self.unschedule()
        self.schedule()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(self.delay, self.show)

    def unschedule(self):
        if self.id:
            self.widget.after_cancel(self.id)
            self.id = None

    def show(self):
        if self.tooltip:
            return

        x, y, _, _ = (
            self.widget.bbox("insert") if hasattr(self.widget, "bbox") else (0, 0, 0, 0)
        )
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 25

        self.tooltip = tk.Toplevel(self.widget)
        self.tooltip.wm_overrideredirect(True)
        self.tooltip.wm_geometry(f"+{x}+{y}")

        label = tk.Label(
            self.tooltip,
            text=self.text,
            background="lightyellow",
            relief="solid",
            borderwidth=1,
            font=("TkDefaultFont", "8"),
            wraplength=300,
            justify="left",
        )
        label.pack()

    def hide(self):
        if self.tooltip:
            self.tooltip.destroy()
            self.tooltip = None


from ..constants.parameters import (
    AnalysisModes,
    AnalysisParameters,
    ParameterRanges,
    ParametersDefault,
)


class GeometryCutoffsDialog:
    """Dialog for configuring comprehensive molecular interaction parameters.

    Provides GUI interface for setting parameters for multiple interaction types:

    - **Hydrogen Bonds:** Classical strong interactions (N/O-H···O/N)
    - **Weak Hydrogen Bonds:** C-H···O interactions (important for binding)
    - **Halogen Bonds:** C-X···A interactions (X = Cl, Br, I) with 150° default
    - **π Interactions:** Multiple subtypes including:

      - Hydrogen-π: C-H···π, N-H···π, O-H···π, S-H···π
      - Halogen-π: C-Cl···π, C-Br···π, C-I···π

    Uses tabbed interface to organize parameters by interaction type.
    """

    def __init__(
        self, parent: tk.Tk, current_params: Optional[AnalysisParameters] = None
    ):
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
        self.dialog.geometry("1200x600")
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
        self.covalent_factor = tk.DoubleVar(
            value=ParametersDefault.COVALENT_CUTOFF_FACTOR
        )
        self.hb_distance = tk.DoubleVar(value=ParametersDefault.HB_DISTANCE_CUTOFF)
        self.hb_angle = tk.DoubleVar(value=ParametersDefault.HB_ANGLE_CUTOFF)
        self.da_distance = tk.DoubleVar(value=ParametersDefault.HB_DA_DISTANCE)
        self.whb_distance = tk.DoubleVar(value=ParametersDefault.WHB_DISTANCE_CUTOFF)
        self.whb_angle = tk.DoubleVar(value=ParametersDefault.WHB_ANGLE_CUTOFF)
        self.whb_da_distance = tk.DoubleVar(value=ParametersDefault.WHB_DA_DISTANCE)
        self.xb_distance = tk.DoubleVar(value=ParametersDefault.XB_DISTANCE_CUTOFF)
        self.xb_angle = tk.DoubleVar(value=ParametersDefault.XB_ANGLE_CUTOFF)
        self.pi_distance = tk.DoubleVar(value=ParametersDefault.PI_DISTANCE_CUTOFF)
        self.pi_angle = tk.DoubleVar(value=ParametersDefault.PI_ANGLE_CUTOFF)

        # Initialize π interaction subtype variables (needed for set_parameters)
        self.pi_ccl_distance = None
        self.pi_ccl_angle = None
        self.pi_cbr_distance = None
        self.pi_cbr_angle = None
        self.pi_ci_distance = None
        self.pi_ci_angle = None
        self.pi_ch_distance = None
        self.pi_ch_angle = None
        self.pi_nh_distance = None
        self.pi_nh_angle = None
        self.pi_oh_distance = None
        self.pi_oh_angle = None
        self.pi_sh_distance = None
        self.pi_sh_angle = None

        # New interaction type variables
        self.pi_pi_distance = tk.DoubleVar(
            value=ParametersDefault.PI_PI_DISTANCE_CUTOFF
        )
        self.pi_pi_parallel_angle = tk.DoubleVar(
            value=ParametersDefault.PI_PI_PARALLEL_ANGLE_CUTOFF
        )
        self.pi_pi_tshaped_angle_min = tk.DoubleVar(
            value=ParametersDefault.PI_PI_TSHAPED_ANGLE_MIN
        )
        self.pi_pi_tshaped_angle_max = tk.DoubleVar(
            value=ParametersDefault.PI_PI_TSHAPED_ANGLE_MAX
        )
        self.pi_pi_offset = tk.DoubleVar(value=ParametersDefault.PI_PI_OFFSET_CUTOFF)

        self.carbonyl_distance = tk.DoubleVar(
            value=ParametersDefault.CARBONYL_DISTANCE_CUTOFF
        )
        self.carbonyl_angle_min = tk.DoubleVar(
            value=ParametersDefault.CARBONYL_ANGLE_MIN
        )
        self.carbonyl_angle_max = tk.DoubleVar(
            value=ParametersDefault.CARBONYL_ANGLE_MAX
        )

        self.n_pi_distance = tk.DoubleVar(value=ParametersDefault.N_PI_DISTANCE_CUTOFF)
        self.n_pi_sulfur_distance = tk.DoubleVar(
            value=ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF
        )
        self.n_pi_angle_min = tk.DoubleVar(value=ParametersDefault.N_PI_ANGLE_MIN)
        self.n_pi_angle_max = tk.DoubleVar(value=ParametersDefault.N_PI_ANGLE_MAX)

        # Store parameter values directly
        self._param_values = {}

    def _store_current_values(self):
        """Store current parameter values before switching categories."""
        try:
            self._param_values.update(
                {
                    "hb_distance": (
                        self.hb_distance.get()
                        if hasattr(self, "hb_distance")
                        else ParametersDefault.HB_DISTANCE_CUTOFF
                    ),
                    "hb_angle": (
                        self.hb_angle.get()
                        if hasattr(self, "hb_angle")
                        else ParametersDefault.HB_ANGLE_CUTOFF
                    ),
                    "da_distance": (
                        self.da_distance.get()
                        if hasattr(self, "da_distance")
                        else ParametersDefault.HB_DA_DISTANCE
                    ),
                    "whb_distance": (
                        self.whb_distance.get()
                        if hasattr(self, "whb_distance")
                        else ParametersDefault.WHB_DISTANCE_CUTOFF
                    ),
                    "whb_angle": (
                        self.whb_angle.get()
                        if hasattr(self, "whb_angle")
                        else ParametersDefault.WHB_ANGLE_CUTOFF
                    ),
                    "whb_da_distance": (
                        self.whb_da_distance.get()
                        if hasattr(self, "whb_da_distance")
                        else ParametersDefault.WHB_DA_DISTANCE
                    ),
                    "xb_distance": (
                        self.xb_distance.get()
                        if hasattr(self, "xb_distance")
                        else ParametersDefault.XB_DISTANCE_CUTOFF
                    ),
                    "xb_angle": (
                        self.xb_angle.get()
                        if hasattr(self, "xb_angle")
                        else ParametersDefault.XB_ANGLE_CUTOFF
                    ),
                    "pi_distance": (
                        self.pi_distance.get()
                        if hasattr(self, "pi_distance")
                        else ParametersDefault.PI_DISTANCE_CUTOFF
                    ),
                    "pi_angle": (
                        self.pi_angle.get()
                        if hasattr(self, "pi_angle")
                        else ParametersDefault.PI_ANGLE_CUTOFF
                    ),
                    "covalent_factor": (
                        self.covalent_factor.get()
                        if hasattr(self, "covalent_factor")
                        else ParametersDefault.COVALENT_CUTOFF_FACTOR
                    ),
                    "analysis_mode": (
                        self.analysis_mode.get()
                        if hasattr(self, "analysis_mode")
                        else ParametersDefault.ANALYSIS_MODE
                    ),
                }
            )

            # Store π interaction subtype values if they exist
            pi_vars = [
                "pi_ccl_distance",
                "pi_ccl_angle",
                "pi_cbr_distance",
                "pi_cbr_angle",
                "pi_ci_distance",
                "pi_ci_angle",
                "pi_ch_distance",
                "pi_ch_angle",
                "pi_nh_distance",
                "pi_nh_angle",
                "pi_oh_distance",
                "pi_oh_angle",
                "pi_sh_distance",
                "pi_sh_angle",
                "pi_pi_distance",
                "pi_pi_parallel_angle",
                "pi_pi_tshaped_angle_min",
                "pi_pi_tshaped_angle_max",
                "pi_pi_offset",
                "carbonyl_distance",
                "carbonyl_angle_min",
                "carbonyl_angle_max",
                "n_pi_distance",
                "n_pi_sulfur_distance",
                "n_pi_angle_min",
                "n_pi_angle_max",
            ]

            for var_name in pi_vars:
                if hasattr(self, var_name):
                    var = getattr(self, var_name)
                    if var:
                        self._param_values[var_name] = var.get()

        except tk.TclError:
            pass  # Variables destroyed, ignore

    def _create_widgets(self) -> None:
        """Create and layout all parameter widgets with list selection interface.

        :returns: None
        :rtype: None
        """
        # Main frame with padding
        main_frame = ttk.Frame(self.dialog, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Create container for content
        content_frame = ttk.Frame(main_frame)
        content_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # Create paned window for list and content
        paned = ttk.PanedWindow(content_frame, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        # Left side - Category list
        list_frame = ttk.Frame(paned, relief=tk.GROOVE, borderwidth=1)
        paned.add(list_frame, weight=1)

        # List label
        ttk.Label(
            list_frame, text="Parameter Categories", font=("TkDefaultFont", 10, "bold")
        ).pack(pady=(10, 5))

        # Create listbox for categories
        self.category_listbox = tk.Listbox(list_frame, selectmode=tk.SINGLE, height=10)
        self.category_listbox.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 10))

        # Add categories
        categories = [
            "General Parameters",
            "Hydrogen Bonds",
            "Weak Hydrogen Bonds",
            "Halogen Bonds",
            "π Interactions",
            "π-π Stacking",
            "Carbonyl Interactions",
            "n→π* Interactions",
        ]
        for cat in categories:
            self.category_listbox.insert(tk.END, cat)

        # Bind selection event
        self.category_listbox.bind("<<ListboxSelect>>", self._on_category_selected)

        # Right side - Content area
        self.content_container = ttk.Frame(paned)
        paned.add(self.content_container, weight=3)

        # Create scrollable area for content
        self.content_canvas = tk.Canvas(self.content_container)
        scrollbar = ttk.Scrollbar(
            self.content_container, orient="vertical", command=self.content_canvas.yview
        )
        self.content_frame = ttk.Frame(self.content_canvas)

        self.content_frame.bind(
            "<Configure>",
            lambda e: self.content_canvas.configure(
                scrollregion=self.content_canvas.bbox("all")
            ),
        )

        self.content_canvas.create_window(
            (0, 0), window=self.content_frame, anchor="nw"
        )
        self.content_canvas.configure(yscrollcommand=scrollbar.set)

        # Enable mouse wheel scrolling
        def _on_mousewheel(event):
            self.content_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        self.content_canvas.bind("<MouseWheel>", _on_mousewheel)  # Windows
        self.content_canvas.bind(
            "<Button-4>", lambda e: self.content_canvas.yview_scroll(-1, "units")
        )  # Linux
        self.content_canvas.bind(
            "<Button-5>", lambda e: self.content_canvas.yview_scroll(1, "units")
        )  # Linux

        # Pack canvas and scrollbar
        self.content_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Store current content widget
        self.current_content = None

        # Select first category by default
        self.category_listbox.selection_set(0)
        self._on_category_selected(None)

        # Buttons at bottom
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, side=tk.BOTTOM)

        ttk.Button(button_frame, text="OK", command=self._ok).pack(
            side=tk.RIGHT, padx=(5, 0)
        )

        ttk.Button(button_frame, text="Cancel", command=self._cancel).pack(
            side=tk.RIGHT
        )

        ttk.Button(
            button_frame, text="Reset to Defaults", command=self._set_defaults
        ).pack(side=tk.LEFT, padx=(0, 5))

        ttk.Button(
            button_frame, text="Manage Presets...", command=self._open_preset_manager
        ).pack(side=tk.LEFT, padx=5)

    def _create_general_parameters(self, parent):
        """Create general analysis parameters."""
        group = ttk.LabelFrame(parent, text="General Parameters", padding=10)
        group.pack(fill=tk.X, padx=10, pady=5)

        # Analysis mode
        ttk.Label(group, text="Analysis Mode:").grid(
            row=0, column=0, sticky=tk.W, pady=2
        )
        stored_mode = self._param_values.get(
            "analysis_mode", ParametersDefault.ANALYSIS_MODE
        )
        self.analysis_mode = tk.StringVar(value=stored_mode)
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
        stored_factor = self._param_values.get(
            "covalent_factor", ParametersDefault.COVALENT_CUTOFF_FACTOR
        )
        self.covalent_factor = tk.DoubleVar(value=stored_factor)
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
            try:
                factor_label.config(text=f"{self.covalent_factor.get():.2f}")
            except tk.TclError:
                pass  # Widget destroyed, ignore

        self.covalent_factor.trace("w", update_factor_label)
        update_factor_label()

    def _on_category_selected(self, event):
        """Handle category selection from list."""
        selection = self.category_listbox.curselection()
        if not selection:
            return

        # Store current values before switching
        self._store_current_values()

        # Clear current content
        if self.current_content:
            self.current_content.destroy()

        # Create new content frame
        self.current_content = ttk.Frame(self.content_frame, padding="20")
        self.current_content.pack(fill=tk.BOTH, expand=True)

        # Show appropriate content based on selection
        category_index = selection[0]
        if category_index == 0:  # General Parameters
            self._create_general_parameters(self.current_content)
        elif category_index == 1:  # Hydrogen Bonds
            self._create_hydrogen_bond_parameters(self.current_content)
        elif category_index == 2:  # Weak Hydrogen Bonds
            self._create_weak_hydrogen_bond_parameters(self.current_content)
        elif category_index == 3:  # Halogen Bonds
            self._create_halogen_bond_parameters(self.current_content)
        elif category_index == 4:  # π Interactions
            self._create_pi_interaction_parameters(self.current_content)
        elif category_index == 5:  # π-π Stacking
            self._create_pi_pi_stacking_parameters(self.current_content)
        elif category_index == 6:  # Carbonyl Interactions
            self._create_carbonyl_interaction_parameters(self.current_content)
        elif category_index == 7:  # n→π* Interactions
            self._create_n_pi_interaction_parameters(self.current_content)

        # Reset scroll position
        self.content_canvas.yview_moveto(0)

    def _create_hydrogen_bond_parameters(self, parent):
        """Create hydrogen bond parameter controls."""
        group = ttk.LabelFrame(parent, text="Hydrogen Bond Parameters", padding=10)
        group.pack(fill=tk.X, padx=10, pady=5)

        # H...A distance
        ttk.Label(group, text="H...A Distance (Å):").grid(
            row=0, column=0, sticky=tk.W, pady=2
        )
        stored_dist = self._param_values.get(
            "hb_distance", ParametersDefault.HB_DISTANCE_CUTOFF
        )
        self.hb_distance = tk.DoubleVar(value=stored_dist)
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
            try:
                hb_dist_label.config(text=f"{self.hb_distance.get():.1f}")
            except tk.TclError:
                pass

        self.hb_distance.trace("w", update_hb_dist)
        update_hb_dist()

        # D-H...A angle
        ttk.Label(group, text="D-H...A Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        stored_angle = self._param_values.get(
            "hb_angle", ParametersDefault.HB_ANGLE_CUTOFF
        )
        self.hb_angle = tk.DoubleVar(value=stored_angle)
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
            try:
                hb_angle_label.config(text=f"{self.hb_angle.get():.0f}")
            except tk.TclError:
                pass

        self.hb_angle.trace("w", update_hb_angle)
        update_hb_angle()

        # D...A distance
        ttk.Label(group, text="D...A Distance (Å):").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )
        stored_da = self._param_values.get(
            "da_distance", ParametersDefault.HB_DA_DISTANCE
        )
        self.da_distance = tk.DoubleVar(value=stored_da)
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
            try:
                da_dist_label.config(text=f"{self.da_distance.get():.1f}")
            except tk.TclError:
                pass

        self.da_distance.trace("w", update_da_dist)
        update_da_dist()

    def _create_weak_hydrogen_bond_parameters(self, parent):
        """Create weak hydrogen bond parameter controls for carbon donors."""
        group = ttk.LabelFrame(
            parent, text="Weak Hydrogen Bond Parameters (Carbon Donors)", padding=10
        )
        group.pack(fill=tk.X, padx=10, pady=5)

        # WHB H...A distance
        ttk.Label(group, text="H...A Distance (Å):").grid(
            row=0, column=0, sticky=tk.W, pady=2
        )
        stored_whb_dist = self._param_values.get(
            "whb_distance", ParametersDefault.WHB_DISTANCE_CUTOFF
        )
        self.whb_distance = tk.DoubleVar(value=stored_whb_dist)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.whb_distance,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        whb_dist_label = ttk.Label(group, text="")
        whb_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_whb_dist(*args):
            try:
                whb_dist_label.config(text=f"{self.whb_distance.get():.1f}")
            except tk.TclError:
                pass

        self.whb_distance.trace("w", update_whb_dist)
        update_whb_dist()

        # WHB D-H...A angle
        ttk.Label(group, text="D-H...A Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        stored_whb_angle = self._param_values.get(
            "whb_angle", ParametersDefault.WHB_ANGLE_CUTOFF
        )
        self.whb_angle = tk.DoubleVar(value=stored_whb_angle)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.whb_angle,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        whb_angle_label = ttk.Label(group, text="")
        whb_angle_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_whb_angle(*args):
            try:
                whb_angle_label.config(text=f"{self.whb_angle.get():.0f}")
            except tk.TclError:
                pass

        self.whb_angle.trace("w", update_whb_angle)
        update_whb_angle()

        # WHB D...A distance
        ttk.Label(group, text="D...A Distance (Å):").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )
        stored_whb_da = self._param_values.get(
            "whb_da_distance", ParametersDefault.WHB_DA_DISTANCE
        )
        self.whb_da_distance = tk.DoubleVar(value=stored_whb_da)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.whb_da_distance,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=2, column=1, sticky=tk.W, padx=10, pady=2)

        whb_da_dist_label = ttk.Label(group, text="")
        whb_da_dist_label.grid(row=2, column=2, sticky=tk.W, padx=5, pady=2)

        def update_whb_da_dist(*args):
            try:
                whb_da_dist_label.config(text=f"{self.whb_da_distance.get():.1f}")
            except tk.TclError:
                pass

        self.whb_da_distance.trace("w", update_whb_da_dist)
        update_whb_da_dist()

    def _create_halogen_bond_parameters(self, parent):
        """Create halogen bond parameter controls."""
        group = ttk.LabelFrame(parent, text="Halogen Bond Parameters", padding=10)
        group.pack(fill=tk.X, padx=10, pady=5)

        # X...A distance
        ttk.Label(group, text="X...A Distance (Å):").grid(
            row=0, column=0, sticky=tk.W, pady=2
        )
        stored_xb_dist = self._param_values.get(
            "xb_distance", ParametersDefault.XB_DISTANCE_CUTOFF
        )
        self.xb_distance = tk.DoubleVar(value=stored_xb_dist)
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
            try:
                xb_dist_label.config(text=f"{self.xb_distance.get():.1f}")
            except tk.TclError:
                pass

        self.xb_distance.trace("w", update_xb_dist)
        update_xb_dist()

        # C-X...A angle
        ttk.Label(group, text="C-X...A Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        stored_xb_angle = self._param_values.get(
            "xb_angle", ParametersDefault.XB_ANGLE_CUTOFF
        )
        self.xb_angle = tk.DoubleVar(value=stored_xb_angle)
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
            try:
                xb_angle_label.config(text=f"{self.xb_angle.get():.0f}")
            except tk.TclError:
                pass

        self.xb_angle.trace("w", update_xb_angle)
        update_xb_angle()

    def _create_pi_interaction_parameters(self, parent):
        """Create π interaction parameter controls."""
        # General π interaction parameters
        general_group = ttk.LabelFrame(
            parent, text="General π Interaction Parameters", padding=10
        )
        general_group.pack(fill=tk.X, padx=10, pady=5)

        # H...π distance
        ttk.Label(general_group, text="H...π Distance (Å):").grid(
            row=0, column=0, sticky=tk.W, pady=2
        )
        stored_pi_dist = self._param_values.get(
            "pi_distance", ParametersDefault.PI_DISTANCE_CUTOFF
        )
        self.pi_distance = tk.DoubleVar(value=stored_pi_dist)
        ttk.Scale(
            general_group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.pi_distance,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        pi_dist_label = ttk.Label(general_group, text="")
        pi_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_dist(*args):
            try:
                pi_dist_label.config(text=f"{self.pi_distance.get():.1f}")
            except tk.TclError:
                pass

        self.pi_distance.trace("w", update_pi_dist)
        update_pi_dist()

        # D-H...π angle
        ttk.Label(general_group, text="D-H...π Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        stored_pi_angle = self._param_values.get(
            "pi_angle", ParametersDefault.PI_ANGLE_CUTOFF
        )
        self.pi_angle = tk.DoubleVar(value=stored_pi_angle)
        ttk.Scale(
            general_group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.pi_angle,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        pi_angle_label = ttk.Label(general_group, text="")
        pi_angle_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_angle(*args):
            try:
                pi_angle_label.config(text=f"{self.pi_angle.get():.0f}")
            except tk.TclError:
                pass

        self.pi_angle.trace("w", update_pi_angle)
        update_pi_angle()

        # π interaction subtype parameters
        subtypes_group = ttk.LabelFrame(
            parent, text="π Interaction Subtype Parameters", padding=10
        )
        subtypes_group.pack(fill=tk.X, padx=10, pady=5)

        # Helper function to create parameter pair
        def create_parameter_pair(
            parent_frame,
            row,
            label_text,
            dist_var_name,
            angle_var_name,
            dist_default,
            angle_default,
        ):
            # Get stored values or use defaults
            stored_dist = self._param_values.get(dist_var_name, dist_default)
            stored_angle = self._param_values.get(angle_var_name, angle_default)

            # Distance parameter
            ttk.Label(parent_frame, text=f"{label_text} Distance (Å):").grid(
                row=row, column=0, sticky=tk.W, pady=2
            )
            dist_var = tk.DoubleVar(value=stored_dist)
            setattr(self, dist_var_name, dist_var)
            ttk.Scale(
                parent_frame,
                from_=ParameterRanges.MIN_DISTANCE,
                to=ParameterRanges.MAX_DISTANCE,
                variable=dist_var,
                orient=tk.HORIZONTAL,
                length=150,
            ).grid(row=row, column=1, sticky=tk.W, padx=5, pady=2)

            dist_label = ttk.Label(parent_frame, text="")
            dist_label.grid(row=row, column=2, sticky=tk.W, padx=5, pady=2)

            # Angle parameter
            ttk.Label(parent_frame, text=f"{label_text} Angle (°):").grid(
                row=row, column=3, sticky=tk.W, pady=2, padx=(20, 0)
            )
            angle_var = tk.DoubleVar(value=stored_angle)
            setattr(self, angle_var_name, angle_var)
            ttk.Scale(
                parent_frame,
                from_=ParameterRanges.MIN_ANGLE,
                to=ParameterRanges.MAX_ANGLE,
                variable=angle_var,
                orient=tk.HORIZONTAL,
                length=150,
            ).grid(row=row, column=4, sticky=tk.W, padx=5, pady=2)

            angle_label = ttk.Label(parent_frame, text="")
            angle_label.grid(row=row, column=5, sticky=tk.W, padx=5, pady=2)

            # Update functions
            def update_dist(*args):
                try:
                    dist_label.config(text=f"{dist_var.get():.1f}")
                except tk.TclError:
                    pass

            def update_angle(*args):
                try:
                    angle_label.config(text=f"{angle_var.get():.0f}")
                except tk.TclError:
                    pass

            dist_var.trace("w", update_dist)
            angle_var.trace("w", update_angle)
            update_dist()
            update_angle()

        # Create all subtype parameters
        create_parameter_pair(
            subtypes_group,
            0,
            "C-Cl...π",
            "pi_ccl_distance",
            "pi_ccl_angle",
            ParametersDefault.PI_CCL_DISTANCE_CUTOFF,
            ParametersDefault.PI_CCL_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            1,
            "C-Br...π",
            "pi_cbr_distance",
            "pi_cbr_angle",
            ParametersDefault.PI_CBR_DISTANCE_CUTOFF,
            ParametersDefault.PI_CBR_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            2,
            "C-I...π",
            "pi_ci_distance",
            "pi_ci_angle",
            ParametersDefault.PI_CI_DISTANCE_CUTOFF,
            ParametersDefault.PI_CI_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            3,
            "C-H...π",
            "pi_ch_distance",
            "pi_ch_angle",
            ParametersDefault.PI_CH_DISTANCE_CUTOFF,
            ParametersDefault.PI_CH_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            4,
            "N-H...π",
            "pi_nh_distance",
            "pi_nh_angle",
            ParametersDefault.PI_NH_DISTANCE_CUTOFF,
            ParametersDefault.PI_NH_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            5,
            "O-H...π",
            "pi_oh_distance",
            "pi_oh_angle",
            ParametersDefault.PI_OH_DISTANCE_CUTOFF,
            ParametersDefault.PI_OH_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            6,
            "S-H...π",
            "pi_sh_distance",
            "pi_sh_angle",
            ParametersDefault.PI_SH_DISTANCE_CUTOFF,
            ParametersDefault.PI_SH_ANGLE_CUTOFF,
        )

    def _create_pi_pi_stacking_parameters(self, parent):
        """Create π-π stacking parameter controls."""
        group = ttk.LabelFrame(parent, text="π-π Stacking Parameters", padding=10)
        group.pack(fill=tk.X, padx=10, pady=5)

        # Distance cutoff
        dist_label = ttk.Label(group, text="Distance Cutoff (Å):")
        dist_label.grid(row=0, column=0, sticky=tk.W, pady=2)
        ToolTip(
            dist_label,
            "Maximum distance between aromatic ring centroids for π-π stacking interactions.\n"
            "Typical range: 3.0-4.5 Å\n"
            "• Parallel stacking: ~3.5-4.0 Å\n"
            "• T-shaped stacking: ~3.5-4.5 Å\n"
            "Based on McGaughey et al. (1998) and crystallographic data.",
        )

        self.pi_pi_distance.set(
            self._param_values.get(
                "pi_pi_distance", ParametersDefault.PI_PI_DISTANCE_CUTOFF
            )
        )
        pi_pi_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.pi_pi_distance,
            orient=tk.HORIZONTAL,
            length=200,
        )
        pi_pi_scale.grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)
        ToolTip(
            pi_pi_scale,
            "Adjust the maximum distance between aromatic ring centroids.\n"
            "Lower values = stricter detection, higher values = more permissive.",
        )

        pi_pi_dist_label = ttk.Label(group, text="")
        pi_pi_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_pi_dist(*args):
            try:
                pi_pi_dist_label.config(text=f"{self.pi_pi_distance.get():.1f}")
            except tk.TclError:
                pass

        self.pi_pi_distance.trace("w", update_pi_pi_dist)
        update_pi_pi_dist()

        # Parallel angle cutoff
        parallel_label = ttk.Label(group, text="Parallel Angle Cutoff (°):")
        parallel_label.grid(row=1, column=0, sticky=tk.W, pady=2)
        ToolTip(
            parallel_label,
            "Maximum angle between aromatic ring planes for parallel π-π stacking.\n"
            "• 0° = perfectly parallel rings\n"
            "• Typical cutoff: 20-30°\n"
            "• Values >30° indicate T-shaped or edge-to-face interactions\n"
            "Based on Hunter & Sanders (1990) π-π interaction classification.",
        )

        self.pi_pi_parallel_angle.set(
            self._param_values.get(
                "pi_pi_parallel_angle", ParametersDefault.PI_PI_PARALLEL_ANGLE_CUTOFF
            )
        )
        parallel_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.pi_pi_parallel_angle,
            orient=tk.HORIZONTAL,
            length=200,
        )
        parallel_scale.grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)
        ToolTip(
            parallel_scale,
            "Adjust angle tolerance for parallel stacking.\n"
            "Lower values = more strict parallel geometry required.",
        )

        pi_pi_parallel_label = ttk.Label(group, text="")
        pi_pi_parallel_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_pi_parallel(*args):
            try:
                pi_pi_parallel_label.config(
                    text=f"{self.pi_pi_parallel_angle.get():.0f}"
                )
            except tk.TclError:
                pass

        self.pi_pi_parallel_angle.trace("w", update_pi_pi_parallel)
        update_pi_pi_parallel()

        # T-shaped angle minimum
        tmin_label = ttk.Label(group, text="T-shaped Angle Min (°):")
        tmin_label.grid(row=2, column=0, sticky=tk.W, pady=2)
        ToolTip(
            tmin_label,
            "Minimum angle between ring planes for T-shaped (edge-to-face) stacking.\n"
            "• Typical range: 60-90°\n"
            "• T-shaped geometry involves one ring edge interacting with another ring face\n"
            "• Complementary to parallel stacking interactions\n"
            "Based on Burley & Petsko (1985) aromatic interaction studies.",
        )

        self.pi_pi_tshaped_angle_min.set(
            self._param_values.get(
                "pi_pi_tshaped_angle_min", ParametersDefault.PI_PI_TSHAPED_ANGLE_MIN
            )
        )
        tmin_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.pi_pi_tshaped_angle_min,
            orient=tk.HORIZONTAL,
            length=200,
        )
        tmin_scale.grid(row=2, column=1, sticky=tk.W, padx=10, pady=2)
        ToolTip(
            tmin_scale,
            "Adjust minimum angle for T-shaped interactions.\n"
            "Higher values = more perpendicular geometry required.",
        )

        pi_pi_tmin_label = ttk.Label(group, text="")
        pi_pi_tmin_label.grid(row=2, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_pi_tmin(*args):
            try:
                pi_pi_tmin_label.config(
                    text=f"{self.pi_pi_tshaped_angle_min.get():.0f}"
                )
            except tk.TclError:
                pass

        self.pi_pi_tshaped_angle_min.trace("w", update_pi_pi_tmin)
        update_pi_pi_tmin()

        # T-shaped angle maximum
        ttk.Label(group, text="T-shaped Angle Max (°):").grid(
            row=3, column=0, sticky=tk.W, pady=2
        )
        self.pi_pi_tshaped_angle_max.set(
            self._param_values.get(
                "pi_pi_tshaped_angle_max", ParametersDefault.PI_PI_TSHAPED_ANGLE_MAX
            )
        )
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.pi_pi_tshaped_angle_max,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=3, column=1, sticky=tk.W, padx=10, pady=2)

        pi_pi_tmax_label = ttk.Label(group, text="")
        pi_pi_tmax_label.grid(row=3, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_pi_tmax(*args):
            try:
                pi_pi_tmax_label.config(
                    text=f"{self.pi_pi_tshaped_angle_max.get():.0f}"
                )
            except tk.TclError:
                pass

        self.pi_pi_tshaped_angle_max.trace("w", update_pi_pi_tmax)
        update_pi_pi_tmax()

        # Offset cutoff
        ttk.Label(group, text="Offset Cutoff (Å):").grid(
            row=4, column=0, sticky=tk.W, pady=2
        )
        self.pi_pi_offset.set(
            self._param_values.get(
                "pi_pi_offset", ParametersDefault.PI_PI_OFFSET_CUTOFF
            )
        )
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.pi_pi_offset,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=4, column=1, sticky=tk.W, padx=10, pady=2)

        pi_pi_offset_label = ttk.Label(group, text="")
        pi_pi_offset_label.grid(row=4, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_pi_offset(*args):
            try:
                pi_pi_offset_label.config(text=f"{self.pi_pi_offset.get():.1f}")
            except tk.TclError:
                pass

        self.pi_pi_offset.trace("w", update_pi_pi_offset)
        update_pi_pi_offset()

    def _create_carbonyl_interaction_parameters(self, parent):
        """Create carbonyl interaction parameter controls."""
        group = ttk.LabelFrame(
            parent, text="Carbonyl n→π* Interaction Parameters", padding=10
        )
        group.pack(fill=tk.X, padx=10, pady=5)

        # Distance cutoff
        carb_dist_label = ttk.Label(group, text="O···C Distance Cutoff (Å):")
        carb_dist_label.grid(row=0, column=0, sticky=tk.W, pady=2)
        ToolTip(
            carb_dist_label,
            "Maximum distance between lone pair donor oxygen and carbonyl carbon.\n"
            "• Typical range: 2.8-3.2 Å\n"
            "• Represents n→π* orbital interaction\n"
            "• Distance based on crystallographic surveys of protein structures\n"
            "• Critical for protein backbone stability and folding",
        )

        self.carbonyl_distance.set(
            self._param_values.get(
                "carbonyl_distance", ParametersDefault.CARBONYL_DISTANCE_CUTOFF
            )
        )
        carb_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.carbonyl_distance,
            orient=tk.HORIZONTAL,
            length=200,
        )
        carb_scale.grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)
        ToolTip(
            carb_scale,
            "Adjust O···C distance cutoff for carbonyl interactions.\n"
            "Based on van der Waals radii and quantum mechanical calculations.",
        )

        carbonyl_dist_label = ttk.Label(group, text="")
        carbonyl_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_carbonyl_dist(*args):
            try:
                carbonyl_dist_label.config(text=f"{self.carbonyl_distance.get():.1f}")
            except tk.TclError:
                pass

        self.carbonyl_distance.trace("w", update_carbonyl_dist)
        update_carbonyl_dist()

        # Angle minimum
        angle_min_label = ttk.Label(group, text="Bürgi-Dunitz Angle Min (°):")
        angle_min_label.grid(row=1, column=0, sticky=tk.W, pady=2)
        ToolTip(
            angle_min_label,
            "Minimum Bürgi-Dunitz approach angle for nucleophilic attack on carbonyl.\n"
            "• Optimal angle: ~107° (tetrahedral trajectory)\n"
            "• Range: 95-125° based on crystal structures\n"
            "• Named after Bürgi & Dunitz (1983) crystallographic studies\n"
            "• Represents stereoelectronically favored approach geometry",
        )

        self.carbonyl_angle_min.set(
            self._param_values.get(
                "carbonyl_angle_min", ParametersDefault.CARBONYL_ANGLE_MIN
            )
        )
        angle_min_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.carbonyl_angle_min,
            orient=tk.HORIZONTAL,
            length=200,
        )
        angle_min_scale.grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)
        ToolTip(
            angle_min_scale,
            "Adjust minimum Bürgi-Dunitz angle.\n"
            "Based on ab initio calculations and crystallographic analysis.",
        )

        carbonyl_min_label = ttk.Label(group, text="")
        carbonyl_min_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_carbonyl_min(*args):
            try:
                carbonyl_min_label.config(text=f"{self.carbonyl_angle_min.get():.0f}")
            except tk.TclError:
                pass

        self.carbonyl_angle_min.trace("w", update_carbonyl_min)
        update_carbonyl_min()

        # Angle maximum
        ttk.Label(group, text="Bürgi-Dunitz Angle Max (°):").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )
        self.carbonyl_angle_max.set(
            self._param_values.get(
                "carbonyl_angle_max", ParametersDefault.CARBONYL_ANGLE_MAX
            )
        )
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.carbonyl_angle_max,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=2, column=1, sticky=tk.W, padx=10, pady=2)

        carbonyl_max_label = ttk.Label(group, text="")
        carbonyl_max_label.grid(row=2, column=2, sticky=tk.W, padx=5, pady=2)

        def update_carbonyl_max(*args):
            try:
                carbonyl_max_label.config(text=f"{self.carbonyl_angle_max.get():.0f}")
            except tk.TclError:
                pass

        self.carbonyl_angle_max.trace("w", update_carbonyl_max)
        update_carbonyl_max()

    def _create_n_pi_interaction_parameters(self, parent):
        """Create n→π* interaction parameter controls."""
        group = ttk.LabelFrame(parent, text="n→π* Interaction Parameters", padding=10)
        group.pack(fill=tk.X, padx=10, pady=5)

        # Distance cutoff
        n_pi_dist_label = ttk.Label(group, text="Distance Cutoff (Å):")
        n_pi_dist_label.grid(row=0, column=0, sticky=tk.W, pady=2)
        ToolTip(
            n_pi_dist_label,
            "Maximum distance between lone pair donor atom and aromatic ring centroid.\n"
            "• Typical range: 3.0-4.0 Å for O/N donors\n"
            "• Represents n→π* orbital interaction\n"
            "• Common in protein-ligand and protein-protein interactions\n"
            "• Important for molecular recognition and binding affinity",
        )

        self.n_pi_distance.set(
            self._param_values.get(
                "n_pi_distance", ParametersDefault.N_PI_DISTANCE_CUTOFF
            )
        )
        n_pi_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.n_pi_distance,
            orient=tk.HORIZONTAL,
            length=200,
        )
        n_pi_scale.grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)
        ToolTip(
            n_pi_scale,
            "Adjust distance cutoff for n→π* interactions.\n"
            "Based on computational studies and structural databases.",
        )

        n_pi_dist_label = ttk.Label(group, text="")
        n_pi_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_n_pi_dist(*args):
            try:
                n_pi_dist_label.config(text=f"{self.n_pi_distance.get():.1f}")
            except tk.TclError:
                pass

        self.n_pi_distance.trace("w", update_n_pi_dist)
        update_n_pi_dist()

        # Sulfur distance cutoff
        sulfur_label = ttk.Label(group, text="Sulfur Distance Cutoff (Å):")
        sulfur_label.grid(row=1, column=0, sticky=tk.W, pady=2)
        ToolTip(
            sulfur_label,
            "Maximum distance for sulfur n→π* interactions (S···aromatic ring).\n"
            "• Typically larger than O/N due to sulfur's larger atomic radius\n"
            "• Range: 3.5-4.5 Å\n"
            "• Important in cysteine-aromatic interactions\n"
            "• Weaker than O/N interactions but geometrically significant",
        )

        self.n_pi_sulfur_distance.set(
            self._param_values.get(
                "n_pi_sulfur_distance", ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF
            )
        )
        sulfur_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self.n_pi_sulfur_distance,
            orient=tk.HORIZONTAL,
            length=200,
        )
        sulfur_scale.grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)
        ToolTip(
            sulfur_scale,
            "Adjust sulfur-aromatic distance cutoff.\n"
            "Larger values account for sulfur's extended electron cloud.",
        )

        n_pi_sulfur_label = ttk.Label(group, text="")
        n_pi_sulfur_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_n_pi_sulfur(*args):
            try:
                n_pi_sulfur_label.config(text=f"{self.n_pi_sulfur_distance.get():.1f}")
            except tk.TclError:
                pass

        self.n_pi_sulfur_distance.trace("w", update_n_pi_sulfur)
        update_n_pi_sulfur()

        # Angle minimum
        ttk.Label(group, text="Angle Min (°):").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )
        self.n_pi_angle_min.set(
            self._param_values.get("n_pi_angle_min", ParametersDefault.N_PI_ANGLE_MIN)
        )
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.n_pi_angle_min,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=2, column=1, sticky=tk.W, padx=10, pady=2)

        n_pi_min_label = ttk.Label(group, text="")
        n_pi_min_label.grid(row=2, column=2, sticky=tk.W, padx=5, pady=2)

        def update_n_pi_min(*args):
            try:
                n_pi_min_label.config(text=f"{self.n_pi_angle_min.get():.0f}")
            except tk.TclError:
                pass

        self.n_pi_angle_min.trace("w", update_n_pi_min)
        update_n_pi_min()

        # Angle maximum
        ttk.Label(group, text="Angle Max (°):").grid(
            row=3, column=0, sticky=tk.W, pady=2
        )
        self.n_pi_angle_max.set(
            self._param_values.get("n_pi_angle_max", ParametersDefault.N_PI_ANGLE_MAX)
        )
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self.n_pi_angle_max,
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=3, column=1, sticky=tk.W, padx=10, pady=2)

        n_pi_max_label = ttk.Label(group, text="")
        n_pi_max_label.grid(row=3, column=2, sticky=tk.W, padx=5, pady=2)

        def update_n_pi_max(*args):
            try:
                n_pi_max_label.config(text=f"{self.n_pi_angle_max.get():.0f}")
            except tk.TclError:
                pass

        self.n_pi_angle_max.trace("w", update_n_pi_max)
        update_n_pi_max()

    def get_parameters(self) -> AnalysisParameters:
        """Get current parameter values.

        :returns: Current analysis parameters
        :rtype: AnalysisParameters
        """
        # Store current values first
        self._store_current_values()

        # Get values from stored parameters or current variables
        def get_value(var_name, default_value):
            if var_name in self._param_values:
                return self._param_values[var_name]
            elif hasattr(self, var_name):
                var = getattr(self, var_name)
                if var:
                    try:
                        return var.get()
                    except tk.TclError:
                        pass
            return default_value

        return AnalysisParameters(
            hb_distance_cutoff=get_value(
                "hb_distance", ParametersDefault.HB_DISTANCE_CUTOFF
            ),
            hb_angle_cutoff=get_value("hb_angle", ParametersDefault.HB_ANGLE_CUTOFF),
            hb_donor_acceptor_cutoff=get_value(
                "da_distance", ParametersDefault.HB_DA_DISTANCE
            ),
            whb_distance_cutoff=get_value(
                "whb_distance", ParametersDefault.WHB_DISTANCE_CUTOFF
            ),
            whb_angle_cutoff=get_value("whb_angle", ParametersDefault.WHB_ANGLE_CUTOFF),
            whb_donor_acceptor_cutoff=get_value(
                "whb_da_distance", ParametersDefault.WHB_DA_DISTANCE
            ),
            xb_distance_cutoff=get_value(
                "xb_distance", ParametersDefault.XB_DISTANCE_CUTOFF
            ),
            xb_angle_cutoff=get_value("xb_angle", ParametersDefault.XB_ANGLE_CUTOFF),
            pi_distance_cutoff=get_value(
                "pi_distance", ParametersDefault.PI_DISTANCE_CUTOFF
            ),
            pi_angle_cutoff=get_value("pi_angle", ParametersDefault.PI_ANGLE_CUTOFF),
            # π interaction subtype parameters
            pi_ccl_distance_cutoff=get_value(
                "pi_ccl_distance", ParametersDefault.PI_CCL_DISTANCE_CUTOFF
            ),
            pi_ccl_angle_cutoff=get_value(
                "pi_ccl_angle", ParametersDefault.PI_CCL_ANGLE_CUTOFF
            ),
            pi_cbr_distance_cutoff=get_value(
                "pi_cbr_distance", ParametersDefault.PI_CBR_DISTANCE_CUTOFF
            ),
            pi_cbr_angle_cutoff=get_value(
                "pi_cbr_angle", ParametersDefault.PI_CBR_ANGLE_CUTOFF
            ),
            pi_ci_distance_cutoff=get_value(
                "pi_ci_distance", ParametersDefault.PI_CI_DISTANCE_CUTOFF
            ),
            pi_ci_angle_cutoff=get_value(
                "pi_ci_angle", ParametersDefault.PI_CI_ANGLE_CUTOFF
            ),
            pi_ch_distance_cutoff=get_value(
                "pi_ch_distance", ParametersDefault.PI_CH_DISTANCE_CUTOFF
            ),
            pi_ch_angle_cutoff=get_value(
                "pi_ch_angle", ParametersDefault.PI_CH_ANGLE_CUTOFF
            ),
            pi_nh_distance_cutoff=get_value(
                "pi_nh_distance", ParametersDefault.PI_NH_DISTANCE_CUTOFF
            ),
            pi_nh_angle_cutoff=get_value(
                "pi_nh_angle", ParametersDefault.PI_NH_ANGLE_CUTOFF
            ),
            pi_oh_distance_cutoff=get_value(
                "pi_oh_distance", ParametersDefault.PI_OH_DISTANCE_CUTOFF
            ),
            pi_oh_angle_cutoff=get_value(
                "pi_oh_angle", ParametersDefault.PI_OH_ANGLE_CUTOFF
            ),
            pi_sh_distance_cutoff=get_value(
                "pi_sh_distance", ParametersDefault.PI_SH_DISTANCE_CUTOFF
            ),
            pi_sh_angle_cutoff=get_value(
                "pi_sh_angle", ParametersDefault.PI_SH_ANGLE_CUTOFF
            ),
            # New interaction parameters
            pi_pi_distance_cutoff=get_value(
                "pi_pi_distance", ParametersDefault.PI_PI_DISTANCE_CUTOFF
            ),
            pi_pi_parallel_angle_cutoff=get_value(
                "pi_pi_parallel_angle", ParametersDefault.PI_PI_PARALLEL_ANGLE_CUTOFF
            ),
            pi_pi_tshaped_angle_min=get_value(
                "pi_pi_tshaped_angle_min", ParametersDefault.PI_PI_TSHAPED_ANGLE_MIN
            ),
            pi_pi_tshaped_angle_max=get_value(
                "pi_pi_tshaped_angle_max", ParametersDefault.PI_PI_TSHAPED_ANGLE_MAX
            ),
            pi_pi_offset_cutoff=get_value(
                "pi_pi_offset", ParametersDefault.PI_PI_OFFSET_CUTOFF
            ),
            carbonyl_distance_cutoff=get_value(
                "carbonyl_distance", ParametersDefault.CARBONYL_DISTANCE_CUTOFF
            ),
            carbonyl_angle_min=get_value(
                "carbonyl_angle_min", ParametersDefault.CARBONYL_ANGLE_MIN
            ),
            carbonyl_angle_max=get_value(
                "carbonyl_angle_max", ParametersDefault.CARBONYL_ANGLE_MAX
            ),
            n_pi_distance_cutoff=get_value(
                "n_pi_distance", ParametersDefault.N_PI_DISTANCE_CUTOFF
            ),
            n_pi_sulfur_distance_cutoff=get_value(
                "n_pi_sulfur_distance", ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF
            ),
            n_pi_angle_min=get_value(
                "n_pi_angle_min", ParametersDefault.N_PI_ANGLE_MIN
            ),
            n_pi_angle_max=get_value(
                "n_pi_angle_max", ParametersDefault.N_PI_ANGLE_MAX
            ),
            covalent_cutoff_factor=get_value(
                "covalent_factor", ParametersDefault.COVALENT_CUTOFF_FACTOR
            ),
            analysis_mode=get_value("analysis_mode", ParametersDefault.ANALYSIS_MODE),
        )

    def set_parameters(self, params: AnalysisParameters) -> None:
        """Set parameter values from AnalysisParameters object.

        :param params: Analysis parameters to set
        :type params: AnalysisParameters
        """
        # Store values directly in our parameter storage
        self._param_values.update(
            {
                "hb_distance": params.hb_distance_cutoff,
                "hb_angle": params.hb_angle_cutoff,
                "da_distance": params.hb_donor_acceptor_cutoff,
                "whb_distance": params.whb_distance_cutoff,
                "whb_angle": params.whb_angle_cutoff,
                "whb_da_distance": params.whb_donor_acceptor_cutoff,
                "xb_distance": params.xb_distance_cutoff,
                "xb_angle": params.xb_angle_cutoff,
                "pi_distance": params.pi_distance_cutoff,
                "pi_angle": params.pi_angle_cutoff,
                "pi_ccl_distance": params.pi_ccl_distance_cutoff,
                "pi_ccl_angle": params.pi_ccl_angle_cutoff,
                "pi_cbr_distance": params.pi_cbr_distance_cutoff,
                "pi_cbr_angle": params.pi_cbr_angle_cutoff,
                "pi_ci_distance": params.pi_ci_distance_cutoff,
                "pi_ci_angle": params.pi_ci_angle_cutoff,
                "pi_ch_distance": params.pi_ch_distance_cutoff,
                "pi_ch_angle": params.pi_ch_angle_cutoff,
                "pi_nh_distance": params.pi_nh_distance_cutoff,
                "pi_nh_angle": params.pi_nh_angle_cutoff,
                "pi_oh_distance": params.pi_oh_distance_cutoff,
                "pi_oh_angle": params.pi_oh_angle_cutoff,
                "pi_sh_distance": params.pi_sh_distance_cutoff,
                "pi_sh_angle": params.pi_sh_angle_cutoff,
                # New interaction parameters
                "pi_pi_distance": params.pi_pi_distance_cutoff,
                "pi_pi_parallel_angle": params.pi_pi_parallel_angle_cutoff,
                "pi_pi_tshaped_angle_min": params.pi_pi_tshaped_angle_min,
                "pi_pi_tshaped_angle_max": params.pi_pi_tshaped_angle_max,
                "pi_pi_offset": params.pi_pi_offset_cutoff,
                "carbonyl_distance": params.carbonyl_distance_cutoff,
                "carbonyl_angle_min": params.carbonyl_angle_min,
                "carbonyl_angle_max": params.carbonyl_angle_max,
                "n_pi_distance": params.n_pi_distance_cutoff,
                "n_pi_sulfur_distance": params.n_pi_sulfur_distance_cutoff,
                "n_pi_angle_min": params.n_pi_angle_min,
                "n_pi_angle_max": params.n_pi_angle_max,
                "covalent_factor": params.covalent_cutoff_factor,
                "analysis_mode": params.analysis_mode,
            }
        )

        # Set values on existing variables if they exist
        def safe_set(var_name, value):
            if hasattr(self, var_name):
                var = getattr(self, var_name)
                if var:
                    try:
                        var.set(value)
                    except tk.TclError:
                        pass  # Variable destroyed, ignore

        safe_set("hb_distance", params.hb_distance_cutoff)
        safe_set("hb_angle", params.hb_angle_cutoff)
        safe_set("da_distance", params.hb_donor_acceptor_cutoff)
        safe_set("whb_distance", params.whb_distance_cutoff)
        safe_set("whb_angle", params.whb_angle_cutoff)
        safe_set("whb_da_distance", params.whb_donor_acceptor_cutoff)
        safe_set("xb_distance", params.xb_distance_cutoff)
        safe_set("xb_angle", params.xb_angle_cutoff)
        safe_set("pi_distance", params.pi_distance_cutoff)
        safe_set("pi_angle", params.pi_angle_cutoff)
        safe_set("pi_ccl_distance", params.pi_ccl_distance_cutoff)
        safe_set("pi_ccl_angle", params.pi_ccl_angle_cutoff)
        safe_set("pi_cbr_distance", params.pi_cbr_distance_cutoff)
        safe_set("pi_cbr_angle", params.pi_cbr_angle_cutoff)
        safe_set("pi_ci_distance", params.pi_ci_distance_cutoff)
        safe_set("pi_ci_angle", params.pi_ci_angle_cutoff)
        safe_set("pi_ch_distance", params.pi_ch_distance_cutoff)
        safe_set("pi_ch_angle", params.pi_ch_angle_cutoff)
        safe_set("pi_nh_distance", params.pi_nh_distance_cutoff)
        safe_set("pi_nh_angle", params.pi_nh_angle_cutoff)
        safe_set("pi_oh_distance", params.pi_oh_distance_cutoff)
        safe_set("pi_oh_angle", params.pi_oh_angle_cutoff)
        safe_set("pi_sh_distance", params.pi_sh_distance_cutoff)
        safe_set("pi_sh_angle", params.pi_sh_angle_cutoff)
        # New interaction parameters
        safe_set("pi_pi_distance", params.pi_pi_distance_cutoff)
        safe_set("pi_pi_parallel_angle", params.pi_pi_parallel_angle_cutoff)
        safe_set("pi_pi_tshaped_angle_min", params.pi_pi_tshaped_angle_min)
        safe_set("pi_pi_tshaped_angle_max", params.pi_pi_tshaped_angle_max)
        safe_set("pi_pi_offset", params.pi_pi_offset_cutoff)
        safe_set("carbonyl_distance", params.carbonyl_distance_cutoff)
        safe_set("carbonyl_angle_min", params.carbonyl_angle_min)
        safe_set("carbonyl_angle_max", params.carbonyl_angle_max)
        safe_set("n_pi_distance", params.n_pi_distance_cutoff)
        safe_set("n_pi_sulfur_distance", params.n_pi_sulfur_distance_cutoff)
        safe_set("n_pi_angle_min", params.n_pi_angle_min)
        safe_set("n_pi_angle_max", params.n_pi_angle_max)
        safe_set("covalent_factor", params.covalent_cutoff_factor)
        safe_set("analysis_mode", params.analysis_mode)

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
            self.hb_distance.set(
                hb.get("h_a_distance_cutoff", ParametersDefault.HB_DISTANCE_CUTOFF)
            )
            self.hb_angle.set(
                hb.get("dha_angle_cutoff", ParametersDefault.HB_ANGLE_CUTOFF)
            )
            self.da_distance.set(
                hb.get("d_a_distance_cutoff", ParametersDefault.HB_DA_DISTANCE)
            )

        # Apply halogen bond parameters
        if "halogen_bonds" in params:
            xb = params["halogen_bonds"]
            self.xb_distance.set(
                xb.get("x_a_distance_cutoff", ParametersDefault.XB_DISTANCE_CUTOFF)
            )
            self.xb_angle.set(
                xb.get("dxa_angle_cutoff", ParametersDefault.XB_ANGLE_CUTOFF)
            )

        # Apply π interaction parameters
        if "pi_interactions" in params:
            pi = params["pi_interactions"]
            self.pi_distance.set(
                pi.get("h_pi_distance_cutoff", ParametersDefault.PI_DISTANCE_CUTOFF)
            )
            self.pi_angle.set(
                pi.get("dh_pi_angle_cutoff", ParametersDefault.PI_ANGLE_CUTOFF)
            )

            # Apply π interaction subtype parameters
            self.pi_ccl_distance.set(
                pi.get(
                    "ccl_pi_distance_cutoff", ParametersDefault.PI_CCL_DISTANCE_CUTOFF
                )
            )
            self.pi_ccl_angle.set(
                pi.get("ccl_pi_angle_cutoff", ParametersDefault.PI_CCL_ANGLE_CUTOFF)
            )
            self.pi_cbr_distance.set(
                pi.get(
                    "cbr_pi_distance_cutoff", ParametersDefault.PI_CBR_DISTANCE_CUTOFF
                )
            )
            self.pi_cbr_angle.set(
                pi.get("cbr_pi_angle_cutoff", ParametersDefault.PI_CBR_ANGLE_CUTOFF)
            )
            self.pi_ci_distance.set(
                pi.get("ci_pi_distance_cutoff", ParametersDefault.PI_CI_DISTANCE_CUTOFF)
            )
            self.pi_ci_angle.set(
                pi.get("ci_pi_angle_cutoff", ParametersDefault.PI_CI_ANGLE_CUTOFF)
            )
            self.pi_ch_distance.set(
                pi.get("ch_pi_distance_cutoff", ParametersDefault.PI_CH_DISTANCE_CUTOFF)
            )
            self.pi_ch_angle.set(
                pi.get("ch_pi_angle_cutoff", ParametersDefault.PI_CH_ANGLE_CUTOFF)
            )
            self.pi_nh_distance.set(
                pi.get("nh_pi_distance_cutoff", ParametersDefault.PI_NH_DISTANCE_CUTOFF)
            )
            self.pi_nh_angle.set(
                pi.get("nh_pi_angle_cutoff", ParametersDefault.PI_NH_ANGLE_CUTOFF)
            )
            self.pi_oh_distance.set(
                pi.get("oh_pi_distance_cutoff", ParametersDefault.PI_OH_DISTANCE_CUTOFF)
            )
            self.pi_oh_angle.set(
                pi.get("oh_pi_angle_cutoff", ParametersDefault.PI_OH_ANGLE_CUTOFF)
            )
            self.pi_sh_distance.set(
                pi.get("sh_pi_distance_cutoff", ParametersDefault.PI_SH_DISTANCE_CUTOFF)
            )
            self.pi_sh_angle.set(
                pi.get("sh_pi_angle_cutoff", ParametersDefault.PI_SH_ANGLE_CUTOFF)
            )

        # Apply general parameters
        if "general" in params:
            gen = params["general"]
            self.covalent_factor.set(
                gen.get(
                    "covalent_cutoff_factor", ParametersDefault.COVALENT_CUTOFF_FACTOR
                )
            )
            self.analysis_mode.set(
                gen.get("analysis_mode", ParametersDefault.ANALYSIS_MODE)
            )

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
