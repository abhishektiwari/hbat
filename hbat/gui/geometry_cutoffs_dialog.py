"""
Geometry cutoffs configuration dialog for HBAT GUI.

This module provides a dialog for configuring molecular interaction
analysis parameters (distances, angles, etc.) without PDB fixing options.
"""

import tkinter as tk
from tkinter import messagebox, ttk
from typing import Any, Dict, Optional

from ..constants.parameters import (
    AnalysisParameters,
    ParameterRanges,
    ParametersDefault,
)
from ..config.parameter_controller import ParameterController


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
        """Initialize parameter controller and variable dictionaries."""
        self.controller = ParameterController()
        self._vars = {}  # field_name -> tk.Var (created lazily per panel)
        self._param_values = {}  # field_name -> value (cross-panel backing store)

    def _store_current_values(self):
        """Store current parameter values before switching categories."""
        for field, var in self._vars.items():
            if var is not None:
                try:
                    self._param_values[field] = var.get()
                except tk.TclError:
                    pass  # Widget destroyed, ignore

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

        # Interaction inclusion mode
        ttk.Label(group, text="Interaction Inclusion:").grid(
            row=0, column=0, sticky=tk.W, pady=2
        )
        stored_mode = self._param_values.get(
            "analysis_mode", ParametersDefault.ANALYSIS_MODE
        )
        self._vars["analysis_mode"] = tk.StringVar(value=stored_mode)
        mode_frame = ttk.Frame(group)
        mode_frame.grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        ttk.Radiobutton(
            mode_frame,
            text="All interactions",
            variable=self._vars["analysis_mode"],
            value="all",
        ).pack(anchor=tk.W)
        ttk.Radiobutton(
            mode_frame,
            text="Inter-residue only",
            variable=self._vars["analysis_mode"],
            value="inter",
        ).pack(anchor=tk.W)

        # Covalent bond cutoff factor
        ttk.Label(group, text="Covalent Bond Factor:").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        stored_factor = self._param_values.get(
            "covalent_cutoff_factor", ParametersDefault.COVALENT_CUTOFF_FACTOR
        )
        self._vars["covalent_cutoff_factor"] = tk.DoubleVar(value=stored_factor)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_COVALENT_FACTOR,
            to=ParameterRanges.MAX_COVALENT_FACTOR,
            variable=self._vars["covalent_cutoff_factor"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        # Value display
        factor_label = ttk.Label(group, text="")
        factor_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_factor_label(*args):
            try:
                factor_label.config(
                    text=f"{self._vars['covalent_cutoff_factor'].get():.2f}"
                )
            except tk.TclError:
                pass  # Widget destroyed, ignore

        self._vars["covalent_cutoff_factor"].trace("w", update_factor_label)
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
            "hb_distance_cutoff", ParametersDefault.HB_DISTANCE_CUTOFF
        )
        self._vars["hb_distance_cutoff"] = tk.DoubleVar(value=stored_dist)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["hb_distance_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        hb_dist_label = ttk.Label(group, text="")
        hb_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_hb_dist(*args):
            try:
                hb_dist_label.config(
                    text=f"{self._vars['hb_distance_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["hb_distance_cutoff"].trace("w", update_hb_dist)
        update_hb_dist()

        # D-H...A angle
        ttk.Label(group, text="D-H...A Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        stored_angle = self._param_values.get(
            "hb_angle_cutoff", ParametersDefault.HB_ANGLE_CUTOFF
        )
        self._vars["hb_angle_cutoff"] = tk.DoubleVar(value=stored_angle)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["hb_angle_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        hb_angle_label = ttk.Label(group, text="")
        hb_angle_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_hb_angle(*args):
            try:
                hb_angle_label.config(text=f"{self._vars['hb_angle_cutoff'].get():.0f}")
            except tk.TclError:
                pass

        self._vars["hb_angle_cutoff"].trace("w", update_hb_angle)
        update_hb_angle()

        # D...A distance
        ttk.Label(group, text="D...A Distance (Å):").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )
        stored_da = self._param_values.get(
            "hb_donor_acceptor_cutoff", ParametersDefault.HB_DA_DISTANCE
        )
        self._vars["hb_donor_acceptor_cutoff"] = tk.DoubleVar(value=stored_da)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["hb_donor_acceptor_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=2, column=1, sticky=tk.W, padx=10, pady=2)

        da_dist_label = ttk.Label(group, text="")
        da_dist_label.grid(row=2, column=2, sticky=tk.W, padx=5, pady=2)

        def update_da_dist(*args):
            try:
                da_dist_label.config(
                    text=f"{self._vars['hb_donor_acceptor_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["hb_donor_acceptor_cutoff"].trace("w", update_da_dist)
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
            "whb_distance_cutoff", ParametersDefault.WHB_DISTANCE_CUTOFF
        )
        self._vars["whb_distance_cutoff"] = tk.DoubleVar(value=stored_whb_dist)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["whb_distance_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        whb_dist_label = ttk.Label(group, text="")
        whb_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_whb_dist(*args):
            try:
                whb_dist_label.config(
                    text=f"{self._vars['whb_distance_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["whb_distance_cutoff"].trace("w", update_whb_dist)
        update_whb_dist()

        # WHB D-H...A angle
        ttk.Label(group, text="D-H...A Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        stored_whb_angle = self._param_values.get(
            "whb_angle_cutoff", ParametersDefault.WHB_ANGLE_CUTOFF
        )
        self._vars["whb_angle_cutoff"] = tk.DoubleVar(value=stored_whb_angle)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["whb_angle_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        whb_angle_label = ttk.Label(group, text="")
        whb_angle_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_whb_angle(*args):
            try:
                whb_angle_label.config(
                    text=f"{self._vars['whb_angle_cutoff'].get():.0f}"
                )
            except tk.TclError:
                pass

        self._vars["whb_angle_cutoff"].trace("w", update_whb_angle)
        update_whb_angle()

        # WHB D...A distance
        ttk.Label(group, text="D...A Distance (Å):").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )
        stored_whb_da = self._param_values.get(
            "whb_donor_acceptor_cutoff", ParametersDefault.WHB_DA_DISTANCE
        )
        self._vars["whb_donor_acceptor_cutoff"] = tk.DoubleVar(value=stored_whb_da)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["whb_donor_acceptor_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=2, column=1, sticky=tk.W, padx=10, pady=2)

        whb_da_dist_label = ttk.Label(group, text="")
        whb_da_dist_label.grid(row=2, column=2, sticky=tk.W, padx=5, pady=2)

        def update_whb_da_dist(*args):
            try:
                whb_da_dist_label.config(
                    text=f"{self._vars['whb_donor_acceptor_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["whb_donor_acceptor_cutoff"].trace("w", update_whb_da_dist)
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
            "xb_distance_cutoff", ParametersDefault.XB_DISTANCE_CUTOFF
        )
        self._vars["xb_distance_cutoff"] = tk.DoubleVar(value=stored_xb_dist)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["xb_distance_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        xb_dist_label = ttk.Label(group, text="")
        xb_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_xb_dist(*args):
            try:
                xb_dist_label.config(
                    text=f"{self._vars['xb_distance_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["xb_distance_cutoff"].trace("w", update_xb_dist)
        update_xb_dist()

        # C-X...A angle
        ttk.Label(group, text="C-X...A Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        stored_xb_angle = self._param_values.get(
            "xb_angle_cutoff", ParametersDefault.XB_ANGLE_CUTOFF
        )
        self._vars["xb_angle_cutoff"] = tk.DoubleVar(value=stored_xb_angle)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["xb_angle_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        xb_angle_label = ttk.Label(group, text="")
        xb_angle_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_xb_angle(*args):
            try:
                xb_angle_label.config(text=f"{self._vars['xb_angle_cutoff'].get():.0f}")
            except tk.TclError:
                pass

        self._vars["xb_angle_cutoff"].trace("w", update_xb_angle)
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
            "pi_distance_cutoff", ParametersDefault.PI_DISTANCE_CUTOFF
        )
        self._vars["pi_distance_cutoff"] = tk.DoubleVar(value=stored_pi_dist)
        ttk.Scale(
            general_group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["pi_distance_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=0, column=1, sticky=tk.W, padx=10, pady=2)

        pi_dist_label = ttk.Label(general_group, text="")
        pi_dist_label.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_dist(*args):
            try:
                pi_dist_label.config(
                    text=f"{self._vars['pi_distance_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["pi_distance_cutoff"].trace("w", update_pi_dist)
        update_pi_dist()

        # D-H...π angle
        ttk.Label(general_group, text="D-H...π Angle (°):").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        stored_pi_angle = self._param_values.get(
            "pi_angle_cutoff", ParametersDefault.PI_ANGLE_CUTOFF
        )
        self._vars["pi_angle_cutoff"] = tk.DoubleVar(value=stored_pi_angle)
        ttk.Scale(
            general_group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["pi_angle_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=1, column=1, sticky=tk.W, padx=10, pady=2)

        pi_angle_label = ttk.Label(general_group, text="")
        pi_angle_label.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_angle(*args):
            try:
                pi_angle_label.config(text=f"{self._vars['pi_angle_cutoff'].get():.0f}")
            except tk.TclError:
                pass

        self._vars["pi_angle_cutoff"].trace("w", update_pi_angle)
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
            dist_field_name,
            angle_field_name,
            dist_default,
            angle_default,
        ):
            # Get stored values or use defaults
            stored_dist = self._param_values.get(dist_field_name, dist_default)
            stored_angle = self._param_values.get(angle_field_name, angle_default)

            # Distance parameter
            ttk.Label(parent_frame, text=f"{label_text} Distance (Å):").grid(
                row=row, column=0, sticky=tk.W, pady=2
            )
            self._vars[dist_field_name] = tk.DoubleVar(value=stored_dist)
            ttk.Scale(
                parent_frame,
                from_=ParameterRanges.MIN_DISTANCE,
                to=ParameterRanges.MAX_DISTANCE,
                variable=self._vars[dist_field_name],
                orient=tk.HORIZONTAL,
                length=150,
            ).grid(row=row, column=1, sticky=tk.W, padx=5, pady=2)

            dist_label = ttk.Label(parent_frame, text="")
            dist_label.grid(row=row, column=2, sticky=tk.W, padx=5, pady=2)

            # Angle parameter
            ttk.Label(parent_frame, text=f"{label_text} Angle (°):").grid(
                row=row, column=3, sticky=tk.W, pady=2, padx=(20, 0)
            )
            self._vars[angle_field_name] = tk.DoubleVar(value=stored_angle)
            ttk.Scale(
                parent_frame,
                from_=ParameterRanges.MIN_ANGLE,
                to=ParameterRanges.MAX_ANGLE,
                variable=self._vars[angle_field_name],
                orient=tk.HORIZONTAL,
                length=150,
            ).grid(row=row, column=4, sticky=tk.W, padx=5, pady=2)

            angle_label = ttk.Label(parent_frame, text="")
            angle_label.grid(row=row, column=5, sticky=tk.W, padx=5, pady=2)

            # Update functions
            def update_dist(*args):
                try:
                    dist_label.config(text=f"{self._vars[dist_field_name].get():.1f}")
                except tk.TclError:
                    pass

            def update_angle(*args):
                try:
                    angle_label.config(text=f"{self._vars[angle_field_name].get():.0f}")
                except tk.TclError:
                    pass

            self._vars[dist_field_name].trace("w", update_dist)
            self._vars[angle_field_name].trace("w", update_angle)
            update_dist()
            update_angle()

        # Create all subtype parameters
        create_parameter_pair(
            subtypes_group,
            0,
            "C-Cl...π",
            "pi_ccl_distance_cutoff",
            "pi_ccl_angle_cutoff",
            ParametersDefault.PI_CCL_DISTANCE_CUTOFF,
            ParametersDefault.PI_CCL_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            1,
            "C-Br...π",
            "pi_cbr_distance_cutoff",
            "pi_cbr_angle_cutoff",
            ParametersDefault.PI_CBR_DISTANCE_CUTOFF,
            ParametersDefault.PI_CBR_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            2,
            "C-I...π",
            "pi_ci_distance_cutoff",
            "pi_ci_angle_cutoff",
            ParametersDefault.PI_CI_DISTANCE_CUTOFF,
            ParametersDefault.PI_CI_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            3,
            "C-H...π",
            "pi_ch_distance_cutoff",
            "pi_ch_angle_cutoff",
            ParametersDefault.PI_CH_DISTANCE_CUTOFF,
            ParametersDefault.PI_CH_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            4,
            "N-H...π",
            "pi_nh_distance_cutoff",
            "pi_nh_angle_cutoff",
            ParametersDefault.PI_NH_DISTANCE_CUTOFF,
            ParametersDefault.PI_NH_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            5,
            "O-H...π",
            "pi_oh_distance_cutoff",
            "pi_oh_angle_cutoff",
            ParametersDefault.PI_OH_DISTANCE_CUTOFF,
            ParametersDefault.PI_OH_ANGLE_CUTOFF,
        )
        create_parameter_pair(
            subtypes_group,
            6,
            "S-H...π",
            "pi_sh_distance_cutoff",
            "pi_sh_angle_cutoff",
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

        stored_pi_pi_dist = self._param_values.get(
            "pi_pi_distance_cutoff", ParametersDefault.PI_PI_DISTANCE_CUTOFF
        )
        self._vars["pi_pi_distance_cutoff"] = tk.DoubleVar(value=stored_pi_pi_dist)
        pi_pi_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["pi_pi_distance_cutoff"],
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
                pi_pi_dist_label.config(
                    text=f"{self._vars['pi_pi_distance_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["pi_pi_distance_cutoff"].trace("w", update_pi_pi_dist)
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

        stored_pi_pi_parallel = self._param_values.get(
            "pi_pi_parallel_angle_cutoff", ParametersDefault.PI_PI_PARALLEL_ANGLE_CUTOFF
        )
        self._vars["pi_pi_parallel_angle_cutoff"] = tk.DoubleVar(
            value=stored_pi_pi_parallel
        )
        parallel_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["pi_pi_parallel_angle_cutoff"],
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
                    text=f"{self._vars['pi_pi_parallel_angle_cutoff'].get():.0f}"
                )
            except tk.TclError:
                pass

        self._vars["pi_pi_parallel_angle_cutoff"].trace("w", update_pi_pi_parallel)
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

        stored_pi_pi_tmin = self._param_values.get(
            "pi_pi_tshaped_angle_min", ParametersDefault.PI_PI_TSHAPED_ANGLE_MIN
        )
        self._vars["pi_pi_tshaped_angle_min"] = tk.DoubleVar(value=stored_pi_pi_tmin)
        tmin_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["pi_pi_tshaped_angle_min"],
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
                    text=f"{self._vars['pi_pi_tshaped_angle_min'].get():.0f}"
                )
            except tk.TclError:
                pass

        self._vars["pi_pi_tshaped_angle_min"].trace("w", update_pi_pi_tmin)
        update_pi_pi_tmin()

        # T-shaped angle maximum
        ttk.Label(group, text="T-shaped Angle Max (°):").grid(
            row=3, column=0, sticky=tk.W, pady=2
        )
        stored_pi_pi_tmax = self._param_values.get(
            "pi_pi_tshaped_angle_max", ParametersDefault.PI_PI_TSHAPED_ANGLE_MAX
        )
        self._vars["pi_pi_tshaped_angle_max"] = tk.DoubleVar(value=stored_pi_pi_tmax)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["pi_pi_tshaped_angle_max"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=3, column=1, sticky=tk.W, padx=10, pady=2)

        pi_pi_tmax_label = ttk.Label(group, text="")
        pi_pi_tmax_label.grid(row=3, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_pi_tmax(*args):
            try:
                pi_pi_tmax_label.config(
                    text=f"{self._vars['pi_pi_tshaped_angle_max'].get():.0f}"
                )
            except tk.TclError:
                pass

        self._vars["pi_pi_tshaped_angle_max"].trace("w", update_pi_pi_tmax)
        update_pi_pi_tmax()

        # Offset cutoff
        ttk.Label(group, text="Offset Cutoff (Å):").grid(
            row=4, column=0, sticky=tk.W, pady=2
        )
        stored_pi_pi_offset = self._param_values.get(
            "pi_pi_offset_cutoff", ParametersDefault.PI_PI_OFFSET_CUTOFF
        )
        self._vars["pi_pi_offset_cutoff"] = tk.DoubleVar(value=stored_pi_pi_offset)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["pi_pi_offset_cutoff"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=4, column=1, sticky=tk.W, padx=10, pady=2)

        pi_pi_offset_label = ttk.Label(group, text="")
        pi_pi_offset_label.grid(row=4, column=2, sticky=tk.W, padx=5, pady=2)

        def update_pi_pi_offset(*args):
            try:
                pi_pi_offset_label.config(
                    text=f"{self._vars['pi_pi_offset_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["pi_pi_offset_cutoff"].trace("w", update_pi_pi_offset)
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

        stored_carb_dist = self._param_values.get(
            "carbonyl_distance_cutoff", ParametersDefault.CARBONYL_DISTANCE_CUTOFF
        )
        self._vars["carbonyl_distance_cutoff"] = tk.DoubleVar(value=stored_carb_dist)
        carb_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["carbonyl_distance_cutoff"],
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
                carbonyl_dist_label.config(
                    text=f"{self._vars['carbonyl_distance_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["carbonyl_distance_cutoff"].trace("w", update_carbonyl_dist)
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

        stored_carb_min = self._param_values.get(
            "carbonyl_angle_min", ParametersDefault.CARBONYL_ANGLE_MIN
        )
        self._vars["carbonyl_angle_min"] = tk.DoubleVar(value=stored_carb_min)
        angle_min_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["carbonyl_angle_min"],
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
                carbonyl_min_label.config(
                    text=f"{self._vars['carbonyl_angle_min'].get():.0f}"
                )
            except tk.TclError:
                pass

        self._vars["carbonyl_angle_min"].trace("w", update_carbonyl_min)
        update_carbonyl_min()

        # Angle maximum
        ttk.Label(group, text="Bürgi-Dunitz Angle Max (°):").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )
        stored_carb_max = self._param_values.get(
            "carbonyl_angle_max", ParametersDefault.CARBONYL_ANGLE_MAX
        )
        self._vars["carbonyl_angle_max"] = tk.DoubleVar(value=stored_carb_max)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["carbonyl_angle_max"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=2, column=1, sticky=tk.W, padx=10, pady=2)

        carbonyl_max_label = ttk.Label(group, text="")
        carbonyl_max_label.grid(row=2, column=2, sticky=tk.W, padx=5, pady=2)

        def update_carbonyl_max(*args):
            try:
                carbonyl_max_label.config(
                    text=f"{self._vars['carbonyl_angle_max'].get():.0f}"
                )
            except tk.TclError:
                pass

        self._vars["carbonyl_angle_max"].trace("w", update_carbonyl_max)
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

        stored_n_pi_dist = self._param_values.get(
            "n_pi_distance_cutoff", ParametersDefault.N_PI_DISTANCE_CUTOFF
        )
        self._vars["n_pi_distance_cutoff"] = tk.DoubleVar(value=stored_n_pi_dist)
        n_pi_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["n_pi_distance_cutoff"],
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
                n_pi_dist_label.config(
                    text=f"{self._vars['n_pi_distance_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["n_pi_distance_cutoff"].trace("w", update_n_pi_dist)
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

        stored_n_pi_sulfur = self._param_values.get(
            "n_pi_sulfur_distance_cutoff", ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF
        )
        self._vars["n_pi_sulfur_distance_cutoff"] = tk.DoubleVar(
            value=stored_n_pi_sulfur
        )
        sulfur_scale = ttk.Scale(
            group,
            from_=ParameterRanges.MIN_DISTANCE,
            to=ParameterRanges.MAX_DISTANCE,
            variable=self._vars["n_pi_sulfur_distance_cutoff"],
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
                n_pi_sulfur_label.config(
                    text=f"{self._vars['n_pi_sulfur_distance_cutoff'].get():.1f}"
                )
            except tk.TclError:
                pass

        self._vars["n_pi_sulfur_distance_cutoff"].trace("w", update_n_pi_sulfur)
        update_n_pi_sulfur()

        # Angle minimum
        ttk.Label(group, text="Angle Min (°):").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )
        stored_n_pi_min = self._param_values.get(
            "n_pi_angle_min", ParametersDefault.N_PI_ANGLE_MIN
        )
        self._vars["n_pi_angle_min"] = tk.DoubleVar(value=stored_n_pi_min)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["n_pi_angle_min"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=2, column=1, sticky=tk.W, padx=10, pady=2)

        n_pi_min_label = ttk.Label(group, text="")
        n_pi_min_label.grid(row=2, column=2, sticky=tk.W, padx=5, pady=2)

        def update_n_pi_min(*args):
            try:
                n_pi_min_label.config(text=f"{self._vars['n_pi_angle_min'].get():.0f}")
            except tk.TclError:
                pass

        self._vars["n_pi_angle_min"].trace("w", update_n_pi_min)
        update_n_pi_min()

        # Angle maximum
        ttk.Label(group, text="Angle Max (°):").grid(
            row=3, column=0, sticky=tk.W, pady=2
        )
        stored_n_pi_max = self._param_values.get(
            "n_pi_angle_max", ParametersDefault.N_PI_ANGLE_MAX
        )
        self._vars["n_pi_angle_max"] = tk.DoubleVar(value=stored_n_pi_max)
        ttk.Scale(
            group,
            from_=ParameterRanges.MIN_ANGLE,
            to=ParameterRanges.MAX_ANGLE,
            variable=self._vars["n_pi_angle_max"],
            orient=tk.HORIZONTAL,
            length=200,
        ).grid(row=3, column=1, sticky=tk.W, padx=10, pady=2)

        n_pi_max_label = ttk.Label(group, text="")
        n_pi_max_label.grid(row=3, column=2, sticky=tk.W, padx=5, pady=2)

        def update_n_pi_max(*args):
            try:
                n_pi_max_label.config(text=f"{self._vars['n_pi_angle_max'].get():.0f}")
            except tk.TclError:
                pass

        self._vars["n_pi_angle_max"].trace("w", update_n_pi_max)
        update_n_pi_max()

    def get_parameters(self) -> AnalysisParameters:
        """Get current parameter values.

        :returns: Current analysis parameters
        :rtype: AnalysisParameters
        """
        # Store current values from tkinter variables
        self._store_current_values()

        # Create new AnalysisParameters with values from backing store
        params = AnalysisParameters()
        for field, value in self._param_values.items():
            if hasattr(params, field):
                setattr(params, field, value)

        return params

    def set_parameters(self, params: AnalysisParameters) -> None:
        """Set parameter values from AnalysisParameters object.

        :param params: Analysis parameters to set
        :type params: AnalysisParameters
        """
        for field, value in vars(params).items():
            self._param_values[field] = value
            var = self._vars.get(field)
            if var is not None:
                try:
                    var.set(value)
                except tk.TclError:
                    pass  # Variable destroyed, ignore

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
            try:
                self._apply_preset_data(result)
                messagebox.showinfo("Success", "Preset loaded successfully")
            except ValueError as e:
                messagebox.showerror("Invalid Preset", str(e))

    def _apply_preset_data(self, data: Dict[str, Any]) -> None:
        """Apply preset data to parameters."""
        if "parameters" not in data:
            raise ValueError("Invalid preset format: missing 'parameters' section")

        params = data["parameters"]
        previous_params = self.get_parameters()

        # Helper to apply value to parameter
        def apply_value(field_name, value):
            self._param_values[field_name] = value
            var = self._vars.get(field_name)
            if var is not None:
                try:
                    var.set(value)
                except tk.TclError:
                    pass

        # Apply hydrogen bond parameters
        if "hydrogen_bonds" in params:
            hb = params["hydrogen_bonds"]
            apply_value(
                "hb_distance_cutoff",
                hb.get("h_a_distance_cutoff", ParametersDefault.HB_DISTANCE_CUTOFF),
            )
            apply_value(
                "hb_angle_cutoff",
                hb.get("dha_angle_cutoff", ParametersDefault.HB_ANGLE_CUTOFF),
            )
            apply_value(
                "hb_donor_acceptor_cutoff",
                hb.get("d_a_distance_cutoff", ParametersDefault.HB_DA_DISTANCE),
            )

        # Apply halogen bond parameters
        if "halogen_bonds" in params:
            xb = params["halogen_bonds"]
            apply_value(
                "xb_distance_cutoff",
                xb.get("x_a_distance_cutoff", ParametersDefault.XB_DISTANCE_CUTOFF),
            )
            apply_value(
                "xb_angle_cutoff",
                xb.get("dxa_angle_cutoff", ParametersDefault.XB_ANGLE_CUTOFF),
            )

        # Apply π interaction parameters
        if "pi_interactions" in params:
            pi = params["pi_interactions"]
            apply_value(
                "pi_distance_cutoff",
                pi.get("h_pi_distance_cutoff", ParametersDefault.PI_DISTANCE_CUTOFF),
            )
            apply_value(
                "pi_angle_cutoff",
                pi.get("dh_pi_angle_cutoff", ParametersDefault.PI_ANGLE_CUTOFF),
            )

            # Apply π interaction subtype parameters
            apply_value(
                "pi_ccl_distance_cutoff",
                pi.get(
                    "ccl_pi_distance_cutoff", ParametersDefault.PI_CCL_DISTANCE_CUTOFF
                ),
            )
            apply_value(
                "pi_ccl_angle_cutoff",
                pi.get("ccl_pi_angle_cutoff", ParametersDefault.PI_CCL_ANGLE_CUTOFF),
            )
            apply_value(
                "pi_cbr_distance_cutoff",
                pi.get(
                    "cbr_pi_distance_cutoff", ParametersDefault.PI_CBR_DISTANCE_CUTOFF
                ),
            )
            apply_value(
                "pi_cbr_angle_cutoff",
                pi.get("cbr_pi_angle_cutoff", ParametersDefault.PI_CBR_ANGLE_CUTOFF),
            )
            apply_value(
                "pi_ci_distance_cutoff",
                pi.get(
                    "ci_pi_distance_cutoff", ParametersDefault.PI_CI_DISTANCE_CUTOFF
                ),
            )
            apply_value(
                "pi_ci_angle_cutoff",
                pi.get("ci_pi_angle_cutoff", ParametersDefault.PI_CI_ANGLE_CUTOFF),
            )
            apply_value(
                "pi_ch_distance_cutoff",
                pi.get(
                    "ch_pi_distance_cutoff", ParametersDefault.PI_CH_DISTANCE_CUTOFF
                ),
            )
            apply_value(
                "pi_ch_angle_cutoff",
                pi.get("ch_pi_angle_cutoff", ParametersDefault.PI_CH_ANGLE_CUTOFF),
            )
            apply_value(
                "pi_nh_distance_cutoff",
                pi.get(
                    "nh_pi_distance_cutoff", ParametersDefault.PI_NH_DISTANCE_CUTOFF
                ),
            )
            apply_value(
                "pi_nh_angle_cutoff",
                pi.get("nh_pi_angle_cutoff", ParametersDefault.PI_NH_ANGLE_CUTOFF),
            )
            apply_value(
                "pi_oh_distance_cutoff",
                pi.get(
                    "oh_pi_distance_cutoff", ParametersDefault.PI_OH_DISTANCE_CUTOFF
                ),
            )
            apply_value(
                "pi_oh_angle_cutoff",
                pi.get("oh_pi_angle_cutoff", ParametersDefault.PI_OH_ANGLE_CUTOFF),
            )
            apply_value(
                "pi_sh_distance_cutoff",
                pi.get(
                    "sh_pi_distance_cutoff", ParametersDefault.PI_SH_DISTANCE_CUTOFF
                ),
            )
            apply_value(
                "pi_sh_angle_cutoff",
                pi.get("sh_pi_angle_cutoff", ParametersDefault.PI_SH_ANGLE_CUTOFF),
            )

        # Apply general parameters
        if "general" in params:
            gen = params["general"]
            apply_value(
                "covalent_cutoff_factor",
                gen.get(
                    "covalent_cutoff_factor", ParametersDefault.COVALENT_CUTOFF_FACTOR
                ),
            )
            apply_value(
                "analysis_mode",
                gen.get("analysis_mode", ParametersDefault.ANALYSIS_MODE),
            )

        try:
            self.get_parameters().validate_or_raise("preset parameters")
        except ValueError:
            self.set_parameters(previous_params)
            raise

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
