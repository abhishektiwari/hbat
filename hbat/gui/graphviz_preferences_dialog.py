"""
GraphViz preferences dialog for HBAT.

This module provides a dialog window for configuring GraphViz visualization
settings including engine selection, rendering options, and export preferences.
"""

import logging
import tkinter as tk
from tkinter import messagebox, ttk
from typing import Dict, Optional

from hbat.core.app_config import HBATConfig
from hbat.utilities.graphviz_utils import GraphVizDetector

# Set up logging
logger = logging.getLogger(__name__)

# Background color options
BACKGROUND_COLORS = [
    ("White", "white"),
    ("Light Gray", "lightgray"),
    ("Light Blue", "lightblue"),
    ("Transparent", "transparent"),
]

# Node shape options
NODE_SHAPES = [
    ("Ellipse", "ellipse"),
    ("Box", "box"),
    ("Circle", "circle"),
    ("Diamond", "diamond"),
    ("Polygon", "polygon"),
    ("Record", "record"),
]

# Graph direction options
RANK_DIRECTIONS = [
    ("Top to Bottom", "TB"),
    ("Bottom to Top", "BT"),
    ("Left to Right", "LR"),
    ("Right to Left", "RL"),
]


class GraphVizPreferencesDialog:
    """Dialog for configuring GraphViz preferences.

    Provides a user interface for setting GraphViz engine preferences,
    rendering options, and export settings.
    """

    def __init__(self, parent: tk.Widget, config: HBATConfig) -> None:
        """Initialize GraphViz preferences dialog.

        :param parent: Parent widget
        :type parent: tk.Widget
        :param config: HBAT configuration instance
        :type config: HBATConfig
        """
        self.parent = parent
        self.config = config
        self.dialog: Optional[tk.Toplevel] = None
        self.result = False

        # Variables for settings
        self.vars: Dict[str, tk.Variable] = {}

        # Check GraphViz availability
        self.graphviz_available = GraphVizDetector.is_graphviz_available()
        self.available_engines = GraphVizDetector.get_available_engines()

        # Create dialog
        self._create_dialog()

    def _create_dialog(self) -> None:
        """Create the preferences dialog window."""
        self.dialog = tk.Toplevel(self.parent)
        self.dialog.title("GraphViz Preferences")
        self.dialog.geometry("550x450")
        self.dialog.resizable(False, False)
        self.dialog.grab_set()  # Make dialog modal

        # Center the dialog
        self.dialog.update_idletasks()
        x = (self.dialog.winfo_screenwidth() - self.dialog.winfo_width()) // 2
        y = (self.dialog.winfo_screenheight() - self.dialog.winfo_height()) // 2
        self.dialog.geometry(f"+{x}+{y}")

        # Create main frame
        main_frame = ttk.Frame(self.dialog, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Create sections
        self._create_status_section(main_frame)
        self._create_engine_section(main_frame)
        self._create_rendering_section(main_frame)
        self._create_export_section(main_frame)
        self._create_buttons(main_frame)

    def _create_status_section(self, parent: ttk.Frame) -> None:
        """Create GraphViz status section.

        :param parent: Parent frame
        :type parent: ttk.Frame
        """
        status_frame = ttk.LabelFrame(parent, text="GraphViz Status", padding="10")
        status_frame.pack(fill=tk.X, pady=(0, 10))

        # GraphViz availability
        status_text = "Available" if self.graphviz_available else "Not Available"
        status_color = "green" if self.graphviz_available else "red"

        status_label = ttk.Label(
            status_frame,
            text=f"GraphViz Status: {status_text}",
            foreground=status_color,
        )
        status_label.pack(anchor=tk.W)

        if self.graphviz_available:
            # Show version if available
            version = GraphVizDetector.get_graphviz_version()
            if version:
                ttk.Label(status_frame, text=f"Version: {version}").pack(anchor=tk.W)

            # Show available engines
            ttk.Label(
                status_frame,
                text=f"Available Engines: {', '.join(self.available_engines)}",
            ).pack(anchor=tk.W)
        else:
            ttk.Label(
                status_frame,
                text="GraphViz is not installed. Please install GraphViz to use these features.",
                foreground="gray",
            ).pack(anchor=tk.W, pady=(5, 0))

        # Enable/Disable GraphViz
        self.vars["enabled"] = tk.BooleanVar(value=self.config.is_graphviz_enabled())
        ttk.Checkbutton(
            status_frame,
            text="Enable GraphViz visualization",
            variable=self.vars["enabled"],
            state="normal" if self.graphviz_available else "disabled",
        ).pack(anchor=tk.W, pady=(10, 0))

    def _create_engine_section(self, parent: ttk.Frame) -> None:
        """Create engine selection section.

        :param parent: Parent frame
        :type parent: ttk.Frame
        """
        engine_frame = ttk.LabelFrame(parent, text="Layout Engine", padding="10")
        engine_frame.pack(fill=tk.X, pady=(0, 10))

        # Current engine
        current_engine = self.config.get_graphviz_engine()
        self.vars["engine"] = tk.StringVar(value=current_engine)

        # Engine selection - ensure "dot" is first
        engines_ordered = self.available_engines.copy()
        if "dot" in engines_ordered:
            engines_ordered = ["dot"] + [
                engine for engine in engines_ordered if engine != "dot"
            ]

        for engine in engines_ordered:
            if engine in ENGINE_DESCRIPTIONS:
                text = f"{engine} - {ENGINE_DESCRIPTIONS[engine]}"
            else:
                text = engine

            ttk.Radiobutton(
                engine_frame,
                text=text,
                variable=self.vars["engine"],
                value=engine,
                state="normal" if self.graphviz_available else "disabled",
            ).pack(anchor=tk.W, pady=2)

    def _create_rendering_section(self, parent: ttk.Frame) -> None:
        """Create rendering options section.

        :param parent: Parent frame
        :type parent: ttk.Frame
        """
        render_frame = ttk.LabelFrame(parent, text="Rendering Options", padding="10")
        render_frame.pack(fill=tk.X, pady=(0, 10))

        # Background color
        bg_frame = ttk.Frame(render_frame)
        bg_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Label(bg_frame, text="Background Color:").pack(side=tk.LEFT)

        current_bg = self.config.get_graphviz_preference("background_color", "white")
        self.vars["background_color"] = tk.StringVar(value=current_bg)

        bg_menu = ttk.Combobox(
            bg_frame,
            textvariable=self.vars["background_color"],
            values=[name for name, _ in BACKGROUND_COLORS],
            state="readonly" if self.graphviz_available else "disabled",
            width=15,
        )
        bg_menu.pack(side=tk.LEFT, padx=(10, 0))

        # Node shape
        shape_frame = ttk.Frame(render_frame)
        shape_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Label(shape_frame, text="Node Shape:").pack(side=tk.LEFT)

        current_shape = self.config.get_graphviz_preference("node_shape", "ellipse")
        self.vars["node_shape"] = tk.StringVar(value=current_shape)

        shape_menu = ttk.Combobox(
            shape_frame,
            textvariable=self.vars["node_shape"],
            values=[name for name, _ in NODE_SHAPES],
            state="readonly" if self.graphviz_available else "disabled",
            width=15,
        )
        shape_menu.pack(side=tk.LEFT, padx=(10, 0))

        # Graph direction
        dir_frame = ttk.Frame(render_frame)
        dir_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Label(dir_frame, text="Graph Direction:").pack(side=tk.LEFT)

        current_dir = self.config.get_graphviz_preference("rankdir", "TB")
        self.vars["rankdir"] = tk.StringVar(value=current_dir)

        dir_menu = ttk.Combobox(
            dir_frame,
            textvariable=self.vars["rankdir"],
            values=[name for name, _ in RANK_DIRECTIONS],
            state="readonly" if self.graphviz_available else "disabled",
            width=15,
        )
        dir_menu.pack(side=tk.LEFT, padx=(10, 0))

    def _create_export_section(self, parent: ttk.Frame) -> None:
        """Create export settings section.

        :param parent: Parent frame
        :type parent: ttk.Frame
        """
        export_frame = ttk.LabelFrame(parent, text="Export Settings", padding="10")
        export_frame.pack(fill=tk.X, pady=(0, 10))

        # Export DPI
        dpi_frame = ttk.Frame(export_frame)
        dpi_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Label(dpi_frame, text="Export DPI:").pack(side=tk.LEFT)

        current_dpi = self.config.get_graphviz_export_dpi()
        self.vars["export_dpi"] = tk.IntVar(value=current_dpi)

        dpi_spinbox = ttk.Spinbox(
            dpi_frame,
            from_=72,
            to=600,
            increment=50,
            textvariable=self.vars["export_dpi"],
            state="normal" if self.graphviz_available else "disabled",
            width=10,
        )
        dpi_spinbox.pack(side=tk.LEFT, padx=(10, 0))

        ttk.Label(dpi_frame, text="(72-600)").pack(side=tk.LEFT, padx=(5, 0))

        # Default format
        format_frame = ttk.Frame(export_frame)
        format_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Label(format_frame, text="Default Format:").pack(side=tk.LEFT)

        current_format = self.config.get_graphviz_render_format()
        self.vars["render_format"] = tk.StringVar(value=current_format)

        format_menu = ttk.Combobox(
            format_frame,
            textvariable=self.vars["render_format"],
            values=["png", "svg", "pdf"],
            state="readonly" if self.graphviz_available else "disabled",
            width=10,
        )
        format_menu.pack(side=tk.LEFT, padx=(10, 0))

    def _create_buttons(self, parent: ttk.Frame) -> None:
        """Create dialog buttons.

        :param parent: Parent frame
        :type parent: ttk.Frame
        """
        button_frame = ttk.Frame(parent)
        button_frame.pack(fill=tk.X, pady=(10, 0))

        # Cancel button
        ttk.Button(button_frame, text="Cancel", command=self._cancel).pack(
            side=tk.RIGHT, padx=(5, 0)
        )

        # OK button
        ttk.Button(button_frame, text="OK", command=self._save_preferences).pack(
            side=tk.RIGHT
        )

        # Reset button
        ttk.Button(
            button_frame, text="Reset to Defaults", command=self._reset_defaults
        ).pack(side=tk.LEFT)

    def _save_preferences(self) -> None:
        """Save preferences and close dialog."""
        try:
            # Save enabled state
            self.config.enable_graphviz(self.vars["enabled"].get())

            # Save engine preference
            self.config.set_graphviz_engine(self.vars["engine"].get())

            # Save rendering options
            # Convert display names back to values
            bg_value = "white"
            for name, value in BACKGROUND_COLORS:
                if name == self.vars["background_color"].get():
                    bg_value = value
                    break
            self.config.set_graphviz_preference("background_color", bg_value)

            shape_value = "ellipse"
            for name, value in NODE_SHAPES:
                if name == self.vars["node_shape"].get():
                    shape_value = value
                    break
            self.config.set_graphviz_preference("node_shape", shape_value)

            dir_value = "TB"
            for name, value in RANK_DIRECTIONS:
                if name == self.vars["rankdir"].get():
                    dir_value = value
                    break
            self.config.set_graphviz_preference("rankdir", dir_value)

            # Save export settings
            dpi = self.vars["export_dpi"].get()
            if 72 <= dpi <= 600:
                self.config.set_graphviz_export_dpi(dpi)
            else:
                raise ValueError("DPI must be between 72 and 600")

            self.config.set_graphviz_render_format(self.vars["render_format"].get())

            # Save config
            self.config.save_config(self.config.load_config())

            self.result = True
            self.dialog.destroy()

            logger.info("GraphViz preferences saved successfully")

        except Exception as e:
            logger.error(f"Failed to save preferences: {e}")
            messagebox.showerror(
                "Error", f"Failed to save preferences: {str(e)}", parent=self.dialog
            )

    def _cancel(self) -> None:
        """Cancel and close dialog."""
        self.result = False
        self.dialog.destroy()

    def _reset_defaults(self) -> None:
        """Reset all settings to defaults."""
        response = messagebox.askyesno(
            "Reset to Defaults",
            "Are you sure you want to reset all GraphViz settings to defaults?",
            parent=self.dialog,
        )

        if response:
            # Reset all variables to defaults
            self.vars["enabled"].set(True)
            self.vars["engine"].set("dot")
            self.vars["background_color"].set("White")
            self.vars["node_shape"].set("Ellipse")
            self.vars["rankdir"].set("Top to Bottom")
            self.vars["export_dpi"].set(300)
            self.vars["render_format"].set("png")

    def show(self) -> bool:
        """Show dialog and wait for result.

        :returns: True if preferences were saved, False if cancelled
        :rtype: bool
        """
        self.dialog.wait_window()
        return self.result


def show_graphviz_preferences(parent: tk.Widget, config: HBATConfig) -> bool:
    """Show GraphViz preferences dialog.

    :param parent: Parent widget
    :type parent: tk.Widget
    :param config: HBAT configuration instance
    :type config: HBATConfig
    :returns: True if preferences were saved
    :rtype: bool
    """
    dialog = GraphVizPreferencesDialog(parent, config)
    return dialog.show()
