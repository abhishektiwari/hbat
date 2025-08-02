"""
Export manager for HBAT visualizations.

This module provides centralized export functionality for different
visualization renderers, handling file dialogs, format selection,
and export operations.
"""

import logging
import os
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, simpledialog, ttk
from typing import Dict, List, Optional, Tuple

from hbat.core.app_config import HBATConfig
from hbat.gui.visualization_renderer import VisualizationRenderer

# Set up logging
logger = logging.getLogger(__name__)

# Export format configurations
EXPORT_FORMATS = {
    "png": {
        "name": "PNG Image",
        "extension": ".png",
        "description": "Portable Network Graphics",
        "filetypes": [("PNG files", "*.png")],
    },
    "svg": {
        "name": "SVG Vector",
        "extension": ".svg",
        "description": "Scalable Vector Graphics",
        "filetypes": [("SVG files", "*.svg")],
    },
    "pdf": {
        "name": "PDF Document",
        "extension": ".pdf",
        "description": "Portable Document Format",
        "filetypes": [("PDF files", "*.pdf")],
    },
    "eps": {
        "name": "EPS Vector",
        "extension": ".eps",
        "description": "Encapsulated PostScript",
        "filetypes": [("EPS files", "*.eps")],
    },
    "dot": {
        "name": "DOT Source",
        "extension": ".dot",
        "description": "GraphViz DOT language source file",
        "filetypes": [("DOT files", "*.dot")],
    },
}


class ExportManager:
    """Manages visualization export functionality.

    Provides a unified interface for exporting visualizations from different
    renderers with format selection, quality settings, and file management.
    """

    def __init__(self, renderer: VisualizationRenderer, config: HBATConfig) -> None:
        """Initialize export manager.

        :param renderer: Visualization renderer instance
        :type renderer: VisualizationRenderer
        :param config: HBAT configuration instance
        :type config: HBATConfig
        """
        self.renderer = renderer
        self.config = config
        self.supported_formats = self._get_supported_formats()
        self.last_export_directory = self.config.get_preference("last_export_directory")

    def export_visualization(
        self,
        filename: Optional[str] = None,
        format: Optional[str] = None,
        resolution: Optional[int] = None,
    ) -> bool:
        """Export visualization to file.

        If filename or format are not provided, shows dialog to get them.

        :param filename: Output filename (optional)
        :type filename: Optional[str]
        :param format: Export format (optional)
        :type format: Optional[str]
        :param resolution: Export resolution/DPI (optional)
        :type resolution: Optional[int]
        :returns: True if export successful
        :rtype: bool
        """
        try:
            # Get filename and format if not provided
            if not filename or not format:
                export_info = self.show_export_dialog()
                if not export_info:
                    return False
                filename, format = export_info

            # Validate format
            if format not in self.supported_formats:
                messagebox.showerror(
                    "Export Error",
                    f"Format '{format}' is not supported by {self.renderer.get_renderer_name()}",
                )
                return False

            # Set resolution if provided
            if resolution:
                if hasattr(self.config, "set_graphviz_export_dpi"):
                    self.config.set_graphviz_export_dpi(resolution)

            # Perform export
            success = self.renderer.export(format, filename)

            if success:
                # Update last export directory
                self._update_last_export_directory(filename)

                # Show success message
                messagebox.showinfo(
                    "Export Successful", f"Visualization exported to:\n{filename}"
                )
                logger.info(f"Successfully exported visualization to {filename}")
            else:
                # Show error message
                messagebox.showerror(
                    "Export Failed", f"Failed to export visualization to {filename}"
                )
                logger.error(f"Failed to export visualization to {filename}")

            return success

        except Exception as e:
            error_msg = f"Export failed: {str(e)}"
            messagebox.showerror("Export Error", error_msg)
            logger.error(error_msg)
            return False

    def show_export_dialog(self) -> Optional[Tuple[str, str]]:
        """Show export dialog to get filename and format.

        :returns: Tuple of (filename, format) or None if cancelled
        :rtype: Optional[Tuple[str, str]]
        """
        try:
            # Create export dialog
            dialog = ExportDialog(
                self.supported_formats, self.config, self.last_export_directory
            )
            result = dialog.show()

            if result:
                filename, format = result
                return filename, format

            return None

        except Exception as e:
            logger.error(f"Error showing export dialog: {e}")
            return None

    def get_supported_formats(self) -> List[str]:
        """Get list of supported export formats.

        :returns: List of supported format names
        :rtype: List[str]
        """
        return self.supported_formats

    def get_format_info(self, format: str) -> Optional[Dict[str, any]]:
        """Get information about an export format.

        :param format: Format name
        :type format: str
        :returns: Format information dictionary
        :rtype: Optional[Dict[str, any]]
        """
        return EXPORT_FORMATS.get(format)

    def _get_supported_formats(self) -> List[str]:
        """Get formats supported by the current renderer.

        :returns: List of supported format names
        :rtype: List[str]
        """
        renderer_formats = self.renderer.get_supported_formats()
        return [fmt for fmt in renderer_formats if fmt in EXPORT_FORMATS]

    def _update_last_export_directory(self, filename: str) -> None:
        """Update the last export directory preference.

        :param filename: Exported filename
        :type filename: str
        """
        try:
            directory = str(Path(filename).parent)
            self.config.set_preference("last_export_directory", directory)
            self.last_export_directory = directory
        except Exception as e:
            logger.warning(f"Failed to update last export directory: {e}")


class ExportDialog:
    """Dialog for selecting export options.

    Provides a GUI for selecting export format, filename, and quality settings.
    """

    def __init__(
        self,
        supported_formats: List[str],
        config: HBATConfig,
        initial_directory: Optional[str] = None,
    ) -> None:
        """Initialize export dialog.

        :param supported_formats: List of supported format names
        :type supported_formats: List[str]
        :param config: HBAT configuration instance
        :type config: HBATConfig
        :param initial_directory: Initial directory for file dialog
        :type initial_directory: Optional[str]
        """
        self.supported_formats = supported_formats
        self.config = config
        self.initial_directory = initial_directory or os.getcwd()
        self.result: Optional[Tuple[str, str]] = None

        # Dialog widgets
        self.dialog: Optional[tk.Toplevel] = None
        self.format_var: Optional[tk.StringVar] = None
        self.resolution_var: Optional[tk.IntVar] = None

    def show(self) -> Optional[Tuple[str, str]]:
        """Show the export dialog.

        :returns: Tuple of (filename, format) or None if cancelled
        :rtype: Optional[Tuple[str, str]]
        """
        # Quick export if only one format supported
        if len(self.supported_formats) == 1:
            format = self.supported_formats[0]
            filename = self._show_file_dialog(format)
            if filename:
                return filename, format
            return None

        # Show full dialog for multiple formats
        return self._show_full_dialog()

    def _show_full_dialog(self) -> Optional[Tuple[str, str]]:
        """Show full export options dialog.

        :returns: Tuple of (filename, format) or None if cancelled
        :rtype: Optional[Tuple[str, str]]
        """
        # Create dialog window
        self.dialog = tk.Toplevel()
        self.dialog.title("Export Visualization")
        self.dialog.geometry("400x300")
        self.dialog.resizable(False, False)
        self.dialog.grab_set()  # Make dialog modal

        # Center the dialog
        self.dialog.update_idletasks()
        x = (self.dialog.winfo_screenwidth() - self.dialog.winfo_width()) // 2
        y = (self.dialog.winfo_screenheight() - self.dialog.winfo_height()) // 2
        self.dialog.geometry(f"+{x}+{y}")

        self._create_dialog_widgets()

        # Wait for dialog to close
        self.dialog.wait_window()

        return self.result

    def _create_dialog_widgets(self) -> None:
        """Create dialog widgets."""
        if not self.dialog:
            return

        # Main frame
        main_frame = ttk.Frame(self.dialog, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Format selection
        format_frame = ttk.LabelFrame(main_frame, text="Export Format", padding="10")
        format_frame.pack(fill=tk.X, pady=(0, 10))

        self.format_var = tk.StringVar(value=self.supported_formats[0])

        for format in self.supported_formats:
            format_info = EXPORT_FORMATS.get(format, {})
            text = f"{format_info.get('name', format.upper())} ({format})"
            desc = format_info.get("description", "")
            if desc:
                text += f" - {desc}"

            ttk.Radiobutton(
                format_frame, text=text, variable=self.format_var, value=format
            ).pack(anchor=tk.W, pady=2)

        # Quality settings
        quality_frame = ttk.LabelFrame(
            main_frame, text="Quality Settings", padding="10"
        )
        quality_frame.pack(fill=tk.X, pady=(0, 10))

        # Resolution setting
        resolution_frame = ttk.Frame(quality_frame)
        resolution_frame.pack(fill=tk.X, pady=2)

        ttk.Label(resolution_frame, text="Resolution (DPI):").pack(side=tk.LEFT)

        self.resolution_var = tk.IntVar(
            value=self.config.get_preference("export_dpi", 300)
        )
        resolution_spinbox = ttk.Spinbox(
            resolution_frame,
            from_=72,
            to=600,
            width=10,
            textvariable=self.resolution_var,
        )
        resolution_spinbox.pack(side=tk.RIGHT)

        # Buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=(10, 0))

        ttk.Button(button_frame, text="Cancel", command=self._cancel).pack(
            side=tk.RIGHT, padx=(5, 0)
        )

        ttk.Button(button_frame, text="Export...", command=self._export).pack(
            side=tk.RIGHT
        )

    def _export(self) -> None:
        """Handle export button click."""
        if not self.format_var:
            return

        format = self.format_var.get()

        # Update resolution preference
        if self.resolution_var:
            resolution = self.resolution_var.get()
            self.config.set_preference("export_dpi", resolution)

        # Show file dialog
        filename = self._show_file_dialog(format)

        if filename:
            self.result = (filename, format)

        if self.dialog:
            self.dialog.destroy()

    def _cancel(self) -> None:
        """Handle cancel button click."""
        self.result = None
        if self.dialog:
            self.dialog.destroy()

    def _show_file_dialog(self, format: str) -> Optional[str]:
        """Show file save dialog for specific format.

        :param format: Export format
        :type format: str
        :returns: Selected filename or None
        :rtype: Optional[str]
        """
        format_info = EXPORT_FORMATS.get(format, {})

        # Prepare file dialog options
        filetypes = format_info.get(
            "filetypes", [(f"{format.upper()} files", f"*.{format}")]
        )
        filetypes.append(("All files", "*.*"))

        # Show save dialog
        filename = filedialog.asksaveasfilename(
            title=f"Export as {format_info.get('name', format.upper())}",
            initialdir=self.initial_directory,
            filetypes=filetypes,
            defaultextension=format_info.get("extension", f".{format}"),
        )

        return filename if filename else None


def show_quick_export_dialog(
    renderer: VisualizationRenderer, config: HBATConfig
) -> bool:
    """Show a quick export dialog for immediate export.

    :param renderer: Visualization renderer
    :type renderer: VisualizationRenderer
    :param config: HBAT configuration instance
    :type config: HBATConfig
    :returns: True if export successful
    :rtype: bool
    """
    export_manager = ExportManager(renderer, config)
    return export_manager.export_visualization()
