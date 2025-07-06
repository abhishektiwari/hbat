"""
Misc. HBAT Constants and Default Parameters
"""

from .app import APP_VERSION


# GUI defaults
class GUIDefaults:
    """Default values for GUI interface."""

    # Window settings
    WINDOW_WIDTH = 1800
    WINDOW_HEIGHT = 900
    MIN_WINDOW_WIDTH = 1200
    MIN_WINDOW_HEIGHT = 800

    # Layout settings
    LEFT_PANEL_WIDTH = 400  # Initial pane position

    # Progress bar settings
    PROGRESS_BAR_INTERVAL = 10  # milliseconds


# Vector mathematics defaults
class VectorDefaults:
    """Default values for vector operations."""

    DEFAULT_X = 0.0
    DEFAULT_Y = 0.0
    DEFAULT_Z = 0.0


# File format constants
class FileFormats:
    """Supported file formats and extensions."""

    PDB_EXTENSIONS = [".pdb"]
    OUTPUT_EXTENSIONS = [".txt", ".csv", ".json"]

    # Export format defaults
    JSON_VERSION = APP_VERSION
