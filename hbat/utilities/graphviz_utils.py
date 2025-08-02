"""
GraphViz detection and utility functions for HBAT.

This module provides functionality to detect GraphViz installation,
check available engines, and validate GraphViz executables on the system.
"""

import logging
import os
import subprocess
import sys
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Set up logging
logger = logging.getLogger(__name__)

# GraphViz executable names to check
GRAPHVIZ_EXECUTABLES = [
    "dot",
    "neato",
    "fdp",
    "sfdp",
    "circo",
    "twopi",
    "osage",
    "patchwork",
]

# Common GraphViz installation paths by platform
PLATFORM_PATHS: Dict[str, List[str]] = {
    "win32": [
        r"C:\Program Files\Graphviz\bin",
        r"C:\Program Files (x86)\Graphviz\bin",
        r"C:\Graphviz\bin",
    ],
    "darwin": [
        "/usr/local/bin",
        "/opt/local/bin",
        "/opt/homebrew/bin",
        "/usr/bin",
    ],
    "linux": [
        "/usr/bin",
        "/usr/local/bin",
        "/opt/bin",
    ],
}


class GraphVizDetector:
    """Detects and validates GraphViz installation.

    This class provides static methods to check for GraphViz availability,
    version information, and available layout engines.
    """

    # Cache for detection results (cleared between sessions)
    _detection_cache: Optional[bool] = None
    _version_cache: Optional[str] = None
    _engines_cache: Optional[List[str]] = None

    @staticmethod
    @lru_cache(maxsize=1)
    def is_graphviz_available() -> bool:
        """Check if GraphViz executables are in PATH.

        :returns: True if GraphViz is available, False otherwise
        :rtype: bool
        """
        if GraphVizDetector._detection_cache is not None:
            return GraphVizDetector._detection_cache

        logger.debug("Detecting GraphViz installation...")

        # First try to find 'dot' in PATH
        result = GraphVizDetector._check_executable_in_path("dot")

        # If not found in PATH, check platform-specific locations
        if not result and sys.platform in PLATFORM_PATHS:
            for path in PLATFORM_PATHS[sys.platform]:
                if GraphVizDetector._check_executable_in_directory("dot", path):
                    result = True
                    break

        GraphVizDetector._detection_cache = result
        logger.info(
            f"GraphViz detection result: {'Available' if result else 'Not Available'}"
        )
        return result

    @staticmethod
    def get_graphviz_version() -> Optional[str]:
        """Get installed GraphViz version.

        :returns: Version string if available, None otherwise
        :rtype: Optional[str]
        """
        if GraphVizDetector._version_cache is not None:
            return GraphVizDetector._version_cache

        if not GraphVizDetector.is_graphviz_available():
            return None

        try:
            # Run 'dot -V' to get version information
            result = subprocess.run(
                ["dot", "-V"],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=5,
            )

            if result.returncode == 0:
                # GraphViz outputs version to stderr
                output = (
                    result.stdout.strip() if result.stdout else result.stderr.strip()
                )
                # Extract version from output like "dot - graphviz version 2.40.1"
                if "version" in output.lower():
                    version_parts = output.split()
                    for i, part in enumerate(version_parts):
                        if part.lower() == "version" and i + 1 < len(version_parts):
                            GraphVizDetector._version_cache = version_parts[i + 1]
                            return GraphVizDetector._version_cache

        except (subprocess.SubprocessError, OSError) as e:
            logger.warning(f"Failed to get GraphViz version: {e}")

        return None

    @staticmethod
    def get_available_engines() -> List[str]:
        """Get list of available GraphViz layout engines.

        :returns: List of available engine names
        :rtype: List[str]
        """
        if GraphVizDetector._engines_cache is not None:
            return GraphVizDetector._engines_cache

        if not GraphVizDetector.is_graphviz_available():
            return []

        available_engines = []

        for engine in GRAPHVIZ_EXECUTABLES:
            if GraphVizDetector._check_executable_in_path(engine):
                available_engines.append(engine)
            else:
                # Check platform-specific paths
                if sys.platform in PLATFORM_PATHS:
                    for path in PLATFORM_PATHS[sys.platform]:
                        if GraphVizDetector._check_executable_in_directory(
                            engine, path
                        ):
                            available_engines.append(engine)
                            break

        GraphVizDetector._engines_cache = available_engines
        logger.debug(f"Available GraphViz engines: {available_engines}")
        return available_engines

    @staticmethod
    def validate_engine(engine: str) -> bool:
        """Validate if a specific engine is available.

        :param engine: Engine name to validate (e.g., 'dot', 'neato')
        :type engine: str
        :returns: True if engine is available, False otherwise
        :rtype: bool
        """
        return engine in GraphVizDetector.get_available_engines()

    @staticmethod
    def get_engine_path(engine: str) -> Optional[str]:
        """Get the full path to a GraphViz engine executable.

        :param engine: Engine name (e.g., 'dot', 'neato')
        :type engine: str
        :returns: Full path to executable if found, None otherwise
        :rtype: Optional[str]
        """
        # Try to find in PATH first
        path = GraphVizDetector._which(engine)
        if path:
            return path

        # Check platform-specific locations
        if sys.platform in PLATFORM_PATHS:
            for directory in PLATFORM_PATHS[sys.platform]:
                full_path = Path(directory) / engine
                if sys.platform == "win32":
                    full_path = full_path.with_suffix(".exe")
                if full_path.exists() and full_path.is_file():
                    return str(full_path)

        return None

    @staticmethod
    def clear_cache() -> None:
        """Clear all cached detection results.

        Useful for re-detecting after GraphViz installation/removal.
        """
        GraphVizDetector._detection_cache = None
        GraphVizDetector._version_cache = None
        GraphVizDetector._engines_cache = None
        GraphVizDetector.is_graphviz_available.cache_clear()
        logger.debug("GraphViz detection cache cleared")

    # Private helper methods

    @staticmethod
    def _check_executable_in_path(executable: str) -> bool:
        """Check if an executable exists in PATH.

        :param executable: Name of the executable
        :type executable: str
        :returns: True if found in PATH
        :rtype: bool
        """
        return GraphVizDetector._which(executable) is not None

    @staticmethod
    def _check_executable_in_directory(executable: str, directory: str) -> bool:
        """Check if an executable exists in a specific directory.

        :param executable: Name of the executable
        :type executable: str
        :param directory: Directory path to check
        :type directory: str
        :returns: True if found in directory
        :rtype: bool
        """
        path = Path(directory) / executable
        if sys.platform == "win32":
            path = path.with_suffix(".exe")
        return path.exists() and path.is_file()

    @staticmethod
    def _which(executable: str) -> Optional[str]:
        """Find executable in PATH (cross-platform which command).

        :param executable: Name of the executable
        :type executable: str
        :returns: Full path if found, None otherwise
        :rtype: Optional[str]
        """
        # Python 3.3+ has shutil.which, but we implement our own for compatibility
        if sys.platform == "win32":
            executable = f"{executable}.exe"

        for path in os.environ.get("PATH", "").split(os.pathsep):
            full_path = Path(path) / executable
            if full_path.exists() and full_path.is_file():
                return str(full_path)

        return None


def get_graphviz_info() -> Dict[str, any]:
    """Get comprehensive GraphViz installation information.

    :returns: Dictionary with GraphViz info
    :rtype: Dict[str, any]
    """
    return {
        "available": GraphVizDetector.is_graphviz_available(),
        "version": GraphVizDetector.get_graphviz_version(),
        "engines": GraphVizDetector.get_available_engines(),
        "platform": sys.platform,
    }
