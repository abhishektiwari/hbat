#!/usr/bin/env python3
"""
HBAT GUI Application Entry Point

Launch the HBAT graphical user interface for interactive analysis
of hydrogen bonds, halogen bonds, and π interactions in PDB files.
"""

import os
import sys

# Add the current directory to the Python path to ensure local imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    """Launch the HBAT GUI application."""
    try:
        from hbat.gui.main_window import MainWindow

        app = MainWindow()
        app.run()

    except ImportError as e:
        print(f"Error: Missing required dependencies: {e}")
        print("Please ensure all required packages are installed.")
        sys.exit(1)
    except Exception as e:
        print(f"Error starting HBAT GUI: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
