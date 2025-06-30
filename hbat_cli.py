#!/usr/bin/env python3
"""
HBAT Command-Line Interface Entry Point

Run HBAT analysis from the command line for batch processing
and automated analysis workflows.
"""

import os
import sys

# Add the current directory to the Python path to ensure local imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    """Launch the HBAT CLI application."""
    try:
        from hbat.cli.main import main as cli_main

        return cli_main()

    except ImportError as e:
        print(f"Error: Missing required dependencies: {e}")
        print("Please ensure all required packages are installed.")
        return 1
    except Exception as e:
        print(f"Error starting HBAT CLI: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
