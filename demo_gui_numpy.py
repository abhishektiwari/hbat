#!/usr/bin/env python3
"""
Demo script showing HBAT GUI with NumPy analyzer integration.

This script demonstrates how to use the HBAT GUI with the new NumPy-optimized
analyzer for enhanced performance on large molecular structures.
"""

import os
import sys

# Add the current directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hbat.gui.main_window import MainWindow


def main():
    """Run the HBAT GUI demo with NumPy integration."""
    print("🧬 HBAT GUI Demo - NumPy Integration")
    print("=" * 50)
    print()
    print("Features demonstrated:")
    print("• NumPy-optimized analyzer for high performance")
    print("• Toggle between standard and NumPy analyzers")
    print("• Performance indicator in toolbar")
    print("• Enhanced status reporting")
    print()
    print("📋 Usage Instructions:")
    print("1. Click 'File' → 'Open PDB File...' to load a structure")
    print("2. Use 'Settings' → 'Use NumPy Analyzer' to toggle analysis engines")
    print("3. Click 'Analysis' → 'Run Analysis' or press F5")
    print("4. View results in the right panel")
    print("5. Export results using 'File' → 'Save Results...'")
    print()
    print("💡 Performance Tips:")
    print("• NumPy analyzer is enabled by default for best performance")
    print("• For large structures (>1000 atoms), expect 10-30x speedup")
    print("• Watch the performance indicator during analysis")
    print()
    print("🚀 Starting GUI...")
    print()

    try:
        # Create and run the main window
        app = MainWindow()
        app.run()

    except KeyboardInterrupt:
        print("\n👋 GUI closed by user")
    except Exception as e:
        print(f"❌ Error starting GUI: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
