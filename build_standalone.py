#!/usr/bin/env python3
"""
Simple build script for HBAT standalone executables.

This script creates standalone executables using PyInstaller.
Run this from the project root directory.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def install_dependencies():
    """Install build dependencies."""
    print("Installing build dependencies...")
    try:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "pyinstaller>=5.0.0",
                "setuptools-scm>=6.2.0",
            ],
            check=True,
        )
        return True
    except subprocess.CalledProcessError:
        print("Failed to install dependencies")
        return False


def clean_build():
    """Clean previous build artifacts."""
    print("Cleaning previous builds...")
    paths_to_clean = ["build", "dist", "__pycache__"]

    for path in paths_to_clean:
        if os.path.exists(path):
            shutil.rmtree(path)
            print(f"  Removed {path}")


def build_gui():
    """Build GUI executable."""
    print("Building GUI executable...")

    cmd = [
        "pyinstaller",
        "--onedir",
        "--windowed",
        "--name",
        "HBAT-GUI",
        "--icon",
        "hbat.ico",
        "--add-data",
        "example_pdb_files:example_pdb_files",
        "--add-data",
        "example_presets:example_presets",
        "--add-data",
        "hbat.png:.",
        "--add-data",
        "README.md:.",
        "--hidden-import",
        "tkinter",
        "--hidden-import",
        "matplotlib.backends.backend_tkagg",
        "--hidden-import",
        "networkx",
        "--exclude-module",
        "PyQt5",
        "--exclude-module",
        "PyQt6",
        "--clean",
        "hbat_gui.py",
    ]

    try:
        subprocess.run(cmd, check=True)
        print("✓ GUI executable built successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ GUI build failed: {e}")
        return False


def main():
    """Main build function."""
    print("HBAT Standalone Build Script")
    print("=" * 40)

    # Check we're in the right directory
    if not os.path.exists("hbat_gui.py"):
        print("Error: Please run this script from the HBAT project root directory")
        return 1

    # Install dependencies
    if not install_dependencies():
        return 1

    # Clean previous builds
    clean_build()

    # Build executables
    gui_success = build_gui()

    if not (gui_success):
        print("All builds failed!")
        return 1

    print("\n" + "=" * 40)
    print("Build Summary")
    print("=" * 40)

    if gui_success:
        print("✓ GUI: dist/HBAT-GUI/")
    else:
        print("✗ GUI build failed")

    print("\nUsage:")
    if gui_success:
        print("  GUI: ./dist/HBAT-GUI/HBAT-GUI")

    print("\nTo distribute:")
    print("  1. Test the executables")
    print("  2. Zip the entire project folder")
    print("  3. Share with end users")

    return 0


if __name__ == "__main__":
    sys.exit(main())
