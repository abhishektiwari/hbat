#!/usr/bin/env python3
"""
Build script for HBAT standalone Windows executables.

This script creates standalone Windows executables using PyInstaller.
Run this from the project root directory on Windows.
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
    """Build GUI executable for Windows."""
    print("Building Windows GUI executable...")

    # Windows-specific PyInstaller options
    cmd = [
        "pyinstaller",
        "--onefile",  # Single file executable for Windows
        "--windowed",
        "--name",
        "HBAT-GUI",
        "--icon",
        "hbat.ico",
        "--add-data",
        "example_pdb_files;example_pdb_files",  # Windows uses semicolon
        "--add-data",
        "example_presets;example_presets",
        "--add-data",
        "hbat.png;.",
        "--add-data",
        "README.md;.",
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
        "--distpath",
        "dist/windows",
        "hbat_gui.py",
    ]

    try:
        subprocess.run(cmd, check=True)
        print("✓ Windows GUI executable built successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ Windows GUI build failed: {e}")
        return False


def build_cli():
    """Build CLI executable for Windows."""
    print("Building Windows CLI executable...")

    cmd = [
        "pyinstaller",
        "--onefile",
        "--console",  # Console application
        "--name",
        "hbat",
        "--icon",
        "hbat.ico",
        "--add-data",
        "example_pdb_files;example_pdb_files",
        "--add-data",
        "example_presets;example_presets",
        "--hidden-import",
        "matplotlib",
        "--hidden-import",
        "networkx",
        "--clean",
        "--distpath",
        "dist/windows",
        "hbat_cli.py",
    ]

    try:
        subprocess.run(cmd, check=True)
        print("✓ Windows CLI executable built successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ Windows CLI build failed: {e}")
        return False


def create_installer():
    """Create Windows installer package."""
    print("\nCreating Windows installer...")

    # Check if NSIS is available
    nsis_path = shutil.which("makensis")
    if not nsis_path:
        print("NSIS not found. Skipping installer creation.")
        print("To create installers, install NSIS from https://nsis.sourceforge.io/")
        return False

    # Create NSIS script
    nsis_script = """
!define PRODUCT_NAME "HBAT"
!define PRODUCT_VERSION "1.0.0"
!define PRODUCT_PUBLISHER "HBAT Team"

Name "${PRODUCT_NAME}"
OutFile "dist\\HBAT-Setup.exe"
InstallDir "$PROGRAMFILES\\${PRODUCT_NAME}"
RequestExecutionLevel admin

Section "Main"
    SetOutPath "$INSTDIR"
    File "dist\\windows\\HBAT-GUI.exe"
    File "dist\\windows\\hbat.exe"
    
    CreateDirectory "$SMPROGRAMS\\${PRODUCT_NAME}"
    CreateShortcut "$SMPROGRAMS\\${PRODUCT_NAME}\\HBAT GUI.lnk" "$INSTDIR\\HBAT-GUI.exe"
    CreateShortcut "$SMPROGRAMS\\${PRODUCT_NAME}\\Uninstall.lnk" "$INSTDIR\\uninstall.exe"
    
    WriteUninstaller "$INSTDIR\\uninstall.exe"
SectionEnd

Section "Uninstall"
    Delete "$INSTDIR\\HBAT-GUI.exe"
    Delete "$INSTDIR\\hbat.exe"
    Delete "$INSTDIR\\uninstall.exe"
    
    Delete "$SMPROGRAMS\\${PRODUCT_NAME}\\HBAT GUI.lnk"
    Delete "$SMPROGRAMS\\${PRODUCT_NAME}\\Uninstall.lnk"
    RMDir "$SMPROGRAMS\\${PRODUCT_NAME}"
    
    RMDir "$INSTDIR"
SectionEnd
"""

    with open("installer.nsi", "w") as f:
        f.write(nsis_script)

    try:
        subprocess.run(["makensis", "installer.nsi"], check=True)
        os.remove("installer.nsi")
        print("✓ Windows installer created successfully")
        return True
    except subprocess.CalledProcessError:
        print("✗ Failed to create Windows installer")
        if os.path.exists("installer.nsi"):
            os.remove("installer.nsi")
        return False


def main():
    """Main build function."""
    print("HBAT Windows Build Script")
    print("=" * 40)

    # Check we're in the right directory
    if not os.path.exists("hbat_gui.py"):
        print("Error: Please run this script from the HBAT project root directory")
        return 1

    # Check if running on Windows
    if sys.platform != "win32":
        print("Warning: This script is designed for Windows.")
        print("Cross-compilation may not work properly.")

    # Install dependencies
    if not install_dependencies():
        return 1

    # Clean previous builds
    clean_build()

    # Build executables
    gui_success = build_gui()
    cli_success = build_cli()

    if not (gui_success or cli_success):
        print("All builds failed!")
        return 1

    # Try to create installer
    installer_success = create_installer()

    print("\n" + "=" * 40)
    print("Build Summary")
    print("=" * 40)

    if gui_success:
        print("✓ GUI: dist/windows/HBAT-GUI.exe")
    else:
        print("✗ GUI build failed")

    if cli_success:
        print("✓ CLI: dist/windows/hbat.exe")
    else:
        print("✗ CLI build failed")

    if installer_success:
        print("✓ Installer: dist/HBAT-Setup.exe")
    else:
        print("✗ Installer creation skipped/failed")

    print("\nUsage:")
    if gui_success:
        print("  GUI: dist\\windows\\HBAT-GUI.exe")
    if cli_success:
        print("  CLI: dist\\windows\\hbat.exe")

    return 0


if __name__ == "__main__":
    sys.exit(main())
