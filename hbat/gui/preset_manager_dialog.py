"""
Preset Manager dialog for HBAT GUI.

This module provides a dedicated dialog for managing analysis parameter presets,
including loading, saving, and organizing preset files.
"""

import json
import os
import tkinter as tk
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Any, Callable, Dict, Optional

from ..constants.parameters import ParametersDefault
from ..core.analysis import AnalysisParameters


class PresetManagerDialog:
    """Dialog for managing analysis parameter presets."""

    def __init__(
        self, parent: tk.Tk, current_params: Optional[AnalysisParameters] = None
    ):
        """Initialize preset manager dialog.

        :param parent: Parent window
        :type parent: tk.Tk
        :param current_params: Current analysis parameters
        :type current_params: Optional[AnalysisParameters]
        """
        self.parent = parent
        self.current_params = current_params or AnalysisParameters()
        self.selected_preset_data = None
        self.result = None
        self.preset_file_paths = {}  # Store mapping of item_id to file path

        # Create dialog window
        self.dialog = tk.Toplevel(parent)
        self.dialog.title("Preset Manager")
        self.dialog.geometry("800x600")
        self.dialog.resizable(True, True)

        # Make dialog modal
        self.dialog.transient(parent)
        self.dialog.grab_set()

        # Create widgets
        self._create_widgets()

        # Center the dialog
        self.dialog.update_idletasks()
        x = (self.dialog.winfo_screenwidth() // 2) - (self.dialog.winfo_width() // 2)
        y = (self.dialog.winfo_screenheight() // 2) - (self.dialog.winfo_height() // 2)
        self.dialog.geometry(f"+{x}+{y}")

        # Handle window closing
        self.dialog.protocol("WM_DELETE_WINDOW", self._cancel)

        # Load available presets
        self._refresh_preset_list()

    def _create_widgets(self):
        """Create dialog widgets."""
        # Main frame with padding
        main_frame = ttk.Frame(self.dialog, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Title
        title_label = ttk.Label(
            main_frame,
            text="HBAT Analysis Parameter Presets",
            font=("TkDefaultFont", 12, "bold"),
        )
        title_label.pack(pady=(0, 20))

        # Create notebook for different sections
        notebook = ttk.Notebook(main_frame)
        notebook.pack(fill=tk.BOTH, expand=True, pady=(0, 20))

        # Load presets tab
        self._create_load_tab(notebook)

        # Save presets tab
        self._create_save_tab(notebook)

        # Manage presets tab
        self._create_manage_tab(notebook)

        # Button frame
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X)

        ttk.Button(
            button_frame, text="Load Selected", command=self._load_selected_preset
        ).pack(side=tk.RIGHT, padx=(5, 0))

        ttk.Button(button_frame, text="Cancel", command=self._cancel).pack(
            side=tk.RIGHT
        )

    def _create_load_tab(self, notebook):
        """Create the load presets tab."""
        load_frame = ttk.Frame(notebook, padding="10")
        notebook.add(load_frame, text="Load Preset")

        # Instructions
        instructions = ttk.Label(
            load_frame,
            text="Select a preset to load parameters:",
            font=("TkDefaultFont", 10),
        )
        instructions.pack(anchor=tk.W, pady=(0, 10))

        # Preset list frame
        list_frame = ttk.Frame(load_frame)
        list_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # Treeview for presets
        columns = ("Name", "Description", "Date", "Location")
        self.preset_tree = ttk.Treeview(
            list_frame, columns=columns, show="tree headings", height=12
        )

        # Configure columns
        self.preset_tree.heading("#0", text="", anchor=tk.W)
        self.preset_tree.column("#0", width=30, minwidth=30)

        for col in columns:
            self.preset_tree.heading(col, text=col, anchor=tk.W)

        self.preset_tree.column("Name", width=150, minwidth=100)
        self.preset_tree.column("Description", width=200, minwidth=150)
        self.preset_tree.column("Date", width=120, minwidth=100)
        self.preset_tree.column("Location", width=100, minwidth=80)

        # Scrollbars for treeview
        v_scrollbar = ttk.Scrollbar(
            list_frame, orient=tk.VERTICAL, command=self.preset_tree.yview
        )
        h_scrollbar = ttk.Scrollbar(
            list_frame, orient=tk.HORIZONTAL, command=self.preset_tree.xview
        )
        self.preset_tree.configure(
            yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set
        )

        # Grid layout for treeview and scrollbars
        self.preset_tree.grid(row=0, column=0, sticky="nsew")
        v_scrollbar.grid(row=0, column=1, sticky="ns")
        h_scrollbar.grid(row=1, column=0, sticky="ew")

        list_frame.grid_rowconfigure(0, weight=1)
        list_frame.grid_columnconfigure(0, weight=1)

        # Bind selection event
        self.preset_tree.bind("<<TreeviewSelect>>", self._on_preset_select)
        self.preset_tree.bind("<Double-1>", lambda e: self._load_selected_preset())

        # Action buttons
        action_frame = ttk.Frame(load_frame)
        action_frame.pack(fill=tk.X)

        ttk.Button(
            action_frame, text="Browse...", command=self._browse_for_preset
        ).pack(side=tk.LEFT, padx=(0, 5))

        ttk.Button(
            action_frame, text="Refresh", command=self._refresh_preset_list
        ).pack(side=tk.LEFT, padx=5)

    def _create_save_tab(self, notebook):
        """Create the save presets tab."""
        save_frame = ttk.Frame(notebook, padding="10")
        notebook.add(save_frame, text="Save Preset")

        # Instructions
        instructions = ttk.Label(
            save_frame,
            text="Save current parameters as a preset:",
            font=("TkDefaultFont", 10),
        )
        instructions.pack(anchor=tk.W, pady=(0, 20))

        # Preset info frame
        info_frame = ttk.LabelFrame(save_frame, text="Preset Information", padding="10")
        info_frame.pack(fill=tk.X, pady=(0, 20))

        # Name entry
        ttk.Label(info_frame, text="Name:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.name_var = tk.StringVar(
            value=f"Custom Preset {datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )
        name_entry = ttk.Entry(info_frame, textvariable=self.name_var, width=40)
        name_entry.grid(row=0, column=1, sticky=tk.EW, padx=(10, 0), pady=2)

        # Description entry
        ttk.Label(info_frame, text="Description:").grid(
            row=1, column=0, sticky=tk.NW, pady=2
        )
        self.description_var = tk.StringVar(value="Custom HBAT Analysis Parameters")
        desc_entry = ttk.Entry(info_frame, textvariable=self.description_var, width=40)
        desc_entry.grid(row=1, column=1, sticky=tk.EW, padx=(10, 0), pady=2)

        info_frame.grid_columnconfigure(1, weight=1)

        # Current parameters preview
        preview_frame = ttk.LabelFrame(
            save_frame, text="Current Parameters", padding="10"
        )
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 20))

        # Create text widget for parameters preview
        self.preview_text = tk.Text(
            preview_frame, height=10, width=60, font=("Courier", 9)
        )
        preview_scrollbar = ttk.Scrollbar(
            preview_frame, orient=tk.VERTICAL, command=self.preview_text.yview
        )
        self.preview_text.configure(yscrollcommand=preview_scrollbar.set)

        self.preview_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        preview_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Update preview
        self._update_parameters_preview()

        # Save button
        ttk.Button(
            save_frame, text="Save Preset", command=self._save_current_preset
        ).pack(pady=(0, 10))

    def _create_manage_tab(self, notebook):
        """Create the manage presets tab."""
        manage_frame = ttk.Frame(notebook, padding="10")
        notebook.add(manage_frame, text="Manage")

        # Instructions
        instructions = ttk.Label(
            manage_frame, text="Manage your preset files:", font=("TkDefaultFont", 10)
        )
        instructions.pack(anchor=tk.W, pady=(0, 20))

        # Directory info
        dir_frame = ttk.LabelFrame(
            manage_frame, text="Preset Directories", padding="10"
        )
        dir_frame.pack(fill=tk.X, pady=(0, 20))

        # User presets directory
        user_dir = self._get_presets_directory()
        ttk.Label(dir_frame, text="User Presets:").pack(anchor=tk.W)
        user_dir_frame = ttk.Frame(dir_frame)
        user_dir_frame.pack(fill=tk.X, pady=(2, 10))

        ttk.Label(user_dir_frame, text=user_dir, font=("Courier", 8)).pack(
            side=tk.LEFT, fill=tk.X, expand=True
        )
        ttk.Button(
            user_dir_frame, text="Open", command=lambda: self._open_directory(user_dir)
        ).pack(side=tk.RIGHT)

        # Example presets directory
        example_dir = self._get_example_presets_directory()
        ttk.Label(dir_frame, text="Example Presets:").pack(anchor=tk.W)
        example_dir_frame = ttk.Frame(dir_frame)
        example_dir_frame.pack(fill=tk.X, pady=2)

        ttk.Label(example_dir_frame, text=example_dir, font=("Courier", 8)).pack(
            side=tk.LEFT, fill=tk.X, expand=True
        )
        ttk.Button(
            example_dir_frame,
            text="Open",
            command=lambda: self._open_directory(example_dir),
        ).pack(side=tk.RIGHT)

        # Actions
        actions_frame = ttk.LabelFrame(manage_frame, text="Actions", padding="10")
        actions_frame.pack(fill=tk.X)

        ttk.Button(
            actions_frame,
            text="Open Preset Directory",
            command=lambda: self._open_directory(user_dir),
        ).pack(anchor=tk.W, pady=2)

        ttk.Button(
            actions_frame, text="Import Preset...", command=self._import_preset
        ).pack(anchor=tk.W, pady=2)

    def _refresh_preset_list(self):
        """Refresh the list of available presets."""
        # Clear existing items and file path mapping
        for item in self.preset_tree.get_children():
            self.preset_tree.delete(item)
        self.preset_file_paths.clear()

        # Load presets from different locations
        self._load_presets_from_directory("User", self._get_presets_directory())
        self._load_presets_from_directory(
            "Examples", self._get_example_presets_directory()
        )

    def _load_presets_from_directory(self, category: str, directory: str):
        """Load presets from a specific directory."""
        if not os.path.exists(directory):
            return

        # Create category node
        category_node = self.preset_tree.insert(
            "", "end", text=category, values=("", "", "", "")
        )

        try:
            for filename in os.listdir(directory):
                if filename.endswith((".hbat", ".json")):
                    filepath = os.path.join(directory, filename)
                    try:
                        with open(filepath, "r") as f:
                            data = json.load(f)

                        # Extract preset info
                        name = data.get("name", os.path.splitext(filename)[0])
                        description = data.get("description", "")
                        created = data.get("created", "")

                        # Format date
                        if created:
                            try:
                                date_obj = datetime.fromisoformat(
                                    created.replace("Z", "+00:00")
                                )
                                formatted_date = date_obj.strftime("%Y-%m-%d %H:%M")
                            except:
                                formatted_date = created
                        else:
                            # Use file modification time
                            mtime = os.path.getmtime(filepath)
                            formatted_date = datetime.fromtimestamp(mtime).strftime(
                                "%Y-%m-%d %H:%M"
                            )

                        # Insert preset item
                        item_id = self.preset_tree.insert(
                            category_node,
                            "end",
                            text="📄",
                            values=(name, description, formatted_date, category),
                        )

                        # Store file path in mapping for later use
                        self.preset_file_paths[item_id] = filepath

                    except Exception as e:
                        print(f"Error loading preset {filename}: {e}")

        except Exception as e:
            print(f"Error reading directory {directory}: {e}")

        # Expand category node
        self.preset_tree.item(category_node, open=True)

    def _on_preset_select(self, event):
        """Handle preset selection."""
        selection = self.preset_tree.selection()
        if not selection:
            return

        item = selection[0]
        filepath = self.preset_file_paths.get(item)

        # Only process if it's a preset file (not a category)
        if filepath and os.path.exists(filepath):
            try:
                with open(filepath, "r") as f:
                    self.selected_preset_data = json.load(f)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load preset: {str(e)}")
                self.selected_preset_data = None

    def _load_selected_preset(self):
        """Load the selected preset."""
        if not self.selected_preset_data:
            messagebox.showwarning("Warning", "Please select a preset to load.")
            return

        self.result = self.selected_preset_data
        self.dialog.destroy()

    def _browse_for_preset(self):
        """Browse for a preset file."""
        filename = filedialog.askopenfilename(
            title="Load Preset File",
            filetypes=[
                ("HBAT Presets", "*.hbat"),
                ("JSON files", "*.json"),
                ("All files", "*.*"),
            ],
        )

        if filename:
            try:
                with open(filename, "r") as f:
                    data = json.load(f)
                self.selected_preset_data = data
                self.result = data
                self.dialog.destroy()
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load preset: {str(e)}")

    def _save_current_preset(self):
        """Save current parameters as a preset."""
        name = self.name_var.get().strip()
        if not name:
            messagebox.showwarning("Warning", "Please enter a name for the preset.")
            return

        # Get save filename
        user_dir = self._get_presets_directory()
        filename = filedialog.asksaveasfilename(
            title="Save Preset",
            initialdir=user_dir,
            initialfile=f"{name}.hbat",
            defaultextension=".hbat",
            filetypes=[("HBAT Presets", "*.hbat"), ("JSON files", "*.json")],
        )

        if filename:
            try:
                preset_data = self._create_preset_data()
                preset_data["name"] = name
                preset_data["description"] = self.description_var.get()

                with open(filename, "w") as f:
                    json.dump(preset_data, f, indent=2)

                messagebox.showinfo(
                    "Success",
                    f"Preset saved successfully to:\n{os.path.basename(filename)}",
                )
                self._refresh_preset_list()

            except Exception as e:
                messagebox.showerror("Error", f"Failed to save preset: {str(e)}")

    def _import_preset(self):
        """Import a preset file to the user directory."""
        filename = filedialog.askopenfilename(
            title="Import Preset File",
            filetypes=[
                ("HBAT Presets", "*.hbat"),
                ("JSON files", "*.json"),
                ("All files", "*.*"),
            ],
        )

        if filename:
            try:
                # Copy to user directory
                user_dir = self._get_presets_directory()
                dest_filename = os.path.join(user_dir, os.path.basename(filename))

                import shutil

                shutil.copy2(filename, dest_filename)

                messagebox.showinfo(
                    "Success", f"Preset imported successfully to:\n{dest_filename}"
                )
                self._refresh_preset_list()

            except Exception as e:
                messagebox.showerror("Error", f"Failed to import preset: {str(e)}")

    def _open_directory(self, directory: str):
        """Open a directory in the file manager."""
        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)

        import subprocess
        import sys

        try:
            if sys.platform == "win32":
                os.startfile(directory)
            elif sys.platform == "darwin":
                subprocess.run(["open", directory])
            else:
                subprocess.run(["xdg-open", directory])
        except Exception as e:
            messagebox.showerror("Error", f"Failed to open directory: {str(e)}")

    def _update_parameters_preview(self):
        """Update the parameters preview text."""
        self.preview_text.delete(1.0, tk.END)

        params = self.current_params

        preview = f"""Hydrogen Bond Parameters:
  H...A Distance: {params.hb_distance_cutoff:.1f} Å
  D-H...A Angle: {params.hb_angle_cutoff:.0f}°
  D...A Distance: {params.hb_donor_acceptor_cutoff:.1f} Å

Weak Hydrogen Bond Parameters (C-H···O):
  H...A Distance: {params.whb_distance_cutoff:.1f} Å
  D-H...A Angle: {params.whb_angle_cutoff:.0f}°
  D...A Distance: {params.whb_donor_acceptor_cutoff:.1f} Å

Halogen Bond Parameters:
  X...A Distance: {params.xb_distance_cutoff:.1f} Å
  C-X...A Angle: {params.xb_angle_cutoff:.0f}°

π Interaction Parameters:
  General H...π Distance: {params.pi_distance_cutoff:.1f} Å
  General D-H...π Angle: {params.pi_angle_cutoff:.0f}°
  
  Hydrogen-π Interactions:
    C-H...π Distance: {params.pi_ch_distance_cutoff:.1f} Å, Angle: {params.pi_ch_angle_cutoff:.0f}°
    N-H...π Distance: {params.pi_nh_distance_cutoff:.1f} Å, Angle: {params.pi_nh_angle_cutoff:.0f}°
    O-H...π Distance: {params.pi_oh_distance_cutoff:.1f} Å, Angle: {params.pi_oh_angle_cutoff:.0f}°
    S-H...π Distance: {params.pi_sh_distance_cutoff:.1f} Å, Angle: {params.pi_sh_angle_cutoff:.0f}°

π-π Stacking Parameters:
  Distance: {params.pi_pi_distance_cutoff:.1f} Å
  Parallel Angle: {params.pi_pi_parallel_angle_cutoff:.0f}°
  T-shaped Angles: {params.pi_pi_tshaped_angle_min:.0f}°-{params.pi_pi_tshaped_angle_max:.0f}°
  Offset: {params.pi_pi_offset_cutoff:.1f} Å

Carbonyl Interaction Parameters:
  Distance: {params.carbonyl_distance_cutoff:.1f} Å
  Angles: {params.carbonyl_angle_min:.0f}°-{params.carbonyl_angle_max:.0f}°

n→π* Interaction Parameters:
  Distance: {params.n_pi_distance_cutoff:.1f} Å
  Sulfur Distance: {params.n_pi_sulfur_distance_cutoff:.1f} Å
  Angles: {params.n_pi_angle_min:.0f}°-{params.n_pi_angle_max:.0f}°
  
  Halogen-π Interactions:
    C-Cl...π Distance: {params.pi_ccl_distance_cutoff:.1f} Å, Angle: {params.pi_ccl_angle_cutoff:.0f}°
    C-Br...π Distance: {params.pi_cbr_distance_cutoff:.1f} Å, Angle: {params.pi_cbr_angle_cutoff:.0f}°
    C-I...π Distance: {params.pi_ci_distance_cutoff:.1f} Å, Angle: {params.pi_ci_angle_cutoff:.0f}°

General Parameters:
  Covalent Bond Factor: {params.covalent_cutoff_factor:.2f}
  Analysis Mode: {params.analysis_mode}
"""

        # Add PDB fixing parameters if they exist
        if hasattr(params, "fix_pdb_enabled"):
            pdb_preview = f"""
PDB Fixing Parameters:
  Enabled: {params.fix_pdb_enabled}
  Method: {params.fix_pdb_method}
  Add Hydrogens: {params.fix_pdb_add_hydrogens}
  Add Heavy Atoms: {params.fix_pdb_add_heavy_atoms}
  Replace Nonstandard: {params.fix_pdb_replace_nonstandard}
  Remove Heterogens: {params.fix_pdb_remove_heterogens}
  Keep Water: {params.fix_pdb_keep_water}
"""
            preview += pdb_preview

        self.preview_text.insert(1.0, preview)
        self.preview_text.config(state=tk.DISABLED)

    def _create_preset_data(self) -> Dict[str, Any]:
        """Create preset data from current parameters."""
        params = self.current_params
        preset_data = {
            "format_version": "1.0",
            "application": "HBAT",
            "created": datetime.now().isoformat(),
            "description": self.description_var.get(),
            "parameters": {
                "hydrogen_bonds": {
                    "h_a_distance_cutoff": params.hb_distance_cutoff,
                    "dha_angle_cutoff": params.hb_angle_cutoff,
                    "d_a_distance_cutoff": params.hb_donor_acceptor_cutoff,
                },
                "weak_hydrogen_bonds": {
                    "h_a_distance_cutoff": params.whb_distance_cutoff,
                    "dha_angle_cutoff": params.whb_angle_cutoff,
                    "d_a_distance_cutoff": params.whb_donor_acceptor_cutoff,
                },
                "halogen_bonds": {
                    "x_a_distance_cutoff": params.xb_distance_cutoff,
                    "dxa_angle_cutoff": params.xb_angle_cutoff,
                },
                "pi_interactions": {
                    "h_pi_distance_cutoff": params.pi_distance_cutoff,
                    "dh_pi_angle_cutoff": params.pi_angle_cutoff,
                    # π interaction subtype parameters
                    "ccl_pi_distance_cutoff": params.pi_ccl_distance_cutoff,
                    "ccl_pi_angle_cutoff": params.pi_ccl_angle_cutoff,
                    "cbr_pi_distance_cutoff": params.pi_cbr_distance_cutoff,
                    "cbr_pi_angle_cutoff": params.pi_cbr_angle_cutoff,
                    "ci_pi_distance_cutoff": params.pi_ci_distance_cutoff,
                    "ci_pi_angle_cutoff": params.pi_ci_angle_cutoff,
                    "ch_pi_distance_cutoff": params.pi_ch_distance_cutoff,
                    "ch_pi_angle_cutoff": params.pi_ch_angle_cutoff,
                    "nh_pi_distance_cutoff": params.pi_nh_distance_cutoff,
                    "nh_pi_angle_cutoff": params.pi_nh_angle_cutoff,
                    "oh_pi_distance_cutoff": params.pi_oh_distance_cutoff,
                    "oh_pi_angle_cutoff": params.pi_oh_angle_cutoff,
                    "sh_pi_distance_cutoff": params.pi_sh_distance_cutoff,
                    "sh_pi_angle_cutoff": params.pi_sh_angle_cutoff,
                },
                "pi_pi_stacking": {
                    "distance_cutoff": params.pi_pi_distance_cutoff,
                    "parallel_angle_cutoff": params.pi_pi_parallel_angle_cutoff,
                    "tshaped_angle_min": params.pi_pi_tshaped_angle_min,
                    "tshaped_angle_max": params.pi_pi_tshaped_angle_max,
                    "offset_cutoff": params.pi_pi_offset_cutoff,
                },
                "carbonyl_interactions": {
                    "distance_cutoff": params.carbonyl_distance_cutoff,
                    "angle_min": params.carbonyl_angle_min,
                    "angle_max": params.carbonyl_angle_max,
                },
                "n_pi_interactions": {
                    "distance_cutoff": params.n_pi_distance_cutoff,
                    "sulfur_distance_cutoff": params.n_pi_sulfur_distance_cutoff,
                    "angle_min": params.n_pi_angle_min,
                    "angle_max": params.n_pi_angle_max,
                },
                "general": {
                    "covalent_cutoff_factor": params.covalent_cutoff_factor,
                    "analysis_mode": params.analysis_mode,
                },
            },
        }

        # Add PDB fixing parameters if they exist
        if hasattr(params, "fix_pdb_enabled"):
            preset_data["parameters"]["pdb_fixing"] = {
                "enabled": params.fix_pdb_enabled,
                "method": params.fix_pdb_method,
                "add_hydrogens": params.fix_pdb_add_hydrogens,
                "add_heavy_atoms": params.fix_pdb_add_heavy_atoms,
                "replace_nonstandard": params.fix_pdb_replace_nonstandard,
                "remove_heterogens": params.fix_pdb_remove_heterogens,
                "keep_water": params.fix_pdb_keep_water,
            }

        return preset_data

    def _get_presets_directory(self) -> str:
        """Get or create the user presets directory."""
        home_dir = os.path.expanduser("~")
        presets_dir = os.path.join(home_dir, ".hbat", "presets")
        os.makedirs(presets_dir, exist_ok=True)
        return presets_dir

    def _get_example_presets_directory(self) -> str:
        """Get the example presets directory relative to the package."""
        current_dir = os.path.dirname(os.path.abspath(__file__))
        package_root = os.path.dirname(os.path.dirname(current_dir))
        example_presets_dir = os.path.join(package_root, "example_presets")
        return example_presets_dir

    def _cancel(self):
        """Cancel and close dialog."""
        self.result = None
        self.dialog.destroy()

    def get_result(self) -> Optional[Dict[str, Any]]:
        """Get the dialog result.

        :returns: Selected preset data or None if cancelled
        :rtype: Optional[Dict[str, Any]]
        """
        self.dialog.wait_window()
        return self.result
