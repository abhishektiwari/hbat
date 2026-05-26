"""
Results display panel for HBAT analysis.

This module provides GUI components for displaying analysis results
including hydrogen bonds, halogen bonds, and π interactions.
"""

import math
import tkinter as tk
from tkinter import messagebox, ttk
from typing import Optional

try:
    from .chain_visualization_window import ChainVisualizationWindow

    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False

from ..core.analysis import MolecularInteractionAnalyzer
from ..config import (
    INTERACTION_CONFIGS,
    extract_interaction_data,
    get_interaction_config,
)


class ResultsPanel:
    """Panel for displaying analysis results.

    This class provides a tabbed interface for viewing different types
    of molecular interaction results including summaries, detailed lists,
    and statistical analysis.

    :param parent: Parent widget to contain this panel
    :type parent: tkinter widget
    """

    def __init__(self, parent) -> None:
        """Initialize the results panel.

        Creates a complete results display interface with multiple tabs
        for different views of analysis results.

        :param parent: Parent widget
        :type parent: tkinter widget
        :returns: None
        :rtype: None
        """
        self.parent = parent
        self.analyzer: Optional[MolecularInteractionAnalyzer] = None
        self._create_widgets()

    def _create_widgets(self):
        """Create result display widgets."""
        # Create main notebook for different result types
        self.notebook = ttk.Notebook(self.parent)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # Summary tab (special case - custom logic)
        self._create_summary_tab()

        # Ligand interactions tab (special case - custom logic)
        self._create_ligand_interactions_tab()

        # Create interaction type tabs dynamically from INTERACTION_CONFIGS
        for interaction_type, config in INTERACTION_CONFIGS.items():
            # Special case: cooperativity chains has visualization button
            if interaction_type == "cooperativity_chains":
                self._create_cooperativity_chains_tab()
            else:
                self._create_interaction_tab(interaction_type, config)

    def _create_summary_tab(self):
        """Create summary results tab."""
        summary_frame = ttk.Frame(self.notebook)
        self.notebook.add(summary_frame, text="Summary")

        # Create text widget with scrollbars
        text_frame = ttk.Frame(summary_frame)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.summary_text = tk.Text(text_frame, wrap=tk.NONE, font=("Courier", 12))
        summary_v_scrollbar = ttk.Scrollbar(
            text_frame, orient=tk.VERTICAL, command=self.summary_text.yview
        )
        summary_h_scrollbar = ttk.Scrollbar(
            text_frame, orient=tk.HORIZONTAL, command=self.summary_text.xview
        )
        self.summary_text.configure(
            yscrollcommand=summary_v_scrollbar.set,
            xscrollcommand=summary_h_scrollbar.set,
        )

        # Use grid layout for proper scrollbar positioning
        self.summary_text.grid(row=0, column=0, sticky="nsew")
        summary_v_scrollbar.grid(row=0, column=1, sticky="ns")
        summary_h_scrollbar.grid(row=1, column=0, sticky="ew")

        text_frame.grid_rowconfigure(0, weight=1)
        text_frame.grid_columnconfigure(0, weight=1)

        # Configure text tags for formatting
        self.summary_text.tag_configure(
            "header", font=("Courier", 12, "bold"), foreground="blue"
        )
        self.summary_text.tag_configure("subheader", font=("Courier", 12, "bold"))
        self.summary_text.tag_configure("highlight", background="cyan")

    def _create_cooperativity_chains_tab(self):
        """Create cooperativity chains results tab."""
        coop_frame = ttk.Frame(self.notebook)
        self.notebook.add(coop_frame, text="Cooperativity Chains")

        # Create treeview for cooperativity chains
        tree_frame = ttk.Frame(coop_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        columns = ("chain_id", "chain_length", "chain_description")

        self.coop_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=15
        )

        # Configure columns
        self.coop_tree.heading("chain_id", text="Chain ID")
        self.coop_tree.heading("chain_length", text="Length")
        self.coop_tree.heading("chain_description", text="Chain Description")

        # Configure column widths
        self.coop_tree.column("chain_id", width=100)
        self.coop_tree.column("chain_length", width=100)
        self.coop_tree.column("chain_description", width=1000)

        # Add scrollbars
        coop_v_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.VERTICAL, command=self.coop_tree.yview
        )
        coop_h_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.HORIZONTAL, command=self.coop_tree.xview
        )
        self.coop_tree.configure(
            yscrollcommand=coop_v_scrollbar.set, xscrollcommand=coop_h_scrollbar.set
        )

        self.coop_tree.grid(row=0, column=0, sticky="nsew")
        coop_v_scrollbar.grid(row=0, column=1, sticky="ns")
        coop_h_scrollbar.grid(row=1, column=0, sticky="ew")

        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        # Bind double-click event to visualize chain
        self.coop_tree.bind("<Double-1>", self._on_chain_double_click)

        # Add info label
        info_frame = ttk.Frame(coop_frame)
        info_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(
            info_frame,
            text="Potential Cooperative Chains: Sequences where acceptors also act as donors",
        ).pack(side=tk.LEFT)

        # Add search functionality
        search_frame = ttk.Frame(coop_frame)
        search_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.coop_search_var = tk.StringVar()
        search_entry = ttk.Entry(
            search_frame, textvariable=self.coop_search_var, width=30
        )
        search_entry.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            search_frame,
            text="Filter",
            command=lambda: self._filter_results(
                self.coop_tree, self.coop_search_var.get()
            ),
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            search_frame,
            text="Clear",
            command=lambda: self._clear_filter(self.coop_tree, self.coop_search_var),
        ).pack(side=tk.LEFT, padx=5)

        # Add visualization button
        if VISUALIZATION_AVAILABLE:
            ttk.Button(
                search_frame,
                text="Visualize Selected Chain",
                command=self._visualize_selected_chain,
            ).pack(side=tk.RIGHT, padx=5)

    def _create_ligand_interactions_tab(self):
        """Create ligand interactions results tab with selector and dual tables."""
        lig_frame = ttk.Frame(self.notebook)
        self.notebook.add(lig_frame, text="Ligand Interactions")

        # Ligand selector at top
        selector_frame = ttk.Frame(lig_frame)
        selector_frame.pack(fill=tk.X, padx=10, pady=10)

        ttk.Label(selector_frame, text="Select Ligand:").pack(side=tk.LEFT, padx=5)
        self.lig_selector_var = tk.StringVar()
        self.lig_selector_combo = ttk.Combobox(
            selector_frame,
            textvariable=self.lig_selector_var,
            width=40,
            state="readonly",
        )
        self.lig_selector_combo.pack(side=tk.LEFT, padx=5)
        self.lig_selector_combo.bind("<<ComboboxSelected>>", self._on_ligand_selected)

        # Create paned window to hold both tables
        paned = ttk.PanedWindow(lig_frame, orient=tk.VERTICAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Regular interactions section
        lig_inter_frame = ttk.LabelFrame(paned, text="Regular Interactions", padding=5)
        paned.add(lig_inter_frame, weight=1)

        tree_frame1 = ttk.Frame(lig_inter_frame)
        tree_frame1.pack(fill=tk.BOTH, expand=True)

        lig_inter_columns = (
            "type",
            "donor_res",
            "donor_atom",
            "acceptor_res",
            "acceptor_atom",
            "distance",
            "properties",
        )

        self.lig_inter_tree = ttk.Treeview(
            tree_frame1, columns=lig_inter_columns, show="headings", height=10
        )

        self.lig_inter_tree.heading("type", text="Type")
        self.lig_inter_tree.heading("donor_res", text="Donor Residue")
        self.lig_inter_tree.heading("donor_atom", text="Donor Atom")
        self.lig_inter_tree.heading("acceptor_res", text="Acceptor Residue")
        self.lig_inter_tree.heading("acceptor_atom", text="Acceptor Atom")
        self.lig_inter_tree.heading("distance", text="Distance (Å)")
        self.lig_inter_tree.heading("properties", text="Properties")

        self.lig_inter_tree.column("type", width=100)
        self.lig_inter_tree.column("donor_res", width=120)
        self.lig_inter_tree.column("donor_atom", width=100)
        self.lig_inter_tree.column("acceptor_res", width=120)
        self.lig_inter_tree.column("acceptor_atom", width=100)
        self.lig_inter_tree.column("distance", width=80)
        self.lig_inter_tree.column("properties", width=80)

        lig_inter_v_scrollbar = ttk.Scrollbar(
            tree_frame1, orient=tk.VERTICAL, command=self.lig_inter_tree.yview
        )
        lig_inter_h_scrollbar = ttk.Scrollbar(
            tree_frame1, orient=tk.HORIZONTAL, command=self.lig_inter_tree.xview
        )
        self.lig_inter_tree.configure(
            yscrollcommand=lig_inter_v_scrollbar.set,
            xscrollcommand=lig_inter_h_scrollbar.set,
        )

        self.lig_inter_tree.grid(row=0, column=0, sticky="nsew")
        lig_inter_v_scrollbar.grid(row=0, column=1, sticky="ns")
        lig_inter_h_scrollbar.grid(row=1, column=0, sticky="ew")

        tree_frame1.grid_rowconfigure(0, weight=1)
        tree_frame1.grid_columnconfigure(0, weight=1)

        # Water bridges section
        lig_wb_frame = ttk.LabelFrame(paned, text="Water Bridges", padding=5)
        paned.add(lig_wb_frame, weight=1)

        tree_frame2 = ttk.Frame(lig_wb_frame)
        tree_frame2.pack(fill=tk.BOTH, expand=True)

        lig_wb_columns = ("start_res", "end_res", "hops", "water_residues", "distance")

        self.lig_wb_tree = ttk.Treeview(
            tree_frame2, columns=lig_wb_columns, show="headings", height=10
        )

        self.lig_wb_tree.heading("start_res", text="Start Residue")
        self.lig_wb_tree.heading("end_res", text="End Residue")
        self.lig_wb_tree.heading("hops", text="Hops")
        self.lig_wb_tree.heading("water_residues", text="Water Residues")
        self.lig_wb_tree.heading("distance", text="Distance (Å)")

        self.lig_wb_tree.column("start_res", width=120)
        self.lig_wb_tree.column("end_res", width=120)
        self.lig_wb_tree.column("hops", width=60)
        self.lig_wb_tree.column("water_residues", width=250)
        self.lig_wb_tree.column("distance", width=80)

        lig_wb_v_scrollbar = ttk.Scrollbar(
            tree_frame2, orient=tk.VERTICAL, command=self.lig_wb_tree.yview
        )
        lig_wb_h_scrollbar = ttk.Scrollbar(
            tree_frame2, orient=tk.HORIZONTAL, command=self.lig_wb_tree.xview
        )
        self.lig_wb_tree.configure(
            yscrollcommand=lig_wb_v_scrollbar.set,
            xscrollcommand=lig_wb_h_scrollbar.set,
        )

        self.lig_wb_tree.grid(row=0, column=0, sticky="nsew")
        lig_wb_v_scrollbar.grid(row=0, column=1, sticky="ns")
        lig_wb_h_scrollbar.grid(row=1, column=0, sticky="ew")

        tree_frame2.grid_rowconfigure(0, weight=1)
        tree_frame2.grid_columnconfigure(0, weight=1)

        # Store original data for filtering
        self.lig_inter_data = {}
        self.lig_wb_data = {}

    def _on_ligand_selected(self, event=None):
        """Handle ligand selection from dropdown."""
        selected_ligand = self.lig_selector_var.get()
        self._update_ligand_tables(selected_ligand)

    def _update_ligand_tables(self, selected_ligand):
        """Update ligand interactions and water bridges tables for selected ligand."""
        # Clear tables
        for item in self.lig_inter_tree.get_children():
            self.lig_inter_tree.delete(item)
        for item in self.lig_wb_tree.get_children():
            self.lig_wb_tree.delete(item)

        # Populate regular interactions
        if selected_ligand in self.lig_inter_data:
            for inter_data in self.lig_inter_data[selected_ligand]:
                self.lig_inter_tree.insert("", tk.END, values=inter_data)

        # Populate water bridges
        if selected_ligand in self.lig_wb_data:
            for wb_data in self.lig_wb_data[selected_ligand]:
                self.lig_wb_tree.insert("", tk.END, values=wb_data)

    def _create_interaction_tab(self, interaction_type: str, config):
        """Create a results tab dynamically from configuration.

        Generates a tab with treeview based on interaction type configuration,
        eliminating code duplication across specific interaction type methods.

        :param interaction_type: Interaction type key (e.g., 'hydrogen_bonds')
        :param config: InteractionConfig object with column definitions
        """
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text=config.label)

        # Create treeview for this interaction type
        tree_frame = ttk.Frame(frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Extract column names and create treeview
        column_names = tuple(col.name for col in config.columns)
        tree = ttk.Treeview(
            tree_frame, columns=column_names, show="headings", height=15
        )

        # Configure columns
        for col in config.columns:
            tree.heading(col.name, text=col.label)
            tree.column(col.name, width=col.width)

        # Add scrollbars
        v_scrollbar = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=tree.yview)
        h_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.HORIZONTAL, command=tree.xview
        )
        tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)

        tree.grid(row=0, column=0, sticky="nsew")
        v_scrollbar.grid(row=0, column=1, sticky="ns")
        h_scrollbar.grid(row=1, column=0, sticky="ew")

        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        # Store tree reference for update_results
        tree_var_name = f"{interaction_type}_tree"
        setattr(self, tree_var_name, tree)

        # Also store old-style variable names for backward compatibility with filter methods
        # Map: hydrogen_bonds -> hb_tree, halogen_bonds -> xb_tree, etc.
        old_name_map = {
            "hydrogen_bonds": ("hb_tree", "hb_search_var"),
            "water_bridges": ("wb_tree", "wb_search_var"),
            "halogen_bonds": ("xb_tree", "xb_search_var"),
            "pi_interactions": ("pi_tree", "pi_search_var"),
            "pi_pi_interactions": ("pi_pi_tree", "pi_pi_search_var"),
            "carbonyl_interactions": ("carbonyl_tree", "carbonyl_search_var"),
            "n_pi_interactions": ("n_pi_tree", "n_pi_search_var"),
        }

        # Add search functionality
        search_frame = ttk.Frame(frame)
        search_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        search_var = tk.StringVar()
        search_entry = ttk.Entry(search_frame, textvariable=search_var, width=30)
        search_entry.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            search_frame,
            text="Filter",
            command=lambda: self._filter_results(tree, search_var.get()),
        ).pack(side=tk.LEFT, padx=5)

        ttk.Button(
            search_frame,
            text="Clear",
            command=lambda: self._clear_filter(tree, search_var),
        ).pack(side=tk.LEFT, padx=5)

        # Store search variable for later access if needed
        setattr(self, f"{interaction_type}_search_var", search_var)

        # Store old-style names for backward compatibility
        if interaction_type in old_name_map:
            tree_old_name, search_var_old_name = old_name_map[interaction_type]
            setattr(self, tree_old_name, tree)
            setattr(self, search_var_old_name, search_var)

    def update_results(self, analyzer: MolecularInteractionAnalyzer) -> None:
        """Update the results panel with new analysis results.

        Refreshes all result displays with data from the provided
        analyzer instance.

        :param analyzer: MolecularInteractionAnalyzer instance with results
        :type analyzer: MolecularInteractionAnalyzer
        :returns: None
        :rtype: None
        """
        self.analyzer = analyzer

        # Update all tabs
        self._update_summary()

        # Update interaction tabs dynamically from config
        for interaction_type, config in INTERACTION_CONFIGS.items():
            self._update_interaction_tab(interaction_type, config)

        # Update special tabs (custom logic)
        self._update_cooperativity_chains()
        self._update_ligand_interactions()

    def _update_summary(self):
        """Update the summary tab."""
        if not self.analyzer:
            return

        self.summary_text.delete(1.0, tk.END)

        # Insert header
        self.summary_text.insert(tk.END, "HBAT Analysis Summary\n", "header")
        self.summary_text.insert(tk.END, "=" * 50 + "\n\n")

        # Get summary
        summary = self.analyzer.get_summary()

        # Timing information
        if "timing" in summary:
            self.summary_text.insert(tk.END, "Analysis Performance:\n", "subheader")
            timing = summary["timing"]
            self.summary_text.insert(
                tk.END,
                f"  Analysis Duration: {timing['analysis_duration_seconds']:.3f} seconds\n\n",
            )

        # PDB fixing information
        if "pdb_fixing" in summary:
            pdb_info = summary["pdb_fixing"]
            self.summary_text.insert(tk.END, "PDB Structure Processing:\n", "subheader")

            if pdb_info.get("applied", False):
                self.summary_text.insert(
                    tk.END, f"  PDB Fixing: Applied using {pdb_info['method']}\n"
                )
                self.summary_text.insert(
                    tk.END, f"  Original Atoms: {pdb_info['original_atoms']}\n"
                )
                self.summary_text.insert(
                    tk.END, f"  Fixed Atoms: {pdb_info['fixed_atoms']}\n"
                )
                if pdb_info.get("added_hydrogens", 0) > 0:
                    self.summary_text.insert(
                        tk.END,
                        f"  Added Hydrogens: {pdb_info['added_hydrogens']} "
                        f"({pdb_info['original_hydrogens']} → {pdb_info['fixed_hydrogens']})\n",
                    )
                self.summary_text.insert(
                    tk.END, f"  Re-detected Bonds: {pdb_info['redetected_bonds']}\n"
                )
            elif "error" in pdb_info:
                self.summary_text.insert(
                    tk.END, f"  PDB Fixing: Failed ({pdb_info['error']})\n"
                )
            else:
                self.summary_text.insert(tk.END, "  PDB Fixing: Not applied\n")
            self.summary_text.insert(tk.END, "\n")

        # Insert summary statistics
        self.summary_text.insert(tk.END, "Interaction Counts:\n", "subheader")
        self.summary_text.insert(
            tk.END, f"  Hydrogen Bonds: {summary['hydrogen_bonds']['count']}\n"
        )
        self.summary_text.insert(
            tk.END, f"  Halogen Bonds: {summary['halogen_bonds']['count']}\n"
        )
        self.summary_text.insert(
            tk.END, f"  π Interactions: {summary['pi_interactions']['count']}\n"
        )

        # Add new interaction types if they exist in summary
        if "pi_pi_stacking" in summary:
            self.summary_text.insert(
                tk.END, f"  π-π Stacking: {summary['pi_pi_stacking']['count']}\n"
            )
        if "carbonyl_interactions" in summary:
            self.summary_text.insert(
                tk.END,
                f"  Carbonyl Interactions: {summary['carbonyl_interactions']['count']}\n",
            )
        if "n_pi_interactions" in summary:
            self.summary_text.insert(
                tk.END,
                f"  n→π* Interactions: {summary['n_pi_interactions']['count']}\n",
            )

        self.summary_text.insert(
            tk.END,
            f"  Cooperativity Chains: {summary['cooperativity_chains']['count']}\n",
        )
        if "water_bridges" in summary:
            self.summary_text.insert(
                tk.END, f"  Water Bridges: {summary['water_bridges']['count']}\n"
            )
        self.summary_text.insert(
            tk.END, f"  Total Interactions: {summary['total_interactions']}\n\n"
        )

        # Bond detection statistics
        if "bond_detection" in summary:
            bond_stats = summary["bond_detection"]
            self.summary_text.insert(tk.END, "Bond Detection:\n", "subheader")
            self.summary_text.insert(
                tk.END, f"  Total Bonds Detected: {bond_stats['total_bonds']}\n"
            )
            if bond_stats["breakdown"]:
                self.summary_text.insert(tk.END, "  Detection Methods:\n")
                for method, stats in bond_stats["breakdown"].items():
                    method_name = method.replace("_", " ").title()
                    self.summary_text.insert(
                        tk.END,
                        f"    {method_name}: {stats['count']} ({stats['percentage']}%)\n",
                    )
            self.summary_text.insert(tk.END, "\n")

        # Detailed interaction statistics
        if summary["hydrogen_bonds"]["count"] > 0:
            self.summary_text.insert(tk.END, "Hydrogen Bond Statistics:\n", "subheader")
            hb_data = summary["hydrogen_bonds"]
            self.summary_text.insert(
                tk.END,
                f"  Average H...A Distance: {hb_data['average_distance']:.2f} Å\n",
            )
            self.summary_text.insert(
                tk.END, f"  Average Angle: {hb_data['average_angle']:.1f}°\n"
            )

            # Bond type distribution
            if "bond_types" in hb_data:
                self.summary_text.insert(tk.END, "  Bond Types:\n")
                for bond_type, count in sorted(hb_data["bond_types"].items()):
                    self.summary_text.insert(tk.END, f"    {bond_type}: {count}\n")
            self.summary_text.insert(tk.END, "\n")

        if summary["halogen_bonds"]["count"] > 0:
            self.summary_text.insert(tk.END, "Halogen Bond Statistics:\n", "subheader")
            xb_data = summary["halogen_bonds"]
            self.summary_text.insert(
                tk.END,
                f"  Average X...A Distance: {xb_data['average_distance']:.2f} Å\n",
            )
            self.summary_text.insert(
                tk.END, f"  Average Angle: {xb_data['average_angle']:.1f}°\n"
            )

            # Bond type distribution
            if "bond_types" in xb_data:
                self.summary_text.insert(tk.END, "  Bond Types:\n")
                for bond_type, count in sorted(xb_data["bond_types"].items()):
                    self.summary_text.insert(tk.END, f"    {bond_type}: {count}\n")
            self.summary_text.insert(tk.END, "\n")

        if summary["pi_interactions"]["count"] > 0:
            self.summary_text.insert(tk.END, "π Interaction Statistics:\n", "subheader")
            pi_data = summary["pi_interactions"]
            self.summary_text.insert(
                tk.END,
                f"  Average H...π Distance: {pi_data['average_distance']:.2f} Å\n",
            )
            self.summary_text.insert(
                tk.END, f"  Average Angle: {pi_data['average_angle']:.1f}°\n\n"
            )

        # Cooperativity chain statistics
        if summary["cooperativity_chains"]["count"] > 0:
            self.summary_text.insert(
                tk.END, "Cooperativity Chain Statistics:\n", "subheader"
            )
            coop_data = summary["cooperativity_chains"]
            self.summary_text.insert(tk.END, f"  Total Chains: {coop_data['count']}\n")

            # Chain types
            if "types" in coop_data and coop_data["types"]:
                type_counts = {}
                for chain_type in coop_data["types"]:
                    type_counts[chain_type] = type_counts.get(chain_type, 0) + 1

                self.summary_text.insert(tk.END, "  Chain Types:\n")
                for chain_type, count in sorted(type_counts.items()):
                    self.summary_text.insert(tk.END, f"    {chain_type}: {count}\n")

            # Chain length distribution
            if "chain_lengths" in coop_data:
                self.summary_text.insert(tk.END, "  Chain Length Distribution:\n")
                for length, count in sorted(coop_data["chain_lengths"].items()):
                    self.summary_text.insert(
                        tk.END, f"    Length {length}: {count} chains\n"
                    )
            self.summary_text.insert(tk.END, "\n")

        # Add some example interactions
        if self.analyzer.hydrogen_bonds:
            self.summary_text.insert(tk.END, "Sample Hydrogen Bonds:\n", "subheader")
            for i, hb in enumerate(self.analyzer.hydrogen_bonds[:5]):
                self.summary_text.insert(tk.END, f"  {i + 1}. {hb}\n")
            if len(self.analyzer.hydrogen_bonds) > 5:
                self.summary_text.insert(
                    tk.END,
                    f"  ... and {len(self.analyzer.hydrogen_bonds) - 5} more\n\n",
                )

        if self.analyzer.halogen_bonds:
            self.summary_text.insert(tk.END, "Sample Halogen Bonds:\n", "subheader")
            for i, xb in enumerate(self.analyzer.halogen_bonds[:3]):
                self.summary_text.insert(tk.END, f"  {i + 1}. {xb}\n")
            if len(self.analyzer.halogen_bonds) > 3:
                self.summary_text.insert(
                    tk.END, f"  ... and {len(self.analyzer.halogen_bonds) - 3} more\n\n"
                )

        if self.analyzer.pi_interactions:
            self.summary_text.insert(tk.END, "Sample π Interactions:\n", "subheader")
            for i, pi in enumerate(self.analyzer.pi_interactions[:3]):
                self.summary_text.insert(tk.END, f"  {i + 1}. {pi}\n")
            if len(self.analyzer.pi_interactions) > 3:
                self.summary_text.insert(
                    tk.END, f"  ... and {len(self.analyzer.pi_interactions) - 3} more\n"
                )

    def _update_interaction_tab(self, interaction_type: str, config):
        """Update interaction tab dynamically from configuration.

        Populates a treeview with data from analyzer based on interaction type
        configuration, eliminating code duplication across specific update methods.

        :param interaction_type: Interaction type key (e.g., 'hydrogen_bonds')
        :param config: InteractionConfig object with column accessors
        """
        if not self.analyzer:
            return

        # Get tree reference
        tree = getattr(self, f"{interaction_type}_tree", None)
        if not tree:
            return

        # Clear existing items
        for item in tree.get_children():
            tree.delete(item)

        # Get interactions from analyzer
        interactions = getattr(self.analyzer, config.analyzer_attr, [])
        if not interactions:
            return

        # Populate tree with data
        for interaction in interactions:
            row_values = []
            for col in config.columns:
                if col.accessor:
                    value = col.accessor(interaction)
                else:
                    value = getattr(interaction, col.name, "")
                row_values.append(value)
            tree.insert("", tk.END, values=tuple(row_values))

    def _update_cooperativity_chains(self):
        """Update the cooperativity chains tab."""
        if not self.analyzer:
            return

        # Clear existing items
        for item in self.coop_tree.get_children():
            self.coop_tree.delete(item)

        # Add cooperativity chains
        for i, chain in enumerate(self.analyzer.cooperativity_chains, 1):
            # Create chain description
            chain_desc = self._format_chain_description(chain)

            self.coop_tree.insert(
                "", tk.END, values=(f"Chain-{i}", chain.chain_length, chain_desc)
            )

    def _format_chain_description(self, chain) -> str:
        """Format a chain description for display."""
        if not chain.interactions:
            return "Empty chain"

        parts = []
        for i, interaction in enumerate(chain.interactions):
            if i == 0:
                # First interaction: show donor
                donor_res = interaction.get_donor_residue()
                donor_atom = interaction.get_donor_atom()
                donor_name = donor_atom.name if donor_atom else "?"
                parts.append(f"{donor_res}({donor_name})")

            # Add interaction symbol and acceptor
            acceptor_res = interaction.get_acceptor_residue()
            if interaction.get_acceptor_atom():
                acceptor_name = interaction.get_acceptor_atom().name
                acceptor_str = f"{acceptor_res}({acceptor_name})"
            else:
                acceptor_str = acceptor_res  # For π interactions

            # Get interaction symbol
            if interaction.interaction_type == "H-Bond":
                symbol = " -> "
            elif interaction.interaction_type == "X-Bond":
                symbol = " =X=> "
            elif interaction.interaction_type == "π–Inter":
                symbol = " ~π~> "
            else:
                symbol = " -> "

            angle_str = f"[{math.degrees(interaction.angle):.1f}°]"
            parts.append(f"{symbol}{acceptor_str} {angle_str}")

        return "".join(parts)

    def _update_ligand_interactions(self):
        """Update the ligand interactions tab with selector and tables."""
        if not self.analyzer:
            return

        # Clear stored data
        self.lig_inter_data = {}
        self.lig_wb_data = {}

        # Check if ligand interactions exist
        if (
            not hasattr(self.analyzer, "ligand_interactions")
            or not self.analyzer.ligand_interactions
        ):
            # No ligands, clear the selector and tables
            self.lig_selector_combo["values"] = []
            self.lig_selector_var.set("")
            for item in self.lig_inter_tree.get_children():
                self.lig_inter_tree.delete(item)
            for item in self.lig_wb_tree.get_children():
                self.lig_wb_tree.delete(item)
            return

        # Get ligand info from ligand_interactions
        ligand_info = self.analyzer.ligand_interactions.ligand_info
        ligand_set = set(ligand_info.keys()) if ligand_info else set()

        # Get all interactions from all types and filter for ligands
        inter_types = {
            "H-Bond": self.analyzer.hydrogen_bonds,
            "Halogen Bond": self.analyzer.halogen_bonds,
            "π-Interaction": self.analyzer.pi_interactions,
            "π-π Stacking": self.analyzer.pi_pi_interactions,
            "Carbonyl": self.analyzer.carbonyl_interactions,
            "n→π*": self.analyzer.n_pi_interactions,
        }

        for inter_type, interactions in inter_types.items():
            for inter in interactions:
                donor_res = inter.get_donor_residue()
                acceptor_res = inter.get_acceptor_residue()

                # Check if this interaction involves any ligand
                for ligand_res in ligand_set:
                    if ligand_res in donor_res or ligand_res in acceptor_res:
                        if ligand_res not in self.lig_inter_data:
                            self.lig_inter_data[ligand_res] = []

                        # Get atom names safely
                        donor_atom_name = "N/A"
                        if hasattr(inter, "get_donor_atom"):
                            donor_atom = inter.get_donor_atom()
                            if donor_atom and hasattr(donor_atom, "name"):
                                donor_atom_name = donor_atom.name
                        elif hasattr(inter, "donor") and inter.donor:
                            donor_atom_name = inter.donor.name

                        acceptor_atom_name = "N/A"
                        if hasattr(inter, "get_acceptor_atom"):
                            acceptor_atom = inter.get_acceptor_atom()
                            if acceptor_atom and hasattr(acceptor_atom, "name"):
                                acceptor_atom_name = acceptor_atom.name
                        elif hasattr(inter, "acceptor") and inter.acceptor:
                            acceptor_atom_name = inter.acceptor.name

                        # Format distance based on interaction type
                        distance = ""
                        try:
                            if hasattr(inter, "distance"):
                                distance = f"{inter.distance:.2f}"
                            elif hasattr(inter, "get_donor_interaction_distance"):
                                distance = (
                                    f"{inter.get_donor_interaction_distance():.2f}"
                                )
                        except:
                            pass

                        # Get properties
                        properties = ""
                        if hasattr(inter, "donor_acceptor_properties"):
                            properties = inter.donor_acceptor_properties

                        inter_data = (
                            inter_type,
                            donor_res,
                            donor_atom_name,
                            acceptor_res,
                            acceptor_atom_name,
                            distance,
                            properties,
                        )
                        if inter_data not in self.lig_inter_data[ligand_res]:
                            self.lig_inter_data[ligand_res].append(inter_data)

        # Get water bridges
        if hasattr(self.analyzer, "water_bridges"):
            for wb in self.analyzer.water_bridges:
                donor_res = wb.get_donor_residue()
                acceptor_res = wb.get_acceptor_residue()

                for ligand_res in ligand_set:
                    if ligand_res in donor_res or ligand_res in acceptor_res:
                        if ligand_res not in self.lig_wb_data:
                            self.lig_wb_data[ligand_res] = []

                        water_res = "; ".join(wb.water_residues)
                        distance = f"{wb.get_donor_acceptor_distance():.2f}"

                        wb_data = (
                            donor_res,
                            acceptor_res,
                            wb.bridge_length,
                            water_res,
                            distance,
                        )
                        if wb_data not in self.lig_wb_data[ligand_res]:
                            self.lig_wb_data[ligand_res].append(wb_data)

        # Update selector dropdown
        ligand_list = sorted(ligand_set)
        self.lig_selector_combo["values"] = ligand_list
        if ligand_list:
            self.lig_selector_combo.current(0)
            self._update_ligand_tables(ligand_list[0])
        else:
            self.lig_selector_var.set("")
            for item in self.lig_inter_tree.get_children():
                self.lig_inter_tree.delete(item)
            for item in self.lig_wb_tree.get_children():
                self.lig_wb_tree.delete(item)

    def _filter_results(self, tree, search_term):
        """Filter tree results based on search term."""
        if not search_term:
            return

        # Hide items that don't match the search term
        for item in tree.get_children():
            values = tree.item(item)["values"]
            match = any(search_term.lower() in str(value).lower() for value in values)
            if not match:
                tree.detach(item)

    def _clear_filter(self, tree, search_var):
        """Clear filter and show all results."""
        search_var.set("")
        # Refresh the tree by updating the corresponding data
        if self.analyzer:
            # Map trees to interaction types
            tree_to_type = {
                getattr(self, "hb_tree", None): "hydrogen_bonds",
                getattr(self, "wb_tree", None): "water_bridges",
                getattr(self, "xb_tree", None): "halogen_bonds",
                getattr(self, "pi_tree", None): "pi_interactions",
                getattr(self, "pi_pi_tree", None): "pi_pi_interactions",
                getattr(self, "carbonyl_tree", None): "carbonyl_interactions",
                getattr(self, "n_pi_tree", None): "n_pi_interactions",
                getattr(self, "coop_tree", None): "cooperativity_chains",
            }

            # Find the interaction type for this tree
            interaction_type = tree_to_type.get(tree)
            if interaction_type and interaction_type in INTERACTION_CONFIGS:
                config = get_interaction_config(interaction_type)
                self._update_interaction_tab(interaction_type, config)
            elif tree == self.coop_tree:
                self._update_cooperativity_chains()
            elif tree == getattr(self, "lig_inter_tree", None) or tree == getattr(
                self, "lig_wb_tree", None
            ):
                self._update_ligand_interactions()

    def clear_results(self) -> None:
        """Clear all results from the panel.

        Removes all displayed results and resets the panel to
        its initial empty state.

        :returns: None
        :rtype: None
        """
        self.analyzer = None

        # Clear text widgets
        self.summary_text.delete(1.0, tk.END)

        # Clear all treeviews dynamically
        trees_to_clear = [
            "hb_tree",
            "xb_tree",
            "pi_tree",
            "pi_pi_tree",
            "carbonyl_tree",
            "n_pi_tree",
            "wb_tree",
            "coop_tree",
            "lig_inter_tree",
            "lig_wb_tree",
            # Also handle new-style names
            "hydrogen_bonds_tree",
            "water_bridges_tree",
            "halogen_bonds_tree",
            "pi_interactions_tree",
            "pi_pi_interactions_tree",
            "carbonyl_interactions_tree",
            "n_pi_interactions_tree",
            "cooperativity_chains_tree",
        ]

        for tree_name in trees_to_clear:
            tree = getattr(self, tree_name, None)
            if tree:
                for item in tree.get_children():
                    tree.delete(item)

        # Clear all search filters dynamically
        search_vars_to_clear = [
            "hb_search_var",
            "xb_search_var",
            "pi_search_var",
            "pi_pi_search_var",
            "carbonyl_search_var",
            "n_pi_search_var",
            "coop_search_var",
            "wb_search_var",
            "lig_selector_var",
            # Also handle new-style names
            "hydrogen_bonds_search_var",
            "water_bridges_search_var",
            "halogen_bonds_search_var",
            "pi_interactions_search_var",
            "pi_pi_interactions_search_var",
            "carbonyl_interactions_search_var",
            "n_pi_interactions_search_var",
            "cooperativity_chains_search_var",
        ]

        for var_name in search_vars_to_clear:
            var = getattr(self, var_name, None)
            if var:
                var.set("")

        # Add placeholder text
        self.summary_text.insert(tk.END, "No analysis results available.\n\n")
        self.summary_text.insert(
            tk.END,
            "Please load a structure file (PDB or CIF) and run analysis to see results.",
        )

    def _visualize_selected_chain(self):
        """Visualize the selected cooperativity chain in a new window."""
        if not VISUALIZATION_AVAILABLE:
            messagebox.showerror(
                "Error",
                "Visualization libraries (networkx, matplotlib) are not available.",
            )
            return

        selection = self.coop_tree.selection()
        if not selection:
            messagebox.showwarning(
                "Warning", "Please select a cooperativity chain to visualize."
            )
            return

        item = selection[0]
        values = self.coop_tree.item(item)["values"]
        chain_id = values[0]  # Chain-1, Chain-2, etc.

        # Get the chain index from the ID
        try:
            chain_index = int(chain_id.split("-")[1]) - 1
            if chain_index < 0 or chain_index >= len(
                self.analyzer.cooperativity_chains
            ):
                raise IndexError
            chain = self.analyzer.cooperativity_chains[chain_index]
        except (ValueError, IndexError):
            messagebox.showerror("Error", "Invalid chain selection.")
            return

        # Create the visualization window using the new module
        ChainVisualizationWindow(self.parent, chain, chain_id)

    def _on_chain_double_click(self, event):
        """Handle double-click on cooperativity chain to open visualization."""
        # Get the item that was double-clicked
        item = self.coop_tree.identify_row(event.y)
        if item:
            # Select the item first
            self.coop_tree.selection_set(item)
            # Then visualize it
            self._visualize_selected_chain()
