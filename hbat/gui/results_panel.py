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
    from .chain_visualization import ChainVisualizationWindow

    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False

from ..core.analysis import HBondAnalyzer


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
        self.analyzer: Optional[HBondAnalyzer] = None
        self._create_widgets()

    def _create_widgets(self):
        """Create result display widgets."""
        # Create main notebook for different result types
        self.notebook = ttk.Notebook(self.parent)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # Summary tab
        self._create_summary_tab()

        # Hydrogen bonds tab
        self._create_hydrogen_bonds_tab()

        # Halogen bonds tab
        self._create_halogen_bonds_tab()

        # Pi interactions tab
        self._create_pi_interactions_tab()

        # Cooperativity chains tab
        self._create_cooperativity_chains_tab()

        # Statistics tab
        self._create_statistics_tab()

    def _create_summary_tab(self):
        """Create summary results tab."""
        summary_frame = ttk.Frame(self.notebook)
        self.notebook.add(summary_frame, text="Summary")

        # Create text widget with scrollbar
        text_frame = ttk.Frame(summary_frame)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.summary_text = tk.Text(text_frame, wrap=tk.WORD, font=("Courier", 12))
        summary_scrollbar = ttk.Scrollbar(
            text_frame, orient=tk.VERTICAL, command=self.summary_text.yview
        )
        self.summary_text.configure(yscrollcommand=summary_scrollbar.set)

        self.summary_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        summary_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Configure text tags for formatting
        self.summary_text.tag_configure(
            "header", font=("Courier", 12, "bold"), foreground="blue"
        )
        self.summary_text.tag_configure("subheader", font=("Courier", 12, "bold"))
        self.summary_text.tag_configure("highlight", background="cyan")

    def _create_hydrogen_bonds_tab(self):
        """Create hydrogen bonds results tab."""
        hb_frame = ttk.Frame(self.notebook)
        self.notebook.add(hb_frame, text="Hydrogen Bonds")

        # Create treeview for hydrogen bonds
        tree_frame = ttk.Frame(hb_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        columns = (
            "donor_res",
            "donor_atom",
            "acceptor_res",
            "acceptor_atom",
            "distance",
            "angle",
            "da_distance",
            "type",
        )

        self.hb_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=15
        )

        # Configure columns
        self.hb_tree.heading("donor_res", text="Donor Residue")
        self.hb_tree.heading("donor_atom", text="Donor Atom")
        self.hb_tree.heading("acceptor_res", text="Acceptor Residue")
        self.hb_tree.heading("acceptor_atom", text="Acceptor Atom")
        self.hb_tree.heading("distance", text="H...A (Å)")
        self.hb_tree.heading("angle", text="Angle (°)")
        self.hb_tree.heading("da_distance", text="D...A (Å)")
        self.hb_tree.heading("type", text="Type")

        # Configure column widths
        self.hb_tree.column("donor_res", width=100)
        self.hb_tree.column("donor_atom", width=80)
        self.hb_tree.column("acceptor_res", width=100)
        self.hb_tree.column("acceptor_atom", width=80)
        self.hb_tree.column("distance", width=80)
        self.hb_tree.column("angle", width=80)
        self.hb_tree.column("da_distance", width=80)
        self.hb_tree.column("type", width=100)

        # Add scrollbars
        hb_v_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.VERTICAL, command=self.hb_tree.yview
        )
        hb_h_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.HORIZONTAL, command=self.hb_tree.xview
        )
        self.hb_tree.configure(
            yscrollcommand=hb_v_scrollbar.set, xscrollcommand=hb_h_scrollbar.set
        )

        self.hb_tree.grid(row=0, column=0, sticky="nsew")
        hb_v_scrollbar.grid(row=0, column=1, sticky="ns")
        hb_h_scrollbar.grid(row=1, column=0, sticky="ew")

        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        # Add search functionality
        search_frame = ttk.Frame(hb_frame)
        search_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.hb_search_var = tk.StringVar()
        search_entry = ttk.Entry(
            search_frame, textvariable=self.hb_search_var, width=30
        )
        search_entry.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            search_frame,
            text="Filter",
            command=lambda: self._filter_results(
                self.hb_tree, self.hb_search_var.get()
            ),
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            search_frame,
            text="Clear",
            command=lambda: self._clear_filter(self.hb_tree, self.hb_search_var),
        ).pack(side=tk.LEFT, padx=5)

    def _create_halogen_bonds_tab(self):
        """Create halogen bonds results tab."""
        xb_frame = ttk.Frame(self.notebook)
        self.notebook.add(xb_frame, text="Halogen Bonds")

        # Create treeview for halogen bonds
        tree_frame = ttk.Frame(xb_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        columns = (
            "halogen_res",
            "halogen_atom",
            "acceptor_res",
            "acceptor_atom",
            "distance",
            "angle",
            "type",
        )

        self.xb_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=15
        )

        # Configure columns
        self.xb_tree.heading("halogen_res", text="Halogen Residue")
        self.xb_tree.heading("halogen_atom", text="Halogen Atom")
        self.xb_tree.heading("acceptor_res", text="Acceptor Residue")
        self.xb_tree.heading("acceptor_atom", text="Acceptor Atom")
        self.xb_tree.heading("distance", text="X...A (Å)")
        self.xb_tree.heading("angle", text="Angle (°)")
        self.xb_tree.heading("type", text="Type")

        # Configure column widths
        self.xb_tree.column("halogen_res", width=120)
        self.xb_tree.column("halogen_atom", width=100)
        self.xb_tree.column("acceptor_res", width=120)
        self.xb_tree.column("acceptor_atom", width=100)
        self.xb_tree.column("distance", width=80)
        self.xb_tree.column("angle", width=80)
        self.xb_tree.column("type", width=100)

        # Add scrollbars
        xb_v_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.VERTICAL, command=self.xb_tree.yview
        )
        xb_h_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.HORIZONTAL, command=self.xb_tree.xview
        )
        self.xb_tree.configure(
            yscrollcommand=xb_v_scrollbar.set, xscrollcommand=xb_h_scrollbar.set
        )

        self.xb_tree.grid(row=0, column=0, sticky="nsew")
        xb_v_scrollbar.grid(row=0, column=1, sticky="ns")
        xb_h_scrollbar.grid(row=1, column=0, sticky="ew")

        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        # Add search functionality
        search_frame = ttk.Frame(xb_frame)
        search_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.xb_search_var = tk.StringVar()
        search_entry = ttk.Entry(
            search_frame, textvariable=self.xb_search_var, width=30
        )
        search_entry.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            search_frame,
            text="Filter",
            command=lambda: self._filter_results(
                self.xb_tree, self.xb_search_var.get()
            ),
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            search_frame,
            text="Clear",
            command=lambda: self._clear_filter(self.xb_tree, self.xb_search_var),
        ).pack(side=tk.LEFT, padx=5)

    def _create_pi_interactions_tab(self):
        """Create π interactions results tab."""
        pi_frame = ttk.Frame(self.notebook)
        self.notebook.add(pi_frame, text="π Interactions")

        # Create treeview for π interactions
        tree_frame = ttk.Frame(pi_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        columns = ("donor_res", "donor_atom", "pi_res", "distance", "angle")

        self.pi_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=15
        )

        # Configure columns
        self.pi_tree.heading("donor_res", text="Donor Residue")
        self.pi_tree.heading("donor_atom", text="Donor Atom")
        self.pi_tree.heading("pi_res", text="π Residue")
        self.pi_tree.heading("distance", text="H...π (Å)")
        self.pi_tree.heading("angle", text="Angle (°)")

        # Configure column widths
        self.pi_tree.column("donor_res", width=120)
        self.pi_tree.column("donor_atom", width=100)
        self.pi_tree.column("pi_res", width=120)
        self.pi_tree.column("distance", width=100)
        self.pi_tree.column("angle", width=100)

        # Add scrollbars
        pi_v_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.VERTICAL, command=self.pi_tree.yview
        )
        pi_h_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.HORIZONTAL, command=self.pi_tree.xview
        )
        self.pi_tree.configure(
            yscrollcommand=pi_v_scrollbar.set, xscrollcommand=pi_h_scrollbar.set
        )

        self.pi_tree.grid(row=0, column=0, sticky="nsew")
        pi_v_scrollbar.grid(row=0, column=1, sticky="ns")
        pi_h_scrollbar.grid(row=1, column=0, sticky="ew")

        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        # Add search functionality
        search_frame = ttk.Frame(pi_frame)
        search_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.pi_search_var = tk.StringVar()
        search_entry = ttk.Entry(
            search_frame, textvariable=self.pi_search_var, width=30
        )
        search_entry.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            search_frame,
            text="Filter",
            command=lambda: self._filter_results(
                self.pi_tree, self.pi_search_var.get()
            ),
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            search_frame,
            text="Clear",
            command=lambda: self._clear_filter(self.pi_tree, self.pi_search_var),
        ).pack(side=tk.LEFT, padx=5)

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
        self.coop_tree.column("chain_id", width=80)
        self.coop_tree.column("chain_length", width=80)
        self.coop_tree.column("chain_description", width=800)

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

    def _create_statistics_tab(self):
        """Create statistics tab."""
        stats_frame = ttk.Frame(self.notebook)
        self.notebook.add(stats_frame, text="Statistics")

        # Create text widget for statistics
        text_frame = ttk.Frame(stats_frame)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.stats_text = tk.Text(text_frame, wrap=tk.WORD, font=("Courier", 12))
        stats_scrollbar = ttk.Scrollbar(
            text_frame, orient=tk.VERTICAL, command=self.stats_text.yview
        )
        self.stats_text.configure(yscrollcommand=stats_scrollbar.set)

        self.stats_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        stats_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Configure text tags
        self.stats_text.tag_configure(
            "header", font=("Courier", 12, "bold"), foreground="blue"
        )
        self.stats_text.tag_configure("subheader", font=("Courier", 12, "bold"))
        self.stats_text.tag_configure(
            "number", foreground="red", font=("Courier", 12, "bold")
        )

    def update_results(self, analyzer: HBondAnalyzer) -> None:
        """Update the results panel with new analysis results.

        Refreshes all result displays with data from the provided
        analyzer instance.

        :param analyzer: HBondAnalyzer instance with results
        :type analyzer: HBondAnalyzer
        :returns: None
        :rtype: None
        """
        self.analyzer = analyzer

        # Update all tabs
        self._update_summary()
        self._update_hydrogen_bonds()
        self._update_halogen_bonds()
        self._update_pi_interactions()
        self._update_cooperativity_chains()
        self._update_statistics()

    def _update_summary(self):
        """Update the summary tab."""
        if not self.analyzer:
            return

        self.summary_text.delete(1.0, tk.END)

        # Insert header
        self.summary_text.insert(tk.END, "HBAT Analysis Summary\n", "header")
        self.summary_text.insert(tk.END, "=" * 50 + "\n\n")

        # Get statistics
        stats = self.analyzer.get_statistics()

        # Insert summary statistics
        self.summary_text.insert(tk.END, "Interaction Counts:\n", "subheader")
        self.summary_text.insert(
            tk.END, f"  Hydrogen Bonds: {stats['hydrogen_bonds']}\n"
        )
        self.summary_text.insert(tk.END, f"  Halogen Bonds: {stats['halogen_bonds']}\n")
        self.summary_text.insert(
            tk.END, f"  π Interactions: {stats['pi_interactions']}\n"
        )
        self.summary_text.insert(
            tk.END, f"  Cooperativity Chains: {stats.get('cooperativity_chains', 0)}\n"
        )
        self.summary_text.insert(
            tk.END, f"  Total Interactions: {stats['total_interactions']}\n\n"
        )

        # Add hydrogen bond details if available
        if stats.get("hb_avg_distance"):
            self.summary_text.insert(tk.END, "Hydrogen Bond Statistics:\n", "subheader")
            self.summary_text.insert(
                tk.END, f"  Average H...A Distance: {stats['hb_avg_distance']:.2f} Å\n"
            )
            self.summary_text.insert(
                tk.END, f"  Average Angle: {stats['hb_avg_angle']:.1f}°\n"
            )
            self.summary_text.insert(
                tk.END,
                f"  Distance Range: {stats['hb_min_distance']:.2f} - {stats['hb_max_distance']:.2f} Å\n\n",
            )

        # Add some example interactions
        if self.analyzer.hydrogen_bonds:
            self.summary_text.insert(tk.END, "Sample Hydrogen Bonds:\n", "subheader")
            for i, hb in enumerate(self.analyzer.hydrogen_bonds[:5]):
                self.summary_text.insert(tk.END, f"  {i+1}. {hb}\n")
            if len(self.analyzer.hydrogen_bonds) > 5:
                self.summary_text.insert(
                    tk.END, f"  ... and {len(self.analyzer.hydrogen_bonds) - 5} more\n"
                )

    def _update_hydrogen_bonds(self):
        """Update the hydrogen bonds tab."""
        if not self.analyzer:
            return

        # Clear existing items
        for item in self.hb_tree.get_children():
            self.hb_tree.delete(item)

        # Add hydrogen bonds
        for hb in self.analyzer.hydrogen_bonds:
            self.hb_tree.insert(
                "",
                tk.END,
                values=(
                    hb.donor_residue,
                    hb.donor.name,
                    hb.acceptor_residue,
                    hb.acceptor.name,
                    f"{hb.distance:.2f}",
                    f"{math.degrees(hb.angle):.1f}",
                    f"{hb.donor_acceptor_distance:.2f}",
                    hb.bond_type,
                ),
            )

    def _update_halogen_bonds(self):
        """Update the halogen bonds tab."""
        if not self.analyzer:
            return

        # Clear existing items
        for item in self.xb_tree.get_children():
            self.xb_tree.delete(item)

        # Add halogen bonds
        for xb in self.analyzer.halogen_bonds:
            self.xb_tree.insert(
                "",
                tk.END,
                values=(
                    xb.halogen_residue,
                    xb.halogen.name,
                    xb.acceptor_residue,
                    xb.acceptor.name,
                    f"{xb.distance:.2f}",
                    f"{math.degrees(xb.angle):.1f}",
                    xb.bond_type,
                ),
            )

    def _update_pi_interactions(self):
        """Update the π interactions tab."""
        if not self.analyzer:
            return

        # Clear existing items
        for item in self.pi_tree.get_children():
            self.pi_tree.delete(item)

        # Add π interactions
        for pi in self.analyzer.pi_interactions:
            self.pi_tree.insert(
                "",
                tk.END,
                values=(
                    pi.donor_residue,
                    pi.donor.name,
                    pi.pi_residue,
                    f"{pi.distance:.2f}",
                    f"{math.degrees(pi.angle):.1f}",
                ),
            )

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
            if interaction.interaction_type == "hydrogen_bond":
                symbol = " -> "
            elif interaction.interaction_type == "halogen_bond":
                symbol = " =X=> "
            elif interaction.interaction_type == "pi_interaction":
                symbol = " ~π~> "
            else:
                symbol = " -> "

            angle_str = f"[{math.degrees(interaction.angle):.1f}°]"
            parts.append(f"{symbol}{acceptor_str} {angle_str}")

        return "".join(parts)

    def _update_statistics(self):
        """Update the statistics tab."""
        if not self.analyzer:
            return

        self.stats_text.delete(1.0, tk.END)

        # Insert header
        self.stats_text.insert(tk.END, "Detailed Statistics\n", "header")
        self.stats_text.insert(tk.END, "=" * 50 + "\n\n")

        stats = self.analyzer.get_statistics()

        # Interaction counts
        self.stats_text.insert(tk.END, "Interaction Counts:\n", "subheader")
        self.stats_text.insert(tk.END, f"  Hydrogen Bonds: ")
        self.stats_text.insert(tk.END, f"{stats['hydrogen_bonds']}\n", "number")
        self.stats_text.insert(tk.END, f"  Halogen Bonds: ")
        self.stats_text.insert(tk.END, f"{stats['halogen_bonds']}\n", "number")
        self.stats_text.insert(tk.END, f"  π Interactions: ")
        self.stats_text.insert(tk.END, f"{stats['pi_interactions']}\n", "number")
        self.stats_text.insert(tk.END, f"  Cooperativity Chains: ")
        self.stats_text.insert(
            tk.END, f"{stats.get('cooperativity_chains', 0)}\n", "number"
        )
        self.stats_text.insert(tk.END, f"  Total: ")
        self.stats_text.insert(tk.END, f"{stats['total_interactions']}\n\n", "number")

        # Hydrogen bond type distribution
        if self.analyzer.hydrogen_bonds:
            hb_types = {}
            for hb in self.analyzer.hydrogen_bonds:
                hb_types[hb.bond_type] = hb_types.get(hb.bond_type, 0) + 1

            self.stats_text.insert(tk.END, "Hydrogen Bond Types:\n", "subheader")
            for bond_type, count in sorted(hb_types.items()):
                self.stats_text.insert(tk.END, f"  {bond_type}: ")
                self.stats_text.insert(tk.END, f"{count}\n", "number")
            self.stats_text.insert(tk.END, "\n")

        # Distance and angle distributions
        if self.analyzer.hydrogen_bonds:
            distances = [hb.distance for hb in self.analyzer.hydrogen_bonds]
            angles = [math.degrees(hb.angle) for hb in self.analyzer.hydrogen_bonds]

            self.stats_text.insert(tk.END, "Hydrogen Bond Geometry:\n", "subheader")
            self.stats_text.insert(tk.END, f"  Distance Statistics (Å):\n")
            self.stats_text.insert(
                tk.END, f"    Mean: {sum(distances)/len(distances):.2f}\n"
            )
            self.stats_text.insert(tk.END, f"    Min: {min(distances):.2f}\n")
            self.stats_text.insert(tk.END, f"    Max: {max(distances):.2f}\n")

            self.stats_text.insert(tk.END, f"  Angle Statistics (°):\n")
            self.stats_text.insert(tk.END, f"    Mean: {sum(angles)/len(angles):.1f}\n")
            self.stats_text.insert(tk.END, f"    Min: {min(angles):.1f}\n")
            self.stats_text.insert(tk.END, f"    Max: {max(angles):.1f}\n")

        # Cooperativity statistics
        if self.analyzer.cooperativity_chains:
            self.stats_text.insert(
                tk.END, f"\nCooperativity Statistics:\n", "subheader"
            )
            self.stats_text.insert(tk.END, f"  Total Chains: ")
            self.stats_text.insert(
                tk.END, f"{len(self.analyzer.cooperativity_chains)}\n", "number"
            )

            chain_lengths = [
                chain.chain_length for chain in self.analyzer.cooperativity_chains
            ]
            self.stats_text.insert(
                tk.END,
                f"  Average Chain Length: {sum(chain_lengths)/len(chain_lengths):.1f}\n",
            )
            self.stats_text.insert(tk.END, f"  Longest Chain: {max(chain_lengths)}\n")

            # Count chain types
            chain_types = {}
            for chain in self.analyzer.cooperativity_chains:
                chain_types[chain.chain_type] = chain_types.get(chain.chain_type, 0) + 1

            self.stats_text.insert(tk.END, f"  Chain Types:\n")
            for chain_type, count in sorted(chain_types.items()):
                self.stats_text.insert(tk.END, f"    {chain_type}: ")
                self.stats_text.insert(tk.END, f"{count}\n", "number")

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
            if tree == self.hb_tree:
                self._update_hydrogen_bonds()
            elif tree == self.xb_tree:
                self._update_halogen_bonds()
            elif tree == self.pi_tree:
                self._update_pi_interactions()
            elif tree == self.coop_tree:
                self._update_cooperativity_chains()

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
        self.stats_text.delete(1.0, tk.END)

        # Clear treeviews
        for item in self.hb_tree.get_children():
            self.hb_tree.delete(item)

        for item in self.xb_tree.get_children():
            self.xb_tree.delete(item)

        for item in self.pi_tree.get_children():
            self.pi_tree.delete(item)

        for item in self.coop_tree.get_children():
            self.coop_tree.delete(item)

        # Clear all search filters
        self.hb_search_var.set("")
        self.xb_search_var.set("")
        self.pi_search_var.set("")
        self.coop_search_var.set("")

        # Add placeholder text
        self.summary_text.insert(tk.END, "No analysis results available.\n\n")
        self.summary_text.insert(
            tk.END, "Please load a PDB file and run analysis to see results."
        )

        self.stats_text.insert(tk.END, "No statistics available.\n\n")
        self.stats_text.insert(
            tk.END, "Please run analysis to see detailed statistics."
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
