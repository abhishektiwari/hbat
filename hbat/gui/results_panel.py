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

from ..core.analysis import MolecularInteractionAnalyzer


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

    def _parse_residue_string(self, residue_str: str) -> str:
        """Parse a residue string like 'A123ALA' into 'A:ALA123'.

        :param residue_str: Residue string in format ChainResSeqResName
        :type residue_str: str
        :returns: Formatted string like 'Chain:ResNameResSeq'
        :rtype: str
        """
        if not residue_str:
            return "Unknown"

        # Handle the format: ChainIdResSeqResName (e.g., "A123ALA")
        # Chain ID is typically 1 character, residue name is typically 3 characters at the end
        if len(residue_str) >= 4:
            chain_id = residue_str[0]
            # The residue name is typically the last 3 characters
            res_name = residue_str[-3:]
            # The residue sequence number is everything in between
            res_seq = residue_str[1:-3]
            return f"{chain_id}:{res_name}{res_seq}"
        else:
            # Fallback for unexpected formats
            return residue_str

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

        # New interaction types tabs
        self._create_pi_pi_stacking_tab()

        self._create_carbonyl_interactions_tab()

        self._create_n_pi_interactions_tab()

        # Cooperativity chains tab
        self._create_cooperativity_chains_tab()

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
            "hydrogen",
            "acceptor_res",
            "acceptor_atom",
            "distance",
            "angle",
            "da_distance",
            "type",
            "da_props",
            "bs_int",
        )

        self.hb_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=15
        )

        # Configure columns
        self.hb_tree.heading("donor_res", text="Donor Residue")
        self.hb_tree.heading("donor_atom", text="Donor Atom")
        self.hb_tree.heading("hydrogen", text="Hydrogen Atom")
        self.hb_tree.heading("acceptor_res", text="Acceptor Residue")
        self.hb_tree.heading("acceptor_atom", text="Acceptor Atom")
        self.hb_tree.heading("distance", text="H...A (Å)")
        self.hb_tree.heading("angle", text="Angle (°)")
        self.hb_tree.heading("da_distance", text="D...A (Å)")
        self.hb_tree.heading("type", text="Type")
        self.hb_tree.heading("da_props", text="D-A Props")
        self.hb_tree.heading("bs_int", text="B/S")

        # Configure column widths
        self.hb_tree.column("donor_res", width=120)
        self.hb_tree.column("hydrogen", width=120)
        self.hb_tree.column("donor_atom", width=100)
        self.hb_tree.column("acceptor_res", width=120)
        self.hb_tree.column("acceptor_atom", width=100)
        self.hb_tree.column("distance", width=90)
        self.hb_tree.column("angle", width=90)
        self.hb_tree.column("da_distance", width=90)
        self.hb_tree.column("type", width=120)
        self.hb_tree.column("da_props", width=100)
        self.hb_tree.column("bs_int", width=70)

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
            "donor_atom",
            "halogen_atom",
            "acceptor_res",
            "acceptor_atom",
            "distance",
            "angle",
            "type",
            "bs_interaction",
            "da_properties",
        )

        self.xb_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=15
        )

        # Configure columns
        self.xb_tree.heading("halogen_res", text="Halogen Residue")
        self.xb_tree.heading("donor_atom", text="Donor Atom")
        self.xb_tree.heading("halogen_atom", text="Halogen Atom")
        self.xb_tree.heading("acceptor_res", text="Acceptor Residue")
        self.xb_tree.heading("acceptor_atom", text="Acceptor Atom")
        self.xb_tree.heading("distance", text="X...A (Å)")
        self.xb_tree.heading("angle", text="Angle (°)")
        self.xb_tree.heading("type", text="Type")
        self.xb_tree.heading("bs_interaction", text="B/S Interaction")
        self.xb_tree.heading("da_properties", text="D-A Properties")

        # Configure column widths
        self.xb_tree.column("halogen_res", width=140)
        self.xb_tree.column("donor_atom", width=100)
        self.xb_tree.column("halogen_atom", width=120)
        self.xb_tree.column("acceptor_res", width=140)
        self.xb_tree.column("acceptor_atom", width=120)
        self.xb_tree.column("distance", width=90)
        self.xb_tree.column("angle", width=90)
        self.xb_tree.column("type", width=120)
        self.xb_tree.column("bs_interaction", width=100)
        self.xb_tree.column("da_properties", width=120)

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

        columns = (
            "donor_res",
            "donor_atom",
            "pi_res",
            "distance",
            "angle",
            "type",
            "da_props",
            "bs_int",
        )

        self.pi_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=15
        )

        # Configure columns
        self.pi_tree.heading("donor_res", text="Donor Residue")
        self.pi_tree.heading("donor_atom", text="Donor Atom")
        self.pi_tree.heading("pi_res", text="π Residue")
        self.pi_tree.heading("distance", text="H...π (Å)")
        self.pi_tree.heading("angle", text="Angle (°)")
        self.pi_tree.heading("type", text="Type")
        self.pi_tree.heading("da_props", text="D-A Props")
        self.pi_tree.heading("bs_int", text="B/S")

        # Configure column widths
        self.pi_tree.column("donor_res", width=140)
        self.pi_tree.column("donor_atom", width=120)
        self.pi_tree.column("pi_res", width=140)
        self.pi_tree.column("distance", width=110)
        self.pi_tree.column("angle", width=110)
        self.pi_tree.column("type", width=90)
        self.pi_tree.column("da_props", width=100)
        self.pi_tree.column("bs_int", width=70)

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

    def _create_pi_pi_stacking_tab(self):
        """Create π-π stacking interactions results tab."""
        pi_pi_frame = ttk.Frame(self.notebook)
        self.notebook.add(pi_pi_frame, text="π-π Stacking")

        # Create treeview for π-π stacking
        tree_frame = ttk.Frame(pi_pi_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        columns = (
            "ring1_res",
            "ring1_atoms",
            "ring2_res",
            "ring2_atoms",
            "distance",
            "plane_angle",
            "offset",
            "stacking_type",
            "bs_int",
        )

        self.pi_pi_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=15
        )

        # Configure columns
        self.pi_pi_tree.heading("ring1_res", text="Ring 1 Residue")
        self.pi_pi_tree.heading("ring1_atoms", text="Ring 1 Atoms")
        self.pi_pi_tree.heading("ring2_res", text="Ring 2 Residue")
        self.pi_pi_tree.heading("ring2_atoms", text="Ring 2 Atoms")
        self.pi_pi_tree.heading("distance", text="Distance (Å)")
        self.pi_pi_tree.heading("plane_angle", text="Plane Angle (°)")
        self.pi_pi_tree.heading("offset", text="Offset (Å)")
        self.pi_pi_tree.heading("stacking_type", text="Stacking Type")
        self.pi_pi_tree.heading("bs_int", text="B/S")

        # Configure column widths
        self.pi_pi_tree.column("ring1_res", width=120)
        self.pi_pi_tree.column("ring1_atoms", width=140)
        self.pi_pi_tree.column("ring2_res", width=120)
        self.pi_pi_tree.column("ring2_atoms", width=140)
        self.pi_pi_tree.column("distance", width=100)
        self.pi_pi_tree.column("plane_angle", width=110)
        self.pi_pi_tree.column("offset", width=100)
        self.pi_pi_tree.column("stacking_type", width=120)
        self.pi_pi_tree.column("bs_int", width=70)

        # Add scrollbars
        pi_pi_v_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.VERTICAL, command=self.pi_pi_tree.yview
        )
        pi_pi_h_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.HORIZONTAL, command=self.pi_pi_tree.xview
        )
        self.pi_pi_tree.configure(
            yscrollcommand=pi_pi_v_scrollbar.set, xscrollcommand=pi_pi_h_scrollbar.set
        )

        self.pi_pi_tree.grid(row=0, column=0, sticky="nsew")
        pi_pi_v_scrollbar.grid(row=0, column=1, sticky="ns")
        pi_pi_h_scrollbar.grid(row=1, column=0, sticky="ew")

        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        # Add info label
        info_frame = ttk.Frame(pi_pi_frame)
        info_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(
            info_frame,
            text="π-π Stacking Interactions: Parallel, T-shaped, and offset aromatic ring interactions",
        ).pack(side=tk.LEFT)

        # Add search functionality
        search_frame = ttk.Frame(pi_pi_frame)
        search_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.pi_pi_search_var = tk.StringVar()
        search_entry = ttk.Entry(
            search_frame, textvariable=self.pi_pi_search_var, width=30
        )
        search_entry.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            search_frame,
            text="Filter",
            command=lambda: self._filter_results(
                self.pi_pi_tree, self.pi_pi_search_var.get()
            ),
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            search_frame,
            text="Clear",
            command=lambda: self._clear_filter(self.pi_pi_tree, self.pi_pi_search_var),
        ).pack(side=tk.LEFT, padx=5)

    def _create_carbonyl_interactions_tab(self):
        """Create carbonyl n→π* interactions results tab."""
        carbonyl_frame = ttk.Frame(self.notebook)
        self.notebook.add(carbonyl_frame, text="Carbonyl Interactions")

        # Create treeview for carbonyl interactions
        tree_frame = ttk.Frame(carbonyl_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        columns = (
            "acceptor_res",
            "acceptor_atom",
            "carbonyl_res",
            "carbonyl_atoms",
            "distance",
            "angle",
            "carbonyl_type",
            "bs_int",
        )

        self.carbonyl_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=15
        )

        # Configure columns
        self.carbonyl_tree.heading("acceptor_res", text="Acceptor Residue")
        self.carbonyl_tree.heading("acceptor_atom", text="Acceptor Atom")
        self.carbonyl_tree.heading("carbonyl_res", text="Carbonyl Residue")
        self.carbonyl_tree.heading("carbonyl_atoms", text="Carbonyl C=O")
        self.carbonyl_tree.heading("distance", text="O···C Distance (Å)")
        self.carbonyl_tree.heading("angle", text="Bürgi-Dunitz Angle (°)")
        self.carbonyl_tree.heading("carbonyl_type", text="Carbonyl Type")
        self.carbonyl_tree.heading("bs_int", text="B/S")

        # Configure column widths
        self.carbonyl_tree.column("acceptor_res", width=130)
        self.carbonyl_tree.column("acceptor_atom", width=120)
        self.carbonyl_tree.column("carbonyl_res", width=130)
        self.carbonyl_tree.column("carbonyl_atoms", width=120)
        self.carbonyl_tree.column("distance", width=140)
        self.carbonyl_tree.column("angle", width=150)
        self.carbonyl_tree.column("carbonyl_type", width=120)
        self.carbonyl_tree.column("bs_int", width=70)

        # Add scrollbars
        carbonyl_v_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.VERTICAL, command=self.carbonyl_tree.yview
        )
        carbonyl_h_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.HORIZONTAL, command=self.carbonyl_tree.xview
        )
        self.carbonyl_tree.configure(
            yscrollcommand=carbonyl_v_scrollbar.set,
            xscrollcommand=carbonyl_h_scrollbar.set,
        )

        self.carbonyl_tree.grid(row=0, column=0, sticky="nsew")
        carbonyl_v_scrollbar.grid(row=0, column=1, sticky="ns")
        carbonyl_h_scrollbar.grid(row=1, column=0, sticky="ew")

        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        # Add info label
        info_frame = ttk.Frame(carbonyl_frame)
        info_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(
            info_frame,
            text="Carbonyl n→π* Interactions: Lone pair electron donation to carbonyl π* orbitals",
        ).pack(side=tk.LEFT)

        # Add search functionality
        search_frame = ttk.Frame(carbonyl_frame)
        search_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.carbonyl_search_var = tk.StringVar()
        search_entry = ttk.Entry(
            search_frame, textvariable=self.carbonyl_search_var, width=30
        )
        search_entry.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            search_frame,
            text="Filter",
            command=lambda: self._filter_results(
                self.carbonyl_tree, self.carbonyl_search_var.get()
            ),
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            search_frame,
            text="Clear",
            command=lambda: self._clear_filter(
                self.carbonyl_tree, self.carbonyl_search_var
            ),
        ).pack(side=tk.LEFT, padx=5)

    def _create_n_pi_interactions_tab(self):
        """Create n→π* interactions results tab."""
        n_pi_frame = ttk.Frame(self.notebook)
        self.notebook.add(n_pi_frame, text="n→π* Interactions")

        # Create treeview for n→π* interactions
        tree_frame = ttk.Frame(n_pi_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        columns = (
            "donor_res",
            "donor_atom",
            "pi_res",
            "pi_atoms",
            "distance",
            "angle",
            "donor_element",
            "bs_int",
        )

        self.n_pi_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=15
        )

        # Configure columns
        self.n_pi_tree.heading("donor_res", text="Donor Residue")
        self.n_pi_tree.heading("donor_atom", text="Donor Atom")
        self.n_pi_tree.heading("pi_res", text="π Residue")
        self.n_pi_tree.heading("pi_atoms", text="π Ring Atoms")
        self.n_pi_tree.heading("distance", text="Distance (Å)")
        self.n_pi_tree.heading("angle", text="Angle (°)")
        self.n_pi_tree.heading("donor_element", text="Donor Element")
        self.n_pi_tree.heading("bs_int", text="B/S")

        # Configure column widths
        self.n_pi_tree.column("donor_res", width=120)
        self.n_pi_tree.column("donor_atom", width=120)
        self.n_pi_tree.column("pi_res", width=120)
        self.n_pi_tree.column("pi_atoms", width=140)
        self.n_pi_tree.column("distance", width=110)
        self.n_pi_tree.column("angle", width=110)
        self.n_pi_tree.column("donor_element", width=120)
        self.n_pi_tree.column("bs_int", width=70)

        # Add scrollbars
        n_pi_v_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.VERTICAL, command=self.n_pi_tree.yview
        )
        n_pi_h_scrollbar = ttk.Scrollbar(
            tree_frame, orient=tk.HORIZONTAL, command=self.n_pi_tree.xview
        )
        self.n_pi_tree.configure(
            yscrollcommand=n_pi_v_scrollbar.set, xscrollcommand=n_pi_h_scrollbar.set
        )

        self.n_pi_tree.grid(row=0, column=0, sticky="nsew")
        n_pi_v_scrollbar.grid(row=0, column=1, sticky="ns")
        n_pi_h_scrollbar.grid(row=1, column=0, sticky="ew")

        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        # Add info label
        info_frame = ttk.Frame(n_pi_frame)
        info_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(
            info_frame,
            text="n→π* Interactions: Lone pair electrons (O, N, S) interacting with aromatic π systems",
        ).pack(side=tk.LEFT)

        # Add search functionality
        search_frame = ttk.Frame(n_pi_frame)
        search_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.n_pi_search_var = tk.StringVar()
        search_entry = ttk.Entry(
            search_frame, textvariable=self.n_pi_search_var, width=30
        )
        search_entry.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            search_frame,
            text="Filter",
            command=lambda: self._filter_results(
                self.n_pi_tree, self.n_pi_search_var.get()
            ),
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            search_frame,
            text="Clear",
            command=lambda: self._clear_filter(self.n_pi_tree, self.n_pi_search_var),
        ).pack(side=tk.LEFT, padx=5)

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
        self._update_hydrogen_bonds()
        self._update_halogen_bonds()
        self._update_pi_interactions()
        self._update_pi_pi_stacking()
        self._update_carbonyl_interactions()
        self._update_n_pi_interactions()
        self._update_cooperativity_chains()

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
                self.summary_text.insert(tk.END, f"  Bond Types:\n")
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
                self.summary_text.insert(tk.END, f"  Bond Types:\n")
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

                self.summary_text.insert(tk.END, f"  Chain Types:\n")
                for chain_type, count in sorted(type_counts.items()):
                    self.summary_text.insert(tk.END, f"    {chain_type}: {count}\n")

            # Chain length distribution
            if "chain_lengths" in coop_data:
                self.summary_text.insert(tk.END, f"  Chain Length Distribution:\n")
                for length, count in sorted(coop_data["chain_lengths"].items()):
                    self.summary_text.insert(
                        tk.END, f"    Length {length}: {count} chains\n"
                    )
            self.summary_text.insert(tk.END, "\n")

        # Add some example interactions
        if self.analyzer.hydrogen_bonds:
            self.summary_text.insert(tk.END, "Sample Hydrogen Bonds:\n", "subheader")
            for i, hb in enumerate(self.analyzer.hydrogen_bonds[:5]):
                self.summary_text.insert(tk.END, f"  {i+1}. {hb}\n")
            if len(self.analyzer.hydrogen_bonds) > 5:
                self.summary_text.insert(
                    tk.END,
                    f"  ... and {len(self.analyzer.hydrogen_bonds) - 5} more\n\n",
                )

        if self.analyzer.halogen_bonds:
            self.summary_text.insert(tk.END, "Sample Halogen Bonds:\n", "subheader")
            for i, xb in enumerate(self.analyzer.halogen_bonds[:3]):
                self.summary_text.insert(tk.END, f"  {i+1}. {xb}\n")
            if len(self.analyzer.halogen_bonds) > 3:
                self.summary_text.insert(
                    tk.END, f"  ... and {len(self.analyzer.halogen_bonds) - 3} more\n\n"
                )

        if self.analyzer.pi_interactions:
            self.summary_text.insert(tk.END, "Sample π Interactions:\n", "subheader")
            for i, pi in enumerate(self.analyzer.pi_interactions[:3]):
                self.summary_text.insert(tk.END, f"  {i+1}. {pi}\n")
            if len(self.analyzer.pi_interactions) > 3:
                self.summary_text.insert(
                    tk.END, f"  ... and {len(self.analyzer.pi_interactions) - 3} more\n"
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
                    hb.hydrogen.name,
                    hb.acceptor_residue,
                    hb.acceptor.name,
                    f"{hb.distance:.2f}",
                    f"{math.degrees(hb.angle):.1f}",
                    f"{hb.donor_acceptor_distance:.2f}",
                    hb.bond_type,
                    hb.donor_acceptor_properties,
                    hb.get_backbone_sidechain_interaction(),
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
                    xb.donor_atom.name,
                    xb.halogen.name,
                    xb.acceptor_residue,
                    xb.acceptor.name,
                    f"{xb.distance:.2f}",
                    f"{math.degrees(xb.angle):.1f}",
                    xb.bond_type,
                    xb.get_backbone_sidechain_interaction(),
                    xb.donor_acceptor_properties,
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
                    pi.get_interaction_type_display(),
                    pi.donor_acceptor_properties,
                    pi.get_backbone_sidechain_interaction(),
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
            elif tree == self.pi_pi_tree:
                self._update_pi_pi_stacking()
            elif tree == self.carbonyl_tree:
                self._update_carbonyl_interactions()
            elif tree == self.n_pi_tree:
                self._update_n_pi_interactions()
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

        # Clear treeviews
        for item in self.hb_tree.get_children():
            self.hb_tree.delete(item)

        for item in self.xb_tree.get_children():
            self.xb_tree.delete(item)

        for item in self.pi_tree.get_children():
            self.pi_tree.delete(item)

        for item in self.coop_tree.get_children():
            self.coop_tree.delete(item)

        # Clear new interaction tree views
        for item in self.pi_pi_tree.get_children():
            self.pi_pi_tree.delete(item)

        for item in self.carbonyl_tree.get_children():
            self.carbonyl_tree.delete(item)

        for item in self.n_pi_tree.get_children():
            self.n_pi_tree.delete(item)

        # Clear all search filters
        self.hb_search_var.set("")
        self.xb_search_var.set("")
        self.pi_search_var.set("")
        self.pi_pi_search_var.set("")
        self.carbonyl_search_var.set("")
        self.n_pi_search_var.set("")
        self.coop_search_var.set("")

        # Add placeholder text
        self.summary_text.insert(tk.END, "No analysis results available.\n\n")
        self.summary_text.insert(
            tk.END, "Please load a PDB file and run analysis to see results."
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

    def _update_pi_pi_stacking(self):
        """Update the π-π stacking tab."""
        if not self.analyzer:
            return

        # Clear existing items
        for item in self.pi_pi_tree.get_children():
            self.pi_pi_tree.delete(item)

        # Add π-π stacking interactions if they exist
        if hasattr(self.analyzer, "pi_pi_interactions"):
            for interaction in self.analyzer.pi_pi_interactions:
                ring1_res = self._parse_residue_string(interaction.ring1_residue)
                ring1_atoms = (
                    ",".join([atom.name for atom in interaction.ring1_atoms[:3]])
                    + "..."
                )
                ring2_res = self._parse_residue_string(interaction.ring2_residue)
                ring2_atoms = (
                    ",".join([atom.name for atom in interaction.ring2_atoms[:3]])
                    + "..."
                )

                bs_int = "B" if interaction.is_between_residues else "S"

                self.pi_pi_tree.insert(
                    "",
                    tk.END,
                    values=(
                        ring1_res,
                        ring1_atoms,
                        ring2_res,
                        ring2_atoms,
                        f"{interaction.distance:.2f}",
                        f"{interaction.plane_angle:.1f}",
                        f"{interaction.offset:.2f}",
                        interaction.stacking_type,
                        bs_int,
                    ),
                )

    def _update_carbonyl_interactions(self):
        """Update the carbonyl interactions tab."""
        if not self.analyzer:
            return

        # Clear existing items
        for item in self.carbonyl_tree.get_children():
            self.carbonyl_tree.delete(item)

        # Add carbonyl interactions if they exist
        if hasattr(self.analyzer, "carbonyl_interactions"):
            for interaction in self.analyzer.carbonyl_interactions:
                acceptor_res = self._parse_residue_string(interaction.acceptor_residue)
                acceptor_atom = interaction.acceptor_carbon.name
                carbonyl_res = self._parse_residue_string(interaction.donor_residue)
                carbonyl_atoms = (
                    f"{interaction.donor_carbon.name}={interaction.donor_oxygen.name}"
                )

                bs_int = "B" if interaction.is_between_residues else "S"

                self.carbonyl_tree.insert(
                    "",
                    tk.END,
                    values=(
                        acceptor_res,
                        acceptor_atom,
                        carbonyl_res,
                        carbonyl_atoms,
                        f"{interaction.distance:.2f}",
                        f"{interaction.burgi_dunitz_angle:.1f}",
                        interaction.carbonyl_type,
                        bs_int,
                    ),
                )

    def _update_n_pi_interactions(self):
        """Update the n→π* interactions tab."""
        if not self.analyzer:
            return

        # Clear existing items
        for item in self.n_pi_tree.get_children():
            self.n_pi_tree.delete(item)

        # Add n→π* interactions if they exist
        if hasattr(self.analyzer, "n_pi_interactions"):
            for interaction in self.analyzer.n_pi_interactions:
                donor_res = self._parse_residue_string(interaction.donor_residue)
                donor_atom = interaction.lone_pair_atom.name
                pi_res = self._parse_residue_string(interaction.acceptor_residue)
                pi_atoms = (
                    ",".join([atom.name for atom in interaction.pi_atoms[:3]]) + "..."
                )

                bs_int = "B" if interaction.is_between_residues else "S"

                self.n_pi_tree.insert(
                    "",
                    tk.END,
                    values=(
                        donor_res,
                        donor_atom,
                        pi_res,
                        pi_atoms,
                        f"{interaction.distance:.2f}",
                        f"{interaction.angle_to_plane:.1f}",
                        interaction.donor_element,
                        bs_int,
                    ),
                )
