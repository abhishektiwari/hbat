"""
Simple ligplot visualization showing ligand with highlighted interaction atoms.
"""

from typing import Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import urllib.request
import json


class LigplotGenerator:
    """Generate simple ligplot SVG showing ligand with highlighted atoms."""

    def __init__(self, ligand_name: str, analyzer):
        """
        Initialize ligplot generator.

        :param ligand_name: Ligand residue name (e.g., "IMP", "ATP")
        :param analyzer: NPMolecularInteractionAnalyzer instance
        """
        self.ligand_name = ligand_name
        self.analyzer = analyzer
        self.mol = None

        # Interaction colors (RGB tuples, values 0-1)
        self.colors = {
            "HydrogenBond": (1,0.5,0.5),  # Blue
            "HalogenBond": (0.7,0.7,0.3),     # Red
            "PiInteraction": (0.6, 0.2, 1), # Purple
            "PiPiInteraction": (0.4, 0, 0.8), # Dark purple
            "WaterBridge": (0, 0.8, 1),     # Cyan
            "Carbonyl": (0.8, 0, 0.8),      # Magenta
            "NPiInteraction": (0, 0.8, 0.4), # Green
        }

    def get_smiles(self) -> Optional[str]:
        """
        Get SMILES for ligand from RCSB PDB API.

        :returns: SMILES string or None
        :rtype: Optional[str]
        """
        rcsb_url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{self.ligand_name}"
        try:
            with urllib.request.urlopen(rcsb_url, timeout=5) as response:
                data = json.loads(response.read())
                if "rcsb_chem_comp_descriptor" in data:
                    descriptor = data["rcsb_chem_comp_descriptor"]
                    return descriptor.get("SMILES")
        except Exception as e:
            print(f"Failed to get SMILES for {self.ligand_name}: {e}")
        return None

    def get_interacting_atoms(self) -> Dict[int, tuple]:
        """
        Get interacting atoms and their colors.

        :returns: Dict mapping RDKit atom index to RGB color tuple
        :rtype: dict
        """
        ligand_res_id = self._find_ligand_residue_id()
        if not ligand_res_id:
            return {}

        interacting = {}

        # H-bonds
        for hbond in self.analyzer.hydrogen_bonds:
            if hbond.donor_residue == ligand_res_id:
                interacting[hbond.donor.serial] = "HydrogenBond"
            elif hbond.acceptor_residue == ligand_res_id:
                interacting[hbond.acceptor.serial] = "HydrogenBond"

        # Halogen bonds
        if hasattr(self.analyzer, "halogen_bonds"):
            for xbond in self.analyzer.halogen_bonds:
                if xbond.donor_residue == ligand_res_id:
                    interacting[xbond.donor.serial] = "HalogenBond"
                elif xbond.acceptor_residue == ligand_res_id:
                    interacting[xbond.acceptor.serial] = "HalogenBond"

        # π-interactions
        if hasattr(self.analyzer, "pi_interactions"):
            for pi_int in self.analyzer.pi_interactions:
                if pi_int.donor_residue == ligand_res_id:
                    interacting[pi_int.donor.serial] = "PiInteraction"
                elif pi_int.acceptor_residue == ligand_res_id:
                    interacting[pi_int.acceptor.serial] = "PiInteraction"

        # π–π stacking interactions
        if hasattr(self.analyzer, "pi_pi_interactions"):
            for pipi in self.analyzer.pi_pi_interactions:
                if pipi.donor_residue == ligand_res_id:
                    interacting[pipi.donor.serial] = "PiPiInteraction"
                elif pipi.acceptor_residue == ligand_res_id:
                    interacting[pipi.acceptor.serial] = "PiPiInteraction"

        # Water bridges
        if hasattr(self.analyzer, "water_bridges"):
            for wb in self.analyzer.water_bridges:
                if wb.donor_residue == ligand_res_id:
                    interacting[wb.donor.serial] = "WaterBridge"
                elif wb.acceptor_residue == ligand_res_id:
                    interacting[wb.acceptor.serial] = "WaterBridge"

        # Map PDB serials to RDKit indices by position
        # Get ligand atoms in PDB order
        ligand_atoms = [
            atom
            for atom in self.analyzer.parser.atoms
            if f"{atom.chain_id}:{atom.res_name}:{atom.res_seq}" == ligand_res_id
        ]

        rdkit_indices = {}
        for i, pdb_atom in enumerate(ligand_atoms):
            if pdb_atom.serial in interacting and i < self.mol.GetNumAtoms():
                rdkit_indices[i] = self.colors[interacting[pdb_atom.serial]]

        return rdkit_indices

    def _find_ligand_residue_id(self) -> Optional[str]:
        """Find ligand residue ID."""
        for atom in self.analyzer.parser.atoms:
            if atom.res_name == self.ligand_name:
                return f"{atom.chain_id}:{atom.res_name}:{atom.res_seq}"
        return None

    def generate_svg(self, width: int = 300, height: int = 250) -> str:
        """
        Generate ligplot SVG with color legend.

        :param width: SVG width in pixels
        :param height: SVG height in pixels
        :returns: HTML string with SVG and legend
        :rtype: str
        """
        # Get SMILES and create molecule
        smiles = self.get_smiles()
        if not smiles:
            return f"<p>Could not find SMILES for {self.ligand_name}</p>"

        self.mol = Chem.MolFromSmiles(smiles)
        if not self.mol:
            return f"<p>Invalid SMILES for {self.ligand_name}</p>"

        # Generate 2D coordinates
        AllChem.Compute2DCoords(self.mol)

        # Get interacting atoms
        atom_highlights = self.get_interacting_atoms()

        # Draw molecule with highlights
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

        # Enable atom indices display
        # drawer.drawOptions().addAtomIndices = True

        # Prepare highlight data
        highlight_atoms = list(atom_highlights.keys())
        highlight_colors = {k: v for k, v in atom_highlights.items()}

        drawer.DrawMolecule(
            self.mol,
            highlightAtoms=highlight_atoms,
            highlightAtomColors=highlight_colors,
        )
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        # Create color legend
        legend_html = self._create_color_legend()

        # Combine SVG and legend in a vertical flex container (diagram center, legend bottom)
        html = f"""
        <div style="display: flex; flex-direction: column; gap: 15px; align-items: center;">
            <div style="display: flex; justify-content: center;">
                {svg}
            </div>
            <div>
                {legend_html}
            </div>
        </div>
        """
        return html

    def _create_color_legend(self) -> str:
        """
        Create HTML color legend showing interaction types from colors dictionary.

        :returns: HTML string with legend
        :rtype: str
        """
        # Map interaction type keys to display labels
        type_to_label = {
            "HydrogenBond": "H-Bond",
            "HalogenBond": "Halogen Bond",
            "PiInteraction": "π–Inter",
            "PiPiInteraction": "π–π Stacking",
            "WaterBridge": "Water Bridge",
            "Carbonyl": "Carbonyl",
            "NPiInteraction": "n-π*",
        }

        legend_html = '<div style="border: 1px solid #ccc; padding: 12px 15px; border-radius: 4px; background: #f9f9f9; display: flex; gap: 20px; flex-wrap: wrap; justify-content: center; align-items: center;">'

        for interaction_type, label in type_to_label.items():
            if interaction_type in self.colors:
                # Convert RGB tuple to hex color
                rgb = self.colors[interaction_type]
                hex_color = f"#{int(rgb[0]*255):02x}{int(rgb[1]*255):02x}{int(rgb[2]*255):02x}"

                legend_html += f'''
                <div style="display: flex; align-items: center; font-size: 11px;">
                    <div style="width: 14px; height: 14px; background: {hex_color}; border-radius: 2px; margin-right: 6px;"></div>
                    <span>{label}</span>
                </div>
                '''

        legend_html += '</div>'
        return legend_html
