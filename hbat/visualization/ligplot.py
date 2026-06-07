"""
Simple ligplot visualization showing ligand with highlighted interaction atoms.
"""

from typing import Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from .minimal_pdb_extractor import _format_atom_as_pdb
import urllib.request
import json
import tempfile
import subprocess


class LigplotGenerator:
    """Generate simple ligplot SVG showing ligand with highlighted atoms."""

    def __init__(self, ligand_name: str, analyzer, residue_id: Optional[str] = None):
        """
        Initialize ligplot generator.

        :param ligand_name: Ligand residue name (e.g., "IMP", "ATP")
        :param analyzer: NPMolecularInteractionAnalyzer instance
        :param residue_id: Full residue ID in format "CHAIN:NAME:SEQ" to identify specific ligand.
                          If not provided, uses ligand_name for backward compatibility.
        """
        self.ligand_name = ligand_name
        self.residue_id = residue_id
        self.analyzer = analyzer
        self.mol = None
        # Map RDKit atom index to PDB info: {idx: {"serial": 1862, "name": "O3V"}}
        self.pdb_ligand_atom_info = {}

        # Interaction colors (RGB tuples, values 0-1)
        self.colors = {
            "HydrogenBond": (1, 0.5, 0.5),  # Red/Pink
            "HalogenBond": (0.7, 0.7, 0.3),  # Yellow
            "PiInteraction": (0.6, 0.2, 1),  # Purple
            "PiPiInteraction": (1, 0.65, 0),  # Orange
            "WaterBridge": (0, 0.8, 1),  # Cyan
            "Carbonyl": (0.8, 0, 0.8),  # Magenta
            "NPiInteraction": (0, 0.8, 0.4),  # Green
        }

    def get_ligand_mol(self) -> Optional[Chem.Mol]:
        """
        Get RDKit molecule with correct bond orders and PDB atom reference.

        Uses SMILES from CCD for correct bond orders while tracking PDB serials
        for interaction matching. Falls back to OpenBabel for bond order inference
        if SMILES template matching fails.

        :returns: RDKit Mol object or None
        :rtype: Optional[Chem.Mol]
        """
        # Get ligand atoms from PDB for serial and name reference
        if self.residue_id:
            # Use specific residue ID if provided (format: "CHAIN:NAME:SEQ")
            ligand_atoms = [
                a for a in self.analyzer.parser.atoms
                if f"{a.chain_id}:{a.res_name}:{a.res_seq}" == self.residue_id
            ]
        else:
            # Fallback: use residue name only (for backward compatibility)
            ligand_atoms = [
                a for a in self.analyzer.parser.atoms
                if a.res_name == self.ligand_name
            ]

        if not ligand_atoms:
            return None

        # Build PDB atom info map: {atom_idx: {"serial": int, "name": str}}
        self.pdb_ligand_atom_info = {
            i: {"serial": a.serial, "name": a.name.strip()}
            for i, a in enumerate(ligand_atoms)
        }

        # Create PDB block from ligand atoms (preserves atom order)
        pdb_lines = [_format_atom_as_pdb(a) for a in ligand_atoms]
        pdb_block = "HEADER    LIGAND\n" + "\n".join(pdb_lines) + "\nEND\n"

        # Load molecule from PDB block (correct atom order)
        mol_pdb = Chem.MolFromPDBBlock(pdb_block, removeHs=True, sanitize=False)
        if not mol_pdb:
            return None

        # Try to get SMILES template for correct bond orders
        smiles = self.get_smiles()
        if smiles:
            try:
                template = Chem.MolFromSmiles(smiles)
                if template and template.GetNumAtoms() == mol_pdb.GetNumAtoms():
                    # Transfer bond orders from template to PDB molecule
                    mol = AllChem.AssignBondOrdersFromTemplate(template, mol_pdb)

                    # Remove stereochemistry for flat 2D
                    Chem.RemoveStereochemistry(mol)
                    for bond in mol.GetBonds():
                        bond.SetStereo(Chem.BondStereo.STEREONONE)

                    # Add hydrogens
                    mol = Chem.AddHs(mol)
                    return mol
            except Exception as e:
                # Template matching failed, try OpenBabel for bond order inference
                print(f"Template matching failed for {self.ligand_name}: {e}")

        # Fallback: use OpenBabel to infer bond orders from coordinates
        try:
            mol = self._get_mol_from_obabel(pdb_block)
            if mol:
                return mol
        except Exception as e:
            print(f"OpenBabel inference failed: {e}")

        # Last resort: return molecule with inferred bonds from PDB
        try:
            Chem.SanitizeMol(mol_pdb)
            Chem.RemoveStereochemistry(mol_pdb)
            for bond in mol_pdb.GetBonds():
                bond.SetStereo(Chem.BondStereo.STEREONONE)
            mol_pdb = Chem.AddHs(mol_pdb)
            return mol_pdb
        except Exception as e:
            print(f"Failed to sanitize molecule from PDB block: {e}")
            return None

    def _get_mol_from_obabel(self, pdb_block: str) -> Optional[Chem.Mol]:
        """
        Use OpenBabel to infer bond orders from PDB coordinates.

        :param pdb_block: PDB format text block
        :returns: RDKit Mol object or None
        :rtype: Optional[Chem.Mol]
        """
        try:
            # Write PDB to temp file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                f.write(pdb_block)
                pdb_path = f.name

            # Use obabel to convert and infer bond orders
            sdf_path = pdb_path.replace('.pdb', '.sdf')
            result = subprocess.run(
                ['obabel', pdb_path, '-O', sdf_path, '-h'],
                capture_output=True,
                timeout=5
            )

            if result.returncode == 0:
                # Load from SDF which has explicit bond orders
                with open(sdf_path) as f:
                    mol = Chem.MolFromMolBlock(f.read())
                if mol:
                    # Clean up temp files
                    import os
                    os.unlink(pdb_path)
                    os.unlink(sdf_path)

                    # Remove stereochemistry for flat 2D
                    Chem.RemoveStereochemistry(mol)
                    for bond in mol.GetBonds():
                        bond.SetStereo(Chem.BondStereo.STEREONONE)

                    # Add hydrogens
                    mol = Chem.AddHs(mol)
                    return mol
        except (subprocess.TimeoutExpired, FileNotFoundError):
            # OpenBabel not available or failed
            pass
        except Exception as e:
            print(f"Error in OpenBabel conversion: {e}")

        return None

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
        ligand_atoms = [
            atom
            for atom in self.analyzer.parser.atoms
            if f"{atom.chain_id}:{atom.res_name}:{atom.res_seq}" == ligand_res_id
        ]

        rdkit_indices = {}
        # Check if molecule is set, otherwise use ligand atoms count
        max_atoms = self.mol.GetNumAtoms() if self.mol else len(ligand_atoms)
        for i, pdb_atom in enumerate(ligand_atoms):
            if pdb_atom.serial in interacting and i < max_atoms:
                rdkit_indices[i] = self.colors[interacting[pdb_atom.serial]]

        return rdkit_indices

    def _find_ligand_residue_id(self) -> Optional[str]:
        """Find ligand residue ID."""
        # Use provided residue_id if available
        if self.residue_id:
            return self.residue_id
        # Otherwise, find the first matching residue by name
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
        svg = self._generate_svg_only(width, height)
        if not svg:
            return f"<p>Could not find structure data for {self.ligand_name}</p>"

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

    def get_svg_only(self, width: int = 300, height: int = 250, include_legend: bool = True) -> Optional[str]:
        """
        Generate SVG without HTML wrapper (for downloading/exporting).

        :param width: SVG width in pixels
        :param height: SVG height in pixels
        :param include_legend: If True, add legend as SVG elements
        :returns: SVG string or None
        :rtype: Optional[str]
        """
        svg = self._generate_svg_only(width, height)
        if not svg or not include_legend:
            return svg

        # Add legend to SVG
        return self._add_legend_to_svg(svg, width)

    def _add_legend_to_svg(self, svg: str, width: int) -> str:
        """
        Add legend as SVG elements to the generated SVG.

        :param svg: SVG string from RDKit
        :param width: SVG width for layout
        :returns: SVG with embedded legend
        :rtype: str
        """
        # Map interaction types to labels
        type_to_label = {
            "HydrogenBond": "H-Bond",
            "HalogenBond": "Halogen Bond",
            "PiInteraction": "π–Inter",
            "PiPiInteraction": "π–π Stacking",
            "WaterBridge": "Water Bridge",
            "Carbonyl": "Carbonyl",
            "NPiInteraction": "n-π*",
        }

        # Create legend SVG elements
        legend_items = []
        x_pos = 10
        y_pos = 10
        item_height = 20

        for interaction_type, label in type_to_label.items():
            if interaction_type in self.colors:
                rgb = self.colors[interaction_type]
                hex_color = f"#{int(rgb[0]*255):02x}{int(rgb[1]*255):02x}{int(rgb[2]*255):02x}"

                # Color box
                legend_items.append(
                    f'<rect x="{x_pos}" y="{y_pos}" width="14" height="14" '
                    f'fill="{hex_color}" stroke="black" stroke-width="0.5"/>'
                )
                # Label text
                legend_items.append(
                    f'<text x="{x_pos + 18}" y="{y_pos + 12}" font-size="11" '
                    f'font-family="Arial,sans-serif">{label}</text>'
                )

                x_pos += 140
                if x_pos > width - 150:
                    x_pos = 10
                    y_pos += item_height + 5

        # Insert legend before closing </svg> tag
        legend_group = f'<g id="ligplot-legend">\n' + '\n'.join(legend_items) + '\n</g>'
        svg = svg.replace('</svg>', f'{legend_group}\n</svg>')

        return svg

    def _generate_svg_only(self, width: int = 300, height: int = 250) -> Optional[str]:
        """
        Generate SVG without legend (for downloading/exporting).

        :param width: SVG width in pixels
        :param height: SVG height in pixels
        :returns: SVG string or None
        :rtype: Optional[str]
        """
        # Get molecule with correct bond orders and stereo removed
        self.mol = self.get_ligand_mol()
        if not self.mol:
            return None

        # Generate 2D coordinates
        AllChem.Compute2DCoords(self.mol)

        # Get interacting atoms
        atom_highlights = self.get_interacting_atoms()

        # Add PDB atom names as notes only for highlighted atoms
        for atom_idx in atom_highlights.keys():
            if atom_idx in self.pdb_ligand_atom_info:
                atom_name = self.pdb_ligand_atom_info[atom_idx]["name"]
                self.mol.GetAtomWithIdx(atom_idx).SetProp("atomNote", atom_name)

        # Draw molecule with highlights
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

        # Prepare highlight data
        highlight_atoms = list(atom_highlights.keys())
        highlight_colors = {k: v for k, v in atom_highlights.items()}

        drawer.DrawMolecule(
            self.mol,
            highlightAtoms=highlight_atoms,
            highlightAtomColors=highlight_colors,
        )
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

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
