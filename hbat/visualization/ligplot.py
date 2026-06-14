"""
Simple ligplot visualization showing ligand with highlighted interaction atoms.
"""

from html import escape
from typing import Dict, List, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from .minimal_pdb_extractor import _format_atom_as_pdb
import urllib.parse
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
        # Track which interaction types are present in current ligand
        self.applicable_interactions = set()

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

    def _interaction_color_key(self, interaction) -> Optional[str]:
        """Return the legend/color key for an interaction object."""
        return {
            "HydrogenBond": "HydrogenBond",
            "HalogenBond": "HalogenBond",
            "PiInteraction": "PiInteraction",
            "PiPiInteraction": "PiPiInteraction",
            "CarbonylInteraction": "Carbonyl",
            "NPiInteraction": "NPiInteraction",
            "WaterBridge": "WaterBridge",
        }.get(interaction.__class__.__name__)

    @staticmethod
    def _residue_id(atom) -> Optional[str]:
        """Return an atom's canonical residue identifier."""
        if not all(hasattr(atom, attr) for attr in ("chain_id", "res_name", "res_seq")):
            return None
        return f"{atom.chain_id}:{atom.res_name}:{atom.res_seq}"

    @staticmethod
    def _coords_array(atom) -> Optional[np.ndarray]:
        """Return an atom-like object's coordinates as a NumPy array."""
        coords = getattr(atom, "coords", None)
        if coords is None:
            return None
        try:
            return np.array([coords.x, coords.y, coords.z], dtype=float)
        except (AttributeError, TypeError, ValueError):
            return None

    def _get_ligand_interactions(self, ligand_res_id: str) -> List:
        """Return interactions involving the selected ligand."""
        container = getattr(self.analyzer, "ligand_interactions", None)
        interactions = getattr(container, "interactions", [])
        return [
            interaction
            for interaction in interactions
            if ligand_res_id
            in (
                interaction.get_donor_residue(),
                interaction.get_acceptor_residue(),
            )
        ]

    def _get_residue_atoms(self, residue_id: str) -> List:
        """Return all parsed atoms belonging to a residue."""
        return [
            atom
            for atom in self.analyzer.parser.atoms
            if self._residue_id(atom) == residue_id
        ]

    def _get_ligand_anchor_indices(self, interaction, ligand_res_id: str) -> List[int]:
        """Return RDKit atom indices representing the ligand interaction endpoint."""
        serials = []
        for endpoint in (interaction.get_donor(), interaction.get_acceptor()):
            if self._residue_id(endpoint) == ligand_res_id and hasattr(endpoint, "serial"):
                serials.append(endpoint.serial)

        for atom in getattr(interaction, "pi_atoms", []):
            if self._residue_id(atom) == ligand_res_id:
                serials.append(atom.serial)

        serial_to_index = {
            info["serial"]: atom_index
            for atom_index, info in self.pdb_ligand_atom_info.items()
        }
        return sorted(
            {
                serial_to_index[serial]
                for serial in serials
                if serial in serial_to_index
            }
        )

    def _get_partner_classification(self, interaction, partner_residue: str) -> str:
        """Return backbone/side-chain classification for a partner endpoint."""
        partner_atoms = []
        for endpoint in (interaction.get_donor(), interaction.get_acceptor()):
            if self._residue_id(endpoint) == partner_residue:
                partner_atoms.append(endpoint)
        partner_atoms.extend(
            atom
            for atom in getattr(interaction, "pi_atoms", [])
            if self._residue_id(atom) == partner_residue
        )

        classifications = {
            getattr(atom, "backbone_sidechain", "")
            for atom in partner_atoms
            if getattr(atom, "backbone_sidechain", "") in {"B", "S"}
        }
        if classifications == {"B"}:
            return "backbone"
        if classifications == {"S"}:
            return "sidechain"
        if classifications == {"B", "S"}:
            return "mixed"
        return "unknown"

    def _build_residue_interactions(self, ligand_res_id: str) -> Dict[str, Dict]:
        """Group direct ligand interactions by the non-ligand residue."""
        grouped = {}
        for interaction in self._get_ligand_interactions(ligand_res_id):
            donor_residue = interaction.get_donor_residue()
            acceptor_residue = interaction.get_acceptor_residue()
            partner_residue = (
                acceptor_residue if donor_residue == ligand_res_id else donor_residue
            )
            if partner_residue in ("Unknown", ligand_res_id):
                continue

            color_key = self._interaction_color_key(interaction)
            if color_key is None:
                continue

            residue_atoms = self._get_residue_atoms(partner_residue)
            residue_coords = [
                coords
                for coords in (self._coords_array(atom) for atom in residue_atoms)
                if coords is not None
            ]
            if not residue_coords:
                continue

            entry = grouped.setdefault(
                partner_residue,
                {
                    "position": np.mean(residue_coords, axis=0),
                    "interactions": [],
                    "classifications": set(),
                },
            )
            entry["classifications"].add(
                self._get_partner_classification(interaction, partner_residue)
            )
            entry["interactions"].append(
                {
                    "type": color_key,
                    "anchor_indices": self._get_ligand_anchor_indices(
                        interaction, ligand_res_id
                    ),
                }
            )
            self.applicable_interactions.add(color_key)
        return grouped

    @staticmethod
    def _aggregate_classification(classifications: set) -> str:
        """Aggregate all interaction classifications for one residue."""
        known = classifications - {"unknown"}
        if known == {"backbone"}:
            return "backbone"
        if known == {"sidechain"}:
            return "sidechain"
        if known:
            return "mixed"
        return "unknown"

    def _layout_residues(
        self,
        grouped: Dict[str, Dict],
        ligand_atoms: List,
        width: int,
        height: int,
        ligand_bounds: Optional[tuple] = None,
    ) -> Dict[str, tuple]:
        """Project residues into deterministic, non-overlapping perimeter slots."""
        if not grouped:
            return {}

        ligand_coords = [
            coords
            for coords in (self._coords_array(atom) for atom in ligand_atoms)
            if coords is not None
        ]
        ligand_center = np.mean(ligand_coords, axis=0)
        residue_ids = list(grouped)
        vectors = np.array(
            [grouped[residue_id]["position"] - ligand_center for residue_id in residue_ids]
        )

        if len(vectors) >= 2:
            _, _, axes = np.linalg.svd(vectors, full_matrices=False)
            projected = vectors @ axes[:2].T
        else:
            projected = vectors[:, :2]

        sides = {"top": [], "right": [], "bottom": [], "left": []}
        for residue_id, position in zip(residue_ids, projected):
            x, y = float(position[0]), float(-position[1])
            if abs(x) >= abs(y):
                side = "right" if x >= 0 else "left"
                sort_value = y
            else:
                side = "bottom" if y >= 0 else "top"
                sort_value = x
            sides[side].append((sort_value, residue_id))

        box_width, box_height, gap = 96, 42, 12
        horizontal_capacity = max(1, int((width - gap) // (box_width + gap)))
        vertical_capacity = max(1, int((height - gap) // (box_height + gap)))
        positions = {}

        for side, entries in sides.items():
            entries.sort()
            capacity = (
                horizontal_capacity if side in {"top", "bottom"} else vertical_capacity
            )
            for index, (_, residue_id) in enumerate(entries):
                lane, slot = divmod(index, capacity)
                lane_entries = min(capacity, len(entries) - lane * capacity)
                if side in {"top", "bottom"}:
                    spacing = width / (lane_entries + 1)
                    x = spacing * (slot + 1)
                    y = (
                        box_height / 2 + gap + lane * (box_height + gap)
                        if side == "top"
                        else height
                        - box_height / 2
                        - gap
                        - lane * (box_height + gap)
                    )
                else:
                    spacing = height / (lane_entries + 1)
                    y = spacing * (slot + 1)
                    x = (
                        width - box_width / 2 - gap - lane * (box_width + gap)
                        if side == "right"
                        else box_width / 2 + gap + lane * (box_width + gap)
                    )
                positions[residue_id] = (x, y)

        return positions

    @staticmethod
    def _hex_color(rgb: tuple) -> str:
        """Convert an RDKit-style RGB tuple to a CSS hex color."""
        return f"#{int(rgb[0] * 255):02x}{int(rgb[1] * 255):02x}{int(rgb[2] * 255):02x}"

    @staticmethod
    def _box_edge_point(
        anchor_x: float, anchor_y: float, box_x: float, box_y: float
    ) -> tuple:
        """Return the intersection between an anchor-to-box line and box edge."""
        half_width, half_height = 48, 21
        delta_x, delta_y = anchor_x - box_x, anchor_y - box_y
        if abs(delta_x) < 1e-8 and abs(delta_y) < 1e-8:
            return box_x, box_y
        scale = min(
            half_width / abs(delta_x) if abs(delta_x) > 1e-8 else float("inf"),
            half_height / abs(delta_y) if abs(delta_y) > 1e-8 else float("inf"),
        )
        return box_x + delta_x * scale, box_y + delta_y * scale

    def _create_interactive_overlay(
        self,
        grouped: Dict[str, Dict],
        positions: Dict[str, tuple],
        anchor_positions: Dict[int, tuple],
        width: int,
        height: int,
    ) -> tuple:
        """Create SVG fragments for hoverable interaction lines and residue nodes."""
        lines = []
        nodes = []
        center_x, center_y = width / 2, height / 2

        for residue_index, (residue_id, entry) in enumerate(grouped.items()):
            x, y = positions[residue_id]
            group_class = f"ligplot-residue-{residue_index}"
            interaction_types = []

            for interaction in entry["interactions"]:
                interaction_type = interaction["type"]
                interaction_types.append(interaction_type)
                color = self._hex_color(self.colors[interaction_type])
                anchors = interaction["anchor_indices"]
                draw_anchors = [
                    anchor_positions[anchor_index]
                    for anchor_index in anchors
                    if anchor_index in anchor_positions
                ]
                if draw_anchors:
                    anchor_x = sum(point[0] for point in draw_anchors) / len(
                        draw_anchors
                    )
                    anchor_y = sum(point[1] for point in draw_anchors) / len(
                        draw_anchors
                    )
                else:
                    anchor_x, anchor_y = center_x, center_y
                box_edge_x, box_edge_y = self._box_edge_point(anchor_x, anchor_y, x, y)
                lines.append(
                    f'<line class="ligplot-line {group_class}" '
                    f'x1="{anchor_x:.1f}" y1="{anchor_y:.1f}" '
                    f'x2="{box_edge_x:.1f}" y2="{box_edge_y:.1f}" '
                    f'stroke="{color}" />'
                )

            chain, residue_name, residue_number = residue_id.split(":", 2)
            classification = self._aggregate_classification(entry["classifications"])
            title = escape(
                f"{residue_id}: {classification}; "
                f"{', '.join(sorted(set(interaction_types)))}"
            )
            nodes.append(
                f'<g class="ligplot-residue-node ligplot-{classification} {group_class}" '
                f'data-target="{group_class}" transform="translate({x:.1f} {y:.1f})">'
                f"<title>{title}</title>"
                '<rect x="-48" y="-21" width="96" height="42" rx="7" />'
                f'<text y="4">{escape(chain)}:{escape(residue_name)}:'
                f"{escape(residue_number)}</text>"
                "</g>"
            )

        style = """
<style>
  .ligplot-line {
    opacity: 0;
    stroke-width: 2.5;
    stroke-dasharray: 6 5;
    stroke-linecap: round;
    pointer-events: none;
    transition: opacity 120ms ease;
  }
  svg:has(.ligplot-ligand-focus:hover) .ligplot-line {
    opacity: 0.72;
  }
  .ligplot-ligand-focus {
    cursor: crosshair;
  }
  .ligplot-residue-node {
    cursor: pointer;
  }
  .ligplot-residue-node rect {
    fill: #ffffff;
    stroke: #374151;
    stroke-width: 2;
    transition: fill 120ms ease, stroke-width 120ms ease;
  }
  .ligplot-residue-node.ligplot-backbone rect {
    fill: #e0f2fe;
    stroke: #0284c7;
  }
  .ligplot-residue-node.ligplot-sidechain rect {
    fill: #dcfce7;
    stroke: #16a34a;
  }
  .ligplot-residue-node.ligplot-mixed rect {
    fill: #f3e8ff;
    stroke: #9333ea;
  }
  .ligplot-residue-node text {
    fill: #111827;
    font-family: Arial, sans-serif;
    font-size: 11px;
    font-weight: 700;
    letter-spacing: 0;
    pointer-events: none;
    text-anchor: middle;
  }
  .ligplot-residue-node:hover rect {
    stroke-width: 3;
  }
"""
        for residue_index in range(len(grouped)):
            group_class = f"ligplot-residue-{residue_index}"
            style += (
                f"  .ligplot-residue-node.{group_class}:hover "
                f"~ .ligplot-unused {{ opacity: 1; }}\n"
                f"  svg:has(.ligplot-residue-node.{group_class}:hover) "
                f".ligplot-line.{group_class} {{ opacity: 1; }}\n"
            )
        style += "</style>"

        interaction_layer = (
            '<g class="ligplot-interaction-layer">' + "".join(lines) + "</g>"
        )
        ligand_focus = (
            f'<rect class="ligplot-ligand-focus" x="{width * 0.19:.1f}" '
            f'y="{height * 0.19:.1f}" width="{width * 0.62:.1f}" '
            f'height="{height * 0.62:.1f}" rx="12" fill="transparent">'
            "<title>Hover to show all ligand interactions</title></rect>"
        )
        return (
            style + interaction_layer,
            ligand_focus + "".join(nodes),
        )

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
                a
                for a in self.analyzer.parser.atoms
                if f"{a.chain_id}:{a.res_name}:{a.res_seq}" == self.residue_id
            ]
        else:
            # Fallback: use residue name only (for backward compatibility)
            ligand_atoms = [
                a for a in self.analyzer.parser.atoms if a.res_name == self.ligand_name
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
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".pdb", delete=False
            ) as f:
                f.write(pdb_block)
                pdb_path = f.name

            # Use obabel to convert and infer bond orders
            sdf_path = pdb_path.replace(".pdb", ".sdf")
            result = subprocess.run(
                ["obabel", pdb_path, "-O", sdf_path, "-h"],
                capture_output=True,
                timeout=5,
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
            parsed = urllib.parse.urlparse(rcsb_url)
            if parsed.scheme != "https" or not parsed.netloc.endswith("rcsb.org"):
                return None
            # URL is validated to be https://data.rcsb.org only
            with urllib.request.urlopen(rcsb_url, timeout=5) as response:  # nosec B310
                data = json.loads(response.read())
                if "rcsb_chem_comp_descriptor" in data:
                    descriptor = data["rcsb_chem_comp_descriptor"]
                    return descriptor.get("SMILES")
        except Exception as e:
            print(f"Failed to get SMILES for {self.ligand_name}: {e}")
        return None

    def _collect_interaction_atoms(
        self, interactions, interaction_type: str, ligand_res_id: str, interacting: Dict
    ) -> None:
        """Collect atoms from a set of interactions."""
        for interaction in interactions:
            if interaction.donor_residue == ligand_res_id:
                interacting[interaction.donor.serial] = interaction_type
                self.applicable_interactions.add(interaction_type)
            elif interaction.acceptor_residue == ligand_res_id:
                interacting[interaction.acceptor.serial] = interaction_type
                self.applicable_interactions.add(interaction_type)

    def get_interacting_atoms(self) -> Dict[int, tuple]:
        """
        Get interacting atoms and their colors.

        :returns: Dict mapping RDKit atom index to RGB color tuple
        :rtype: dict
        """
        ligand_res_id = self._find_ligand_residue_id()
        if not ligand_res_id:
            self.applicable_interactions.clear()
            return {}

        interacting = {}
        self.applicable_interactions.clear()

        # Collect all interaction types
        self._collect_interaction_atoms(
            self.analyzer.hydrogen_bonds, "HydrogenBond", ligand_res_id, interacting
        )
        if hasattr(self.analyzer, "halogen_bonds"):
            self._collect_interaction_atoms(
                self.analyzer.halogen_bonds, "HalogenBond", ligand_res_id, interacting
            )
        if hasattr(self.analyzer, "pi_interactions"):
            self._collect_interaction_atoms(
                self.analyzer.pi_interactions,
                "PiInteraction",
                ligand_res_id,
                interacting,
            )
        if hasattr(self.analyzer, "pi_pi_interactions"):
            self._collect_interaction_atoms(
                self.analyzer.pi_pi_interactions,
                "PiPiInteraction",
                ligand_res_id,
                interacting,
            )
        if hasattr(self.analyzer, "water_bridges"):
            self._collect_interaction_atoms(
                self.analyzer.water_bridges, "WaterBridge", ligand_res_id, interacting
            )

        # Map PDB serials to RDKit indices by position
        ligand_atoms = [
            atom
            for atom in self.analyzer.parser.atoms
            if f"{atom.chain_id}:{atom.res_name}:{atom.res_seq}" == ligand_res_id
        ]

        rdkit_indices = {}
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

        ligand_res_id = self._find_ligand_residue_id()
        if not ligand_res_id:
            return None

        ligand_atoms = self._get_residue_atoms(ligand_res_id)
        self.applicable_interactions.clear()
        grouped = self._build_residue_interactions(ligand_res_id)
        atom_highlights = {}
        for entry in grouped.values():
            for interaction in entry["interactions"]:
                for atom_index in interaction["anchor_indices"]:
                    colors = atom_highlights.setdefault(atom_index, [])
                    color = self.colors[interaction["type"]]
                    if color not in colors:
                        colors.append(color)

        # Add PDB atom names as notes only for highlighted atoms
        for atom_idx in atom_highlights.keys():
            if atom_idx in self.pdb_ligand_atom_info:
                atom_name = self.pdb_ligand_atom_info[atom_idx]["name"]
                self.mol.GetAtomWithIdx(atom_idx).SetProp("atomNote", atom_name)

        # Draw molecule with highlights
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.drawOptions().fixedBondLength = 20
        drawer.drawOptions().clearBackground = False

        # Prepare highlight data
        drawer.DrawMoleculeWithHighlights(
            self.mol,
            "",
            atom_highlights,
            {},
            {},
            {},
        )

        anchor_positions = {}
        for atom_index in atom_highlights:
            point = drawer.GetDrawCoords(atom_index)
            anchor_positions[atom_index] = (point.x, point.y)

        ligand_draw_positions = [
            drawer.GetDrawCoords(atom_index)
            for atom_index in range(self.mol.GetNumAtoms())
        ]
        ligand_bounds = (
            min(point.x for point in ligand_draw_positions),
            min(point.y for point in ligand_draw_positions),
            max(point.x for point in ligand_draw_positions),
            max(point.y for point in ligand_draw_positions),
        )
        residue_positions = self._layout_residues(
            grouped, ligand_atoms, width, height, ligand_bounds
        )

        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        interaction_layer, residue_nodes = self._create_interactive_overlay(
            grouped, residue_positions, anchor_positions, width, height
        )
        svg_start = svg.find("<svg")
        opening_end = svg.find(">", svg_start)
        closing_start = svg.rfind("</svg>")
        if svg_start < 0 or opening_end < 0 or closing_start < 0:
            return svg

        return (
            svg[: opening_end + 1]
            + interaction_layer
            + svg[opening_end + 1 : closing_start]
            + residue_nodes
            + svg[closing_start:]
        ).replace("<svg ", '<svg style="max-width: 100%; height: auto;" ', 1)

    def _create_color_legend(self) -> str:
        """
        Create HTML color legend showing only applicable interaction types.

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
            # Only show interactions applicable to this ligand
            if (
                interaction_type in self.applicable_interactions
                and interaction_type in self.colors
            ):
                # Convert RGB tuple to hex color
                rgb = self.colors[interaction_type]
                hex_color = f"#{int(rgb[0] * 255):02x}{int(rgb[1] * 255):02x}{int(rgb[2] * 255):02x}"

                legend_html += f"""
                <div style="display: flex; align-items: center; font-size: 11px;">
                    <div style="width: 14px; height: 14px; background: {hex_color}; border-radius: 2px; margin-right: 6px;"></div>
                    <span>{label}</span>
                </div>
                """

        for label, fill, border in (
            ("Backbone", "#e0f2fe", "#0284c7"),
            ("Side-chain", "#dcfce7", "#16a34a"),
            ("Mixed", "#f3e8ff", "#9333ea"),
        ):
            legend_html += f"""
            <div style="display: flex; align-items: center; font-size: 11px;">
                <div style="width: 20px; height: 12px; background: {fill}; border: 2px solid {border}; border-radius: 3px; margin-right: 6px;"></div>
                <span>{label}</span>
            </div>
            """

        legend_html += "</div>"
        return legend_html
