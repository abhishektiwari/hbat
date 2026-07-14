"""Unit tests for interactive and static LigPlot SVG generation."""

import math
import re
from types import SimpleNamespace

import pytest
from rdkit import Chem

from hbat.core.interactions import HydrogenBond, LigandInteraction, WaterBridge
from hbat.core.np_vector import NPVec3D
from hbat.core.structure import Atom
from hbat.visualization.ligplot import LigplotGenerator


def make_atom(
    serial,
    name,
    residue,
    sequence,
    coords,
    *,
    element="O",
    record_type="ATOM",
    backbone_sidechain="S",
):
    """Create a compact test atom."""
    return Atom(
        serial=serial,
        name=name,
        alt_loc="",
        res_name=residue,
        chain_id="A",
        res_seq=sequence,
        i_code="",
        coords=NPVec3D(*coords),
        occupancy=1.0,
        temp_factor=20.0,
        element=element,
        charge="",
        record_type=record_type,
        backbone_sidechain=backbone_sidechain,
    )


@pytest.fixture
def ligplot_generator(monkeypatch):
    """Create a LigPlot generator with direct and water-bridge interactions."""
    ligand_c = make_atom(
        1, "C1", "LIG", 100, (0, 0, 0), element="C", record_type="HETATM"
    )
    ligand_o = make_atom(
        2, "O1", "LIG", 100, (1.4, 0, 0), record_type="HETATM"
    )
    backbone_n = make_atom(
        3,
        "N",
        "ASN",
        10,
        (-3, 2, 0),
        element="N",
        backbone_sidechain="B",
    )
    backbone_h = make_atom(
        4,
        "H",
        "ASN",
        10,
        (-2, 1, 0),
        element="H",
        backbone_sidechain="B",
    )
    sidechain_o = make_atom(
        5, "OD1", "ASN", 10, (-3, 3, 0), backbone_sidechain="S"
    )
    ser_o = make_atom(6, "OG", "SER", 20, (4, -2, 0), backbone_sidechain="S")

    hydrogen_bond = HydrogenBond(
        _donor=backbone_n,
        hydrogen=backbone_h,
        _acceptor=ligand_o,
        distance=1.8,
        angle=math.radians(165),
        _donor_acceptor_distance=2.7,
        bond_type="N-H...O",
    )
    water_bridge = WaterBridge(
        donor_atom=ligand_o,
        acceptor_atom=sidechain_o,
        bridge_path=[],
        water_residues=["A:HOH:500"],
        bridge_length=1,
        total_distance=4.0,
    )
    sidechain_hbond = HydrogenBond(
        _donor=ligand_o,
        hydrogen=backbone_h,
        _acceptor=ser_o,
        distance=1.9,
        angle=math.radians(160),
        _donor_acceptor_distance=2.8,
        bond_type="O-H...O",
    )

    analyzer = SimpleNamespace(
        parser=SimpleNamespace(
            atoms=[
                ligand_c,
                ligand_o,
                backbone_n,
                backbone_h,
                sidechain_o,
                ser_o,
            ]
        ),
        ligand_interactions=LigandInteraction(
            [hydrogen_bond, water_bridge, sidechain_hbond]
        ),
    )
    generator = LigplotGenerator("LIG", analyzer, residue_id="A:LIG:100")

    def get_test_molecule():
        generator.pdb_ligand_atom_info = {
            0: {"serial": ligand_c.serial, "name": ligand_c.name},
            1: {"serial": ligand_o.serial, "name": ligand_o.name},
        }
        return Chem.MolFromSmiles("CO")

    monkeypatch.setattr(generator, "get_ligand_mol", get_test_molecule)
    return generator


@pytest.mark.unit
class TestLigplotSvg:
    """Test generated LigPlot SVG behavior."""

    def test_interactive_and_static_line_visibility(self, ligplot_generator):
        interactive = ligplot_generator.generate_interactive_svg()
        static = ligplot_generator.generate_static_svg()

        assert ".ligplot-line {" in interactive
        assert "opacity: 0;" in interactive
        assert "opacity: 1 !important;" not in interactive
        assert ".ligplot-line { opacity: 1 !important; }" in static

    def test_water_bridge_renders_cyan_line(self, ligplot_generator):
        svg = ligplot_generator.generate_interactive_svg()

        assert "WaterBridge" in ligplot_generator.applicable_interactions
        assert 'stroke="#00ccff"' in svg

    def test_multi_interaction_atom_has_split_highlight_colors(self, ligplot_generator):
        svg = ligplot_generator.generate_interactive_svg()

        assert "#FF7F7F" in svg
        assert "#00CCFF" in svg

    def test_residue_boxes_include_mixed_and_sidechain_classes(self, ligplot_generator):
        svg = ligplot_generator.generate_interactive_svg()

        assert "ligplot-residue-node ligplot-mixed" in svg
        assert "ligplot-residue-node ligplot-sidechain" in svg
        assert "<rect x=\"-48\" y=\"-21\" width=\"96\" height=\"42\" rx=\"7\"" in svg

    def test_perimeter_boxes_do_not_overlap(self, ligplot_generator):
        svg = ligplot_generator.generate_interactive_svg()
        positions = [
            (float(x), float(y))
            for x, y in re.findall(
                r"ligplot-residue-node[^>]+translate\(([\d.]+) ([\d.]+)\)", svg
            )
        ]

        for index, (x1, y1) in enumerate(positions):
            for x2, y2 in positions[index + 1 :]:
                assert abs(x1 - x2) >= 96 or abs(y1 - y2) >= 42

    def test_classification_aggregation(self, ligplot_generator):
        aggregate = ligplot_generator._aggregate_classification

        assert aggregate({"backbone"}) == "backbone"
        assert aggregate({"sidechain"}) == "sidechain"
        assert aggregate({"backbone", "sidechain"}) == "mixed"
        assert aggregate({"unknown"}) == "unknown"
