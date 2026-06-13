"""Unit tests for interaction inclusion mode semantics."""

import math

import pytest

from hbat.constants import RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS
from hbat.constants.parameters import AnalysisParameters
from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
from hbat.core.np_vector import NPVec3D
from hbat.core.structure import Atom, Bond, Residue


def make_atom(
    serial: int,
    name: str,
    element: str,
    coords: tuple[float, float, float],
    res_name: str = "PHE",
    res_seq: int = 1,
) -> Atom:
    """Create a minimal atom for detector tests."""
    return Atom(
        serial=serial,
        name=name,
        alt_loc="",
        res_name=res_name,
        chain_id="A",
        res_seq=res_seq,
        i_code="",
        coords=NPVec3D(*coords),
        occupancy=1.0,
        temp_factor=20.0,
        element=element,
        charge="",
        record_type="ATOM",
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    ("mode", "same_residue", "expected"),
    [
        ("inter", True, True),
        ("inter", False, False),
        ("all", True, False),
        ("all", False, False),
    ],
)
def test_same_residue_skip_policy(mode, same_residue, expected):
    """The shared policy skips same-residue candidates only in inter mode."""
    analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters(analysis_mode=mode))
    residue1 = ("A", 1, "PHE")
    residue2 = residue1 if same_residue else ("A", 2, "TYR")

    assert analyzer._should_skip_same_residue(residue1, residue2) is expected


@pytest.mark.unit
@pytest.mark.parametrize(("mode", "expected_count"), [("inter", 0), ("all", 1)])
def test_hydrogen_bond_detector_respects_analysis_mode(
    monkeypatch, mode, expected_count
):
    """Hydrogen bond detection applies the mode policy to same-residue pairs."""
    analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters(analysis_mode=mode))
    donor = make_atom(1, "N", "N", (0.0, 0.0, 0.0))
    hydrogen = make_atom(2, "H", "H", (1.0, 0.0, 0.0))
    acceptor = make_atom(3, "O", "O", (2.5, 0.0, 0.0))
    analyzer.parser.atoms = [donor, hydrogen, acceptor]
    analyzer._prepare_vectorized_data()
    monkeypatch.setattr(
        analyzer,
        "_get_hydrogen_bond_donors",
        lambda: [(donor, hydrogen, 0, 1)],
    )

    analyzer._find_hydrogen_bonds_vectorized()

    assert len(analyzer.hydrogen_bonds) == expected_count


@pytest.mark.unit
@pytest.mark.parametrize(("mode", "expected_count"), [("inter", 0), ("all", 1)])
def test_halogen_bond_detector_respects_analysis_mode(mode, expected_count):
    """Halogen bond detection applies the mode policy to same-residue pairs."""
    analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters(analysis_mode=mode))
    carbon = make_atom(1, "C", "C", (0.0, 0.0, 0.0))
    halogen = make_atom(2, "CL", "CL", (1.0, 0.0, 0.0))
    acceptor = make_atom(3, "O", "O", (3.0, 0.0, 0.0))
    analyzer.parser.atoms = [carbon, halogen, acceptor]
    analyzer.parser.bonds = [Bond(carbon.serial, halogen.serial)]
    analyzer._prepare_vectorized_data()

    analyzer._find_halogen_bonds_vectorized()

    assert len(analyzer.halogen_bonds) == expected_count


@pytest.mark.unit
@pytest.mark.parametrize(("mode", "expected_count"), [("inter", 0), ("all", 1)])
def test_pi_detector_respects_analysis_mode(monkeypatch, mode, expected_count):
    """Pi detection applies the mode policy to same-residue candidates."""
    analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters(analysis_mode=mode))
    ring_names = RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS["PHE"]
    ring_coords = [
        (1.0, 0.0, 0.0),
        (0.5, 0.866, 0.0),
        (-0.5, 0.866, 0.0),
        (-1.0, 0.0, 0.0),
        (-0.5, -0.866, 0.0),
        (0.5, -0.866, 0.0),
    ]
    ring_atoms = [
        make_atom(serial, name, "C", coords)
        for serial, (name, coords) in enumerate(
            zip(ring_names[:6], ring_coords), start=1
        )
    ]
    donor = make_atom(20, "N", "N", (0.0, 0.0, 4.0))
    hydrogen = make_atom(21, "H", "H", (0.0, 0.0, 3.0))
    residue = Residue("PHE", "A", 1, "", ring_atoms + [donor, hydrogen])
    analyzer.parser.atoms = residue.atoms
    analyzer.parser.residues = {"A:PHE:1": residue}
    analyzer._prepare_vectorized_data()
    monkeypatch.setattr(
        analyzer, "_get_pi_interaction_pairs", lambda: [(donor, hydrogen)]
    )

    analyzer._find_pi_interactions_vectorized()

    assert len(analyzer.pi_interactions) == expected_count


@pytest.mark.unit
@pytest.mark.parametrize(
    ("legacy_mode", "canonical_mode", "skips_same_residue"),
    [("local", "inter", True), ("complete", "all", False)],
)
def test_canonical_modes_preserve_legacy_semantics(
    legacy_mode, canonical_mode, skips_same_residue
):
    """Record the intended semantic mapping for the breaking mode rename."""
    analyzer = NPMolecularInteractionAnalyzer(
        AnalysisParameters(analysis_mode=canonical_mode)
    )
    residue = ("A", 1, "PHE")

    assert analyzer._should_skip_same_residue(residue, residue) is skips_same_residue, (
        f"{canonical_mode} must preserve the previous {legacy_mode} behavior"
    )


@pytest.mark.unit
@pytest.mark.parametrize(("mode", "expected_count"), [("inter", 0), ("all", 1)])
def test_pi_pi_detector_respects_analysis_mode(mode, expected_count):
    """Pi-pi detection applies the mode policy to same-residue ring pairs."""
    analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters(analysis_mode=mode))
    ring_names = RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS["PHE"]
    ring_coords = [
        (1.0, 0.0, 0.0),
        (0.5, 0.866, 0.0),
        (-0.5, 0.866, 0.0),
        (-1.0, 0.0, 0.0),
        (-0.5, -0.866, 0.0),
        (0.5, -0.866, 0.0),
    ]
    ring_atoms = [
        make_atom(serial, name, "C", coords)
        for serial, (name, coords) in enumerate(
            zip(ring_names[:6], ring_coords), start=1
        )
    ]
    residue = Residue("PHE", "A", 1, "", ring_atoms)
    # Duplicate values exercise the same-residue detector branch directly.
    analyzer.parser.residues = {"first": residue, "second": residue}

    analyzer._find_pi_pi_interactions_vectorized()

    assert len(analyzer.pi_pi_interactions) == expected_count


@pytest.mark.unit
@pytest.mark.parametrize(("mode", "expected_count"), [("inter", 0), ("all", 2)])
def test_carbonyl_detector_respects_analysis_mode(monkeypatch, mode, expected_count):
    """Carbonyl detection includes same-residue candidates only in all mode."""
    analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters(analysis_mode=mode))
    analyzer.parser.atoms = [
        make_atom(1, "C", "C", (0.0, 0.0, 0.0), "ASP"),
        make_atom(2, "O", "O", (1.2, 0.0, 0.0), "ASP"),
        make_atom(3, "CG", "C", (4.2, 0.0, 0.0), "ASP"),
        make_atom(4, "OD1", "O", (3.0, 0.0, 0.0), "ASP"),
    ]
    monkeypatch.setattr(
        analyzer,
        "_identify_carbonyl_groups",
        lambda: [(0, 1, True, "ASP1"), (2, 3, False, "ASP1")],
    )
    monkeypatch.setattr(
        "hbat.core.np_analyzer.batch_angle_between",
        lambda *_args: math.radians(107.0),
    )

    analyzer._find_carbonyl_interactions_vectorized()

    assert len(analyzer.carbonyl_interactions) == expected_count


@pytest.mark.unit
@pytest.mark.parametrize(("mode", "expected_count"), [("inter", 0), ("all", 1)])
def test_n_pi_detector_respects_analysis_mode(mode, expected_count):
    """n-pi detection includes same-residue candidates only in all mode."""
    analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters(analysis_mode=mode))
    ring_names = RING_ATOMS_FOR_RESIDUES_WITH_AROMATIC_RINGS["PHE"]
    ring_coords = [
        (1.0, 0.0, 0.0),
        (0.5, 0.866, 0.0),
        (-0.5, 0.866, 0.0),
        (-1.0, 0.0, 0.0),
        (-0.5, -0.866, 0.0),
        (0.5, -0.866, 0.0),
    ]
    ring_atoms = [
        make_atom(serial, name, "C", coords)
        for serial, (name, coords) in enumerate(
            zip(ring_names[:6], ring_coords), start=1
        )
    ]
    donor = make_atom(20, "O", "O", (3.2, 0.0, 0.0))
    residue = Residue("PHE", "A", 1, "", ring_atoms + [donor])
    analyzer.parser.atoms = residue.atoms
    analyzer.parser.residues = {"A:PHE:1": residue}

    analyzer._find_n_pi_interactions_vectorized()

    assert len(analyzer.n_pi_interactions) == expected_count
