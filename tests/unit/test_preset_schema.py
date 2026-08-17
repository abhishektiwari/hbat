"""Tests for the shared preset schema and parameter persistence paths."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from hbat.config.parameter_controller import ParameterController
from hbat.config.preset_schema import (
    PRESET_PARAMETER_FIELDS,
    parameters_to_preset,
    preset_to_parameters,
)
from hbat.constants.parameters import AnalysisParameters
from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
from hbat.gui.main_window import MainWindow
from hbat.gui.preset_manager_dialog import PresetManagerDialog


def customized_parameters() -> AnalysisParameters:
    """Return parameters with representative values from every preset section."""
    return AnalysisParameters(
        hb_distance_cutoff=3.1,
        hb_angle_cutoff=130.0,
        hb_donor_acceptor_cutoff=3.8,
        whb_distance_cutoff=3.9,
        whb_angle_cutoff=145.0,
        whb_donor_acceptor_cutoff=3.7,
        xb_distance_cutoff=3.8,
        xb_angle_cutoff=155.0,
        pi_distance_cutoff=3.9,
        pi_angle_cutoff=115.0,
        pi_ccl_distance_cutoff=3.6,
        pi_ccl_angle_cutoff=150.0,
        pi_cbr_distance_cutoff=3.7,
        pi_cbr_angle_cutoff=155.0,
        pi_ci_distance_cutoff=3.8,
        pi_ci_angle_cutoff=165.0,
        pi_ch_distance_cutoff=3.6,
        pi_ch_angle_cutoff=120.0,
        pi_nh_distance_cutoff=3.3,
        pi_nh_angle_cutoff=125.0,
        pi_oh_distance_cutoff=3.1,
        pi_oh_angle_cutoff=125.0,
        pi_sh_distance_cutoff=3.9,
        pi_sh_angle_cutoff=110.0,
        pi_pi_distance_cutoff=4.0,
        pi_pi_parallel_angle_cutoff=25.0,
        pi_pi_tshaped_angle_min=65.0,
        pi_pi_tshaped_angle_max=88.0,
        pi_pi_offset_cutoff=1.8,
        carbonyl_distance_cutoff=3.1,
        carbonyl_angle_min=100.0,
        carbonyl_angle_max=120.0,
        n_pi_distance_cutoff=3.5,
        n_pi_sulfur_distance_cutoff=3.9,
        n_pi_angle_min=5.0,
        n_pi_angle_max=40.0,
        covalent_cutoff_factor=0.9,
        analysis_mode="all",
        fix_pdb_enabled=False,
        fix_pdb_method="pdbfixer",
        fix_pdb_add_hydrogens=False,
        fix_pdb_add_heavy_atoms=True,
        fix_pdb_replace_nonstandard=True,
        fix_pdb_remove_heterogens=True,
        fix_pdb_keep_water=False,
    )


@pytest.mark.unit
def test_preset_schema_covers_every_analysis_parameter():
    """The shared preset mapping must cover the complete core parameter set."""
    mapped_fields = {
        parameter_name
        for fields in PRESET_PARAMETER_FIELDS.values()
        for parameter_name in fields.values()
    }
    assert mapped_fields == set(AnalysisParameters().to_dict())


@pytest.mark.unit
def test_parameters_round_trip_through_preset_schema():
    """All parameter groups survive serialization and deserialization."""
    original = customized_parameters()

    serialized = parameters_to_preset(original)
    loaded = preset_to_parameters(serialized)

    assert loaded.to_dict() == original.to_dict()


@pytest.mark.unit
def test_legacy_example_presets_remain_loadable():
    """Built-in presets using legacy prefixed keys remain compatible."""
    for preset_path in sorted(Path("example_presets").glob("*.hbat")):
        data = json.loads(preset_path.read_text())
        params = preset_to_parameters(data)
        assert isinstance(params, AnalysisParameters), preset_path.name


@pytest.mark.unit
def test_partial_preset_preserves_current_parameters():
    """GUI-style partial presets update only the supplied fields."""
    original = customized_parameters()
    data = {"parameters": {"general": {"analysis_mode": "inter"}}}

    loaded = preset_to_parameters(data, base_params=original)

    assert loaded.analysis_mode == "inter"
    assert loaded.hb_distance_cutoff == original.hb_distance_cutoff
    assert loaded.pi_pi_distance_cutoff == original.pi_pi_distance_cutoff
    assert loaded.fix_pdb_method == original.fix_pdb_method


@pytest.mark.unit
def test_geometry_dialog_applies_all_preset_sections():
    """The geometry dialog applies every shared-schema parameter group."""
    dialog = GeometryCutoffsDialog.__new__(GeometryCutoffsDialog)
    dialog._vars = {}
    dialog._param_values = {}
    original = AnalysisParameters()

    dialog._apply_preset_data(parameters_to_preset(customized_parameters()))

    assert dialog.get_parameters().to_dict() == customized_parameters().to_dict()
    assert dialog.get_parameters() != original


@pytest.mark.unit
def test_main_window_applies_all_preset_sections():
    """The desktop session loader uses the complete shared schema."""
    window = MainWindow.__new__(MainWindow)
    window.session_parameters = AnalysisParameters()

    applied = window._apply_preset_to_session(
        parameters_to_preset(customized_parameters())
    )

    assert applied is True
    assert window.session_parameters.to_dict() == customized_parameters().to_dict()


@pytest.mark.unit
def test_preset_manager_uses_shared_schema():
    """The GUI writer emits the same schema consumed by the readers."""
    dialog = PresetManagerDialog.__new__(PresetManagerDialog)
    dialog.current_params = customized_parameters()
    dialog.description_var = SimpleNamespace(get=lambda: "test preset")

    data = dialog._create_preset_data()

    assert preset_to_parameters(data).to_dict() == dialog.current_params.to_dict()


@pytest.mark.unit
def test_parameter_controller_to_dict_uses_analysis_parameters_serializer():
    """ParameterController serialization returns the core parameter dictionary."""
    params = customized_parameters()

    assert ParameterController(params).to_dict() == params.to_dict()


@pytest.mark.unit
def test_preset_schema_rejects_unknown_fields():
    """Typos in preset fields fail instead of being silently ignored."""
    data = {"parameters": {"hydrogen_bonds": {"h_a_distance_cuttoff": 3.0}}}

    with pytest.raises(ValueError, match="unknown fields"):
        preset_to_parameters(data)
