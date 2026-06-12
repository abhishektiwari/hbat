"""Unit tests for preset validation boundaries."""

import json

import pytest

from hbat.cli.main import load_preset_file
from hbat.constants.parameters import AnalysisParameters
from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
from hbat.gui.main_window import MainWindow


def preset_with_mode(mode: str) -> dict:
    """Create a minimal preset using the requested analysis mode."""
    return {"parameters": {"general": {"analysis_mode": mode}}}


@pytest.mark.unit
@pytest.mark.parametrize("legacy_mode", ["local", "complete"])
def test_cli_rejects_legacy_preset_modes(tmp_path, capsys, legacy_mode):
    """CLI preset loading rejects removed analysis modes with a useful error."""
    preset_path = tmp_path / "legacy.hbat"
    preset_path.write_text(json.dumps(preset_with_mode(legacy_mode)))

    with pytest.raises(SystemExit):
        load_preset_file(str(preset_path))

    error = capsys.readouterr().err
    assert "Invalid preset parameters" in error
    assert "Analysis mode must be one of: inter, all" in error


@pytest.mark.unit
@pytest.mark.parametrize("legacy_mode", ["local", "complete"])
def test_main_window_rejects_legacy_preset_modes(monkeypatch, legacy_mode):
    """Desktop preset loading rejects old modes without replacing the session."""
    window = MainWindow.__new__(MainWindow)
    original_params = AnalysisParameters(analysis_mode="all")
    window.session_parameters = original_params
    errors = []
    monkeypatch.setattr(
        "hbat.gui.main_window.messagebox.showerror",
        lambda title, message: errors.append((title, message)),
    )

    applied = window._apply_preset_to_session(preset_with_mode(legacy_mode))

    assert applied is False
    assert window.session_parameters is original_params
    assert errors[0][0] == "Invalid Preset"
    assert "Analysis mode must be one of: inter, all" in errors[0][1]


@pytest.mark.unit
def test_main_window_applies_valid_preset_mode():
    """Desktop preset loading applies valid canonical analysis modes."""
    window = MainWindow.__new__(MainWindow)
    window.session_parameters = AnalysisParameters(analysis_mode="inter")

    applied = window._apply_preset_to_session(preset_with_mode("all"))

    assert applied is True
    assert window.session_parameters.analysis_mode == "all"


@pytest.mark.unit
@pytest.mark.parametrize("legacy_mode", ["local", "complete"])
def test_geometry_dialog_rejects_legacy_preset_modes(legacy_mode):
    """Geometry preset loading rejects old modes and restores prior values."""
    dialog = GeometryCutoffsDialog.__new__(GeometryCutoffsDialog)
    dialog._vars = {}
    dialog._param_values = {"analysis_mode": "all"}

    with pytest.raises(ValueError, match="Analysis mode must be one of: inter, all"):
        dialog._apply_preset_data(preset_with_mode(legacy_mode))

    assert dialog._param_values["analysis_mode"] == "all"
