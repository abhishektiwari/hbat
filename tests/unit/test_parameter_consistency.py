"""Consistency checks for the CLI, UI, and core parameter definitions."""

import argparse
import inspect

import pytest

from hbat.cli.main import (
    CLI_PARAMETER_OVERRIDES,
    create_parser,
    load_parameters_from_args,
)
from hbat.config.ui_config import PARAMETER_CONFIGS
from hbat.constants.parameters import AnalysisParameters


CLI_NON_PARAMETER_DESTS = {
    "help",
    "version",
    "input",
    "output",
    "json",
    "csv",
    "preset",
    "list_presets",
    "verbose",
    "quiet",
    "summary_only",
    "no_hydrogen_bonds",
    "no_halogen_bonds",
    "no_pi_interactions",
    "no_pi_pi_stacking",
    "no_carbonyl_interactions",
    "no_n_pi_interactions",
}


def analysis_parameter_names():
    """Return the public constructor fields of AnalysisParameters."""
    return set(inspect.signature(AnalysisParameters).parameters) - {"kwargs"}


@pytest.mark.unit
def test_ui_parameters_match_analysis_parameters():
    """Every core parameter has exactly one UI definition with its default."""
    core_names = analysis_parameter_names()
    ui_by_name = {config.name: config for config in PARAMETER_CONFIGS}

    assert set(ui_by_name) == core_names

    defaults = AnalysisParameters()
    assert {name: config.default for name, config in ui_by_name.items()} == {
        name: getattr(defaults, name) for name in core_names
    }


@pytest.mark.unit
def test_cli_parameters_match_analysis_parameters():
    """Every parameter CLI option maps to a core parameter and propagates."""
    parser = create_parser()
    parser_dests = {action.dest for action in parser._actions}

    assert parser_dests - set(CLI_PARAMETER_OVERRIDES) == CLI_NON_PARAMETER_DESTS
    assert set(CLI_PARAMETER_OVERRIDES.values()) == analysis_parameter_names()

    # Exercise every CLI parameter using its parser default. This catches a
    # mapping that exists by name but is not actually wired into the loader.
    args_list = ["dummy.pdb"]
    defaults = AnalysisParameters()
    expected = defaults.to_dict()
    actions = {action.dest: action for action in parser._actions}
    for cli_name, parameter_name in CLI_PARAMETER_OVERRIDES.items():
        action = actions[cli_name]
        option = next(
            option for option in action.option_strings if option.startswith("--")
        )
        if isinstance(action, argparse._StoreTrueAction):
            args_list.append(option)
            expected[parameter_name] = True
        else:
            args_list.extend([option, str(getattr(defaults, parameter_name))])

    params = load_parameters_from_args(parser.parse_args(args_list))
    assert params.to_dict() == expected
