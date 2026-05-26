"""Configuration system for HBAT analysis results.

This module provides centralized configuration for interaction types,
parameters, and data extraction patterns used across GUI, Server, and
Export modules.
"""

from hbat.config.ui_config import (
    # Dataclasses
    ColumnConfig,
    InteractionConfig,
    ParameterConfig,
    # Interaction Configs
    INTERACTION_CONFIGS,
    HYDROGEN_BONDS_CONFIG,
    WATER_BRIDGES_CONFIG,
    HALOGEN_BONDS_CONFIG,
    PI_INTERACTIONS_CONFIG,
    PI_PI_STACKING_CONFIG,
    CARBONYL_INTERACTIONS_CONFIG,
    N_PI_INTERACTIONS_CONFIG,
    COOPERATIVITY_CHAINS_CONFIG,
    # Parameter Configs
    PARAMETER_CONFIGS,
    PARAMETER_LOOKUP,
    # Helper Functions
    get_interaction_config,
    get_parameter_config,
    extract_interaction_data,
    get_parameters_by_category,
)

from hbat.config.parameter_controller import ParameterController

__all__ = [
    # Dataclasses
    "ColumnConfig",
    "InteractionConfig",
    "ParameterConfig",
    # Interaction Configs
    "INTERACTION_CONFIGS",
    "HYDROGEN_BONDS_CONFIG",
    "WATER_BRIDGES_CONFIG",
    "HALOGEN_BONDS_CONFIG",
    "PI_INTERACTIONS_CONFIG",
    "PI_PI_STACKING_CONFIG",
    "CARBONYL_INTERACTIONS_CONFIG",
    "N_PI_INTERACTIONS_CONFIG",
    "COOPERATIVITY_CHAINS_CONFIG",
    # Parameter Configs
    "PARAMETER_CONFIGS",
    "PARAMETER_LOOKUP",
    # Helper Functions
    "get_interaction_config",
    "get_parameter_config",
    "extract_interaction_data",
    "get_parameters_by_category",
    # Controller
    "ParameterController",
]
