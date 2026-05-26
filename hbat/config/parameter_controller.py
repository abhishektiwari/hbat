"""Parameter controller for HBAT analysis parameters.

This module provides a framework-agnostic controller for managing
analysis parameters, independent of any UI framework.
"""

from dataclasses import asdict, fields
from typing import Any, Dict, List, Optional

from ..constants.parameters import AnalysisParameters
from .ui_config import PARAMETER_CONFIGS, PARAMETER_LOOKUP, ParameterConfig


class ParameterController:
    """Manage analysis parameters independent of UI framework.

    Provides get/set operations, validation, persistence, and
    parameter grouping by category.
    """

    def __init__(self, params: Optional[AnalysisParameters] = None):
        """Initialize parameter controller.

        Args:
            params: Initial AnalysisParameters instance. If None, creates default.
        """
        self.params = params or AnalysisParameters()

    def get_parameter(self, name: str) -> Any:
        """Get a parameter value by name.

        Args:
            name: Parameter name (e.g., 'hb_distance_cutoff')

        Returns:
            Parameter value

        Raises:
            AttributeError: If parameter doesn't exist
        """
        return getattr(self.params, name)

    def set_parameter(self, name: str, value: Any) -> None:
        """Set a parameter value by name.

        Args:
            name: Parameter name
            value: New parameter value

        Raises:
            AttributeError: If parameter doesn't exist
            ValueError: If value fails validation
        """
        if not self.validate_parameter(name, value):
            config = PARAMETER_LOOKUP.get(name)
            if config:
                raise ValueError(
                    f"{name} must be between {config.min_val} and {config.max_val}"
                )
            raise ValueError(f"Invalid value for {name}")
        setattr(self.params, name, value)

    def validate_parameter(self, name: str, value: Any) -> bool:
        """Validate a parameter value.

        Args:
            name: Parameter name
            value: Value to validate

        Returns:
            True if valid, False otherwise
        """
        config = PARAMETER_LOOKUP.get(name)
        if not config:
            return False

        # Check type
        if config.param_type == "float":
            try:
                value = float(value)
            except (TypeError, ValueError):
                return False
        elif config.param_type == "int":
            try:
                value = int(value)
            except (TypeError, ValueError):
                return False
        elif config.param_type == "bool":
            if not isinstance(value, bool):
                return False
        elif config.param_type == "str":
            if not isinstance(value, str):
                return False

        # Check range (if applicable)
        if config.min_val is not None and value < config.min_val:
            return False
        if config.max_val is not None and value > config.max_val:
            return False

        return True

    def reset_to_defaults(self) -> None:
        """Reset all parameters to their default values."""
        self.params = AnalysisParameters()

    def to_dict(self) -> Dict[str, Any]:
        """Convert parameters to dictionary.

        Returns:
            Dictionary mapping parameter names to values
        """
        return asdict(self.params)

    def from_dict(self, params_dict: Dict[str, Any]) -> None:
        """Update parameters from dictionary.

        Args:
            params_dict: Dictionary mapping parameter names to values

        Raises:
            ValueError: If any value fails validation
        """
        for name, value in params_dict.items():
            if hasattr(self.params, name):
                if not self.validate_parameter(name, value):
                    raise ValueError(f"Invalid value for {name}: {value}")
                setattr(self.params, name, value)

    def get_parameters_by_category(self) -> Dict[str, Dict[str, ParameterConfig]]:
        """Get all parameters grouped by category.

        Returns:
            Dictionary mapping categories to parameter configs.
            Example structure::

                {
                    "Hydrogen Bonds": {
                        "hb_distance_cutoff": ParameterConfig(...),
                        ...
                    },
                    ...
                }
        """
        grouped = {}
        for config in PARAMETER_CONFIGS:
            category = config.category
            if category not in grouped:
                grouped[category] = {}
            grouped[category][config.name] = config
        return grouped

    def get_parameters_list(self) -> List[ParameterConfig]:
        """Get all parameter configurations.

        Returns:
            List of ParameterConfig objects
        """
        return PARAMETER_CONFIGS

    def get_category_parameters(self, category: str) -> Dict[str, ParameterConfig]:
        """Get all parameters for a specific category.

        Args:
            category: Category name

        Returns:
            Dictionary mapping parameter names to configs for the category
        """
        by_category = self.get_parameters_by_category()
        return by_category.get(category, {})
