"""Shared serialization and validation for HBAT parameter presets."""

from typing import Any, Dict, Mapping, Optional

from ..constants.parameters import AnalysisParameters


PRESET_FORMAT_VERSION = "1.0"

# The public .hbat names are intentionally kept stable. This is the single
# mapping used by the GUI writer, GUI reader, and CLI reader.
PRESET_PARAMETER_FIELDS: Dict[str, Dict[str, str]] = {
    "hydrogen_bonds": {
        "h_a_distance_cutoff": "hb_distance_cutoff",
        "dha_angle_cutoff": "hb_angle_cutoff",
        "d_a_distance_cutoff": "hb_donor_acceptor_cutoff",
    },
    "weak_hydrogen_bonds": {
        "h_a_distance_cutoff": "whb_distance_cutoff",
        "dha_angle_cutoff": "whb_angle_cutoff",
        "d_a_distance_cutoff": "whb_donor_acceptor_cutoff",
    },
    "halogen_bonds": {
        "x_a_distance_cutoff": "xb_distance_cutoff",
        "dxa_angle_cutoff": "xb_angle_cutoff",
    },
    "pi_interactions": {
        "h_pi_distance_cutoff": "pi_distance_cutoff",
        "dh_pi_angle_cutoff": "pi_angle_cutoff",
        "ccl_pi_distance_cutoff": "pi_ccl_distance_cutoff",
        "ccl_pi_angle_cutoff": "pi_ccl_angle_cutoff",
        "cbr_pi_distance_cutoff": "pi_cbr_distance_cutoff",
        "cbr_pi_angle_cutoff": "pi_cbr_angle_cutoff",
        "ci_pi_distance_cutoff": "pi_ci_distance_cutoff",
        "ci_pi_angle_cutoff": "pi_ci_angle_cutoff",
        "ch_pi_distance_cutoff": "pi_ch_distance_cutoff",
        "ch_pi_angle_cutoff": "pi_ch_angle_cutoff",
        "nh_pi_distance_cutoff": "pi_nh_distance_cutoff",
        "nh_pi_angle_cutoff": "pi_nh_angle_cutoff",
        "oh_pi_distance_cutoff": "pi_oh_distance_cutoff",
        "oh_pi_angle_cutoff": "pi_oh_angle_cutoff",
        "sh_pi_distance_cutoff": "pi_sh_distance_cutoff",
        "sh_pi_angle_cutoff": "pi_sh_angle_cutoff",
    },
    "pi_pi_stacking": {
        "distance_cutoff": "pi_pi_distance_cutoff",
        "parallel_angle_cutoff": "pi_pi_parallel_angle_cutoff",
        "tshaped_angle_min": "pi_pi_tshaped_angle_min",
        "tshaped_angle_max": "pi_pi_tshaped_angle_max",
        "offset_cutoff": "pi_pi_offset_cutoff",
    },
    "carbonyl_interactions": {
        "distance_cutoff": "carbonyl_distance_cutoff",
        "angle_min": "carbonyl_angle_min",
        "angle_max": "carbonyl_angle_max",
    },
    "n_pi_interactions": {
        "distance_cutoff": "n_pi_distance_cutoff",
        "sulfur_distance_cutoff": "n_pi_sulfur_distance_cutoff",
        "angle_min": "n_pi_angle_min",
        "angle_max": "n_pi_angle_max",
    },
    "general": {
        "covalent_cutoff_factor": "covalent_cutoff_factor",
        "analysis_mode": "analysis_mode",
    },
    "pdb_fixing": {
        "enabled": "fix_pdb_enabled",
        "method": "fix_pdb_method",
        "add_hydrogens": "fix_pdb_add_hydrogens",
        "add_heavy_atoms": "fix_pdb_add_heavy_atoms",
        "replace_nonstandard": "fix_pdb_replace_nonstandard",
        "remove_heterogens": "fix_pdb_remove_heterogens",
        "keep_water": "fix_pdb_keep_water",
    },
}

# Older built-in/user presets used core-style names for these sections. Keep
# them readable while emitting only the canonical names above.
PRESET_PARAMETER_ALIASES: Dict[str, Dict[str, str]] = {
    "pi_pi_stacking": {
        "pi_pi_distance_cutoff": "distance_cutoff",
        "pi_pi_parallel_angle_cutoff": "parallel_angle_cutoff",
        "pi_pi_tshaped_angle_min": "tshaped_angle_min",
        "pi_pi_tshaped_angle_max": "tshaped_angle_max",
        "pi_pi_offset_cutoff": "offset_cutoff",
    },
    "carbonyl_interactions": {
        "carbonyl_distance_cutoff": "distance_cutoff",
        "carbonyl_angle_min": "angle_min",
        "carbonyl_angle_max": "angle_max",
    },
    "n_pi_interactions": {
        "n_pi_distance_cutoff": "distance_cutoff",
        "n_pi_sulfur_distance_cutoff": "sulfur_distance_cutoff",
        "n_pi_angle_min": "angle_min",
        "n_pi_angle_max": "angle_max",
    },
}


def parameters_to_preset(params: AnalysisParameters) -> Dict[str, Any]:
    """Serialize analysis parameters into the public nested preset schema."""
    values = params.to_dict()
    return {
        "parameters": {
            section: {
                preset_name: values[parameter_name]
                for preset_name, parameter_name in fields.items()
            }
            for section, fields in PRESET_PARAMETER_FIELDS.items()
        }
    }


def preset_to_parameters(
    data: Mapping[str, Any],
    base_params: Optional[AnalysisParameters] = None,
) -> AnalysisParameters:
    """Deserialize a nested preset into validated analysis parameters.

    Missing sections or fields retain values from ``base_params`` when given,
    otherwise they use ``AnalysisParameters`` defaults.
    """
    if not isinstance(data, Mapping):
        raise ValueError("Invalid preset format: expected an object")

    format_version = data.get("format_version")
    if format_version is not None and str(format_version) != PRESET_FORMAT_VERSION:
        raise ValueError(
            f"Unsupported preset format version: {format_version}. "
            f"Expected {PRESET_FORMAT_VERSION}."
        )

    raw_parameters = data.get("parameters")
    if not isinstance(raw_parameters, Mapping):
        raise ValueError("Invalid preset format: missing 'parameters' section")

    unknown_sections = set(raw_parameters) - set(PRESET_PARAMETER_FIELDS)
    if unknown_sections:
        raise ValueError(
            "Invalid preset parameters: unknown sections: "
            + ", ".join(sorted(unknown_sections))
        )

    values = (base_params or AnalysisParameters()).to_dict()
    for section, fields in PRESET_PARAMETER_FIELDS.items():
        section_values = raw_parameters.get(section, {})
        if not isinstance(section_values, Mapping):
            raise ValueError(
                f"Invalid preset parameters: '{section}' must be an object"
            )

        normalized_section = dict(section_values)
        aliases = PRESET_PARAMETER_ALIASES.get(section, {})
        for legacy_name, canonical_name in aliases.items():
            if legacy_name in normalized_section:
                if canonical_name in normalized_section:
                    raise ValueError(
                        f"Invalid preset parameters: both '{legacy_name}' and "
                        f"'{canonical_name}' are set in '{section}'"
                    )
                normalized_section[canonical_name] = normalized_section.pop(legacy_name)

        unknown_fields = set(normalized_section) - set(fields)
        if unknown_fields:
            raise ValueError(
                f"Invalid preset parameters: unknown fields in '{section}': "
                + ", ".join(sorted(unknown_fields))
            )

        for preset_name, parameter_name in fields.items():
            if preset_name in normalized_section:
                values[parameter_name] = normalized_section[preset_name]

    params = AnalysisParameters.from_dict(values)
    params.validate_or_raise("preset parameters")
    return params
