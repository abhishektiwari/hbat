"""
Parameter configuration component for HBAT web interface.

This module provides UI controls for configuring analysis parameters.
"""

from nicegui import ui

from ...constants.parameters import ParametersDefault
from ...core.analysis import AnalysisParameters


class ParameterPanel:
    """Parameter configuration panel component."""

    def __init__(self, param_drawer, param_drawer_content):
        """Initialize the parameter panel.

        :param param_drawer: Right drawer UI element for parameter editing
        :param param_drawer_content: Container for drawer content
        """
        self.param_drawer = param_drawer
        self.drawer_content = param_drawer_content

        # PDB fixing parameters
        self.fix_pdb_enabled = ParametersDefault.FIX_PDB_ENABLED
        self.fix_pdb_method = ParametersDefault.FIX_PDB_METHOD
        self.fix_pdb_add_hydrogens = ParametersDefault.FIX_PDB_ADD_HYDROGENS
        self.fix_pdb_add_heavy_atoms = ParametersDefault.FIX_PDB_ADD_HEAVY_ATOMS
        self.fix_pdb_replace_nonstandard = ParametersDefault.FIX_PDB_REPLACE_NONSTANDARD
        self.fix_pdb_remove_heterogens = ParametersDefault.FIX_PDB_REMOVE_HETEROGENS
        self.fix_pdb_keep_water = ParametersDefault.FIX_PDB_KEEP_WATER

        # Hydrogen bond parameters
        self.hb_distance = ParametersDefault.HB_DISTANCE_CUTOFF
        self.hb_angle = ParametersDefault.HB_ANGLE_CUTOFF
        self.hb_da_distance = ParametersDefault.HB_DA_DISTANCE

        # Weak hydrogen bond parameters
        self.whb_distance = ParametersDefault.WHB_DISTANCE_CUTOFF
        self.whb_angle = ParametersDefault.WHB_ANGLE_CUTOFF
        self.whb_da_distance = ParametersDefault.WHB_DA_DISTANCE

        # Halogen bond parameters
        self.xb_distance = ParametersDefault.XB_DISTANCE_CUTOFF
        self.xb_angle = ParametersDefault.XB_ANGLE_CUTOFF

        # π interaction subtype parameters
        self.pi_ccl_distance = ParametersDefault.PI_CCL_DISTANCE_CUTOFF
        self.pi_ccl_angle = ParametersDefault.PI_CCL_ANGLE_CUTOFF
        self.pi_cbr_distance = ParametersDefault.PI_CBR_DISTANCE_CUTOFF
        self.pi_cbr_angle = ParametersDefault.PI_CBR_ANGLE_CUTOFF
        self.pi_ci_distance = ParametersDefault.PI_CI_DISTANCE_CUTOFF
        self.pi_ci_angle = ParametersDefault.PI_CI_ANGLE_CUTOFF
        self.pi_ch_distance = ParametersDefault.PI_CH_DISTANCE_CUTOFF
        self.pi_ch_angle = ParametersDefault.PI_CH_ANGLE_CUTOFF
        self.pi_nh_distance = ParametersDefault.PI_NH_DISTANCE_CUTOFF
        self.pi_nh_angle = ParametersDefault.PI_NH_ANGLE_CUTOFF
        self.pi_oh_distance = ParametersDefault.PI_OH_DISTANCE_CUTOFF
        self.pi_oh_angle = ParametersDefault.PI_OH_ANGLE_CUTOFF
        self.pi_sh_distance = ParametersDefault.PI_SH_DISTANCE_CUTOFF
        self.pi_sh_angle = ParametersDefault.PI_SH_ANGLE_CUTOFF

        # π-π stacking parameters
        self.pi_pi_distance = ParametersDefault.PI_PI_DISTANCE_CUTOFF
        self.pi_pi_parallel_angle = ParametersDefault.PI_PI_PARALLEL_ANGLE_CUTOFF
        self.pi_pi_tshaped_angle_min = ParametersDefault.PI_PI_TSHAPED_ANGLE_MIN
        self.pi_pi_tshaped_angle_max = ParametersDefault.PI_PI_TSHAPED_ANGLE_MAX
        self.pi_pi_offset = ParametersDefault.PI_PI_OFFSET_CUTOFF

        # Carbonyl interaction parameters
        self.carbonyl_distance = ParametersDefault.CARBONYL_DISTANCE_CUTOFF
        self.carbonyl_angle_min = ParametersDefault.CARBONYL_ANGLE_MIN
        self.carbonyl_angle_max = ParametersDefault.CARBONYL_ANGLE_MAX

        # n→π* interaction parameters
        self.n_pi_distance = ParametersDefault.N_PI_DISTANCE_CUTOFF
        self.n_pi_sulfur_distance = ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF
        self.n_pi_angle_min = ParametersDefault.N_PI_ANGLE_MIN
        self.n_pi_angle_max = ParametersDefault.N_PI_ANGLE_MAX

        # General parameters
        self.covalent_cutoff_factor = ParametersDefault.COVALENT_CUTOFF_FACTOR
        self.analysis_mode = ParametersDefault.ANALYSIS_MODE

    def _should_enable_keep_water(self) -> bool:
        """Check if Keep Water switch should be enabled.

        :returns: True if both PDBFixer is selected AND Remove Heterogens is checked
        :rtype: bool
        """
        return self.fix_pdb_method == "pdbfixer" and self.fix_pdb_remove_heterogens

    def _update_keep_water_state(self):
        """Update the enabled state of Keep Water switch based on method and remove_heterogens."""
        if hasattr(self, "keep_water_switch"):
            should_enable = self._should_enable_keep_water()
            self.keep_water_switch.enabled = should_enable

    def _reset_pdb_fixing_options(self):
        """Reset PDB fixing options to defaults."""
        self.fix_pdb_add_hydrogens = ParametersDefault.FIX_PDB_ADD_HYDROGENS
        self.fix_pdb_add_heavy_atoms = ParametersDefault.FIX_PDB_ADD_HEAVY_ATOMS
        self.fix_pdb_replace_nonstandard = ParametersDefault.FIX_PDB_REPLACE_NONSTANDARD
        self.fix_pdb_remove_heterogens = ParametersDefault.FIX_PDB_REMOVE_HETEROGENS
        self.fix_pdb_keep_water = ParametersDefault.FIX_PDB_KEEP_WATER

    def create_ui(self):
        """Create the parameter configuration UI."""
        # Parameter summary cards
        with ui.row().classes("w-full gap-4 mt-5"):
            # PDB Fixing card
            with ui.card().classes("flex-1"):
                with ui.card_section():
                    with ui.row().classes("items-center justify-between"):
                        ui.label("PDB Fixing").classes("text-h6")
                        ui.button(
                            icon="edit", on_click=lambda: self._edit_pdb_fixing()
                        ).props("flat round dense color=primary")
                    ui.label().classes("text-caption text-grey").bind_text_from(
                        self, "fix_pdb_method", backward=lambda v: f"Method: {v}"
                    )
                    ui.label().classes("text-caption text-grey").bind_text_from(
                        self, "fix_pdb_enabled", backward=lambda v: f"Enabled: {v}"
                    )

            # Hydrogen Bonds card
            with ui.card().classes("flex-1"):
                with ui.card_section():
                    with ui.row().classes("items-center justify-between"):
                        ui.label("Hydrogen Bonds").classes("text-h6")
                        ui.button(
                            icon="edit", on_click=lambda: self._edit_hydrogen_bonds()
                        ).props("flat round dense color=primary")
                    ui.label().classes("text-caption text-grey").bind_text_from(
                        self, "hb_distance", backward=lambda v: f"Distance: {v}Å"
                    )
                    ui.label().classes("text-caption text-grey").bind_text_from(
                        self, "hb_angle", backward=lambda v: f"Angle: {v}°"
                    )

        with ui.row().classes("w-full gap-4 q-mt-md"):
            # Halogen Bonds card
            with ui.card().classes("flex-1"):
                with ui.card_section():
                    with ui.row().classes("items-center justify-between"):
                        ui.label("Halogen Bonds").classes("text-h6")
                        ui.button(
                            icon="edit", on_click=lambda: self._edit_halogen_bonds()
                        ).props("flat round dense color=primary")
                    ui.label().classes("text-caption text-grey").bind_text_from(
                        self, "xb_distance", backward=lambda v: f"Distance: {v}Å"
                    )
                    ui.label().classes("text-caption text-grey").bind_text_from(
                        self, "xb_angle", backward=lambda v: f"Angle: {v}°"
                    )

            # π Interactions card
            with ui.card().classes("flex-1"):
                with ui.card_section():
                    with ui.row().classes("items-center justify-between"):
                        ui.label("π Interactions").classes("text-h6")
                        ui.button(
                            icon="edit", on_click=lambda: self._edit_pi_interactions()
                        ).props("flat round dense color=primary")
                    ui.label("8 subtypes configured").classes("text-caption text-grey")

        with ui.row().classes("w-full gap-4 q-mt-md"):
            # π-π Stacking card
            with ui.card().classes("flex-1"):
                with ui.card_section():
                    with ui.row().classes("items-center justify-between"):
                        ui.label("π-π Stacking").classes("text-h6")
                        ui.button(
                            icon="edit", on_click=lambda: self._edit_pi_pi_stacking()
                        ).props("flat round dense color=primary")
                    ui.label().classes("text-caption text-grey").bind_text_from(
                        self, "pi_pi_distance", backward=lambda v: f"Distance: {v}Å"
                    )

            # Carbonyl n→π* card
            with ui.card().classes("flex-1"):
                with ui.card_section():
                    with ui.row().classes("items-center justify-between"):
                        ui.label("Carbonyl n→π*").classes("text-h6")
                        ui.button(
                            icon="edit", on_click=lambda: self._edit_carbonyl()
                        ).props("flat round dense color=primary")
                    ui.label().classes("text-caption text-grey").bind_text_from(
                        self, "carbonyl_distance", backward=lambda v: f"Distance: {v}Å"
                    )

        with ui.row().classes("w-full gap-4 q-mt-md"):
            # n→π* Interactions card
            with ui.card().classes("flex-1"):
                with ui.card_section():
                    with ui.row().classes("items-center justify-between"):
                        ui.label("n→π* Interactions").classes("text-h6")
                        ui.button(
                            icon="edit", on_click=lambda: self._edit_n_pi()
                        ).props("flat round dense color=primary")
                    ui.label().classes("text-caption text-grey").bind_text_from(
                        self, "n_pi_distance", backward=lambda v: f"Distance: {v}Å"
                    )

            # General Settings card
            with ui.card().classes("flex-1"):
                with ui.card_section():
                    with ui.row().classes("items-center justify-between"):
                        ui.label("General Settings").classes("text-h6")
                        ui.button(
                            icon="edit", on_click=lambda: self._edit_general()
                        ).props("flat round dense color=primary")
                    ui.label().classes("text-caption text-grey").bind_text_from(
                        self, "analysis_mode", backward=lambda v: f"Mode: {v}"
                    )

    def _edit_pdb_fixing(self):
        """Open right drawer to edit PDB fixing parameters."""
        self.drawer_content.clear()

        with self.drawer_content:
            # Save button at top
            with ui.row().classes("w-full justify-between items-center q-mb-md"):
                ui.label("PDB Fixing Parameters").classes("text-h6")
                ui.button(
                    "Save", icon="save", on_click=lambda: self._save_and_close()
                ).props("color=primary")

            ui.separator()

            with ui.column().classes("w-full q-mt-md gap-4"):
                ui.switch("Enable PDB Fixing", value=self.fix_pdb_enabled).bind_value(
                    self, "fix_pdb_enabled"
                )

                method_select = ui.select(
                    options=["openbabel", "pdbfixer"],
                    label="Fixing Method",
                    value=self.fix_pdb_method,
                ).bind_value(self, "fix_pdb_method")
                method_select.on_value_change(lambda: self._reset_pdb_fixing_options())

                # Add Hydrogens - disabled for OpenBabel (always true), enabled for PDBFixer
                add_hydrogens_switch = ui.switch(
                    "Add Hydrogens", value=self.fix_pdb_add_hydrogens
                ).bind_value(self, "fix_pdb_add_hydrogens")
                add_hydrogens_switch.bind_enabled_from(
                    self, "fix_pdb_method", lambda m: m == "pdbfixer"
                )

                # PDBFixer-only options (disabled when OpenBabel is selected)
                add_heavy_switch = ui.switch(
                    "Add Heavy Atoms (PDBFixer only)",
                    value=self.fix_pdb_add_heavy_atoms,
                ).bind_value(self, "fix_pdb_add_heavy_atoms")
                add_heavy_switch.bind_enabled_from(
                    self, "fix_pdb_method", lambda m: m == "pdbfixer"
                )

                replace_nonstandard_switch = ui.switch(
                    "Replace Nonstandard Residues (PDBFixer only)",
                    value=self.fix_pdb_replace_nonstandard,
                ).bind_value(self, "fix_pdb_replace_nonstandard")
                replace_nonstandard_switch.bind_enabled_from(
                    self, "fix_pdb_method", lambda m: m == "pdbfixer"
                )

                self.remove_heterogens_switch = ui.switch(
                    "Remove Heterogens (PDBFixer only)",
                    value=self.fix_pdb_remove_heterogens,
                ).bind_value(self, "fix_pdb_remove_heterogens")
                self.remove_heterogens_switch.bind_enabled_from(
                    self, "fix_pdb_method", lambda m: m == "pdbfixer"
                )

                # Keep Water switch - only enabled when Remove Heterogens is checked AND PDBFixer is selected
                self.keep_water_switch = ui.switch(
                    "Keep Water (PDBFixer only)", value=self.fix_pdb_keep_water
                ).bind_value(self, "fix_pdb_keep_water")

                # Update enable state when method changes
                ui.timer(0.1, lambda: self._update_keep_water_state(), once=False)

        self.param_drawer.show()

    def _edit_hydrogen_bonds(self):
        """Open right drawer to edit hydrogen bond parameters."""
        self.drawer_content.clear()

        with self.drawer_content:
            with ui.row().classes("w-full justify-between items-center q-mb-md"):
                ui.label("Hydrogen Bond Parameters").classes("text-h6")
                ui.button(
                    "Save", icon="save", on_click=lambda: self._save_and_close()
                ).props("color=primary")

            ui.separator()

            with ui.column().classes("w-full q-mt-md gap-4"):
                ui.label("Strong H-Bonds (N-H, O-H)").classes("text-subtitle2")
                ui.number(
                    label="H...A Distance (Å)",
                    value=self.hb_distance,
                    min=1.0,
                    max=5.0,
                    step=0.1,
                    precision=1,
                ).bind_value(self, "hb_distance")

                ui.number(
                    label="D-H...A Angle (°)",
                    value=self.hb_angle,
                    min=90,
                    max=180,
                    step=1,
                    precision=0,
                ).bind_value(self, "hb_angle")

                ui.number(
                    label="D...A Distance (Å)",
                    value=self.hb_da_distance,
                    min=1.0,
                    max=5.0,
                    step=0.1,
                    precision=1,
                ).bind_value(self, "hb_da_distance")

                ui.separator()

                ui.label("Weak H-Bonds (C-H)").classes("text-subtitle2")
                ui.number(
                    label="H...A Distance (Å)",
                    value=self.whb_distance,
                    min=1.0,
                    max=5.0,
                    step=0.1,
                    precision=1,
                ).bind_value(self, "whb_distance")

                ui.number(
                    label="D-H...A Angle (°)",
                    value=self.whb_angle,
                    min=90,
                    max=180,
                    step=1,
                    precision=0,
                ).bind_value(self, "whb_angle")

                ui.number(
                    label="D...A Distance (Å)",
                    value=self.whb_da_distance,
                    min=1.0,
                    max=5.0,
                    step=0.1,
                    precision=1,
                ).bind_value(self, "whb_da_distance")

        self.param_drawer.show()

    def _edit_halogen_bonds(self):
        """Open right drawer to edit halogen bond parameters."""
        self.drawer_content.clear()

        with self.drawer_content:
            with ui.row().classes("w-full justify-between items-center q-mb-md"):
                ui.label("Halogen Bond Parameters").classes("text-h6")
                ui.button(
                    "Save", icon="save", on_click=lambda: self._save_and_close()
                ).props("color=primary")

            ui.separator()

            with ui.column().classes("w-full q-mt-md gap-4"):
                ui.number(
                    label="X...A Distance (Å)",
                    value=self.xb_distance,
                    min=2.0,
                    max=5.0,
                    step=0.1,
                    precision=1,
                ).bind_value(self, "xb_distance")

                ui.number(
                    label="C-X...A Angle (°)",
                    value=self.xb_angle,
                    min=120,
                    max=180,
                    step=1,
                    precision=0,
                ).bind_value(self, "xb_angle")

        self.param_drawer.show()

    def _edit_pi_interactions(self):
        """Open right drawer to edit π interaction parameters."""
        self.drawer_content.clear()

        with self.drawer_content:
            with ui.row().classes("w-full justify-between items-center q-mb-md"):
                ui.label("π Interaction Parameters").classes("text-h6")
                ui.button(
                    "Save", icon="save", on_click=lambda: self._save_and_close()
                ).props("color=primary")

            ui.separator()

            with ui.column().classes("w-full q-mt-md gap-4"):
                ui.label("Halogen-π Interactions").classes("text-subtitle2")

                with ui.expansion("C-Cl...π", icon="whatshot").classes("w-full"):
                    ui.number(
                        "Distance (Å)",
                        value=self.pi_ccl_distance,
                        min=2.0,
                        max=6.0,
                        step=0.1,
                    ).bind_value(self, "pi_ccl_distance")
                    ui.number(
                        "Angle (°)", value=self.pi_ccl_angle, min=0, max=180, step=1
                    ).bind_value(self, "pi_ccl_angle")

                with ui.expansion("C-Br...π", icon="whatshot").classes("w-full"):
                    ui.number(
                        "Distance (Å)",
                        value=self.pi_cbr_distance,
                        min=2.0,
                        max=6.0,
                        step=0.1,
                    ).bind_value(self, "pi_cbr_distance")
                    ui.number(
                        "Angle (°)", value=self.pi_cbr_angle, min=0, max=180, step=1
                    ).bind_value(self, "pi_cbr_angle")

                with ui.expansion("C-I...π", icon="whatshot").classes("w-full"):
                    ui.number(
                        "Distance (Å)",
                        value=self.pi_ci_distance,
                        min=2.0,
                        max=6.0,
                        step=0.1,
                    ).bind_value(self, "pi_ci_distance")
                    ui.number(
                        "Angle (°)", value=self.pi_ci_angle, min=0, max=180, step=1
                    ).bind_value(self, "pi_ci_angle")

                ui.separator()
                ui.label("Hydrogen-π Interactions").classes("text-subtitle2")

                with ui.expansion("C-H...π", icon="link").classes("w-full"):
                    ui.number(
                        "Distance (Å)",
                        value=self.pi_ch_distance,
                        min=2.0,
                        max=6.0,
                        step=0.1,
                    ).bind_value(self, "pi_ch_distance")
                    ui.number(
                        "Angle (°)", value=self.pi_ch_angle, min=0, max=180, step=1
                    ).bind_value(self, "pi_ch_angle")

                with ui.expansion("N-H...π", icon="link").classes("w-full"):
                    ui.number(
                        "Distance (Å)",
                        value=self.pi_nh_distance,
                        min=2.0,
                        max=6.0,
                        step=0.1,
                    ).bind_value(self, "pi_nh_distance")
                    ui.number(
                        "Angle (°)", value=self.pi_nh_angle, min=0, max=180, step=1
                    ).bind_value(self, "pi_nh_angle")

                with ui.expansion("O-H...π", icon="link").classes("w-full"):
                    ui.number(
                        "Distance (Å)",
                        value=self.pi_oh_distance,
                        min=2.0,
                        max=6.0,
                        step=0.1,
                    ).bind_value(self, "pi_oh_distance")
                    ui.number(
                        "Angle (°)", value=self.pi_oh_angle, min=0, max=180, step=1
                    ).bind_value(self, "pi_oh_angle")

                with ui.expansion("S-H...π", icon="link").classes("w-full"):
                    ui.number(
                        "Distance (Å)",
                        value=self.pi_sh_distance,
                        min=2.0,
                        max=6.0,
                        step=0.1,
                    ).bind_value(self, "pi_sh_distance")
                    ui.number(
                        "Angle (°)", value=self.pi_sh_angle, min=0, max=180, step=1
                    ).bind_value(self, "pi_sh_angle")

        self.param_drawer.show()

    def _edit_pi_pi_stacking(self):
        """Open right drawer to edit π-π stacking parameters."""
        self.drawer_content.clear()

        with self.drawer_content:
            with ui.row().classes("w-full justify-between items-center q-mb-md"):
                ui.label("π-π Stacking Parameters").classes("text-h6")
                ui.button(
                    "Save", icon="save", on_click=lambda: self._save_and_close()
                ).props("color=primary")

            ui.separator()

            with ui.column().classes("w-full q-mt-md gap-4"):
                ui.number(
                    label="Centroid Distance (Å)",
                    value=self.pi_pi_distance,
                    min=2.0,
                    max=6.0,
                    step=0.1,
                    precision=1,
                ).bind_value(self, "pi_pi_distance")

                ui.number(
                    label="Parallel Angle Cutoff (°)",
                    value=self.pi_pi_parallel_angle,
                    min=0,
                    max=90,
                    step=1,
                    precision=0,
                ).bind_value(self, "pi_pi_parallel_angle")

                ui.number(
                    label="T-shaped Angle Min (°)",
                    value=self.pi_pi_tshaped_angle_min,
                    min=0,
                    max=90,
                    step=1,
                    precision=0,
                ).bind_value(self, "pi_pi_tshaped_angle_min")

                ui.number(
                    label="T-shaped Angle Max (°)",
                    value=self.pi_pi_tshaped_angle_max,
                    min=0,
                    max=90,
                    step=1,
                    precision=0,
                ).bind_value(self, "pi_pi_tshaped_angle_max")

                ui.number(
                    label="Offset Cutoff (Å)",
                    value=self.pi_pi_offset,
                    min=0.0,
                    max=5.0,
                    step=0.1,
                    precision=1,
                ).bind_value(self, "pi_pi_offset")

        self.param_drawer.show()

    def _edit_carbonyl(self):
        """Open right drawer to edit carbonyl n→π* parameters."""
        self.drawer_content.clear()

        with self.drawer_content:
            with ui.row().classes("w-full justify-between items-center q-mb-md"):
                ui.label("Carbonyl n→π* Parameters").classes("text-h6")
                ui.button(
                    "Save", icon="save", on_click=lambda: self._save_and_close()
                ).props("color=primary")

            ui.separator()

            with ui.column().classes("w-full q-mt-md gap-4"):
                ui.number(
                    label="O...C Distance (Å)",
                    value=self.carbonyl_distance,
                    min=2.0,
                    max=5.0,
                    step=0.1,
                    precision=1,
                ).bind_value(self, "carbonyl_distance")

                ui.number(
                    label="Min O...C=O Angle (°)",
                    value=self.carbonyl_angle_min,
                    min=0,
                    max=180,
                    step=1,
                    precision=0,
                ).bind_value(self, "carbonyl_angle_min")

                ui.number(
                    label="Max O...C=O Angle (°)",
                    value=self.carbonyl_angle_max,
                    min=0,
                    max=180,
                    step=1,
                    precision=0,
                ).bind_value(self, "carbonyl_angle_max")

        self.param_drawer.show()

    def _edit_n_pi(self):
        """Open right drawer to edit n→π* parameters."""
        self.drawer_content.clear()

        with self.drawer_content:
            with ui.row().classes("w-full justify-between items-center q-mb-md"):
                ui.label("n→π* Interaction Parameters").classes("text-h6")
                ui.button(
                    "Save", icon="save", on_click=lambda: self._save_and_close()
                ).props("color=primary")

            ui.separator()

            with ui.column().classes("w-full q-mt-md gap-4"):
                ui.number(
                    label="Distance Cutoff (Å)",
                    value=self.n_pi_distance,
                    min=2.0,
                    max=6.0,
                    step=0.1,
                    precision=1,
                ).bind_value(self, "n_pi_distance")

                ui.number(
                    label="Sulfur Distance Cutoff (Å)",
                    value=self.n_pi_sulfur_distance,
                    min=2.0,
                    max=6.0,
                    step=0.1,
                    precision=1,
                ).bind_value(self, "n_pi_sulfur_distance")

                ui.number(
                    label="Min Angle to π Plane (°)",
                    value=self.n_pi_angle_min,
                    min=0,
                    max=90,
                    step=1,
                    precision=0,
                ).bind_value(self, "n_pi_angle_min")

                ui.number(
                    label="Max Angle to π Plane (°)",
                    value=self.n_pi_angle_max,
                    min=0,
                    max=90,
                    step=1,
                    precision=0,
                ).bind_value(self, "n_pi_angle_max")

        self.param_drawer.show()

    def _edit_general(self):
        """Open right drawer to edit general parameters."""
        self.drawer_content.clear()

        with self.drawer_content:
            with ui.row().classes("w-full justify-between items-center q-mb-md"):
                ui.label("General Parameters").classes("text-h6")
                ui.button(
                    "Save", icon="save", on_click=lambda: self._save_and_close()
                ).props("color=primary")

            ui.separator()

            with ui.column().classes("w-full q-mt-md gap-4"):
                ui.select(
                    options={
                        "inter": "Inter-residue only",
                        "all": "All interactions",
                    },
                    label="Interaction Inclusion",
                    value=self.analysis_mode,
                    with_input=True,
                ).bind_value(self, "analysis_mode").tooltip(
                    "inter: Exclude interactions within the same residue; "
                    "all: Include inter-residue and intra-residue interactions"
                )

                ui.number(
                    label="Covalent Cutoff Factor",
                    value=self.covalent_cutoff_factor,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                    precision=2,
                ).bind_value(self, "covalent_cutoff_factor")

        self.param_drawer.show()

    def _save_and_close(self):
        """Save parameters and close drawer."""
        # ui.notify("Parameters updated", type="positive", position="top-left")
        self.param_drawer.hide()

    def get_parameters(self) -> AnalysisParameters:
        """Get current parameters as AnalysisParameters object.

        :returns: Analysis parameters
        :rtype: AnalysisParameters
        """
        return AnalysisParameters(
            # PDB fixing
            fix_pdb_enabled=self.fix_pdb_enabled,
            fix_pdb_method=self.fix_pdb_method,
            fix_pdb_add_hydrogens=self.fix_pdb_add_hydrogens,
            fix_pdb_add_heavy_atoms=self.fix_pdb_add_heavy_atoms,
            fix_pdb_replace_nonstandard=self.fix_pdb_replace_nonstandard,
            fix_pdb_remove_heterogens=self.fix_pdb_remove_heterogens,
            fix_pdb_keep_water=self.fix_pdb_keep_water,
            # Hydrogen bonds
            hb_distance_cutoff=self.hb_distance,
            hb_angle_cutoff=self.hb_angle,
            hb_donor_acceptor_cutoff=self.hb_da_distance,
            # Weak hydrogen bonds
            whb_distance_cutoff=self.whb_distance,
            whb_angle_cutoff=self.whb_angle,
            whb_donor_acceptor_cutoff=self.whb_da_distance,
            # Halogen bonds
            xb_distance_cutoff=self.xb_distance,
            xb_angle_cutoff=self.xb_angle,
            # π interactions
            pi_ccl_distance_cutoff=self.pi_ccl_distance,
            pi_ccl_angle_cutoff=self.pi_ccl_angle,
            pi_cbr_distance_cutoff=self.pi_cbr_distance,
            pi_cbr_angle_cutoff=self.pi_cbr_angle,
            pi_ci_distance_cutoff=self.pi_ci_distance,
            pi_ci_angle_cutoff=self.pi_ci_angle,
            pi_ch_distance_cutoff=self.pi_ch_distance,
            pi_ch_angle_cutoff=self.pi_ch_angle,
            pi_nh_distance_cutoff=self.pi_nh_distance,
            pi_nh_angle_cutoff=self.pi_nh_angle,
            pi_oh_distance_cutoff=self.pi_oh_distance,
            pi_oh_angle_cutoff=self.pi_oh_angle,
            pi_sh_distance_cutoff=self.pi_sh_distance,
            pi_sh_angle_cutoff=self.pi_sh_angle,
            # π-π stacking
            pi_pi_distance_cutoff=self.pi_pi_distance,
            pi_pi_parallel_angle_cutoff=self.pi_pi_parallel_angle,
            pi_pi_tshaped_angle_min=self.pi_pi_tshaped_angle_min,
            pi_pi_tshaped_angle_max=self.pi_pi_tshaped_angle_max,
            pi_pi_offset_cutoff=self.pi_pi_offset,
            # Carbonyl
            carbonyl_distance_cutoff=self.carbonyl_distance,
            carbonyl_angle_min=self.carbonyl_angle_min,
            carbonyl_angle_max=self.carbonyl_angle_max,
            # n→π*
            n_pi_distance_cutoff=self.n_pi_distance,
            n_pi_sulfur_distance_cutoff=self.n_pi_sulfur_distance,
            n_pi_angle_min=self.n_pi_angle_min,
            n_pi_angle_max=self.n_pi_angle_max,
            # General
            covalent_cutoff_factor=self.covalent_cutoff_factor,
            analysis_mode=self.analysis_mode,
        )

    def reset(self):
        """Reset all parameters to default values."""
        # PDB fixing parameters
        self.fix_pdb_enabled = ParametersDefault.FIX_PDB_ENABLED
        self.fix_pdb_method = ParametersDefault.FIX_PDB_METHOD
        self.fix_pdb_add_hydrogens = ParametersDefault.FIX_PDB_ADD_HYDROGENS
        self.fix_pdb_add_heavy_atoms = ParametersDefault.FIX_PDB_ADD_HEAVY_ATOMS
        self.fix_pdb_replace_nonstandard = ParametersDefault.FIX_PDB_REPLACE_NONSTANDARD
        self.fix_pdb_remove_heterogens = ParametersDefault.FIX_PDB_REMOVE_HETEROGENS
        self.fix_pdb_keep_water = ParametersDefault.FIX_PDB_KEEP_WATER

        # Hydrogen bond parameters
        self.hb_distance = ParametersDefault.HB_DISTANCE_CUTOFF
        self.hb_angle = ParametersDefault.HB_ANGLE_CUTOFF
        self.hb_da_distance = ParametersDefault.HB_DA_DISTANCE

        # Weak hydrogen bond parameters
        self.whb_distance = ParametersDefault.WHB_DISTANCE_CUTOFF
        self.whb_angle = ParametersDefault.WHB_ANGLE_CUTOFF
        self.whb_da_distance = ParametersDefault.WHB_DA_DISTANCE

        # Halogen bond parameters
        self.xb_distance = ParametersDefault.XB_DISTANCE_CUTOFF
        self.xb_angle = ParametersDefault.XB_ANGLE_CUTOFF

        # π interaction subtype parameters
        self.pi_ccl_distance = ParametersDefault.PI_CCL_DISTANCE_CUTOFF
        self.pi_ccl_angle = ParametersDefault.PI_CCL_ANGLE_CUTOFF
        self.pi_cbr_distance = ParametersDefault.PI_CBR_DISTANCE_CUTOFF
        self.pi_cbr_angle = ParametersDefault.PI_CBR_ANGLE_CUTOFF
        self.pi_ci_distance = ParametersDefault.PI_CI_DISTANCE_CUTOFF
        self.pi_ci_angle = ParametersDefault.PI_CI_ANGLE_CUTOFF
        self.pi_ch_distance = ParametersDefault.PI_CH_DISTANCE_CUTOFF
        self.pi_ch_angle = ParametersDefault.PI_CH_ANGLE_CUTOFF
        self.pi_nh_distance = ParametersDefault.PI_NH_DISTANCE_CUTOFF
        self.pi_nh_angle = ParametersDefault.PI_NH_ANGLE_CUTOFF
        self.pi_oh_distance = ParametersDefault.PI_OH_DISTANCE_CUTOFF
        self.pi_oh_angle = ParametersDefault.PI_OH_ANGLE_CUTOFF
        self.pi_sh_distance = ParametersDefault.PI_SH_DISTANCE_CUTOFF
        self.pi_sh_angle = ParametersDefault.PI_SH_ANGLE_CUTOFF

        # π-π stacking parameters
        self.pi_pi_distance = ParametersDefault.PI_PI_DISTANCE_CUTOFF
        self.pi_pi_parallel_angle = ParametersDefault.PI_PI_PARALLEL_ANGLE_CUTOFF
        self.pi_pi_tshaped_angle_min = ParametersDefault.PI_PI_TSHAPED_ANGLE_MIN
        self.pi_pi_tshaped_angle_max = ParametersDefault.PI_PI_TSHAPED_ANGLE_MAX
        self.pi_pi_offset = ParametersDefault.PI_PI_OFFSET_CUTOFF

        # Carbonyl interaction parameters
        self.carbonyl_distance = ParametersDefault.CARBONYL_DISTANCE_CUTOFF
        self.carbonyl_angle_min = ParametersDefault.CARBONYL_ANGLE_MIN
        self.carbonyl_angle_max = ParametersDefault.CARBONYL_ANGLE_MAX

        # n→π* interaction parameters
        self.n_pi_distance = ParametersDefault.N_PI_DISTANCE_CUTOFF
        self.n_pi_sulfur_distance = ParametersDefault.N_PI_SULFUR_DISTANCE_CUTOFF
        self.n_pi_angle_min = ParametersDefault.N_PI_ANGLE_MIN
        self.n_pi_angle_max = ParametersDefault.N_PI_ANGLE_MAX

        # General parameters
        self.covalent_cutoff_factor = ParametersDefault.COVALENT_CUTOFF_FACTOR
        self.analysis_mode = ParametersDefault.ANALYSIS_MODE
