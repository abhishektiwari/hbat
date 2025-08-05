"""
Unit tests for GUI component functionality.

These tests focus on pure logic and data validation without requiring
full GUI setup, tkinter dependencies, or complex integrations.
"""

import pytest
import json
from unittest.mock import Mock


@pytest.mark.unit
class TestGUIImports:
    """Test that GUI modules can be imported."""
    
    def test_gui_module_imports(self):
        """Test importing GUI modules."""
        try:
            from hbat.gui.main_window import MainWindow
            from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
            from hbat.gui.pdb_fixing_dialog import PDBFixingDialog
            from hbat.gui.preset_manager_dialog import PresetManagerDialog
            from hbat.gui.results_panel import ResultsPanel
            assert True, "GUI modules imported successfully"
        except ImportError as e:
            pytest.skip(f"GUI modules not available: {e}")
    
    def test_chain_visualization_import(self):
        """Test importing chain visualization module."""
        try:
            from hbat.gui.chain_visualization import ChainVisualizationWindow
            assert ChainVisualizationWindow is not None
        except ImportError as e:
            pytest.skip(f"Chain visualization not available: {e}")
    
    def test_renderer_imports(self):
        """Test importing visualization renderers."""
        try:
            from hbat.gui.graphviz_renderer import GraphVizRenderer
            from hbat.gui.matplotlib_renderer import MatplotlibRenderer
            assert GraphVizRenderer is not None
            assert MatplotlibRenderer is not None
        except ImportError as e:
            pytest.skip(f"Renderer modules not available: {e}")


@pytest.mark.unit
class TestDialogLogic:
    """Test dialog logic without GUI components."""
    
    def test_analysis_parameters_integration(self):
        """Test that analysis parameters work with dialog architecture."""
        try:
            from hbat.core.analysis import AnalysisParameters
            
            # Test creating parameters that would be used by dialogs
            params = AnalysisParameters(
                hb_distance_cutoff=3.2,
                hb_angle_cutoff=125.0,
                fix_pdb_enabled=True,
                fix_pdb_method="pdbfixer"
            )
            
            # Verify parameters have expected attributes
            assert params.hb_distance_cutoff == 3.2
            assert params.hb_angle_cutoff == 125.0
            assert params.fix_pdb_enabled is True
            assert params.fix_pdb_method == "pdbfixer"
            
        except ImportError as e:
            pytest.skip(f"Analysis parameters not available: {e}")
    
    def test_preset_data_structure(self):
        """Test preset data structure for dialog compatibility."""
        # Test preset data format that would be used by PresetManagerDialog
        preset_data = {
            "format_version": "1.0",
            "parameters": {
                "hydrogen_bonds": {
                    "h_a_distance_cutoff": 3.5,
                    "dha_angle_cutoff": 120.0,
                    "d_a_distance_cutoff": 4.0
                },
                "halogen_bonds": {
                    "x_a_distance_cutoff": 4.0,
                    "cxa_angle_cutoff": 120.0
                },
                "pi_interactions": {
                    "h_pi_distance_cutoff": 4.5,
                    "dh_pi_angle_cutoff": 90.0
                },
                "general": {
                    "covalent_cutoff_factor": 0.85,
                    "analysis_mode": "complete"
                },
                "pdb_fixing": {
                    "enabled": False,
                    "method": "openbabel",
                    "add_hydrogens": True,
                    "add_heavy_atoms": False,
                    "replace_nonstandard": False,
                    "remove_heterogens": False,
                    "keep_water": True
                }
            }
        }
        
        # Verify preset structure
        assert "format_version" in preset_data
        assert "parameters" in preset_data
        assert "hydrogen_bonds" in preset_data["parameters"]
        assert "pdb_fixing" in preset_data["parameters"]
        
        # Test JSON serialization (important for preset saving)
        json_string = json.dumps(preset_data)
        loaded_data = json.loads(json_string)
        assert loaded_data == preset_data


@pytest.mark.unit  
class TestParameterValidation:
    """Test parameter validation logic."""
    
    def test_parameter_ranges(self):
        """Test parameter range validation."""
        try:
            from hbat.core.analysis import AnalysisParameters
            
            # Test valid parameters
            valid_params = AnalysisParameters(
                hb_distance_cutoff=3.5,
                hb_angle_cutoff=120.0,
                hb_donor_acceptor_cutoff=4.0
            )
            
            assert valid_params.hb_distance_cutoff == 3.5
            assert valid_params.hb_angle_cutoff == 120.0
            assert valid_params.hb_donor_acceptor_cutoff == 4.0
            
        except ImportError as e:
            pytest.skip(f"Analysis parameters not available: {e}")
    
    def test_pdb_fixing_parameters(self):
        """Test PDB fixing parameter validation."""
        try:
            from hbat.core.analysis import AnalysisParameters
            
            # Test PDB fixing parameters
            pdb_params = AnalysisParameters(
                fix_pdb_enabled=True,
                fix_pdb_method="pdbfixer",
                fix_pdb_add_hydrogens=True,
                fix_pdb_add_heavy_atoms=True,
                fix_pdb_replace_nonstandard=False
            )
            
            assert pdb_params.fix_pdb_enabled is True
            assert pdb_params.fix_pdb_method == "pdbfixer"
            assert pdb_params.fix_pdb_add_hydrogens is True
            assert pdb_params.fix_pdb_add_heavy_atoms is True
            assert pdb_params.fix_pdb_replace_nonstandard is False
            
        except ImportError as e:
            pytest.skip(f"Analysis parameters not available: {e}")


@pytest.mark.unit
class TestConfigIntegration:
    """Test configuration integration with GUI components."""
    
    def test_app_config_import(self):
        """Test that app config can be imported."""
        try:
            from hbat.core.app_config import HBATConfig
            
            config = HBATConfig()
            assert config is not None
            
        except ImportError as e:
            pytest.skip(f"App config not available: {e}")
    
    def test_parameter_defaults(self):
        """Test parameter defaults integration."""
        try:
            from hbat.constants.parameters import ParametersDefault
            
            # Test that defaults exist for dialog use
            assert hasattr(ParametersDefault, 'HB_DISTANCE_CUTOFF')
            assert hasattr(ParametersDefault, 'HB_ANGLE_CUTOFF')
            assert hasattr(ParametersDefault, 'ANALYSIS_MODE')
            
        except ImportError as e:
            pytest.skip(f"Parameter constants not available: {e}")


@pytest.mark.unit
class TestAnalysisIntegration:
    """Test analysis component integration."""
    
    def test_analyzer_import(self):
        """Test that analyzer can be imported."""
        try:
            from hbat.core.analysis import NPMolecularInteractionAnalyzer
            
            assert NPMolecularInteractionAnalyzer is not None
            
        except ImportError as e:
            pytest.skip(f"Analyzer not available: {e}")
    
    def test_analyzer_parameter_integration(self):
        """Test analyzer works with GUI parameters."""
        try:
            from hbat.core.analysis import AnalysisParameters, NPMolecularInteractionAnalyzer
            
            # Test creating analyzer with GUI parameters
            params = AnalysisParameters(
                hb_distance_cutoff=3.2,
                analysis_mode="complete"
            )
            
            analyzer = NPMolecularInteractionAnalyzer(params)
            assert analyzer is not None
            
        except ImportError as e:
            pytest.skip(f"Analysis components not available: {e}")


@pytest.mark.unit
class TestDialogDataStructures:
    """Test data structures used by dialogs."""
    
    def test_geometry_parameters_structure(self):
        """Test geometry parameters data structure."""
        try:
            from hbat.core.analysis import AnalysisParameters
            
            # Test parameters used by GeometryCutoffsDialog
            geo_params = AnalysisParameters(
                hb_distance_cutoff=3.5,
                hb_angle_cutoff=120.0,
                hb_donor_acceptor_cutoff=4.0,
                xb_distance_cutoff=4.0,
                xb_angle_cutoff=120.0,
                pi_distance_cutoff=4.5,
                pi_angle_cutoff=90.0,
                analysis_mode="complete"
            )
            
            # Verify all geometry parameters exist
            assert hasattr(geo_params, 'hb_distance_cutoff')
            assert hasattr(geo_params, 'hb_angle_cutoff')
            assert hasattr(geo_params, 'hb_donor_acceptor_cutoff')
            assert hasattr(geo_params, 'analysis_mode')
            
        except ImportError as e:
            pytest.skip(f"Geometry parameters not available: {e}")
    
    def test_pdb_fixing_structure(self):
        """Test PDB fixing parameters structure."""
        try:
            from hbat.core.analysis import AnalysisParameters
            
            # Test parameters used by PDBFixingDialog
            pdb_params = AnalysisParameters(
                fix_pdb_enabled=True,
                fix_pdb_method="pdbfixer",
                fix_pdb_add_hydrogens=True,
                fix_pdb_add_heavy_atoms=True,
                fix_pdb_replace_nonstandard=False,
                fix_pdb_remove_heterogens=False,
                fix_pdb_keep_water=True
            )
            
            # Verify all PDB fixing parameters exist
            assert hasattr(pdb_params, 'fix_pdb_enabled')
            assert hasattr(pdb_params, 'fix_pdb_method')
            assert hasattr(pdb_params, 'fix_pdb_add_hydrogens')
            assert hasattr(pdb_params, 'fix_pdb_add_heavy_atoms')
            
        except ImportError as e:
            pytest.skip(f"PDB fixing parameters not available: {e}")


@pytest.mark.unit
class TestUtilityIntegration:
    """Test utility integration with GUI components."""
    
    def test_graphviz_utils_import(self):
        """Test GraphViz utilities can be imported."""
        try:
            from hbat.utilities.graphviz_utils import GraphVizDetector
            
            assert GraphVizDetector is not None
            assert hasattr(GraphVizDetector, 'is_graphviz_available')
            
        except ImportError as e:
            pytest.skip(f"GraphViz utilities not available: {e}")
    
    def test_export_manager_import(self):
        """Test export manager can be imported."""
        try:
            from hbat.gui.export_manager import ExportManager
            
            assert ExportManager is not None
            
        except ImportError as e:
            pytest.skip(f"Export manager not available: {e}")