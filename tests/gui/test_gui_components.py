"""
GUI component tests for HBAT.

Tests the new modal dialog architecture including GeometryCutoffsDialog,
PDBFixingDialog, PresetManagerDialog, and MainWindow integration.
"""

import pytest
import tkinter as tk
import tempfile
import os
from unittest.mock import Mock, patch
from pathlib import Path


@pytest.mark.gui
class TestGeometryCutoffsDialog:
    """Test geometry cutoffs dialog functionality."""
    
    def setup_method(self):
        """Set up test environment."""
        self.root = tk.Tk()
        self.root.withdraw()  # Hide window during testing
        
    def teardown_method(self):
        """Clean up test environment."""
        try:
            self.root.quit()
            self.root.destroy()
        except Exception:
            pass
    
    def test_dialog_creation(self):
        """Test that geometry cutoffs dialog can be created."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        test_params = AnalysisParameters(hb_distance_cutoff=3.2)
        dialog = GeometryCutoffsDialog(self.root, test_params)
        
        try:
            assert dialog is not None
            assert hasattr(dialog, 'dialog')
            assert hasattr(dialog, 'result')
            assert hasattr(dialog, 'current_params')
            assert dialog.current_params.hb_distance_cutoff == 3.2
        finally:
            dialog.dialog.destroy()
    
    def test_dialog_modal_properties(self):
        """Test that dialog is properly configured as modal."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Check modal properties
            assert dialog.dialog.transient() is not None
            assert dialog.dialog.winfo_class() == 'Toplevel'
        finally:
            dialog.dialog.destroy()
    
    def test_preset_manager_integration(self):
        """Test preset manager can be opened from geometry dialog."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            assert hasattr(dialog, '_open_preset_manager')
            
            # Mock preset manager to avoid creating another modal dialog
            with patch('hbat.gui.preset_manager_dialog.PresetManagerDialog') as mock_preset:
                mock_instance = Mock()
                mock_instance.get_result.return_value = None  # User cancelled
                mock_preset.return_value = mock_instance
                
                # Should not raise exceptions
                dialog._open_preset_manager()
                mock_preset.assert_called_once()
                
        finally:
            dialog.dialog.destroy()
    
    def test_parameter_storage_initialization(self):
        """Test parameter storage is initialized correctly."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        test_params = AnalysisParameters(hb_distance_cutoff=3.2)
        dialog = GeometryCutoffsDialog(self.root, test_params)
        
        try:
            # Test _param_values dictionary exists
            assert hasattr(dialog, '_param_values')
            assert isinstance(dialog._param_values, dict)
            
            # Test _store_current_values method exists
            assert hasattr(dialog, '_store_current_values')
            
        finally:
            dialog.dialog.destroy()
    
    def test_parameter_persistence_across_categories(self):
        """Test parameters persist when switching categories."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        # Create dialog with specific parameters
        test_params = AnalysisParameters(
            hb_distance_cutoff=3.2,
            hb_angle_cutoff=125.0,
            whb_distance_cutoff=3.8,
            xb_distance_cutoff=4.2
        )
        dialog = GeometryCutoffsDialog(self.root, test_params)
        
        try:
            # Initially should be on first category (General Parameters)
            assert dialog.category_listbox.curselection() == (0,)
            
            # Switch to Hydrogen Bonds category
            dialog.category_listbox.selection_clear(0)
            dialog.category_listbox.selection_set(1)
            dialog._on_category_selected(None)
            
            # Verify parameters are stored
            dialog._store_current_values()
            
            # Check that stored values contain expected parameters
            stored_values = dialog._param_values
            if 'hb_distance' in stored_values:
                assert stored_values['hb_distance'] == 3.2
            if 'hb_angle' in stored_values:
                assert stored_values['hb_angle'] == 125.0
                
        finally:
            dialog.dialog.destroy()
    
    def test_get_parameters_with_storage(self):
        """Test get_parameters method uses stored values correctly."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Set parameters through the proper method (which stores them)
            test_params = AnalysisParameters(
                hb_distance_cutoff=3.2,
                hb_angle_cutoff=125.0,
                whb_distance_cutoff=3.8,
                analysis_mode='local'
            )
            dialog.set_parameters(test_params)
            
            # Get parameters should use stored values
            result_params = dialog.get_parameters()
            
            assert result_params.hb_distance_cutoff == 3.2
            assert result_params.hb_angle_cutoff == 125.0
            assert result_params.whb_distance_cutoff == 3.8
            assert result_params.analysis_mode == 'local'
            
        finally:
            dialog.dialog.destroy()
    
    def test_set_parameters_stores_values(self):
        """Test set_parameters method stores values in _param_values."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Set new parameters
            new_params = AnalysisParameters(
                hb_distance_cutoff=3.1,
                hb_angle_cutoff=130.0,
                whb_distance_cutoff=3.7,
                analysis_mode='complete'
            )
            
            dialog.set_parameters(new_params)
            
            # Check that values were stored
            assert dialog._param_values['hb_distance'] == 3.1
            assert dialog._param_values['hb_angle'] == 130.0
            assert dialog._param_values['whb_distance'] == 3.7
            assert dialog._param_values['analysis_mode'] == 'complete'
            
        finally:
            dialog.dialog.destroy()
    
    def test_reset_to_defaults_no_errors(self):
        """Test reset to defaults doesn't cause trace callback errors."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        from hbat.constants.parameters import ParametersDefault
        
        # Create dialog with modified parameters
        test_params = AnalysisParameters(
            hb_distance_cutoff=3.2,
            hb_angle_cutoff=125.0
        )
        dialog = GeometryCutoffsDialog(self.root, test_params)
        
        try:
            # Reset to defaults should not raise exceptions
            dialog._set_defaults()
            
            # Check that parameters were reset
            result_params = dialog.get_parameters()
            assert result_params.hb_distance_cutoff == ParametersDefault.HB_DISTANCE_CUTOFF
            assert result_params.hb_angle_cutoff == ParametersDefault.HB_ANGLE_CUTOFF
            
        finally:
            dialog.dialog.destroy()
    
    def test_category_switching_preserves_changes(self):
        """Test that parameter changes are preserved when switching categories."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Set parameters first so they exist in variables
            initial_params = AnalysisParameters(
                hb_distance_cutoff=3.3,
                hb_angle_cutoff=135.0
            )
            dialog.set_parameters(initial_params)
            
            # Switch to different category and back
            dialog.category_listbox.selection_clear(0)
            dialog.category_listbox.selection_set(2)  # Weak Hydrogen Bonds
            dialog._on_category_selected(None)
            
            dialog.category_listbox.selection_clear(2)
            dialog.category_listbox.selection_set(1)  # Back to Hydrogen Bonds
            dialog._on_category_selected(None)
            
            # Changes should be preserved
            final_params = dialog.get_parameters()
            assert final_params.hb_distance_cutoff == 3.3
            assert final_params.hb_angle_cutoff == 135.0
            
        finally:
            dialog.dialog.destroy()
    
    def test_pi_interaction_subtype_persistence(self):
        """Test π interaction subtype parameters are persisted correctly."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        from hbat.constants.parameters import ParametersDefault
        
        # Create parameters with specific π subtype values
        test_params = AnalysisParameters(
            pi_ccl_distance_cutoff=3.9,
            pi_ccl_angle_cutoff=125.0,
            pi_ch_distance_cutoff=4.1,
            pi_ch_angle_cutoff=115.0
        )
        
        dialog = GeometryCutoffsDialog(self.root, test_params)
        
        try:
            # Set the parameters (this should store them)
            dialog.set_parameters(test_params)
            
            # Switch to π interactions category to create the subtype widgets
            dialog.category_listbox.selection_clear(0)
            dialog.category_listbox.selection_set(4)  # π Interactions
            dialog._on_category_selected(None)
            
            # Get final parameters
            result_params = dialog.get_parameters()
            
            # Verify π subtype parameters were preserved
            assert result_params.pi_ccl_distance_cutoff == 3.9
            assert result_params.pi_ccl_angle_cutoff == 125.0
            assert result_params.pi_ch_distance_cutoff == 4.1
            assert result_params.pi_ch_angle_cutoff == 115.0
            
        finally:
            dialog.dialog.destroy()
    
    def test_error_handling_destroyed_widgets(self):
        """Test error handling when widgets are destroyed."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Set parameters properly first
            test_params = AnalysisParameters(hb_distance_cutoff=3.4)
            dialog.set_parameters(test_params)
            
            # Destroy current content (simulating category switch)
            if dialog.current_content:
                dialog.current_content.destroy()
                dialog.current_content = None
            
            # get_parameters should still work using stored values
            result_params = dialog.get_parameters()
            assert result_params.hb_distance_cutoff == 3.4
            
        finally:
            dialog.dialog.destroy()


@pytest.mark.gui
class TestGeometryCutoffsDialogAdvanced:
    """Advanced tests for GeometryCutoffsDialog parameter persistence."""
    
    def setup_method(self):
        """Set up test environment."""
        self.root = tk.Tk()
        self.root.withdraw()  # Hide window during testing
        
    def teardown_method(self):
        """Clean up test environment."""
        try:
            self.root.quit()
            self.root.destroy()
        except Exception:
            pass
    
    def test_complete_parameter_workflow(self):
        """Test complete workflow: set parameters → modify → switch categories → get final."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        # Step 1: Create dialog with initial parameters
        initial_params = AnalysisParameters(
            hb_distance_cutoff=3.1,
            hb_angle_cutoff=115.0,
            whb_distance_cutoff=3.7,
            xb_distance_cutoff=4.1,
            analysis_mode='local'
        )
        
        dialog = GeometryCutoffsDialog(self.root, initial_params)
        
        try:
            # Step 2: Verify initial parameters are set
            result1 = dialog.get_parameters()
            assert result1.hb_distance_cutoff == 3.1
            assert result1.analysis_mode == 'local'
            
            # Step 3: Simulate parameter modifications by setting new parameters
            modified_params = AnalysisParameters(
                hb_distance_cutoff=3.2,  # Modified
                hb_angle_cutoff=115.0,   # Original
                whb_distance_cutoff=3.8, # Modified
                xb_distance_cutoff=4.2,  # Modified
                analysis_mode='local'    # Original
            )
            dialog.set_parameters(modified_params)
            
            # Step 4: Switch between categories
            for category_idx in [1, 2, 3, 0]:  # HB → WHB → XB → General
                dialog.category_listbox.selection_clear(0, tk.END)
                dialog.category_listbox.selection_set(category_idx)
                dialog._on_category_selected(None)
            
            # Step 5: Get final parameters
            final_params = dialog.get_parameters()
            
            # Verify modifications were preserved
            assert final_params.hb_distance_cutoff == 3.2      # Modified
            assert final_params.hb_angle_cutoff == 115.0       # Original
            assert final_params.whb_distance_cutoff == 3.8     # Modified
            assert final_params.xb_distance_cutoff == 4.2      # Modified
            assert final_params.analysis_mode == 'local'       # Original
            
        finally:
            dialog.dialog.destroy()
    
    def test_preset_loading_with_persistence(self):
        """Test that preset loading works with parameter persistence."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Simulate loading a preset by calling _apply_preset_data
            preset_data = {
                "parameters": {
                    "hydrogen_bonds": {
                        "h_a_distance_cutoff": 3.3,
                        "dha_angle_cutoff": 125.0,
                        "d_a_distance_cutoff": 4.1
                    },
                    "halogen_bonds": {
                        "x_a_distance_cutoff": 4.2,
                        "dxa_angle_cutoff": 130.0
                    },
                    "general": {
                        "covalent_cutoff_factor": 0.9,
                        "analysis_mode": "complete"
                    }
                }
            }
            
            # Apply preset data
            dialog._apply_preset_data(preset_data)
            
            # Verify parameters were loaded
            result_params = dialog.get_parameters()
            assert result_params.hb_distance_cutoff == 3.3
            assert result_params.hb_angle_cutoff == 125.0
            assert result_params.xb_distance_cutoff == 4.2
            assert result_params.analysis_mode == "complete"
            
        finally:
            dialog.dialog.destroy()
    
    def test_parameter_validation_edge_cases(self):
        """Test parameter handling with edge cases and invalid values."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        from hbat.constants.parameters import ParametersDefault
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Test that default parameters are returned when dialog is first created
            result = dialog.get_parameters()
            
            # Should return default values (the dialog initializes with defaults)
            assert result.hb_distance_cutoff == ParametersDefault.HB_DISTANCE_CUTOFF
            assert result.hb_angle_cutoff == ParametersDefault.HB_ANGLE_CUTOFF
            
            # Test setting specific parameter and ensuring others remain defaults
            partial_params = AnalysisParameters(hb_distance_cutoff=3.5)
            dialog.set_parameters(partial_params)
            result = dialog.get_parameters()
            
            # Should use set value for hb_distance, defaults for others
            assert result.hb_distance_cutoff == 3.5
            assert result.hb_angle_cutoff == ParametersDefault.HB_ANGLE_CUTOFF
            
        finally:
            dialog.dialog.destroy()
    
    def test_rapid_category_switching(self):
        """Test rapid category switching doesn't cause errors."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Rapidly switch between all categories multiple times
            categories = [0, 1, 2, 3, 4, 0, 4, 1, 3, 2]
            
            for category_idx in categories:
                dialog.category_listbox.selection_clear(0, tk.END)
                dialog.category_listbox.selection_set(category_idx)
                # Should not raise exceptions
                dialog._on_category_selected(None)
                
                # Verify dialog is still functional
                assert dialog.current_content is not None
                assert hasattr(dialog, '_param_values')
            
            # Should still be able to get parameters
            final_params = dialog.get_parameters()
            assert final_params is not None
            
        finally:
            dialog.dialog.destroy()
    
    def test_tcl_error_handling_in_traces(self):
        """Test that TclError exceptions in trace callbacks are handled gracefully."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        import tkinter as tk
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Go to Hydrogen Bonds category to create widgets with traces
            dialog.category_listbox.selection_clear(0)
            dialog.category_listbox.selection_set(1)
            dialog._on_category_selected(None)
            
            # Get reference to a widget that should have traces
            hb_distance_var = getattr(dialog, 'hb_distance', None)
            
            if hb_distance_var is not None:
                # Destroy the current content while traces might still be active
                if dialog.current_content:
                    dialog.current_content.destroy()
                    dialog.current_content = None
                
                # Trying to trigger trace callback on destroyed widget should not crash
                try:
                    hb_distance_var.set(3.7)  # This might trigger TclError
                except tk.TclError:
                    pass  # Expected - widget destroyed
                
                # Dialog should still be functional
                result = dialog.get_parameters()
                assert result is not None
            
        finally:
            dialog.dialog.destroy()
    
    def test_widget_destruction_during_parameter_retrieval(self):
        """Test parameter retrieval when widgets are being destroyed."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Set parameters properly first
            test_params = AnalysisParameters(
                hb_distance_cutoff=3.4,
                hb_angle_cutoff=128.0,
                analysis_mode='local'
            )
            dialog.set_parameters(test_params)
            
            # Create widgets for hydrogen bonds
            dialog.category_listbox.selection_clear(0)
            dialog.category_listbox.selection_set(1)
            dialog._on_category_selected(None)
            
            # Destroy content while getting parameters
            if dialog.current_content:
                dialog.current_content.destroy()
                dialog.current_content = None
            
            # get_parameters should still work using stored values
            result = dialog.get_parameters()
            assert result.hb_distance_cutoff == 3.4
            assert result.hb_angle_cutoff == 128.0
            assert result.analysis_mode == 'local'
            
        finally:
            dialog.dialog.destroy()
    
    def test_safe_parameter_variable_access(self):
        """Test safe access to parameter variables when they may not exist."""
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        from hbat.constants.parameters import ParametersDefault
        
        dialog = GeometryCutoffsDialog(self.root, AnalysisParameters())
        
        try:
            # Test the get_value logic used in get_parameters
            def get_value(var_name, default_value):
                if var_name in dialog._param_values:
                    return dialog._param_values[var_name]
                elif hasattr(dialog, var_name):
                    var = getattr(dialog, var_name)
                    if var:
                        try:
                            return var.get()
                        except tk.TclError:
                            pass
                return default_value
            
            # Test with no stored values and no variables
            result = get_value('hb_distance', ParametersDefault.HB_DISTANCE_CUTOFF)
            assert result == ParametersDefault.HB_DISTANCE_CUTOFF
            
            # Test with stored values
            dialog._param_values['hb_distance'] = 3.6
            result = get_value('hb_distance', ParametersDefault.HB_DISTANCE_CUTOFF)
            assert result == 3.6
            
        finally:
            dialog.dialog.destroy()


@pytest.mark.gui  
class TestPDBFixingDialog:
    """Test PDB fixing dialog functionality."""
    
    def setup_method(self):
        """Set up test environment."""
        self.root = tk.Tk()
        self.root.withdraw()
        
    def teardown_method(self):
        """Clean up test environment."""
        try:
            self.root.quit()
            self.root.destroy()
        except Exception:
            pass
    
    def test_dialog_creation(self):
        """Test PDB fixing dialog creation."""
        from hbat.gui.pdb_fixing_dialog import PDBFixingDialog
        
        dialog = PDBFixingDialog(self.root)
        
        try:
            assert dialog is not None
            assert hasattr(dialog, 'dialog')
            assert hasattr(dialog, 'result')
        finally:
            dialog.dialog.destroy()
    
    def test_preset_manager_integration(self):
        """Test preset manager integration in PDB fixing dialog."""
        from hbat.gui.pdb_fixing_dialog import PDBFixingDialog
        
        dialog = PDBFixingDialog(self.root)
        
        try:
            assert hasattr(dialog, '_open_preset_manager')
            
            # Mock preset manager
            with patch('hbat.gui.preset_manager_dialog.PresetManagerDialog') as mock_preset:
                mock_instance = Mock()
                mock_instance.get_result.return_value = None
                mock_preset.return_value = mock_instance
                
                dialog._open_preset_manager()
                mock_preset.assert_called_once()
                
        finally:
            dialog.dialog.destroy()


@pytest.mark.gui
class TestPresetManagerDialog:
    """Test preset manager dialog functionality."""
    
    def setup_method(self):
        """Set up test environment."""
        self.root = tk.Tk()
        self.root.withdraw()
        
    def teardown_method(self):
        """Clean up test environment."""
        try:
            self.root.quit()
            self.root.destroy()
        except Exception:
            pass
    
    def test_dialog_creation(self):
        """Test preset manager dialog creation."""
        from hbat.gui.preset_manager_dialog import PresetManagerDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = PresetManagerDialog(self.root, AnalysisParameters())
        
        try:
            assert dialog is not None
            assert hasattr(dialog, 'dialog')
            assert hasattr(dialog, 'result')
            assert hasattr(dialog, 'preset_tree')
            assert hasattr(dialog, 'preset_file_paths')
        finally:
            dialog.dialog.destroy()
    
    def test_preset_tree_functionality(self):
        """Test preset tree view functionality."""
        from hbat.gui.preset_manager_dialog import PresetManagerDialog
        from hbat.core.analysis import AnalysisParameters
        
        dialog = PresetManagerDialog(self.root, AnalysisParameters())
        
        try:
            # Test tree refresh doesn't crash
            dialog._refresh_preset_list()
            
            # Verify tree structure
            assert hasattr(dialog.preset_tree, 'get_children')
            
        finally:
            dialog.dialog.destroy()


@pytest.mark.gui
class TestMainWindow:
    """Test main window functionality with new dialog architecture."""
    
    def setup_method(self):
        """Set up test environment."""
        # Stop any existing async executor to prevent event loop conflicts  
        try:
            import tk_async_execute as tae
            tae.stop()
        except Exception:
            pass
        
    def teardown_method(self):
        """Clean up test environment."""
        # Stop async executor to prevent event loop conflicts
        try:
            import tk_async_execute as tae
            tae.stop()
        except Exception:
            pass
    
    def test_main_window_import(self):
        """Test importing main window class."""
        from hbat.gui.main_window import MainWindow
        assert MainWindow is not None
    
    def test_main_window_creation(self):
        """Test main window creation and cleanup."""
        from hbat.gui.main_window import MainWindow
        import tk_async_execute as tae
        
        main_window = None
        try:
            main_window = MainWindow()
            main_window.root.withdraw()  # Hide during testing
            
            # Test essential attributes for new architecture
            assert hasattr(main_window, 'root'), "MainWindow should have root attribute"
            assert hasattr(main_window, 'results_panel'), "MainWindow should have results_panel attribute"
            assert hasattr(main_window, 'analyzer'), "MainWindow should have analyzer attribute"
            assert main_window.session_parameters is None, "session_parameters should be initialized as None"
            
            # Test that menu actions exist for new dialogs
            assert hasattr(main_window, '_open_parameters_window'), "Should have _open_parameters_window method"
            assert hasattr(main_window, '_open_pdb_fixing_window'), "Should have _open_pdb_fixing_window method"
            
        finally:
            if main_window and hasattr(main_window, 'root'):
                try:
                    tae.stop()  # Stop async executor first
                    main_window.root.quit()
                    main_window.root.destroy()
                except Exception:
                    pass
    
    def test_dialog_integration(self):
        """Test main window integration with new dialogs."""
        from hbat.gui.main_window import MainWindow
        import tk_async_execute as tae
        
        main_window = None
        try:
            main_window = MainWindow()
            main_window.root.withdraw()
            
            # Mock dialogs to avoid creating modal windows
            with patch('hbat.gui.main_window.GeometryCutoffsDialog') as mock_geo:
                mock_instance = Mock()
                mock_instance.get_result.return_value = None
                mock_geo.return_value = mock_instance
                
                # Test opening geometry cutoffs dialog
                main_window._open_parameters_window()
                mock_geo.assert_called_once()
            
            with patch('hbat.gui.main_window.PDBFixingDialog') as mock_pdb:
                mock_instance = Mock()
                mock_instance.get_result.return_value = None
                mock_pdb.return_value = mock_instance
                
                # Test opening PDB fixing dialog
                main_window._open_pdb_fixing_window()
                mock_pdb.assert_called_once()
                
        finally:
            if main_window and hasattr(main_window, 'root'):
                try:
                    tae.stop()  # Stop async executor first
                    main_window.root.quit()
                    main_window.root.destroy()
                except Exception:
                    pass


@pytest.mark.gui
class TestResultsPanel:
    """Test results panel functionality (unchanged from before)."""
    
    def setup_method(self):
        """Set up test environment."""
        self.root = tk.Tk()
        self.root.withdraw()
        
    def teardown_method(self):
        """Clean up test environment."""
        try:
            self.root.quit()
            self.root.destroy()
        except Exception:
            pass
    
    def test_results_panel_creation(self):
        """Test results panel creation."""
        from hbat.gui.results_panel import ResultsPanel
        
        panel = ResultsPanel(self.root)
        assert panel is not None
        assert hasattr(panel, 'notebook')
    
    def test_results_display_methods(self):
        """Test methods for displaying results."""
        from hbat.gui.results_panel import ResultsPanel
        
        panel = ResultsPanel(self.root)
        
        # Test that display methods exist
        assert hasattr(panel, 'update_results'), "Should have update_results method"
        assert hasattr(panel, 'clear_results'), "Should have clear_results method"
        
        # Test calling clear_results doesn't raise errors
        panel.clear_results()


@pytest.mark.gui
class TestChainVisualization:
    """Test chain visualization functionality."""
    
    def test_chain_visualization_import(self):
        """Test importing chain visualization components."""
        try:
            from hbat.gui.chain_visualization import ChainVisualizationWindow
            assert ChainVisualizationWindow is not None
        except ImportError:
            pytest.skip("Chain visualization module not available")
    
    def test_chain_visualization_creation(self):
        """Test chain visualization window creation."""
        try:
            import tkinter as tk
            from hbat.gui.chain_visualization import ChainVisualizationWindow
            from hbat.core.app_config import HBATConfig
            from unittest.mock import Mock
            
            # Skip if dependencies not available
            try:
                import networkx as nx
                import matplotlib.pyplot as plt
            except ImportError:
                pytest.skip("Visualization dependencies not available")
            
            root = tk.Tk()
            root.withdraw()
            
            try:
                # Create mock chain
                mock_chain = Mock()
                mock_chain.interactions = []
                mock_chain.chain_length = 0
                mock_chain.chain_type = "test"
                
                config = HBATConfig()
                
                # Mock the Toplevel creation to avoid display issues
                with patch('tkinter.Toplevel'):
                    viz_window = ChainVisualizationWindow(root, mock_chain, "test", config)
                    assert viz_window is not None
                    assert hasattr(viz_window, 'G')  # NetworkX graph
                    
            finally:
                root.quit()
                root.destroy()
                
        except ImportError as e:
            pytest.skip(f"Chain visualization dependencies not available: {e}")


@pytest.mark.gui
class TestGUIImports:
    """Test that all GUI modules can be imported."""
    
    def test_gui_module_imports(self):
        """Test importing all GUI modules."""
        try:
            from hbat.gui.main_window import MainWindow
            from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
            from hbat.gui.pdb_fixing_dialog import PDBFixingDialog
            from hbat.gui.preset_manager_dialog import PresetManagerDialog
            from hbat.gui.results_panel import ResultsPanel
            
            # All imports successful
            assert True
            
        except ImportError as e:
            pytest.skip(f"GUI modules not available: {e}")
    
    def test_renderer_imports(self):
        """Test importing visualization renderers."""
        try:
            from hbat.gui.graphviz_renderer import GraphVizRenderer
            from hbat.gui.matplotlib_renderer import MatplotlibRenderer
            
            assert True
            
        except ImportError as e:
            pytest.skip(f"Renderer modules not available: {e}")