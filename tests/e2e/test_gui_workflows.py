"""
End-to-end GUI workflow tests for HBAT.

Tests complete user workflows using the new modal dialog architecture.
These tests focus on realistic user interactions with the GUI.
"""

import pytest
import tkinter as tk
import tempfile
import os
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path


class BaseGUIWorkflowTest:
    """Base class for GUI workflow tests with common setup."""
    
    def setup_method(self):
        """Set up test environment."""
        try:
            import tkinter as tk
            self.gui_available = True
            self.created_windows = []  # Track windows for cleanup
            
            # Stop any existing async executor to prevent event loop conflicts
            try:
                import tk_async_execute as tae
                tae.stop()
            except Exception:
                pass
        except ImportError:
            self.gui_available = False
            
    def teardown_method(self):
        """Clean up GUI resources."""
        # Stop async executor first to prevent event loop conflicts
        try:
            import tk_async_execute as tae
            tae.stop()
        except Exception:
            pass
            
        # Clean up any tracked windows
        for window in self.created_windows:
            try:
                if hasattr(window, 'root'):
                    window.root.quit()
                    window.root.destroy()
                elif hasattr(window, 'destroy'):
                    window.destroy()
                elif hasattr(window, 'quit'):
                    window.quit()
            except Exception:
                pass
        
        # Force garbage collection
        import gc
        gc.collect()
        
    def track_window(self, window):
        """Track a window for cleanup."""
        self.created_windows.append(window)
        return window


@pytest.mark.gui
class TestGUIDialogWorkflows(BaseGUIWorkflowTest):
    """Test workflows involving the new modal dialogs."""
    
    def test_geometry_cutoffs_dialog_workflow(self):
        """Test complete geometry cutoffs dialog workflow."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        root = self.track_window(tk.Tk())
        root.withdraw()
        
        # Test dialog creation with custom parameters
        test_params = AnalysisParameters(
            hb_distance_cutoff=3.2,
            hb_angle_cutoff=125.0
        )
        
        dialog = GeometryCutoffsDialog(root, test_params)
        
        try:
            # Verify dialog was created with correct parameters
            assert dialog.current_params.hb_distance_cutoff == 3.2
            assert dialog.current_params.hb_angle_cutoff == 125.0
            
            # Test modal behavior
            assert dialog.dialog.transient() is not None
            
            # Test that result is initially None
            assert dialog.result is None
            
        finally:
            dialog.dialog.destroy()
    
    def test_pdb_fixing_dialog_workflow(self):
        """Test PDB fixing dialog workflow."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        from hbat.gui.pdb_fixing_dialog import PDBFixingDialog
        
        root = self.track_window(tk.Tk())
        root.withdraw()
        
        dialog = PDBFixingDialog(root)
        
        try:
            # Test dialog creation
            assert dialog is not None
            assert hasattr(dialog, 'result')
            
            # Test modal behavior
            assert dialog.dialog.transient() is not None
            
        finally:
            dialog.dialog.destroy()
    
    def test_preset_manager_workflow(self):
        """Test preset manager dialog workflow."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        from hbat.gui.preset_manager_dialog import PresetManagerDialog
        from hbat.core.analysis import AnalysisParameters
        
        root = self.track_window(tk.Tk())
        root.withdraw()
        
        test_params = AnalysisParameters(hb_distance_cutoff=3.5)
        dialog = PresetManagerDialog(root, test_params)
        
        try:
            # Test dialog creation
            assert dialog is not None
            assert hasattr(dialog, 'preset_tree')
            assert hasattr(dialog, 'preset_file_paths')
            
            # Test tree functionality
            dialog._refresh_preset_list()
            
        finally:
            dialog.dialog.destroy()


@pytest.mark.gui
class TestMainWindowWorkflows(BaseGUIWorkflowTest):
    """Test main window workflows with new dialog architecture."""
    
    def test_main_window_dialog_integration(self):
        """Test main window integration with new dialogs."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        from hbat.gui.main_window import MainWindow
        
        main_window = None
        try:
            main_window = self.track_window(MainWindow())
            main_window.root.withdraw()
            
            # Test that dialog methods exist
            assert hasattr(main_window, '_open_parameters_window')
            assert hasattr(main_window, '_open_pdb_fixing_window')
            
            # Mock dialogs to test integration without creating modal windows
            with patch('hbat.gui.main_window.GeometryCutoffsDialog') as mock_geo:
                mock_instance = Mock()
                mock_instance.get_result.return_value = None
                mock_geo.return_value = mock_instance
                
                # Test geometry cutoffs dialog integration
                main_window._open_parameters_window()
                mock_geo.assert_called_once()
            
            with patch('hbat.gui.main_window.PDBFixingDialog') as mock_pdb:
                mock_instance = Mock()
                mock_instance.get_result.return_value = None
                mock_pdb.return_value = mock_instance
                
                # Test PDB fixing dialog integration
                main_window._open_pdb_fixing_window()
                mock_pdb.assert_called_once()
                
        finally:
            if main_window and hasattr(main_window, 'root'):
                try:
                    import tk_async_execute as tae
                    tae.stop()  # Stop async executor first
                    main_window.root.quit()
                    main_window.root.destroy()
                except Exception:
                    pass
    
    def test_main_window_session_parameters(self):
        """Test main window session parameter management."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        from hbat.gui.main_window import MainWindow
        from hbat.core.analysis import AnalysisParameters
        
        main_window = None
        try:
            main_window = self.track_window(MainWindow())
            main_window.root.withdraw()
            
            # Test session parameter management
            assert hasattr(main_window, 'session_parameters')
            
            # Test setting session parameters
            test_params = AnalysisParameters(hb_distance_cutoff=3.3)
            main_window.session_parameters = test_params
            
            assert main_window.session_parameters.hb_distance_cutoff == 3.3
            
        finally:
            if main_window and hasattr(main_window, 'root'):
                try:
                    import tk_async_execute as tae
                    tae.stop()  # Stop async executor first
                    main_window.root.quit()
                    main_window.root.destroy()
                except Exception:
                    pass


@pytest.mark.gui
class TestAnalysisWorkflows(BaseGUIWorkflowTest):
    """Test analysis workflows through the GUI."""
    
    def test_analysis_parameter_workflow(self, sample_pdb_file):
        """Test analysis with GUI-configured parameters."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        from hbat.core.analysis import AnalysisParameters, NPMolecularInteractionAnalyzer
        
        # Test that parameters can be created and used for analysis
        gui_params = AnalysisParameters(
            hb_distance_cutoff=3.2,
            hb_angle_cutoff=125.0,
            analysis_mode="complete"
        )
        
        # Test analysis with GUI parameters
        analyzer = NPMolecularInteractionAnalyzer(gui_params)
        success = analyzer.analyze_file(sample_pdb_file)
        
        assert success, "Analysis with GUI parameters should succeed"
        
        # Verify results
        summary = analyzer.get_summary()
        assert summary['hydrogen_bonds']['count'] >= 0
        assert summary['total_interactions'] >= 0
    
    def test_pdb_fixing_analysis_workflow(self, pdb_fixing_test_file):
        """Test analysis workflow with PDB fixing enabled."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        from hbat.core.analysis import AnalysisParameters, NPMolecularInteractionAnalyzer
        
        # Test PDB fixing parameters (would come from PDBFixingDialog)
        fixing_params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True
        )
        
        # Test analysis with PDB fixing
        analyzer = NPMolecularInteractionAnalyzer(fixing_params)
        success = analyzer.analyze_file(pdb_fixing_test_file)
        
        assert success, "Analysis with PDB fixing should succeed"


@pytest.mark.gui
class TestPresetWorkflows(BaseGUIWorkflowTest):
    """Test preset management workflows."""
    
    def test_preset_data_structure(self):
        """Test preset data structure compatibility."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        from hbat.core.analysis import AnalysisParameters
        
        # Test creating parameters that would be saved as presets
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
    
    def test_preset_file_workflow(self):
        """Test preset file save/load workflow simulation."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        import json
        from hbat.core.analysis import AnalysisParameters
        
        # Simulate preset save/load workflow
        original_params = AnalysisParameters(
            hb_distance_cutoff=3.1,
            hb_angle_cutoff=135.0,
            analysis_mode="local"
        )
        
        # Create temporary preset file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            preset_data = {
                "format_version": "1.0",
                "parameters": {
                    "hydrogen_bonds": {
                        "h_a_distance_cutoff": original_params.hb_distance_cutoff,
                        "dha_angle_cutoff": original_params.hb_angle_cutoff,
                        "d_a_distance_cutoff": original_params.hb_donor_acceptor_cutoff
                    },
                    "general": {
                        "analysis_mode": original_params.analysis_mode
                    }
                }
            }
            json.dump(preset_data, f)
            temp_path = f.name
        
        try:
            # Test loading preset data
            with open(temp_path, 'r') as f:
                loaded_data = json.load(f)
            
            # Verify preset structure
            assert "format_version" in loaded_data
            assert "parameters" in loaded_data
            assert loaded_data["parameters"]["hydrogen_bonds"]["h_a_distance_cutoff"] == 3.1
            
        finally:
            os.unlink(temp_path)


@pytest.mark.gui
class TestVisualizationWorkflows(BaseGUIWorkflowTest):
    """Test visualization component workflows."""
    
    def test_results_panel_workflow(self):
        """Test results panel integration workflow."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        from hbat.gui.results_panel import ResultsPanel
        
        root = self.track_window(tk.Tk())
        root.withdraw()
        
        panel = ResultsPanel(root)
        
        # Test basic panel functionality
        assert panel is not None
        assert hasattr(panel, 'notebook')
        assert hasattr(panel, 'clear_results')
        assert hasattr(panel, 'update_results')
        
        # Test clearing results
        panel.clear_results()
    
    def test_chain_visualization_integration(self):
        """Test chain visualization integration."""
        try:
            from hbat.gui.chain_visualization import ChainVisualizationWindow
            from hbat.core.app_config import HBATConfig
            
            # Skip if visualization dependencies not available
            try:
                import networkx as nx
                import matplotlib.pyplot as plt
            except ImportError:
                pytest.skip("Visualization dependencies not available")
            
            # Test can create visualization config
            config = HBATConfig()
            assert config is not None
            
        except ImportError:
            pytest.skip("Chain visualization not available")


@pytest.mark.gui
class TestErrorHandlingWorkflows(BaseGUIWorkflowTest):
    """Test error handling in GUI workflows."""
    
    def test_invalid_file_handling(self):
        """Test handling of invalid files in GUI workflow."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        from hbat.core.analysis import NPMolecularInteractionAnalyzer, AnalysisParameters
        
        # Test analysis with non-existent file
        analyzer = NPMolecularInteractionAnalyzer(AnalysisParameters())
        
        # Should handle invalid file gracefully
        success = analyzer.analyze_file("nonexistent_file.pdb")
        assert not success, "Analysis should fail gracefully for invalid file"
    
    def test_dialog_error_handling(self):
        """Test dialog error handling."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        # Test that dialog creation handles errors gracefully
        root = self.track_window(tk.Tk())
        root.withdraw()
        
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        # Test dialog creation with valid parameters
        dialog = GeometryCutoffsDialog(root, AnalysisParameters())
        
        try:
            assert dialog is not None
            assert dialog.result is None  # Should start with no result
        finally:
            dialog.dialog.destroy()


@pytest.mark.gui
class TestPerformanceWorkflows(BaseGUIWorkflowTest):
    """Test GUI performance characteristics."""
    
    def test_dialog_creation_performance(self):
        """Test that dialogs can be created quickly."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        import time
        from hbat.gui.geometry_cutoffs_dialog import GeometryCutoffsDialog
        from hbat.core.analysis import AnalysisParameters
        
        root = self.track_window(tk.Tk())
        root.withdraw()
        
        # Test dialog creation time
        start_time = time.time()
        dialog = GeometryCutoffsDialog(root, AnalysisParameters())
        creation_time = time.time() - start_time
        
        try:
            # Dialog should be created quickly (under 2 seconds)
            assert creation_time < 2.0, f"Dialog creation took {creation_time:.2f}s"
        finally:
            dialog.dialog.destroy()
    
    def test_window_cleanup_performance(self):
        """Test that windows are cleaned up efficiently."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        import time
        
        # Test creating and destroying multiple windows
        windows = []
        
        start_time = time.time()
        for i in range(5):
            root = tk.Tk()
            root.withdraw()
            windows.append(root)
        
        # Clean up all windows
        for window in windows:
            window.quit()
            window.destroy()
        
        cleanup_time = time.time() - start_time
        
        # Cleanup should be fast
        assert cleanup_time < 1.0, f"Window cleanup took {cleanup_time:.2f}s"