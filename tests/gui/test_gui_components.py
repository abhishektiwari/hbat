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