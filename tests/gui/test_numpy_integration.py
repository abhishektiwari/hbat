"""
Test NumPy analyzer integration with GUI.
"""

import pytest
import tkinter as tk
from unittest.mock import MagicMock, patch

from hbat.gui.main_window import MainWindow
from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer
from hbat.core.analyzer import MolecularInteractionAnalyzer


class TestNumpyGUIIntegration:
    """Test NumPy analyzer integration with GUI components."""
    
    def test_main_window_has_numpy_option(self):
        """Test that main window has NumPy analyzer option."""
        # Create main window without starting mainloop
        with patch('tkinter.Tk'):
            window = MainWindow()
            
            # Check that numpy analyzer is enabled by default
            assert window.use_numpy_analyzer is True
            assert hasattr(window, 'numpy_analyzer_var')
    
    def test_analyzer_toggle_functionality(self):
        """Test that analyzer type can be toggled."""
        with patch('tkinter.Tk'), patch('tkinter.messagebox.showinfo'):
            window = MainWindow()
            
            # Initially should be NumPy
            assert window.use_numpy_analyzer is True
            
            # Toggle to standard
            window.numpy_analyzer_var.set(False)
            window._toggle_analyzer_type()
            assert window.use_numpy_analyzer is False
            
            # Toggle back to NumPy
            window.numpy_analyzer_var.set(True)
            window._toggle_analyzer_type()
            assert window.use_numpy_analyzer is True
    
    def test_analyzer_creation_with_numpy(self):
        """Test that correct analyzer type is created."""
        with patch('tkinter.Tk'):
            window = MainWindow()
            window.current_file = "test.pdb"
            window.use_numpy_analyzer = True
            
            # Mock the analyzer classes
            with patch('hbat.gui.main_window.NPMolecularInteractionAnalyzer') as mock_np_analyzer, \
                 patch('hbat.gui.main_window.MolecularInteractionAnalyzer') as mock_std_analyzer:
                
                mock_np_instance = MagicMock()
                mock_np_instance.analyze_file.return_value = True
                mock_np_analyzer.return_value = mock_np_instance
                
                # Run analysis
                from hbat.constants.parameters import AnalysisParameters
                params = AnalysisParameters()
                window._perform_analysis(params)
                
                # Check that NumPy analyzer was created
                mock_np_analyzer.assert_called_once_with(params)
                mock_std_analyzer.assert_not_called()
    
    def test_analyzer_creation_with_standard(self):
        """Test that standard analyzer is created when NumPy is disabled."""
        with patch('tkinter.Tk'):
            window = MainWindow()
            window.current_file = "test.pdb"
            window.use_numpy_analyzer = False
            
            # Mock the analyzer classes
            with patch('hbat.gui.main_window.NPMolecularInteractionAnalyzer') as mock_np_analyzer, \
                 patch('hbat.gui.main_window.MolecularInteractionAnalyzer') as mock_std_analyzer:
                
                mock_std_instance = MagicMock()
                mock_std_instance.analyze_file.return_value = True
                mock_std_analyzer.return_value = mock_std_instance
                
                # Run analysis
                from hbat.constants.parameters import AnalysisParameters
                params = AnalysisParameters()
                window._perform_analysis(params)
                
                # Check that standard analyzer was created
                mock_std_analyzer.assert_called_once_with(params)
                mock_np_analyzer.assert_not_called()
    
    def test_performance_indicator_updates(self):
        """Test that performance indicator updates with analyzer type."""
        with patch('tkinter.Tk'), patch('tkinter.messagebox.showinfo'):
            window = MainWindow()
            
            # Mock the performance label
            window.performance_label = MagicMock()
            
            # Test NumPy mode
            window.numpy_analyzer_var.set(True)
            window._toggle_analyzer_type()
            
            window.performance_label.config.assert_called_with(
                text="⚡ NumPy High-Performance Mode",
                foreground="green"
            )
            
            # Test standard mode
            window.numpy_analyzer_var.set(False)
            window._toggle_analyzer_type()
            
            window.performance_label.config.assert_called_with(
                text="🐌 Standard Mode",
                foreground="orange"
            )
    
    def test_export_includes_analyzer_type(self):
        """Test that exported results include analyzer type information."""
        with patch('tkinter.Tk'):
            window = MainWindow()
            window.current_file = "test.pdb"
            window.use_numpy_analyzer = True
            
            # Mock analyzer with some results
            window.analyzer = MagicMock()
            window.analyzer.get_statistics.return_value = {
                'hydrogen_bonds': 5,
                'halogen_bonds': 2,
                'pi_interactions': 1,
                'total_interactions': 8
            }
            window.analyzer.hydrogen_bonds = []
            window.analyzer.halogen_bonds = []
            window.analyzer.pi_interactions = []
            
            # Test export
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
                temp_filename = f.name
            
            try:
                window._export_results_to_file(temp_filename)
                
                # Check that file contains analyzer type
                with open(temp_filename, 'r') as f:
                    content = f.read()
                    assert "Analysis engine: NumPy-optimized" in content
                    
            finally:
                import os
                os.unlink(temp_filename)
    
    def test_status_updates_show_analyzer_type(self):
        """Test that status updates show which analyzer was used."""
        with patch('tkinter.Tk'):
            window = MainWindow()
            window.status_var = MagicMock()
            window.use_numpy_analyzer = True
            
            # Mock analyzer with statistics
            window.analyzer = MagicMock()
            window.analyzer.get_statistics.return_value = {
                'hydrogen_bonds': 5,
                'halogen_bonds': 2,
                'pi_interactions': 1
            }
            
            # Mock results panel
            window.results_panel = MagicMock()
            
            # Test analysis completion
            with patch('tkinter.messagebox.showinfo'):
                window._analysis_complete()
            
            # Check that status includes "NumPy-optimized"
            status_call = window.status_var.set.call_args[0][0]
            assert "NumPy-optimized analysis complete" in status_call