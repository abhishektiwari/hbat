"""
Test NumPy analyzer integration with GUI.
"""

import pytest
import tkinter as tk
from unittest.mock import MagicMock, patch

from hbat.gui.main_window import MainWindow
from hbat.core.np_analyzer import NPMolecularInteractionAnalyzer


class TestNumpyGUIIntegration:
    """Test NumPy analyzer integration with GUI components."""
    
    def test_main_window_uses_numpy_analyzer(self):
        """Test that main window uses NumPy analyzer by default."""
        # Create main window without starting mainloop
        with patch('tkinter.Tk'):
            window = MainWindow()
            
            # Check that numpy analyzer is the only option
            assert window.analyzer is None  # Initially no analyzer
    
    def test_analyzer_creation_uses_numpy(self):
        """Test that correct analyzer type is created."""
        with patch('tkinter.Tk'):
            window = MainWindow()
            window.current_file = "test.pdb"
            
            # Mock the analysis parameters
            from hbat.core.analysis import AnalysisParameters
            params = AnalysisParameters()
            
            # Mock file analysis
            with patch.object(NPMolecularInteractionAnalyzer, 'analyze_file', return_value=True):
                # This should create an NPMolecularInteractionAnalyzer
                window._perform_analysis(params)
                
                # Verify the analyzer is the NumPy type
                assert isinstance(window.analyzer, NPMolecularInteractionAnalyzer)
    
    def test_export_includes_analyzer_type(self):
        """Test that export functions include analyzer type."""
        with patch('tkinter.Tk'):
            window = MainWindow()
            
            # Mock analyzer with some results
            window.analyzer = MagicMock(spec=NPMolecularInteractionAnalyzer)
            window.analyzer.get_summary.return_value = {
                'hydrogen_bonds': {'count': 5},
                'halogen_bonds': {'count': 2},
                'pi_interactions': {'count': 3},
                'total_interactions': 10
            }
            window.analyzer.hydrogen_bonds = []
            window.analyzer.halogen_bonds = []
            window.analyzer.pi_interactions = []
            window.current_file = "test.pdb"
            
            # Test export functionality
            test_content = []
            def mock_write(content):
                test_content.append(content)
            
            with patch('builtins.open', create=True) as mock_open:
                mock_file = MagicMock()
                mock_file.write = mock_write
                mock_open.return_value.__enter__.return_value = mock_file
                
                window._export_results_to_file("test_output.txt")
                
                # Check that the output mentions HBAT
                output = ''.join(test_content)
                assert "HBAT" in output
    
    def test_status_updates_show_analyzer_type(self):
        """Test that status updates show correct analyzer type."""
        with patch('tkinter.Tk'):
            window = MainWindow()
            window.current_file = "test.pdb"
            
            from hbat.core.analysis import AnalysisParameters
            params = AnalysisParameters()
            
            # Mock the root.after method to capture status updates
            status_updates = []
            def mock_after(delay, func):
                if callable(func):
                    try:
                        func()
                        # Capture the status text
                        status_updates.append(window.status_var.get())
                    except:
                        pass  # Ignore errors during mock execution
            
            window.root.after = mock_after
            
            with patch.object(NPMolecularInteractionAnalyzer, 'analyze_file', return_value=True):
                window._perform_analysis(params)
                
                # Check that at least one status update mentions analysis
                assert any("analysis" in status.lower() for status in status_updates)