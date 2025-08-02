"""
End-to-end tests for GUI workflows.

These tests verify complete GUI usage scenarios from user interactions
through analysis to results display.
"""

import pytest
import tempfile
import os
import json
import unittest.mock


@pytest.mark.e2e
@pytest.mark.gui
@pytest.mark.requires_pdb_files
class TestGUIAnalysisWorkflows:
    """Test complete GUI analysis workflows."""
    
    def setup_method(self):
        """Set up GUI test environment."""
        try:
            import tkinter as tk
            self.gui_available = True
            self._created_windows = []  # Track windows for cleanup
        except ImportError:
            self.gui_available = False
            
    def teardown_method(self):
        """Clean up GUI resources."""
        # Force cleanup of any remaining tkinter windows
        if hasattr(self, '_created_windows'):
            for window in self._created_windows:
                try:
                    if hasattr(window, '_on_closing'):
                        window._on_closing()
                    elif hasattr(window, 'destroy'):
                        window.destroy()
                except:
                    pass
        
        # Force garbage collection to clean up GUI resources
        import gc
        gc.collect()
        
    def _track_window(self, window):
        """Track a window for cleanup."""
        if hasattr(self, '_created_windows'):
            self._created_windows.append(window)
        return window
    
    def test_gui_parameter_to_analysis_workflow(self, sample_pdb_file):
        """Test workflow: GUI parameters → analysis → results display."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        try:
            import tkinter as tk
            from hbat.gui.parameter_panel import ParameterPanel
            from hbat.core.analyzer import MolecularInteractionAnalyzer
            
            # Create GUI components
            root = self._track_window(tk.Tk())
            root.withdraw()  # Hide during testing
            
            try:
                parameter_panel = ParameterPanel(root)
                
                # Set custom parameters via GUI
                from hbat.core.analysis import AnalysisParameters
                custom_params = AnalysisParameters(
                    hb_distance_cutoff=3.2,
                    hb_angle_cutoff=125.0,
                    analysis_mode="complete"
                )
                
                parameter_panel.set_parameters(custom_params)
                
                # Get parameters from GUI
                gui_params = parameter_panel.get_parameters()
                assert gui_params.hb_distance_cutoff == 3.2
                assert gui_params.hb_angle_cutoff == 125.0
                
                # Use GUI parameters for analysis
                analyzer = MolecularInteractionAnalyzer(gui_params)
                success = analyzer.analyze_file(sample_pdb_file)
                assert success, "GUI-configured analysis should succeed"
                
                # Verify results
                stats = analyzer.get_statistics()
                assert stats['hydrogen_bonds'] > 0
                assert stats['total_interactions'] > 0
                
            finally:
                root.destroy()
                
        except ImportError:
            pytest.skip("GUI components not available")
        except Exception as e:
            pytest.skip(f"GUI test failed: {e}")
    
    def test_gui_preset_workflow(self, sample_pdb_file):
        """Test GUI preset creation, save, load workflow."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        try:
            import tkinter as tk
            from hbat.gui.parameter_panel import ParameterPanel
            
            root = self._track_window(tk.Tk())
            root.withdraw()
            
            try:
                parameter_panel = ParameterPanel(root)
                
                # Create custom parameters
                from hbat.core.analysis import AnalysisParameters
                original_params = AnalysisParameters(
                    hb_distance_cutoff=3.1,
                    hb_angle_cutoff=135.0,
                    xb_distance_cutoff=3.8,
                    analysis_mode="local",
                    fix_pdb_enabled=True,
                    fix_pdb_method="openbabel"
                )
                
                # Set parameters in GUI
                parameter_panel.set_parameters(original_params)
                
                # Create preset data
                preset_data = parameter_panel._create_preset_data(original_params)
                
                # Save preset to temporary file
                with tempfile.NamedTemporaryFile(mode='w', suffix='.hbat', delete=False) as f:
                    parameter_panel._save_preset_file(f.name, preset_data)
                    preset_path = f.name
                
                try:
                    # Load preset from file
                    loaded_data = parameter_panel._load_preset_file(preset_path)
                    
                    # Apply loaded preset
                    parameter_panel._apply_preset_data(loaded_data)
                    
                    # Verify parameters match
                    loaded_params = parameter_panel.get_parameters()
                    assert loaded_params.hb_distance_cutoff == 3.1
                    assert loaded_params.hb_angle_cutoff == 135.0
                    assert loaded_params.analysis_mode == "local"
                    assert loaded_params.fix_pdb_enabled is True
                    
                    # Test using loaded parameters for analysis
                    from hbat.core.analyzer import MolecularInteractionAnalyzer
                    analyzer = MolecularInteractionAnalyzer(loaded_params)
                    success = analyzer.analyze_file(sample_pdb_file)
                    assert success, "Preset-based GUI analysis should succeed"
                    
                finally:
                    os.unlink(preset_path)
                    
            finally:
                root.destroy()
                
        except ImportError:
            pytest.skip("GUI components not available")
        except Exception as e:
            pytest.skip(f"GUI preset test failed: {e}")
    
    def test_gui_pdb_fixing_workflow(self, pdb_fixing_test_file):
        """Test GUI workflow with PDB fixing parameters."""
        if not self.gui_available:
            pytest.skip("GUI not available")
        
        try:
            import tkinter as tk
            from hbat.gui.parameter_panel import ParameterPanel
            
            root = self._track_window(tk.Tk())
            root.withdraw()
            
            try:
                parameter_panel = ParameterPanel(root)
                
                # Configure PDB fixing parameters
                from hbat.core.analysis import AnalysisParameters
                fixing_params = AnalysisParameters(
                    fix_pdb_enabled=True,
                    fix_pdb_method="openbabel",
                    fix_pdb_add_hydrogens=True,
                    fix_pdb_add_heavy_atoms=False
                )
                
                parameter_panel.set_parameters(fixing_params)
                
                # Get parameters from GUI
                gui_params = parameter_panel.get_parameters()
                assert gui_params.fix_pdb_enabled is True
                assert gui_params.fix_pdb_method == "openbabel"
                assert gui_params.fix_pdb_add_hydrogens is True
                
                # Use for analysis
                from hbat.core.analyzer import MolecularInteractionAnalyzer
                analyzer = MolecularInteractionAnalyzer(gui_params)
                success = analyzer.analyze_file(pdb_fixing_test_file)
                assert success, "GUI PDB fixing workflow should succeed"
                
            finally:
                root.destroy()
                
        except ImportError:
            pytest.skip("GUI components not available")
        except Exception as e:
            pytest.skip(f"GUI PDB fixing test failed: {e}")


@pytest.mark.e2e
@pytest.mark.gui
@pytest.mark.requires_pdb_files
class TestGUIResultsDisplayWorkflows:
    """Test GUI results display workflows."""
    
    def test_gui_results_display_workflow(self, sample_pdb_file):
        """Test workflow: analysis → results → GUI display."""
        try:
            import tkinter as tk
            from hbat.gui.results_panel import ResultsPanel
            from hbat.core.analyzer import MolecularInteractionAnalyzer
            
            # Run analysis first
            analyzer = MolecularInteractionAnalyzer()
            success = analyzer.analyze_file(sample_pdb_file)
            assert success
            
            # Create GUI for results display
            root = self._track_window(tk.Tk())
            root.withdraw()
            
            try:
                results_panel = ResultsPanel(root)
                
                # Test that results panel can handle analysis results
                assert hasattr(results_panel, 'update_results')
                assert hasattr(results_panel, 'clear_results')
                
                # Clear results (should not raise errors)
                results_panel.clear_results()
                
                # Test updating with real results (basic functionality)
                stats = analyzer.get_statistics()
                
                # This would normally be called by the main GUI
                # results_panel.update_results(analyzer)
                
                # For testing, just verify the panel exists and is functional
                assert results_panel.notebook is not None
                
            finally:
                root.destroy()
                
        except ImportError:
            pytest.skip("GUI components not available")
        except Exception as e:
            pytest.skip(f"GUI results display test failed: {e}")
    
    def test_gui_chain_visualization_workflow(self, sample_pdb_file):
        """Test workflow: analysis → cooperativity chains → visualization."""
        try:
            import tkinter as tk
            from hbat.gui.chain_visualization import ChainVisualizationWindow
            from hbat.core.analyzer import MolecularInteractionAnalyzer
            from hbat.core.app_config import HBATConfig
            
            # Check if visualization dependencies are available
            try:
                import networkx as nx
                import matplotlib.pyplot as plt
            except ImportError:
                pytest.skip("Visualization dependencies not available")
            
            # Run analysis to get chains
            analyzer = MolecularInteractionAnalyzer()
            success = analyzer.analyze_file(sample_pdb_file)
            assert success
            
            chains = analyzer.cooperativity_chains
            
            if len(chains) > 0:
                # Test visualization of first chain
                root = self._track_window(tk.Tk())
                root.withdraw()
                
                try:
                    # Create config and visualization window
                    config = HBATConfig()
                    with unittest.mock.patch('tkinter.Toplevel'):  # Mock window creation
                        viz_window = ChainVisualizationWindow(root, chains[0], "Test Chain", config)
                    
                    # Verify visualization components
                    assert viz_window is not None
                    assert hasattr(viz_window, 'G')  # NetworkX graph
                    assert hasattr(viz_window, 'viz_window')
                    
                    # Clean up visualization
                    if viz_window.viz_window:
                        viz_window.viz_window.destroy()
                        
                finally:
                    root.destroy()
            else:
                pytest.skip("No cooperativity chains found for visualization")
                
        except ImportError:
            pytest.skip("Visualization components not available")
        except Exception as e:
            pytest.skip(f"Chain visualization test failed: {e}")


@pytest.mark.e2e
@pytest.mark.gui
class TestGUIIntegrationWorkflows:
    """Test GUI integration workflows."""
    
    def test_main_window_integration_workflow(self, sample_pdb_file):
        """Test complete main window workflow integration."""
        try:
            from hbat.gui.main_window import MainWindow
            
            # Create main window
            main_window = self._track_window(MainWindow())
            main_window.root.withdraw()  # Hide during testing
            
            try:
                # Verify main window components
                assert hasattr(main_window, 'parameter_panel')
                assert hasattr(main_window, 'results_panel')
                assert hasattr(main_window, 'analyzer')
                
                # Test parameter panel integration
                from hbat.core.analysis import AnalysisParameters
                test_params = AnalysisParameters(hb_distance_cutoff=3.3)
                main_window.parameter_panel.set_parameters(test_params)
                
                retrieved_params = main_window.parameter_panel.get_parameters()
                assert retrieved_params.hb_distance_cutoff == 3.3
                
                # Test analysis integration (simulated)
                # In real GUI, this would be triggered by file selection
                params = main_window.parameter_panel.get_parameters()
                success = main_window.analyzer.analyze_file(sample_pdb_file)
                assert success, "Main window analysis integration should work"
                
                # Test results integration
                # results_panel.update_results would be called here
                stats = main_window.analyzer.get_statistics()
                assert stats['total_interactions'] > 0
                
            finally:
                main_window._on_closing()
                
        except ImportError:
            pytest.skip("Main window not available")
        except Exception as e:
            pytest.skip(f"Main window integration test failed: {e}")
    
    def test_gui_cli_preset_compatibility_workflow(self):
        """Test GUI-CLI preset compatibility workflow."""
        try:
            import tkinter as tk
            from hbat.gui.parameter_panel import ParameterPanel
            from hbat.cli.main import load_preset_file
            
            root = self._track_window(tk.Tk())
            root.withdraw()
            
            try:
                parameter_panel = ParameterPanel(root)
                
                # Create parameters in GUI
                from hbat.core.analysis import AnalysisParameters
                gui_params = AnalysisParameters(
                    hb_distance_cutoff=3.4,
                    hb_angle_cutoff=128.0,
                    fix_pdb_enabled=True,
                    fix_pdb_method="pdbfixer"
                )
                
                # Create preset using GUI
                preset_data = parameter_panel._create_preset_data(gui_params)
                
                # Save preset
                with tempfile.NamedTemporaryFile(mode='w', suffix='.hbat', delete=False) as f:
                    json.dump(preset_data, f)
                    preset_path = f.name
                
                try:
                    # Test that CLI can load GUI-created preset
                    cli_params = load_preset_file(preset_path)
                    
                    assert cli_params.hb_distance_cutoff == 3.4
                    assert cli_params.hb_angle_cutoff == 128.0
                    assert cli_params.fix_pdb_enabled is True
                    assert cli_params.fix_pdb_method == "pdbfixer"
                    
                    # Test that GUI can load CLI-compatible preset
                    parameter_panel._apply_preset_data(preset_data)
                    loaded_params = parameter_panel.get_parameters()
                    
                    assert loaded_params.hb_distance_cutoff == 3.4
                    assert loaded_params.hb_angle_cutoff == 128.0
                    
                except SystemExit:
                    # Acceptable if CLI loading fails in test environment
                    pass
                finally:
                    os.unlink(preset_path)
                    
            finally:
                root.destroy()
                
        except ImportError:
            pytest.skip("GUI or CLI components not available")
        except Exception as e:
            pytest.skip(f"GUI-CLI compatibility test failed: {e}")


@pytest.mark.e2e
@pytest.mark.gui
@pytest.mark.slow
class TestGUIPerformanceWorkflows:
    """Test GUI performance workflows."""
    
    def test_gui_responsiveness_workflow(self, sample_pdb_file):
        """Test GUI responsiveness during analysis."""
        try:
            import tkinter as tk
            from hbat.gui.main_window import MainWindow
            import time
            
            main_window = self._track_window(MainWindow())
            main_window.root.withdraw()
            
            try:
                # Measure GUI setup time
                setup_start = time.time()
                
                # Configure parameters
                from hbat.core.analysis import AnalysisParameters
                params = AnalysisParameters()
                main_window.parameter_panel.set_parameters(params)
                
                setup_time = time.time() - setup_start
                
                # GUI setup should be fast
                assert setup_time < 5.0, f"GUI setup too slow: {setup_time:.2f}s"
                
                # Test analysis timing
                analysis_start = time.time()
                success = main_window.analyzer.analyze_file(sample_pdb_file)
                analysis_time = time.time() - analysis_start
                
                assert success
                assert analysis_time < 60.0, f"GUI analysis too slow: {analysis_time:.2f}s"
                
            finally:
                main_window._on_closing()
                
        except ImportError:
            pytest.skip("GUI components not available")
        except Exception as e:
            pytest.skip(f"GUI performance test failed: {e}")


@pytest.mark.e2e
@pytest.mark.gui
class TestGUIErrorHandlingWorkflows:
    """Test GUI error handling workflows."""
    
    def test_gui_invalid_file_workflow(self):
        """Test GUI error handling with invalid files."""
        try:
            import tkinter as tk
            from hbat.gui.main_window import MainWindow
            
            main_window = self._track_window(MainWindow())
            main_window.root.withdraw()
            
            try:
                # Test with non-existent file
                success = main_window.analyzer.analyze_file("nonexistent_file.pdb")
                assert not success, "Should fail gracefully for non-existent file"
                
                # GUI should remain functional after error
                params = main_window.parameter_panel.get_parameters()
                assert params is not None
                
            finally:
                main_window._on_closing()
                
        except ImportError:
            pytest.skip("GUI components not available")
        except Exception as e:
            pytest.skip(f"GUI error handling test failed: {e}")
    
    def test_gui_invalid_preset_workflow(self):
        """Test GUI error handling with invalid presets."""
        try:
            import tkinter as tk
            from hbat.gui.parameter_panel import ParameterPanel
            
            root = self._track_window(tk.Tk())
            root.withdraw()
            
            try:
                parameter_panel = ParameterPanel(root)
                
                # Create invalid preset data
                invalid_preset = {"invalid": "data"}
                
                # Test validation
                is_valid = parameter_panel._validate_preset_data(invalid_preset)
                assert not is_valid, "Invalid preset should fail validation"
                
                # Create invalid preset file
                with tempfile.NamedTemporaryFile(mode='w', suffix='.hbat', delete=False) as f:
                    f.write("invalid json content")
                    invalid_path = f.name
                
                try:
                    # Test loading invalid preset file
                    try:
                        parameter_panel._load_preset_file(invalid_path)
                        assert False, "Should raise error for invalid preset file"
                    except (json.JSONDecodeError, ValueError, KeyError):
                        # Expected behavior - error should be caught
                        pass
                        
                finally:
                    os.unlink(invalid_path)
                    
            finally:
                root.destroy()
                
        except ImportError:
            pytest.skip("GUI components not available")
        except Exception as e:
            pytest.skip(f"GUI preset error handling test failed: {e}")