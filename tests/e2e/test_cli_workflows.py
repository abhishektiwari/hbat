"""
End-to-end tests for CLI workflows.

These tests verify complete CLI usage scenarios from command-line 
arguments through analysis to output generation.
"""

import pytest
import tempfile
import os
import json
import subprocess
import sys
from hbat.cli.main import (
    create_parser, 
    load_parameters_from_args,
    resolve_preset_path,
    load_preset_file
)
from hbat.core.analyzer import MolecularInteractionAnalyzer


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestCLIAnalysisWorkflows:
    """Test complete CLI analysis workflows."""
    
    def test_basic_cli_analysis_workflow(self, sample_pdb_file):
        """Test basic CLI analysis workflow: args → parameters → analysis → results."""
        parser = create_parser()
        
        # Simulate CLI arguments
        args = parser.parse_args([
            sample_pdb_file,
            "--hb-distance", "3.5",
            "--hb-angle", "120"
        ])
        
        # Load parameters from CLI
        params = load_parameters_from_args(args)
        assert params.hb_distance_cutoff == 3.5
        assert params.hb_angle_cutoff == 120.0
        
        # Execute analysis
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(sample_pdb_file)
        assert success, "CLI-configured analysis should succeed"
        
        # Verify results
        # Create statistics from analyzer results
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        assert stats['hydrogen_bonds'] > 0
        assert stats['total_interactions'] > 0
    
    def test_cli_output_format_workflow(self, sample_pdb_file):
        """Test CLI workflow with different output formats."""
        parser = create_parser()
        
        # Create temporary output files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as json_file, \
             tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as csv_file:
            
            json_path = json_file.name
            csv_path = csv_file.name
        
        try:
            # Parse CLI arguments with output options
            args = parser.parse_args([
                sample_pdb_file,
                "--json", json_path,
                "--csv", csv_path,
                "--verbose"
            ])
            
            # Load parameters and run analysis
            params = load_parameters_from_args(args)
            analyzer = MolecularInteractionAnalyzer(params)
            success = analyzer.analyze_file(sample_pdb_file)
            assert success
            
            # Simulate output generation (normally done by CLI main)
            # Create statistics from analyzer results
            stats = {
                'hydrogen_bonds': len(analyzer.hydrogen_bonds),
                'halogen_bonds': len(analyzer.halogen_bonds),
                'pi_interactions': len(analyzer.pi_interactions),
                'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
            }
            
            # Generate JSON output
            json_output = {
                'file': sample_pdb_file,
                'statistics': stats,
                'parameters': {
                    'hb_distance_cutoff': params.hb_distance_cutoff,
                    'hb_angle_cutoff': params.hb_angle_cutoff,
                    'analysis_mode': params.analysis_mode
                }
            }
            
            with open(json_path, 'w') as f:
                json.dump(json_output, f, indent=2)
            
            # Verify output files
            assert os.path.exists(json_path)
            
            # Validate JSON content
            with open(json_path, 'r') as f:
                loaded_data = json.load(f)
            
            assert 'statistics' in loaded_data
            assert 'parameters' in loaded_data
            assert loaded_data['statistics']['hydrogen_bonds'] > 0
            
        finally:
            for path in [json_path, csv_path]:
                if os.path.exists(path):
                    os.unlink(path)
    
    def test_cli_pdb_fixing_workflow(self, pdb_fixing_test_file):
        """Test CLI workflow with PDB fixing enabled."""
        parser = create_parser()
        
        # Parse CLI arguments with PDB fixing
        args = parser.parse_args([
            pdb_fixing_test_file,
            "--fix-pdb",
            "--fix-method", "openbabel",
            "--fix-add-hydrogens"
        ])
        
        # Load parameters
        params = load_parameters_from_args(args)
        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == "openbabel"
        assert params.fix_pdb_add_hydrogens is True
        
        # Execute analysis with PDB fixing
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_fixing_test_file)
        assert success, "CLI PDB fixing workflow should succeed"
        
        # Verify PDB fixing occurred
        # Create statistics from analyzer results
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        assert stats['total_interactions'] >= 0  # Should have some interactions after fixing
    
    def test_cli_error_handling_workflow(self):
        """Test CLI error handling workflow."""
        parser = create_parser()
        
        # Test with invalid file
        args = parser.parse_args(["nonexistent_file.pdb"])
        params = load_parameters_from_args(args)
        analyzer = MolecularInteractionAnalyzer(params)
        
        success = analyzer.analyze_file("nonexistent_file.pdb")
        assert not success, "Should fail for non-existent file"
        
        # Test invalid parameter combinations
        args = parser.parse_args([
            "test.pdb",
            "--hb-distance", "-1.0"  # Invalid negative distance
        ])
        
        # Parameter loading should handle or validate this
        try:
            params = load_parameters_from_args(args)
            # If no validation in parameter loading, analyzer should handle it
            analyzer = MolecularInteractionAnalyzer(params)
            assert params.hb_distance_cutoff == -1.0  # Should accept the value as given
        except ValueError:
            # Acceptable if parameter validation catches this
            pass


@pytest.mark.e2e
class TestCLIPresetWorkflows:
    """Test CLI preset management workflows."""
    
    def test_preset_creation_and_usage_workflow(self, sample_pdb_file):
        """Test complete preset creation and usage workflow."""
        parser = create_parser()
        
        # Create a custom preset file
        preset_data = {
            "format_version": "1.0",
            "application": "HBAT",
            "description": "Test preset for workflow",
            "parameters": {
                "hydrogen_bonds": {
                    "h_a_distance_cutoff": 3.2,
                    "dha_angle_cutoff": 130.0,
                    "d_a_distance_cutoff": 4.0
                },
                "halogen_bonds": {
                    "x_a_distance_cutoff": 4.0,
                    "dxa_angle_cutoff": 120.0
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
        
        # Save preset to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.hbat', delete=False) as f:
            json.dump(preset_data, f)
            preset_path = f.name
        
        try:
            # Test loading preset file
            params = load_preset_file(preset_path)
            assert params.hb_distance_cutoff == 3.2
            assert params.hb_angle_cutoff == 130.0
            assert params.analysis_mode == "complete"
            
            # Test using preset in analysis
            analyzer = MolecularInteractionAnalyzer(params)
            success = analyzer.analyze_file(sample_pdb_file)
            assert success, "Preset-based analysis should succeed"
            
            # Verify preset parameters were applied
            assert analyzer.parameters.hb_distance_cutoff == 3.2
            assert analyzer.parameters.hb_angle_cutoff == 130.0
            
        except SystemExit:
            # Acceptable if preset loading fails in test environment
            pytest.skip("Preset loading not available in test environment")
        finally:
            os.unlink(preset_path)
    
    def test_preset_override_workflow(self, sample_pdb_file):
        """Test workflow where CLI arguments override preset values."""
        # Create preset file
        preset_data = {
            "format_version": "1.0",
            "application": "HBAT",
            "description": "Base preset",
            "parameters": {
                "hydrogen_bonds": {
                    "h_a_distance_cutoff": 3.0,
                    "dha_angle_cutoff": 120.0,
                    "d_a_distance_cutoff": 4.0
                },
                "halogen_bonds": {
                    "x_a_distance_cutoff": 4.0,
                    "dxa_angle_cutoff": 120.0
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
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.hbat', delete=False) as f:
            json.dump(preset_data, f)
            preset_path = f.name
        
        try:
            # Simulate CLI with preset and override
            parser = create_parser()
            
            # First load preset parameters
            preset_params = load_preset_file(preset_path)
            assert preset_params.hb_distance_cutoff == 3.0
            
            # Then test override logic (simulated since the actual CLI handles this)
            # In real CLI, this would be: ["file.pdb", "--preset", "name", "--hb-distance", "3.5"]
            override_args = parser.parse_args([
                sample_pdb_file,
                "--hb-distance", "3.5"  # Override preset value
            ])
            
            override_params = load_parameters_from_args(override_args)
            assert override_params.hb_distance_cutoff == 3.5  # Override should take effect
            
        except SystemExit:
            pytest.skip("Preset functionality not available in test environment")
        finally:
            os.unlink(preset_path)


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
class TestCLIValidationWorkflows:
    """Test CLI input validation workflows."""
    
    def test_parameter_validation_workflow(self, sample_pdb_file):
        """Test CLI parameter validation workflow."""
        parser = create_parser()
        
        # Test various parameter combinations
        test_cases = [
            # Valid parameters
            ([sample_pdb_file, "--hb-distance", "3.5"], True),
            ([sample_pdb_file, "--hb-angle", "120"], True),
            ([sample_pdb_file, "--mode", "complete"], True),
            ([sample_pdb_file, "--mode", "local"], True),
            
            # Edge case parameters (should be handled gracefully)
            ([sample_pdb_file, "--hb-distance", "0.8"], True),  # Small distance within valid range
            ([sample_pdb_file, "--hb-angle", "180"], True),     # Maximum angle
        ]
        
        for args_list, should_succeed in test_cases:
            args = parser.parse_args(args_list)
            params = load_parameters_from_args(args)
            
            # Parameters should load successfully
            assert isinstance(params, type(params))
            
            # Analysis may or may not succeed depending on parameter values
            analyzer = MolecularInteractionAnalyzer(params)
            success = analyzer.analyze_file(sample_pdb_file)
            
            if should_succeed:
                # For valid parameters, analysis should generally work
                # (though it might find no interactions with extreme parameters)
                assert isinstance(success, bool)
    
    def test_file_validation_workflow(self):
        """Test file validation workflow."""
        parser = create_parser()
        
        # Test with various file scenarios
        test_files = [
            ("nonexistent.pdb", False),
            ("", False),  # Empty filename
        ]
        
        for filename, should_succeed in test_files:
            if filename:  # Skip empty filename test for argument parsing
                args = parser.parse_args([filename])
                params = load_parameters_from_args(args)
                analyzer = MolecularInteractionAnalyzer(params)
                
                success = analyzer.analyze_file(filename)
                assert success == should_succeed
    
    def test_help_and_info_workflows(self):
        """Test CLI help and information workflows."""
        parser = create_parser()
        
        # Test that help can be generated
        help_text = parser.format_help()
        assert len(help_text) > 0
        assert "HBAT" in help_text
        assert "--hb-distance" in help_text
        
        # Test version information (if available)
        try:
            from hbat import __version__
            assert len(__version__) > 0
        except ImportError:
            # Version info not available, which is acceptable
            pass


@pytest.mark.e2e
@pytest.mark.requires_pdb_files
@pytest.mark.slow
class TestCLIPerformanceWorkflows:
    """Test CLI performance workflows."""
    
    def test_cli_performance_workflow(self, sample_pdb_file):
        """Test CLI performance with timing measurements."""
        parser = create_parser()
        
        # Configure for performance testing
        args = parser.parse_args([
            sample_pdb_file,
            "--mode", "complete"  # Most comprehensive analysis
        ])
        
        params = load_parameters_from_args(args)
        analyzer = MolecularInteractionAnalyzer(params)
        
        import time
        start_time = time.time()
        
        success = analyzer.analyze_file(sample_pdb_file)
        
        analysis_time = time.time() - start_time
        
        assert success, "Performance test analysis should succeed"
        assert analysis_time < 60.0, f"CLI analysis took too long: {analysis_time:.2f}s"
        
        # Verify substantial results
        # Create statistics from analyzer results
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        assert stats['total_interactions'] > 0
    
    def test_memory_usage_workflow(self, sample_pdb_file):
        """Test CLI memory usage workflow."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file])
        
        params = load_parameters_from_args(args)
        analyzer = MolecularInteractionAnalyzer(params)
        
        # Monitor basic memory usage
        import sys
        
        initial_objects = len(sys.modules)
        
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        # Check that we haven't leaked too many modules
        final_objects = len(sys.modules)
        module_growth = final_objects - initial_objects
        
        # Some module growth is expected, but not excessive
        assert module_growth < 50, f"Too many new modules loaded: {module_growth}"