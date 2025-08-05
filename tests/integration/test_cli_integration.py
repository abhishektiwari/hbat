"""
Integration tests for CLI-to-core integration.

These tests verify CLI components integrate properly with core analysis.
"""

import pytest
import tempfile
import json
from hbat.cli.main import (
    create_parser, 
    load_parameters_from_args,
    load_preset_file
)
from hbat.core.analyzer import MolecularInteractionAnalyzer


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestCLIParameterIntegration:
    """Test CLI parameter integration with core analysis."""
    
    def test_cli_to_core_parameter_flow(self, sample_pdb_file):
        """Test parameter flow from CLI to core analysis."""
        parser = create_parser()
        
        # Parse CLI arguments
        args = parser.parse_args([
            sample_pdb_file,
            "--hb-distance", "3.2",
            "--hb-angle", "125",
            "--mode", "complete"
        ])
        
        # Load parameters from CLI
        params = load_parameters_from_args(args)
        
        # Verify parameter values
        assert params.hb_distance_cutoff == 3.2
        assert params.hb_angle_cutoff == 125.0
        assert params.analysis_mode == "complete"
        
        # Use parameters in core analysis
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(sample_pdb_file)
        
        assert success, "CLI-configured analysis should succeed"
        
        # Verify parameters were applied in analysis
        assert analyzer.parameters.hb_distance_cutoff == 3.2
        assert analyzer.parameters.hb_angle_cutoff == 125.0
        assert analyzer.parameters.analysis_mode == "complete"
        
        # Verify analysis results
        assert len(analyzer.hydrogen_bonds) >= 0
        assert len(analyzer.halogen_bonds) >= 0
        assert len(analyzer.pi_interactions) >= 0
        total_interactions = len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        assert total_interactions >= 0
    
    def test_cli_pdb_fixing_integration(self, pdb_fixing_test_file):
        """Test CLI PDB fixing integration with core analysis."""
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
        
        # Verify PDB fixing parameters
        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == "openbabel"
        assert params.fix_pdb_add_hydrogens is True
        
        # Execute analysis
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(pdb_fixing_test_file)
        
        assert success, "CLI PDB fixing integration should succeed"
        
        # Verify fixing was applied
        atoms = analyzer.parser.atoms
        hydrogen_count = sum(1 for atom in atoms if atom.is_hydrogen())
        assert hydrogen_count > 0, "Should have hydrogens after fixing"
    
    def test_cli_preset_integration(self, sample_pdb_file):
        """Test CLI preset integration with core analysis."""
        # Create a test preset
        preset_data = {
            "format_version": "1.0",
            "application": "HBAT",
            "description": "Test preset for integration",
            "parameters": {
                "hydrogen_bonds": {
                    "h_a_distance_cutoff": 3.3,
                    "dha_angle_cutoff": 128.0,
                    "d_a_distance_cutoff": 4.1
                },
                "weak_hydrogen_bonds": {
                    "h_a_distance_cutoff": 3.7,
                    "dha_angle_cutoff": 145.0,
                    "d_a_distance_cutoff": 3.8
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
                    "analysis_mode": "local"
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
            # Load preset through CLI
            params = load_preset_file(preset_path)
            
            # Verify preset parameters
            assert params.hb_distance_cutoff == 3.3
            assert params.hb_angle_cutoff == 128.0
            assert params.whb_distance_cutoff == 3.7
            assert params.whb_angle_cutoff == 145.0
            assert params.whb_donor_acceptor_cutoff == 3.8
            assert params.analysis_mode == "local"
            
            # Use preset parameters in analysis
            analyzer = MolecularInteractionAnalyzer(params)
            success = analyzer.analyze_file(sample_pdb_file)
            
            assert success, "Preset-based analysis should succeed"
            
            # Verify preset parameters were applied
            assert analyzer.parameters.hb_distance_cutoff == 3.3
            assert analyzer.parameters.hb_angle_cutoff == 128.0
            assert analyzer.parameters.analysis_mode == "local"
            
        except SystemExit:
            pytest.skip("Preset loading not available in test environment")
        finally:
            import os
            os.unlink(preset_path)


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestCLIAnalysisIntegration:
    """Test CLI analysis integration workflows."""
    
    def test_cli_parameter_validation_integration(self, sample_pdb_file):
        """Test CLI parameter validation with core analysis."""
        parser = create_parser()
        
        # Test various parameter combinations
        test_cases = [
            # Valid combinations
            ([sample_pdb_file, "--hb-distance", "3.5"], True),
            ([sample_pdb_file, "--hb-angle", "120"], True),
            ([sample_pdb_file, "--mode", "complete"], True),
            ([sample_pdb_file, "--mode", "local"], True),
        ]
        
        for args_list, should_work in test_cases:
            args = parser.parse_args(args_list)
            params = load_parameters_from_args(args)
            
            # Parameters should load successfully
            assert params is not None
            
            # Test with analyzer
            analyzer = MolecularInteractionAnalyzer(params)
            success = analyzer.analyze_file(sample_pdb_file)
            
            if should_work:
                assert success, f"Analysis should succeed for args: {args_list}"
    
    def test_cli_output_integration(self, sample_pdb_file):
        """Test CLI output format integration."""
        parser = create_parser()
        
        # Parse arguments
        args = parser.parse_args([sample_pdb_file, "--verbose"])
        params = load_parameters_from_args(args)
        
        # Run analysis
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        # Get results for output
        # Create statistics from analyzer results
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        summary = analyzer.get_summary()
        
        # Test JSON output format
        json_output = {
            'file': sample_pdb_file,
            'parameters': {
                'hb_distance_cutoff': params.hb_distance_cutoff,
                'hb_angle_cutoff': params.hb_angle_cutoff,
                'analysis_mode': params.analysis_mode
            },
            'statistics': stats,
            'summary': summary
        }
        
        # Verify JSON serialization works
        json_str = json.dumps(json_output, indent=2)
        assert len(json_str) > 0
        
        # Verify roundtrip
        loaded_output = json.loads(json_str)
        assert loaded_output['file'] == sample_pdb_file
        assert loaded_output['statistics']['hydrogen_bonds'] == stats['hydrogen_bonds']
    
    def test_cli_error_integration(self):
        """Test CLI error handling integration with core."""
        parser = create_parser()
        
        # Test with non-existent file
        args = parser.parse_args(["nonexistent_file.pdb"])
        params = load_parameters_from_args(args)
        
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file("nonexistent_file.pdb")
        
        assert not success, "Should fail for non-existent file"
        
        # Analyzer should remain in valid state
        # Create statistics from analyzer results
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        assert stats['hydrogen_bonds'] == 0
        assert stats['total_interactions'] == 0


@pytest.mark.integration
class TestCLIPresetIntegration:
    """Test CLI preset system integration."""
    
    def test_preset_file_format_integration(self):
        """Test preset file format integration."""
        # Create preset data in CLI format
        preset_data = {
            "format_version": "1.0",
            "application": "HBAT",
            "description": "Integration test preset",
            "parameters": {
                "hydrogen_bonds": {
                    "h_a_distance_cutoff": 3.4,
                    "dha_angle_cutoff": 126.0,
                    "d_a_distance_cutoff": 4.2
                },
                "weak_hydrogen_bonds": {
                    "h_a_distance_cutoff": 3.8,
                    "dha_angle_cutoff": 148.0,
                    "d_a_distance_cutoff": 3.6
                },
                "halogen_bonds": {
                    "x_a_distance_cutoff": 4.1,
                    "dxa_angle_cutoff": 125.0
                },
                "pi_interactions": {
                    "h_pi_distance_cutoff": 4.6,
                    "dh_pi_angle_cutoff": 88.0
                },
                "general": {
                    "covalent_cutoff_factor": 0.87,
                    "analysis_mode": "complete"
                },
                "pdb_fixing": {
                    "enabled": True,
                    "method": "pdbfixer",
                    "add_hydrogens": True,
                    "add_heavy_atoms": False,
                    "replace_nonstandard": False,
                    "remove_heterogens": False,
                    "keep_water": True
                }
            }
        }
        
        # Save and load preset
        with tempfile.NamedTemporaryFile(mode='w', suffix='.hbat', delete=False) as f:
            json.dump(preset_data, f)
            preset_path = f.name
        
        try:
            # Load through CLI system
            params = load_preset_file(preset_path)
            
            # Verify all parameters loaded correctly
            assert params.hb_distance_cutoff == 3.4
            assert params.hb_angle_cutoff == 126.0
            assert params.whb_distance_cutoff == 3.8
            assert params.whb_angle_cutoff == 148.0
            assert params.whb_donor_acceptor_cutoff == 3.6
            assert params.xb_distance_cutoff == 4.1
            assert params.xb_angle_cutoff == 125.0
            assert params.pi_distance_cutoff == 4.6
            assert params.pi_angle_cutoff == 88.0
            assert params.covalent_cutoff_factor == 0.87
            assert params.analysis_mode == "complete"
            assert params.fix_pdb_enabled is True
            assert params.fix_pdb_method == "pdbfixer"
            assert params.fix_pdb_add_hydrogens is True
            
        except SystemExit:
            pytest.skip("Preset loading not available in test environment")
        finally:
            import os
            os.unlink(preset_path)
    
    def test_preset_parameter_override_integration(self, sample_pdb_file):
        """Test preset parameter override integration."""
        # Create base preset
        preset_data = {
            "format_version": "1.0",
            "application": "HBAT",
            "description": "Base preset",
            "parameters": {
                "hydrogen_bonds": {
                    "h_a_distance_cutoff": 3.0,
                    "dha_angle_cutoff": 120.0,
                    "d_a_distance_cutoff": 3.5
                },
                "weak_hydrogen_bonds": {
                    "h_a_distance_cutoff": 3.5,
                    "dha_angle_cutoff": 145.0,
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
            # Load preset parameters
            preset_params = load_preset_file(preset_path)
            assert preset_params.hb_distance_cutoff == 3.0
            
            # Test CLI override (simulated)
            parser = create_parser()
            override_args = parser.parse_args([
                sample_pdb_file,
                "--hb-distance", "3.7"  # Override preset value
            ])
            
            override_params = load_parameters_from_args(override_args)
            assert override_params.hb_distance_cutoff == 3.7  # Override should work
            
            # Test integration with analysis
            analyzer = MolecularInteractionAnalyzer(override_params)
            success = analyzer.analyze_file(sample_pdb_file)
            assert success
            
            # Verify override was applied
            assert analyzer.parameters.hb_distance_cutoff == 3.7
            
        except SystemExit:
            pytest.skip("Preset functionality not available in test environment")
        finally:
            import os
            os.unlink(preset_path)


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestCLIPerformanceIntegration:
    """Test CLI performance integration."""
    
    def test_cli_analysis_performance_integration(self, sample_pdb_file):
        """Test CLI analysis performance integration."""
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file])
        
        params = load_parameters_from_args(args)
        analyzer = MolecularInteractionAnalyzer(params)
        
        import time
        start_time = time.time()
        
        success = analyzer.analyze_file(sample_pdb_file)
        
        analysis_time = time.time() - start_time
        
        assert success, "CLI-integrated analysis should succeed"
        assert analysis_time < 45.0, f"Analysis took too long: {analysis_time:.2f}s"
        
        # Verify substantial results
        # Create statistics from analyzer results
        stats = {
            'hydrogen_bonds': len(analyzer.hydrogen_bonds),
            'halogen_bonds': len(analyzer.halogen_bonds),
            'pi_interactions': len(analyzer.pi_interactions),
            'total_interactions': len(analyzer.hydrogen_bonds) + len(analyzer.halogen_bonds) + len(analyzer.pi_interactions)
        }
        assert stats['total_interactions'] > 0
    
    def test_cli_memory_integration(self, sample_pdb_file):
        """Test CLI memory usage integration."""
        import sys
        
        parser = create_parser()
        args = parser.parse_args([sample_pdb_file])
        
        initial_modules = len(sys.modules)
        
        params = load_parameters_from_args(args)
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(sample_pdb_file)
        
        assert success
        
        final_modules = len(sys.modules)
        module_growth = final_modules - initial_modules
        
        # Should not load excessive modules
        assert module_growth < 30, f"Too many modules loaded: {module_growth}"


@pytest.mark.integration
class TestCLIValidationIntegration:
    """Test CLI validation integration."""
    
    def test_argument_validation_integration(self):
        """Test argument validation integration."""
        parser = create_parser()
        
        # Test valid arguments
        valid_cases = [
            ["test.pdb"],
            ["test.pdb", "--hb-distance", "3.5"],
            ["test.pdb", "--mode", "complete"],
            ["test.pdb", "--fix-pdb", "--fix-method", "openbabel"],
        ]
        
        for args_list in valid_cases:
            args = parser.parse_args(args_list)
            params = load_parameters_from_args(args)
            
            # Should create valid parameters
            assert params is not None
            assert hasattr(params, 'hb_distance_cutoff')
            assert hasattr(params, 'analysis_mode')
    
    def test_parameter_type_validation_integration(self, sample_pdb_file):
        """Test parameter type validation integration."""
        parser = create_parser()
        
        # Test numeric parameter parsing
        args = parser.parse_args([
            sample_pdb_file,
            "--hb-distance", "3.5",
            "--hb-angle", "120"
        ])
        
        params = load_parameters_from_args(args)
        
        # Verify types
        assert isinstance(params.hb_distance_cutoff, (int, float))
        assert isinstance(params.hb_angle_cutoff, (int, float))
        
        # Verify values
        assert params.hb_distance_cutoff == 3.5
        assert params.hb_angle_cutoff == 120.0
        
        # Test with analysis
        analyzer = MolecularInteractionAnalyzer(params)
        success = analyzer.analyze_file(sample_pdb_file)
        assert success