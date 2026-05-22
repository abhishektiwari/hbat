"""End-to-end CLI workflow tests for CIF format support.

Tests for complete CLI workflows with CIF files, including:
- Basic CIF parsing via CLI
- PDB fixing with CIF input (--fix-pdb flag)
- Export to multiple output formats (JSON, CSV, TXT)
- Format equivalence verification (CIF vs PDB)
- Integration with PDBFixer and OpenBabel methods
"""

import pytest
import tempfile
import os
import json
import shutil
from pathlib import Path

from hbat.cli.main import create_parser, load_parameters_from_args, run_analysis
from hbat.core.analysis import AnalysisParameters, NPMolecularInteractionAnalyzer


@pytest.mark.e2e
class TestCLICIFWorkflows:
    """Test CLI workflows with CIF format files."""

    def test_cli_cif_basic_parsing(self):
        """Test basic CLI workflow with CIF file: parse and analyze without fixing."""
        parser = create_parser()

        # Simulate CLI arguments with CIF input
        args = parser.parse_args([
            'example_pdb_files/6RSA.cif'
        ])

        # Verify arguments were parsed correctly
        assert args.input == 'example_pdb_files/6RSA.cif'
        assert args.fix_pdb is False

    def test_cli_cif_analysis_without_fixing(self):
        """Test CLI workflow with CIF input without any PDB fixing.

        Verifies:
        - CIF file can be analyzed without PDB fixing
        - Results are valid and consistent
        - No temporary fixed files are created
        """
        parser = create_parser()

        # CLI arguments without fixing
        args = parser.parse_args([
            'example_pdb_files/6RSA.cif'
        ])

        # Load parameters from CLI
        params = load_parameters_from_args(args)
        assert params.fix_pdb_enabled is False

        # Execute analysis without fixing
        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file('example_pdb_files/6RSA.cif')
        assert result is True

        # Verify results
        assert len(analyzer.hydrogen_bonds) > 0
        assert analyzer._pdb_fixing_info.get('applied') is False

        # Verify no temporary fixed file was created
        fixed_file = analyzer._pdb_fixing_info.get('fixed_file_path')
        assert fixed_file is None, "No fixed file should be created when fixing is disabled"

    def test_cli_cif_with_pdbfixer(self):
        """Test CLI workflow with CIF input and PDBFixer fixing.

        Verifies:
        - CIF file can be processed with --fix-pdb
        - PDBFixer method is applied
        - Results are valid
        """
        parser = create_parser()

        # CLI arguments with PDBFixer
        args = parser.parse_args([
            'example_pdb_files/6RSA.cif',
            '--fix-pdb',
            '--fix-method', 'pdbfixer'
        ])

        # Load parameters from CLI
        params = load_parameters_from_args(args)
        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == 'pdbfixer'
        assert params.fix_pdb_add_hydrogens is True

        # Execute analysis
        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file('example_pdb_files/6RSA.cif')
        assert result is True

        # Verify results
        assert len(analyzer.hydrogen_bonds) > 0
        assert analyzer._pdb_fixing_info.get('applied') is True

        # Clean up
        fixed_file = analyzer._pdb_fixing_info.get('fixed_file_path')
        if fixed_file and os.path.exists(fixed_file):
            os.unlink(fixed_file)

    def test_cli_cif_with_openbabel(self):
        """Test CLI workflow with CIF input and OpenBabel fixing.

        Verifies:
        - CIF file can be processed with OpenBabel
        - Output is PDB format (as per implementation)
        - Results are valid
        """
        parser = create_parser()

        # CLI arguments with OpenBabel
        args = parser.parse_args([
            'example_pdb_files/6RSA.cif',
            '--fix-pdb',
            '--fix-method', 'openbabel'
        ])

        # Load parameters
        params = load_parameters_from_args(args)
        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == 'openbabel'

        # Execute analysis
        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file('example_pdb_files/6RSA.cif')
        assert result is True

        # Verify results
        assert len(analyzer.hydrogen_bonds) > 0
        assert analyzer._pdb_fixing_info.get('applied') is True

        # OpenBabel outputs PDB even for CIF input
        fixed_file = analyzer._pdb_fixing_info.get('fixed_file_path')
        assert fixed_file.endswith('.pdb'), f"Expected .pdb output, got {fixed_file}"

        # Clean up
        if os.path.exists(fixed_file):
            os.unlink(fixed_file)

    def test_cli_cif_export_json(self):
        """Test CLI workflow with CIF input and JSON export."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, 'results.json')

            parser = create_parser()
            args = parser.parse_args([
                'example_pdb_files/6RSA.cif',
                '-o', output_file
            ])

            # Load parameters and run analysis
            params = load_parameters_from_args(args)
            analyzer = NPMolecularInteractionAnalyzer(params)
            result = analyzer.analyze_file('example_pdb_files/6RSA.cif')
            assert result is True

            # Export results (simulate CLI behavior)
            from hbat.export.results import export_to_json_single_file
            export_to_json_single_file(analyzer, output_file)

            # Verify output file exists and is valid JSON
            assert os.path.exists(output_file), f"Output file not created: {output_file}"

            with open(output_file, 'r') as f:
                data = json.load(f)

            assert 'hydrogen_bonds' in data
            assert len(data['hydrogen_bonds']) > 0

    def test_cli_cif_export_csv(self):
        """Test CLI workflow with CIF input and CSV export."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base_filename = os.path.join(tmpdir, 'results.csv')

            parser = create_parser()
            args = parser.parse_args([
                'example_pdb_files/6RSA.cif',
                '--csv', base_filename
            ])

            # Load parameters and run analysis
            params = load_parameters_from_args(args)
            analyzer = NPMolecularInteractionAnalyzer(params)
            result = analyzer.analyze_file('example_pdb_files/6RSA.cif')
            assert result is True

            # Export results to CSV (simulate CLI behavior)
            from hbat.export.results import export_to_csv_files
            export_to_csv_files(analyzer, base_filename)

            # Verify output directory exists
            assert os.path.exists(tmpdir), f"Output directory not created: {tmpdir}"

            # Check for expected CSV files
            files = os.listdir(tmpdir)
            assert len(files) > 0, "No CSV files exported"

    def test_cli_cif_vs_pdb_equivalence(self):
        """Test CLI workflow equivalence between CIF and PDB inputs."""
        # Analyze CIF
        parser = create_parser()
        args_cif = parser.parse_args([
            'example_pdb_files/6RSA.cif'
        ])
        params_cif = load_parameters_from_args(args_cif)
        analyzer_cif = NPMolecularInteractionAnalyzer(params_cif)
        analyzer_cif.analyze_file('example_pdb_files/6RSA.cif')

        # Analyze PDB
        args_pdb = parser.parse_args([
            'example_pdb_files/6rsa.pdb'
        ])
        params_pdb = load_parameters_from_args(args_pdb)
        analyzer_pdb = NPMolecularInteractionAnalyzer(params_pdb)
        analyzer_pdb.analyze_file('example_pdb_files/6rsa.pdb')

        # Compare results
        hb_cif = len(analyzer_cif.hydrogen_bonds)
        hb_pdb = len(analyzer_pdb.hydrogen_bonds)

        assert hb_cif == hb_pdb, \
            f"H-bond counts differ: CIF={hb_cif}, PDB={hb_pdb}"

    def test_cli_cif_with_custom_parameters(self):
        """Test CLI workflow with CIF input and custom analysis parameters."""
        parser = create_parser()

        args = parser.parse_args([
            'example_pdb_files/6RSA.cif',
            '--hb-distance', '3.5',
            '--hb-angle', '120',
            '--pi-distance', '5.0'
        ])

        # Load and verify custom parameters
        params = load_parameters_from_args(args)
        assert params.hb_distance_cutoff == 3.5
        assert params.hb_angle_cutoff == 120.0
        assert params.pi_distance_cutoff == 5.0

        # Run analysis with custom parameters
        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file('example_pdb_files/6RSA.cif')
        assert result is True
        assert len(analyzer.hydrogen_bonds) > 0

    def test_cli_cif_complete_workflow_with_fixing(self):
        """Test complete end-to-end CLI workflow: CIF → fixing → analysis → export."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, 'results.json')

            # Complete CLI workflow
            parser = create_parser()
            args = parser.parse_args([
                'example_pdb_files/6RSA.cif',
                '--fix-pdb',
                '--fix-method', 'pdbfixer',
                '--hb-distance', '3.5',
                '-o', output_file
            ])

            # Load parameters
            params = load_parameters_from_args(args)
            assert params.fix_pdb_enabled is True
            assert params.fix_pdb_method == 'pdbfixer'
            assert params.hb_distance_cutoff == 3.5

            # Run analysis
            analyzer = NPMolecularInteractionAnalyzer(params)
            result = analyzer.analyze_file('example_pdb_files/6RSA.cif')
            assert result is True

            # Export results
            from hbat.export.results import export_to_json_single_file
            export_to_json_single_file(analyzer, output_file)

            # Verify complete workflow
            assert os.path.exists(output_file), "Output file not created"

            with open(output_file, 'r') as f:
                data = json.load(f)

            assert 'hydrogen_bonds' in data
            assert len(data['hydrogen_bonds']) > 0, "No H-bonds detected in complete workflow"

            # Clean up fixed file
            fixed_file = analyzer._pdb_fixing_info.get('fixed_file_path')
            if fixed_file and os.path.exists(fixed_file):
                os.unlink(fixed_file)

    def test_cli_cif_multi_file_comparison_with_pdbfixer(self):
        """Test CLI workflow comparing multiple PDB files with PDBFixer.

        Analyzes 3 different protein structures (6rsa, 1crn, 4x21) in both
        CIF and PDB formats using PDBFixer method to ensure:
        - CLI workflow handles multiple protein sizes/complexities
        - Results are equivalent between CIF and PDB formats
        - PDBFixer method produces consistent results across different files
        """
        test_files = [
            ('6RSA', 'example_pdb_files/6RSA.cif', 'example_pdb_files/6rsa.pdb'),
            ('1CRN', 'example_pdb_files/1CRN.cif', 'example_pdb_files/1crn.pdb'),
            ('4X21', 'example_pdb_files/4X21.cif', 'example_pdb_files/4x21.pdb'),
        ]

        parser = create_parser()
        results = {}

        for file_id, cif_file, pdb_file in test_files:
            # Test CIF with PDBFixer
            args_cif = parser.parse_args([
                cif_file,
                '--fix-pdb',
                '--fix-method', 'pdbfixer'
            ])
            params_cif = load_parameters_from_args(args_cif)
            analyzer_cif = NPMolecularInteractionAnalyzer(params_cif)
            result_cif = analyzer_cif.analyze_file(cif_file)
            assert result_cif is True, f"CIF analysis failed for {file_id}"

            # Test PDB with PDBFixer
            args_pdb = parser.parse_args([
                pdb_file,
                '--fix-pdb',
                '--fix-method', 'pdbfixer'
            ])
            params_pdb = load_parameters_from_args(args_pdb)
            analyzer_pdb = NPMolecularInteractionAnalyzer(params_pdb)
            result_pdb = analyzer_pdb.analyze_file(pdb_file)
            assert result_pdb is True, f"PDB analysis failed for {file_id}"

            # Store results for comparison
            results[file_id] = {
                'cif_hbonds': len(analyzer_cif.hydrogen_bonds),
                'pdb_hbonds': len(analyzer_pdb.hydrogen_bonds),
                'cif_xbonds': len(analyzer_cif.halogen_bonds),
                'pdb_xbonds': len(analyzer_pdb.halogen_bonds),
                'cif_pi': len(analyzer_cif.pi_interactions),
                'pdb_pi': len(analyzer_pdb.pi_interactions),
            }

            # Verify CIF and PDB results are similar (within tolerance of ±20% for multi-file comparison)
            # Note: PDBFixer uses stochastic hydrogen placement, resulting in variation between runs
            cif_count = results[file_id]['cif_hbonds']
            pdb_count = results[file_id]['pdb_hbonds']
            tolerance = max(cif_count, pdb_count) * 0.20  # 20% tolerance for stochastic variation
            diff = abs(cif_count - pdb_count)
            assert diff <= tolerance, \
                f"H-bond count difference exceeds tolerance for {file_id}: CIF={cif_count}, PDB={pdb_count}, diff={diff}, tolerance={tolerance}"

            # Clean up fixed files
            for analyzer in [analyzer_cif, analyzer_pdb]:
                fixed_file = analyzer._pdb_fixing_info.get('fixed_file_path')
                if fixed_file and os.path.exists(fixed_file):
                    os.unlink(fixed_file)

        # Verify all files produced valid results
        for file_id, data in results.items():
            assert data['cif_hbonds'] > 0, f"No H-bonds detected for {file_id} (CIF)"
            assert data['pdb_hbonds'] > 0, f"No H-bonds detected for {file_id} (PDB)"

        # Log results as formatted table
        import sys
        print("\n" + "="*100, file=sys.stderr)
        print("Multi-file PDBFixer Comparison Results (6rsa, 1crn, 4x21)", file=sys.stderr)
        print("="*100, file=sys.stderr)
        print(f"{'File':<10} {'Metric':<15} {'CIF':<12} {'PDB':<12} {'Diff':<12} {'Status':<10}", file=sys.stderr)
        print("-"*100, file=sys.stderr)

        for file_id, data in results.items():
            cif_hb = data['cif_hbonds']
            pdb_hb = data['pdb_hbonds']
            diff_hb = abs(cif_hb - pdb_hb)
            cif_xb = data['cif_xbonds']
            pdb_xb = data['pdb_xbonds']
            diff_xb = abs(cif_xb - pdb_xb)
            cif_pi = data['cif_pi']
            pdb_pi = data['pdb_pi']
            diff_pi = abs(cif_pi - pdb_pi)

            hb_status = "✓ PASS" if diff_hb <= max(cif_hb, pdb_hb) * 0.20 else "✗ FAIL"

            print(f"{file_id:<10} {'H-bonds':<15} {cif_hb:<12} {pdb_hb:<12} {diff_hb:<12} {hb_status:<10}", file=sys.stderr)
            print(f"{'':10} {'X-bonds':<15} {cif_xb:<12} {pdb_xb:<12} {diff_xb:<12}", file=sys.stderr)
            print(f"{'':10} {'π-interact':<15} {cif_pi:<12} {pdb_pi:<12} {diff_pi:<12}", file=sys.stderr)

        print("="*100, file=sys.stderr)


@pytest.mark.e2e
class TestCLICIFIntegration:
    """Integration tests for CLI with CIF format."""

    def test_cli_cif_with_heteroatoms(self):
        """Test CLI workflow properly handles heteroatoms in CIF files."""
        parser = create_parser()
        args = parser.parse_args([
            'example_pdb_files/6RSA.cif'
        ])

        params = load_parameters_from_args(args)
        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file('example_pdb_files/6RSA.cif')

        # Verify analysis completes despite heteroatoms
        assert result is True
        assert len(analyzer.hydrogen_bonds) > 0

    def test_cli_cif_method_comparison(self):
        """Test CLI workflow to compare PDBFixer vs OpenBabel methods with CIF."""
        # Test with PDBFixer
        parser = create_parser()
        args_pdbfixer = parser.parse_args([
            'example_pdb_files/6RSA.cif',
            '--fix-pdb',
            '--fix-method', 'pdbfixer'
        ])
        params_pdbfixer = load_parameters_from_args(args_pdbfixer)
        analyzer_pdbfixer = NPMolecularInteractionAnalyzer(params_pdbfixer)
        analyzer_pdbfixer.analyze_file('example_pdb_files/6RSA.cif')
        hb_pdbfixer = len(analyzer_pdbfixer.hydrogen_bonds)

        # Test with OpenBabel
        args_openbabel = parser.parse_args([
            'example_pdb_files/6RSA.cif',
            '--fix-pdb',
            '--fix-method', 'openbabel'
        ])
        params_openbabel = load_parameters_from_args(args_openbabel)
        analyzer_openbabel = NPMolecularInteractionAnalyzer(params_openbabel)
        analyzer_openbabel.analyze_file('example_pdb_files/6RSA.cif')
        hb_openbabel = len(analyzer_openbabel.hydrogen_bonds)

        # Results should differ (different algorithms)
        # but both should be > 0
        assert hb_pdbfixer > 0
        assert hb_openbabel > 0
        # OpenBabel typically finds more H-bonds
        assert hb_openbabel > hb_pdbfixer

        # Clean up
        for analyzer in [analyzer_pdbfixer, analyzer_openbabel]:
            fixed_file = analyzer._pdb_fixing_info.get('fixed_file_path')
            if fixed_file and os.path.exists(fixed_file):
                os.unlink(fixed_file)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
