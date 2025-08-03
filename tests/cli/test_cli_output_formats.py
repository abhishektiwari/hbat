"""
Tests for CLI output format functionality.

Tests the new format detection and multiple file export features.
"""

import pytest
import tempfile
import os
import json
import csv
from pathlib import Path
from hbat.cli.main import (
    detect_output_format,
    export_to_csv,
    export_to_json,
    export_to_csv_files,
    export_to_json_files,
    create_parser,
    run_analysis
)
from hbat.core.analysis import NPMolecularInteractionAnalyzer, AnalysisParameters
from unittest.mock import Mock, patch, MagicMock


@pytest.mark.cli
class TestOutputFormatDetection:
    """Test output format detection functionality."""
    
    def test_detect_txt_format(self):
        """Test detection of text format."""
        assert detect_output_format("results.txt") == "text"
        assert detect_output_format("Results.TXT") == "text"
        assert detect_output_format("/path/to/file.txt") == "text"
    
    def test_detect_csv_format(self):
        """Test detection of CSV format."""
        assert detect_output_format("results.csv") == "csv"
        assert detect_output_format("Results.CSV") == "csv"
        assert detect_output_format("/path/to/file.csv") == "csv"
    
    def test_detect_json_format(self):
        """Test detection of JSON format."""
        assert detect_output_format("results.json") == "json"
        assert detect_output_format("Results.JSON") == "json"
        assert detect_output_format("/path/to/file.json") == "json"
    
    def test_unsupported_format_raises_error(self):
        """Test that unsupported formats raise appropriate error."""
        with pytest.raises(ValueError, match="Unsupported output format"):
            detect_output_format("results.pdf")
        
        with pytest.raises(ValueError, match="Unsupported output format"):
            detect_output_format("results.xml")
        
        with pytest.raises(ValueError, match="Unsupported output format"):
            detect_output_format("results")  # No extension


@pytest.mark.cli
class TestSingleFileExports:
    """Test single file export functionality."""
    
    @pytest.fixture
    def mock_analyzer(self):
        """Create a mock analyzer with sample results."""
        analyzer = Mock(spec=NPMolecularInteractionAnalyzer)
        
        # Mock hydrogen bonds
        hb1 = Mock()
        hb1.donor_residue = "A123GLY"
        hb1.donor.name = "N"
        hb1.hydrogen.name = "H"
        hb1.acceptor_residue = "A124ALA"
        hb1.acceptor.name = "O"
        hb1.distance = 2.8
        hb1.angle = 2.79  # radians (~160 degrees)
        hb1.donor_acceptor_distance = 3.2
        hb1.bond_type = "N-H...O"
        hb1.donor_acceptor_properties = "PBN-PBN"
        hb1.get_backbone_sidechain_interaction = Mock(return_value="B-B")
        
        # Mock halogen bonds
        xb1 = Mock()
        xb1.halogen_residue = "A125TYR"
        xb1.donor_residue = "A125TYR"  # For compatibility
        xb1.halogen.name = "CL"
        xb1.acceptor_residue = "A126ASP"
        xb1.acceptor.name = "OD1"
        xb1.distance = 3.5
        xb1.angle = 2.62  # radians (~150 degrees)
        xb1.bond_type = "C-Cl...O"
        xb1.donor_acceptor_properties = "PSN-PSN"
        xb1.get_backbone_sidechain_interaction = Mock(return_value="S-S")
        
        # Mock pi interactions
        pi1 = Mock()
        pi1.donor_residue = "A127LYS"
        pi1.donor.name = "NZ"
        pi1.hydrogen.name = "HZ1"
        pi1.pi_residue = "A128PHE"
        pi1.distance = 3.8
        pi1.angle = 2.44  # radians (~140 degrees)
        pi1.donor_acceptor_properties = "PSN-PSN"
        pi1.get_backbone_sidechain_interaction = Mock(return_value="S-S")
        pi1.get_interaction_type_display = Mock(return_value="NH-π")
        
        # Mock cooperativity chains
        chain1 = Mock()
        chain1.chain_length = 3
        chain1.chain_type = "H-bond chain"
        
        interaction1 = Mock()
        interaction1.get_donor_residue = Mock(return_value="A123GLY")
        interaction1.get_donor_atom = Mock(return_value=Mock(name="N"))
        chain1.interactions = [interaction1, interaction1, interaction1]
        
        analyzer.hydrogen_bonds = [hb1]
        analyzer.halogen_bonds = [xb1]
        analyzer.pi_interactions = [pi1]
        analyzer.cooperativity_chains = [chain1]
        
        return analyzer
    
    def test_export_to_csv_single_file(self, mock_analyzer):
        """Test exporting to single CSV file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            output_path = f.name
        
        try:
            export_to_csv(mock_analyzer, output_path)
            
            # Verify file exists
            assert os.path.exists(output_path)
            
            # Read and verify content
            with open(output_path, 'r') as f:
                content = f.read()
            
            # Check sections
            assert "# Hydrogen Bonds" in content
            assert "# Halogen Bonds" in content
            assert "# Pi Interactions" in content
            assert "# Cooperativity Chains" in content
            
            # Check headers
            assert "D-A_Properties" in content
            assert "B/S" in content
            
            # Check data
            assert "A123GLY" in content
            assert "A125TYR" in content
            assert "A127LYS" in content
            
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)
    
    def test_export_to_json_single_file(self, mock_analyzer):
        """Test exporting to single JSON file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            output_path = f.name
        
        try:
            # Mock get_summary method
            mock_analyzer.get_summary = Mock(return_value={
                'hydrogen_bonds': {'count': 1},
                'halogen_bonds': {'count': 1},
                'pi_interactions': {'count': 1},
                'cooperativity_chains': {'count': 1},
                'total_interactions': 3
            })
            
            export_to_json(mock_analyzer, "test.pdb", output_path)
            
            # Verify file exists
            assert os.path.exists(output_path)
            
            # Read and verify JSON structure
            with open(output_path, 'r') as f:
                data = json.load(f)
            
            assert 'metadata' in data
            assert 'summary' in data
            assert 'hydrogen_bonds' in data
            assert 'halogen_bonds' in data
            assert 'pi_interactions' in data
            assert 'cooperativity_chains' in data
            
            # Check data content
            assert len(data['hydrogen_bonds']) == 1
            assert data['hydrogen_bonds'][0]['donor_residue'] == "A123GLY"
            
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)


@pytest.mark.cli
class TestMultipleFileExports:
    """Test multiple file export functionality."""
    
    @pytest.fixture
    def mock_analyzer(self):
        """Create a mock analyzer with sample results."""
        analyzer = Mock(spec=NPMolecularInteractionAnalyzer)
        
        # Mock hydrogen bonds
        hb1 = Mock()
        hb1.donor_residue = "A123GLY"
        hb1.donor.name = "N"
        hb1.hydrogen.name = "H"
        hb1.acceptor_residue = "A124ALA"
        hb1.acceptor.name = "O"
        hb1.distance = 2.8
        hb1.angle = 2.79  # radians
        hb1.donor_acceptor_distance = 3.2
        hb1.bond_type = "N-H...O"
        hb1.donor_acceptor_properties = "PBN-PBN"
        hb1.get_backbone_sidechain_interaction = Mock(return_value="B-B")
        
        # Mock halogen bonds with all required attributes
        xb1 = Mock()
        xb1.halogen_residue = "A125TYR"
        xb1.donor_residue = "A125TYR"  # For compatibility
        xb1.halogen.name = "CL"
        xb1.acceptor_residue = "A126ASP"
        xb1.acceptor.name = "OD1"
        xb1.distance = 3.5
        xb1.angle = 2.62  # radians
        xb1.bond_type = "C-Cl...O"
        xb1.donor_acceptor_properties = "PSN-PSN"
        xb1.get_backbone_sidechain_interaction = Mock(return_value="S-S")
        
        # Mock cooperativity chains
        chain1 = Mock()
        chain1.chain_length = 3
        chain1.chain_type = "H-bond chain"
        
        interaction1 = Mock()
        interaction1.get_donor_residue = Mock(return_value="A123GLY")
        interaction1.get_donor_atom = Mock(return_value=Mock(name="N"))
        chain1.interactions = [interaction1]
        
        analyzer.hydrogen_bonds = [hb1]
        analyzer.halogen_bonds = [xb1]
        analyzer.pi_interactions = []
        analyzer.cooperativity_chains = [chain1]
        
        return analyzer
    
    def test_export_to_csv_files(self, mock_analyzer):
        """Test exporting to multiple CSV files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base_path = os.path.join(tmpdir, "results")
            
            export_to_csv_files(mock_analyzer, base_path)
            
            # Check that files were created
            assert os.path.exists(os.path.join(tmpdir, "results_h_bonds.csv"))
            assert os.path.exists(os.path.join(tmpdir, "results_x_bonds.csv"))
            assert os.path.exists(os.path.join(tmpdir, "results_cooperativity_chains.csv"))
            # No pi interactions file should be created since list is empty
            assert not os.path.exists(os.path.join(tmpdir, "results_pi_interactions.csv"))
            
            # Verify hydrogen bonds CSV content
            with open(os.path.join(tmpdir, "results_h_bonds.csv"), 'r') as f:
                reader = csv.reader(f)
                rows = list(reader)
                
                # Check header
                assert rows[0] == [
                    "Donor_Residue", "Donor_Atom", "Hydrogen_Atom",
                    "Acceptor_Residue", "Acceptor_Atom", "Distance_Angstrom",
                    "Angle_Degrees", "Donor_Acceptor_Distance_Angstrom",
                    "Bond_Type", "B/S_Interaction", "D-A_Properties"
                ]
                
                # Check data
                assert len(rows) == 2  # Header + 1 data row
                assert rows[1][0] == "A123GLY"
                assert rows[1][9] == "B-B"
                assert rows[1][10] == "PBN-PBN"
            
            # Verify halogen bonds CSV has D-A Properties and B/S columns
            with open(os.path.join(tmpdir, "results_x_bonds.csv"), 'r') as f:
                reader = csv.reader(f)
                rows = list(reader)
                
                # Check header includes new columns
                assert "D-A_Properties" in rows[0]
                assert "B/S_Interaction" in rows[0]
                
                # Check data
                assert len(rows) == 2  # Header + 1 data row
                assert "PSN-PSN" in rows[1]
                assert "S-S" in rows[1]
    
    def test_export_to_json_files(self, mock_analyzer):
        """Test exporting to multiple JSON files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base_path = os.path.join(tmpdir, "results")
            
            export_to_json_files(mock_analyzer, base_path, "test.pdb")
            
            # Check that files were created
            assert os.path.exists(os.path.join(tmpdir, "results_h_bonds.json"))
            assert os.path.exists(os.path.join(tmpdir, "results_x_bonds.json"))
            assert os.path.exists(os.path.join(tmpdir, "results_cooperativity_chains.json"))
            
            # Verify hydrogen bonds JSON content
            with open(os.path.join(tmpdir, "results_h_bonds.json"), 'r') as f:
                data = json.load(f)
                
                assert 'metadata' in data
                assert data['metadata']['interaction_type'] == "Hydrogen Bonds"
                assert 'interactions' in data
                assert len(data['interactions']) == 1
                assert data['interactions'][0]['donor_residue'] == "A123GLY"
                assert 'backbone_sidechain_interaction' in data['interactions'][0]
                assert 'donor_acceptor_properties' in data['interactions'][0]
            
            # Verify cooperativity chains JSON
            with open(os.path.join(tmpdir, "results_cooperativity_chains.json"), 'r') as f:
                data = json.load(f)
                
                assert 'chains' in data
                assert len(data['chains']) == 1
                assert data['chains'][0]['chain_length'] == 3


@pytest.mark.cli
@pytest.mark.integration
class TestCLIOutputIntegration:
    """Test CLI output integration with run_analysis function."""
    
    def test_output_format_detection_in_cli(self, sample_pdb_file):
        """Test that -o option respects file extensions."""
        parser = create_parser()
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Test text output
            txt_path = os.path.join(tmpdir, "results.txt")
            args = parser.parse_args([sample_pdb_file, "-o", txt_path, "-q"])
            
            with patch('hbat.cli.main.NPMolecularInteractionAnalyzer') as mock_analyzer_class:
                mock_analyzer = MagicMock()
                mock_analyzer.analyze_file.return_value = True
                mock_analyzer.hydrogen_bonds = []
                mock_analyzer.halogen_bonds = []
                mock_analyzer.pi_interactions = []
                mock_analyzer.cooperativity_chains = []
                mock_analyzer.get_summary.return_value = {
                    'hydrogen_bonds': {'count': 0},
                    'halogen_bonds': {'count': 0},
                    'pi_interactions': {'count': 0},
                    'cooperativity_chains': {'count': 0},
                    'total_interactions': 0
                }
                mock_analyzer_class.return_value = mock_analyzer
                
                result = run_analysis(args)
                assert result == 0
                assert os.path.exists(txt_path)
    
    def test_unsupported_format_error(self, sample_pdb_file):
        """Test that unsupported formats raise appropriate errors."""
        parser = create_parser()
        
        with tempfile.TemporaryDirectory() as tmpdir:
            pdf_path = os.path.join(tmpdir, "results.pdf")
            args = parser.parse_args([sample_pdb_file, "-o", pdf_path, "-q"])
            
            with patch('hbat.cli.main.NPMolecularInteractionAnalyzer') as mock_analyzer_class:
                mock_analyzer = MagicMock()
                mock_analyzer.analyze_file.return_value = True
                mock_analyzer.get_summary.return_value = {
                    'hydrogen_bonds': {'count': 0},
                    'halogen_bonds': {'count': 0},
                    'pi_interactions': {'count': 0},
                    'cooperativity_chains': {'count': 0},
                    'total_interactions': 0
                }
                mock_analyzer_class.return_value = mock_analyzer
                
                result = run_analysis(args)
                assert result == 1  # Should fail with unsupported format