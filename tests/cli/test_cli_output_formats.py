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
    create_parser,
    run_analysis
)
from hbat.export.results import (
    export_to_csv_files,
    export_to_json_files,
    export_to_json_single_file,
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
    
    def test_detect_csv_format_raises_error(self):
        """Test that CSV format raises error (not supported for single file)."""
        with pytest.raises(ValueError, match="Single CSV file output is not supported"):
            detect_output_format("results.csv")
        with pytest.raises(ValueError, match="Single CSV file output is not supported"):
            detect_output_format("Results.CSV")
        with pytest.raises(ValueError, match="Single CSV file output is not supported"):
            detect_output_format("/path/to/file.csv")
    
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
        
        # Create mock coordinate objects that support to_list()
        def mock_coords(x, y, z):
            coords = Mock()
            coords.to_list.return_value = [x, y, z]
            return coords
        
        # Create mock atom objects with proper attributes
        def mock_atom(name, x=0.0, y=0.0, z=0.0):
            atom = Mock()
            atom.name = name
            atom.coords = mock_coords(x, y, z)
            return atom
        
        # Create hydrogen bond object using a simple namespace class
        class SimpleHBond:
            def __init__(self):
                self.donor_residue = "A123GLY"
                self.donor = mock_atom("N", 1.0, 2.0, 3.0)
                self.hydrogen = mock_atom("H", 1.1, 2.1, 3.1)
                self.acceptor_residue = "A124ALA"
                self.acceptor = mock_atom("O", 4.0, 5.0, 6.0)
                self.distance = 2.8
                self.angle = 2.79  # radians (~160 degrees)
                self.donor_acceptor_distance = 3.2
                self.bond_type = "N-H...O"
                self.donor_acceptor_properties = "PBN-PBN"
            
            def get_backbone_sidechain_interaction(self):
                return "B-B"
        
        hb1 = SimpleHBond()
        
        # Create halogen bond object using a simple namespace class
        class SimpleXBond:
            def __init__(self):
                self.halogen_residue = "A125TYR"
                self.donor_residue = "A125TYR"  # For compatibility
                self.halogen = mock_atom("CL", 7.0, 8.0, 9.0)
                self.donor = mock_atom("C", 6.5, 7.5, 8.5)  # Carbon bonded to halogen
                self.acceptor_residue = "A126ASP"
                self.acceptor = mock_atom("OD1", 10.0, 11.0, 12.0)
                self.distance = 3.5
                self.angle = 2.62  # radians (~150 degrees)
                self.bond_type = "C-Cl...O"
                self.donor_acceptor_properties = "PSN-PSN"

            def get_backbone_sidechain_interaction(self):
                return "S-S"
        
        xb1 = SimpleXBond()
        
        # Create pi interaction object using a simple namespace class
        class SimplePiBond:
            def __init__(self):
                self.donor_residue = "A127LYS"
                self.donor = mock_atom("NZ", 13.0, 14.0, 15.0)
                self.hydrogen = mock_atom("HZ1", 13.1, 14.1, 15.1)
                self.pi_residue = "A128PHE"
                self.pi_center = mock_coords(16.0, 17.0, 18.0)
                self.distance = 3.8
                self.angle = 2.44  # radians (~140 degrees)
                self.donor_acceptor_properties = "PSN-PSN"
            
            def get_backbone_sidechain_interaction(self):
                return "S-S"
            
            def get_interaction_type_display(self):
                return "NH-π"
        
        pi1 = SimplePiBond()
        
        # Create cooperativity chain with simple namespace classes
        class SimpleChainInteraction:
            def __init__(self):
                self.interaction_type = "H-Bond"
                self.distance = 2.8
                self.angle = 2.79
                self.bond_type = "N-H...O"
            
            def get_donor_residue(self):
                return "A123GLY"
            
            def get_acceptor_residue(self):
                return "A124ALA"
            
            def get_donor_atom(self):
                atom = Mock()
                atom.name = "N"
                return atom
            
            def get_acceptor_atom(self):
                atom = Mock()
                atom.name = "O"
                return atom
            
            def get_interaction_type(self):
                return self.interaction_type
        
        class SimpleChain:
            def __init__(self):
                self.chain_length = 3
                self.chain_type = "H-bond chain"
                self.interactions = [SimpleChainInteraction(), SimpleChainInteraction(), SimpleChainInteraction()]
        
        chain1 = SimpleChain()
        
        analyzer.hydrogen_bonds = [hb1]
        analyzer.halogen_bonds = [xb1]
        analyzer.pi_interactions = [pi1]
        analyzer.pi_pi_stacking = []
        analyzer.carbonyl_interactions = []
        analyzer.n_pi_interactions = []
        analyzer.cooperativity_chains = [chain1]

        return analyzer
    
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

            export_to_json_single_file(mock_analyzer, output_path, "test.pdb")

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
        
        # Create mock coordinate objects that support to_list()
        def mock_coords(x, y, z):
            coords = Mock()
            coords.to_list.return_value = [x, y, z]
            return coords
        
        # Create mock atom objects with proper attributes
        def mock_atom(name, x=0.0, y=0.0, z=0.0):
            atom = Mock()
            atom.name = name
            atom.coords = mock_coords(x, y, z)
            return atom
        
        # Create hydrogen bond object using a simple namespace class
        class SimpleHBond:
            def __init__(self):
                self.donor_residue = "A123GLY"
                self.donor = mock_atom("N", 1.0, 2.0, 3.0)
                self.hydrogen = mock_atom("H", 1.1, 2.1, 3.1)
                self.acceptor_residue = "A124ALA"
                self.acceptor = mock_atom("O", 4.0, 5.0, 6.0)
                self.distance = 2.8
                self.angle = 2.79  # radians
                self.donor_acceptor_distance = 3.2
                self.bond_type = "N-H...O"
                self.donor_acceptor_properties = "PBN-PBN"
            
            def get_backbone_sidechain_interaction(self):
                return "B-B"
        
        hb1 = SimpleHBond()
        
        # Create halogen bond object using a simple namespace class
        class SimpleXBond:
            def __init__(self):
                self.halogen_residue = "A125TYR"
                self.donor_residue = "A125TYR"  # For compatibility
                self.halogen = mock_atom("CL", 7.0, 8.0, 9.0)
                self.donor = mock_atom("C", 6.5, 7.5, 8.5)  # Carbon bonded to halogen
                self.acceptor_residue = "A126ASP"
                self.acceptor = mock_atom("OD1", 10.0, 11.0, 12.0)
                self.distance = 3.5
                self.angle = 2.62  # radians
                self.bond_type = "C-Cl...O"
                self.donor_acceptor_properties = "PSN-PSN"

            def get_backbone_sidechain_interaction(self):
                return "S-S"
        
        xb1 = SimpleXBond()
        
        # Create cooperativity chain with simple namespace classes
        class SimpleChainInteraction:
            def __init__(self):
                self.interaction_type = "H-Bond"
                self.distance = 2.8
                self.angle = 2.79
                self.bond_type = "N-H...O"
            
            def get_donor_residue(self):
                return "A123GLY"
            
            def get_acceptor_residue(self):
                return "A124ALA"
            
            def get_donor_atom(self):
                atom = Mock()
                atom.name = "N"
                return atom
            
            def get_acceptor_atom(self):
                atom = Mock()
                atom.name = "O"
                return atom
            
            def get_interaction_type(self):
                return self.interaction_type
        
        class SimpleChain:
            def __init__(self):
                self.chain_length = 3
                self.chain_type = "H-bond chain"
                self.interactions = [SimpleChainInteraction()]
        
        chain1 = SimpleChain()
        
        analyzer.hydrogen_bonds = [hb1]
        analyzer.halogen_bonds = [xb1]
        analyzer.pi_interactions = []
        analyzer.pi_pi_stacking = []
        analyzer.carbonyl_interactions = []
        analyzer.n_pi_interactions = []
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