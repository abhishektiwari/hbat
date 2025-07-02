"""
Tests for PDB fixer functionality using real PDB files.
"""

import pytest
import os
import tempfile

from hbat.core.pdb_fixer import PDBFixer, PDBFixerError
from hbat.core.pdb_parser import PDBParser


def has_openbabel():
    """Check if OpenBabel is available."""
    try:
        import openbabel
        return True
    except ImportError:
        return False


def has_pdbfixer():
    """Check if PDBFixer is available."""
    try:
        import pdbfixer
        try:
            from openmm.app import PDBFile
        except ImportError:
            from simtk.openmm.app import PDBFile
        return True
    except ImportError:
        return False


@pytest.fixture
def ubi_pdb_file():
    """Provide path to 1ubi.pdb file."""
    file_path = os.path.join(os.path.dirname(__file__), "..", "..", "example_pdb_files", "1ubi.pdb")
    if not os.path.exists(file_path):
        pytest.skip("1ubi.pdb file not found")
    return file_path


@pytest.fixture
def ubi_atoms(ubi_pdb_file):
    """Load atoms from 1ubi.pdb file."""
    parser = PDBParser()
    success = parser.parse_file(ubi_pdb_file)
    assert success, "Failed to parse 1ubi.pdb"
    return parser.atoms


class TestPDBFixerBasics:
    """Basic functionality tests."""
    
    def test_fixer_creation(self):
        """Test PDB fixer can be created."""
        fixer = PDBFixer()
        assert fixer is not None
        assert fixer.supported_methods == ["openbabel", "pdbfixer"]
        assert len(fixer.standard_residues) > 0


class TestAddMissingHydrogens:
    """Test adding missing hydrogen atoms."""
    
    def test_ubi_structure_initial_state(self, ubi_atoms):
        """Verify 1ubi.pdb initial state (no hydrogens)."""
        # Should have 683 atoms total (602 protein + 81 water)
        assert len(ubi_atoms) == 683
        
        # Count ATOM vs HETATM records
        atom_records = [atom for atom in ubi_atoms if atom.record_type == "ATOM"]
        hetatm_records = [atom for atom in ubi_atoms if atom.record_type == "HETATM"]
        assert len(atom_records) == 602  # Protein atoms
        assert len(hetatm_records) == 81  # Water molecules
        
        # Should have no hydrogen atoms initially
        hydrogen_atoms = [atom for atom in ubi_atoms if atom.is_hydrogen()]
        assert len(hydrogen_atoms) == 0
        
        # All atoms should be heavy atoms
        heavy_atoms = [atom for atom in ubi_atoms if not atom.is_hydrogen()]
        assert len(heavy_atoms) == 683
    
    @pytest.mark.skipif(not has_openbabel(), reason="OpenBabel not available")
    def test_add_hydrogens_openbabel(self, ubi_atoms):
        """Test adding hydrogens with OpenBabel."""
        fixer = PDBFixer()
        
        # Add hydrogens using OpenBabel
        result_atoms = fixer.add_missing_hydrogens(ubi_atoms, method="openbabel")
        
        # Should have more atoms than original (hydrogens added)
        assert len(result_atoms) > len(ubi_atoms)
        
        # Should now contain hydrogen atoms
        hydrogen_atoms = [atom for atom in result_atoms if atom.is_hydrogen()]
        assert len(hydrogen_atoms) > 0
        
        # Heavy atom count should remain the same or similar
        heavy_atoms = [atom for atom in result_atoms if not atom.is_hydrogen()]
        assert len(heavy_atoms) >= len(ubi_atoms) * 0.9  # Allow some variation
        
        # Total atoms should be heavy atoms + hydrogens
        assert len(result_atoms) == len(heavy_atoms) + len(hydrogen_atoms)
        
        # For 1ubi.pdb, OpenBabel adds ~787 hydrogens (683 → 1470 total atoms)
        assert 750 <= len(hydrogen_atoms) <= 820, f"Expected 750-820 hydrogens for OpenBabel, got {len(hydrogen_atoms)}"
    
    @pytest.mark.skipif(not has_pdbfixer(), reason="PDBFixer not available")
    def test_add_hydrogens_pdbfixer(self, ubi_atoms):
        """Test adding hydrogens with PDBFixer."""
        fixer = PDBFixer()
        
        # Add hydrogens using PDBFixer
        result_atoms = fixer.add_missing_hydrogens(ubi_atoms, method="pdbfixer", pH=7.0)
        
        # Should have more atoms than original (hydrogens added)
        assert len(result_atoms) > len(ubi_atoms)
        
        # Should now contain hydrogen atoms
        hydrogen_atoms = [atom for atom in result_atoms if atom.is_hydrogen()]
        assert len(hydrogen_atoms) > 0
        
        # Heavy atom count should remain the same or similar  
        heavy_atoms = [atom for atom in result_atoms if not atom.is_hydrogen()]
        assert len(heavy_atoms) >= len(ubi_atoms) * 0.9  # Allow some variation
        
        # Total atoms should be heavy atoms + hydrogens
        assert len(result_atoms) == len(heavy_atoms) + len(hydrogen_atoms)
        
        # For 1ubi.pdb, PDBFixer adds ~791 hydrogens (683 → 1474 total atoms)
        assert 750 <= len(hydrogen_atoms) <= 820, f"Expected 750-820 hydrogens for PDBFixer, got {len(hydrogen_atoms)}"


class TestFileOperations:
    """Test file-based operations."""
    
    @pytest.mark.skipif(not (has_openbabel() or has_pdbfixer()), reason="No PDB fixer available")
    def test_fix_structure_file_output_location(self, ubi_pdb_file):
        """Test that fix_structure_file writes to the correct location."""
        fixer = PDBFixer()
        
        # Use a specific output path
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as output_file:
            output_path = output_file.name
        
        try:
            # Remove the temp file so fixer can create it
            os.unlink(output_path)
            
            # Choose method based on availability
            method = "openbabel" if has_openbabel() else "pdbfixer"
            
            # Fix the structure
            result_path = fixer.fix_structure_file(
                ubi_pdb_file, 
                output_path, 
                method=method, 
                overwrite=True
            )
            
            # Should return the path we specified
            assert result_path == output_path
            
            # File should exist at the specified location
            assert os.path.exists(output_path)
            
            # File should not be empty
            assert os.path.getsize(output_path) > 0
            
            # Should be able to parse the output file
            parser = PDBParser()
            success = parser.parse_file(output_path)
            assert success
            assert len(parser.atoms) > 0
            
        finally:
            # Clean up
            if os.path.exists(output_path):
                os.unlink(output_path)
    
    @pytest.mark.skipif(not (has_openbabel() or has_pdbfixer()), reason="No PDB fixer available")
    def test_fix_structure_file_auto_output_path(self, ubi_pdb_file):
        """Test automatic output path generation."""
        fixer = PDBFixer()
        
        # Choose method based on availability
        method = "openbabel" if has_openbabel() else "pdbfixer"
        
        # Fix the structure without specifying output path
        result_path = fixer.fix_structure_file(ubi_pdb_file, method=method, overwrite=True)
        
        try:
            # Should generate an automatic path
            assert result_path is not None
            assert result_path.endswith('.pdb')
            
            # File should exist
            assert os.path.exists(result_path)
            
            # Should be able to parse the output
            parser = PDBParser()
            success = parser.parse_file(result_path)
            assert success
            assert len(parser.atoms) > 0
            
        finally:
            # Clean up
            if os.path.exists(result_path):
                os.unlink(result_path)
    
    def test_fix_structure_file_overwrite_protection(self, ubi_pdb_file):
        """Test that existing files are protected when overwrite=False."""
        fixer = PDBFixer()
        
        # Create a temporary file that already exists
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as existing_file:
            existing_file.write("EXISTING CONTENT")
            existing_path = existing_file.name
        
        try:
            # Should raise error when trying to overwrite without permission
            with pytest.raises(PDBFixerError, match="already exists"):
                fixer.fix_structure_file(ubi_pdb_file, existing_path, overwrite=False)
                
        finally:
            # Clean up
            if os.path.exists(existing_path):
                os.unlink(existing_path)


class TestOtherMethods:
    """Test other PDB fixer methods."""
    
    @pytest.mark.skipif(not has_pdbfixer(), reason="PDBFixer not available")
    def test_add_missing_heavy_atoms(self, ubi_atoms):
        """Test adding missing heavy atoms."""
        fixer = PDBFixer()
        
        # Use only backbone atoms to simulate incomplete structure
        backbone_atoms = [atom for atom in ubi_atoms if atom.name in ["N", "CA", "C", "O"]]
        original_count = len(backbone_atoms)
        assert original_count > 0
        
        # Add missing heavy atoms
        result_atoms = fixer.add_missing_heavy_atoms(backbone_atoms, method="pdbfixer")
        
        # Should have same or more atoms (side chains added)
        assert len(result_atoms) >= original_count
        
        # Should contain more diverse atom names
        original_names = {atom.name for atom in backbone_atoms}
        result_names = {atom.name for atom in result_atoms}
        assert len(result_names) >= len(original_names)
    
    @pytest.mark.skipif(not has_pdbfixer(), reason="PDBFixer not available")
    def test_convert_nonstandard_residues(self, ubi_atoms):
        """Test converting non-standard residues."""
        fixer = PDBFixer()
        
        # This should work even if there are no non-standard residues
        result_atoms = fixer.convert_nonstandard_residues(ubi_atoms)
        
        # Should have similar number of atoms
        assert len(result_atoms) > 0
        assert len(result_atoms) >= len(ubi_atoms) * 0.9  # Allow some variation
        
        # Should contain recognizable residues
        residue_names = {atom.res_name for atom in result_atoms}
        standard_aa = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
                      "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", 
                      "THR", "TRP", "TYR", "VAL"}
        
        # Most residues should be standard amino acids
        common_residues = residue_names.intersection(standard_aa)
        assert len(common_residues) > 0
    
    @pytest.mark.skipif(not has_pdbfixer(), reason="PDBFixer not available")
    def test_remove_heterogens(self, ubi_atoms):
        """Test removing heterogens."""
        fixer = PDBFixer()
        
        # Count original ATOM vs HETATM records
        atom_records = [atom for atom in ubi_atoms if atom.record_type == "ATOM"]
        hetatm_records = [atom for atom in ubi_atoms if atom.record_type == "HETATM"]
        
        # Remove heterogens (keep water)
        result_keep_water = fixer.remove_heterogens(ubi_atoms, keep_water=True)
        
        # Remove heterogens (remove water)
        result_no_water = fixer.remove_heterogens(ubi_atoms, keep_water=False)
        
        # Should keep all ATOM records
        result_atoms_keep = [atom for atom in result_keep_water if atom.record_type == "ATOM"]
        result_atoms_no = [atom for atom in result_no_water if atom.record_type == "ATOM"]
        
        # ATOM records should be preserved
        assert len(result_atoms_keep) >= len(atom_records) * 0.9
        assert len(result_atoms_no) >= len(atom_records) * 0.9
        
        # If there were heterogens, removing them should reduce total count
        if len(hetatm_records) > 0:
            assert len(result_no_water) <= len(result_keep_water)


class TestErrorHandling:
    """Test error handling."""
    
    def test_unsupported_method(self, ubi_atoms):
        """Test error for unsupported method."""
        fixer = PDBFixer()
        
        with pytest.raises(PDBFixerError, match="Unsupported method"):
            fixer.add_missing_hydrogens(ubi_atoms, method="invalid_method")
    
    def test_heavy_atoms_method_restriction(self, ubi_atoms):
        """Test that heavy atoms only works with pdbfixer."""
        fixer = PDBFixer()
        
        with pytest.raises(PDBFixerError, match="only supported with 'pdbfixer'"):
            fixer.add_missing_heavy_atoms(ubi_atoms, method="openbabel")
    
    def test_nonexistent_file(self):
        """Test error for non-existent input file."""
        fixer = PDBFixer()
        
        with pytest.raises(PDBFixerError, match="does not exist"):
            fixer.fix_structure_file("/nonexistent/file.pdb")


class TestHydrogenAnalysis:
    """Test hydrogen analysis functionality."""
    
    def test_hydrogen_info_no_hydrogens(self, ubi_atoms):
        """Test hydrogen analysis on structure without hydrogens."""
        fixer = PDBFixer()
        
        info = fixer.get_missing_hydrogen_info(ubi_atoms)
        
        assert info["total_atoms"] == 683
        assert info["hydrogen_atoms"] == 0
        assert info["heavy_atoms"] == 683
        assert info["hydrogen_percentage"] == 0.0
        assert info["estimated_missing_hydrogens"] > 0
        assert info["has_sufficient_hydrogens"] == False
    
    @pytest.mark.skipif(not (has_openbabel() or has_pdbfixer()), reason="No PDB fixer available")
    def test_hydrogen_info_with_hydrogens(self, ubi_atoms):
        """Test hydrogen analysis after adding hydrogens."""
        fixer = PDBFixer()
        
        # Add hydrogens
        method = "openbabel" if has_openbabel() else "pdbfixer"
        atoms_with_h = fixer.add_missing_hydrogens(ubi_atoms, method=method)
        
        # Analyze the structure with hydrogens
        info = fixer.get_missing_hydrogen_info(atoms_with_h)
        
        # Should have ~1470-1474 total atoms after adding hydrogens
        assert 1450 <= info["total_atoms"] <= 1500, f"Expected 1450-1500 total atoms, got {info['total_atoms']}"
        assert 750 <= info["hydrogen_atoms"] <= 820, f"Expected 750-820 hydrogens, got {info['hydrogen_atoms']}"
        assert info["heavy_atoms"] >= 683 * 0.9  # Allow some variation
        assert info["hydrogen_percentage"] > 50  # Should be >50% hydrogens
        assert info["has_sufficient_hydrogens"] == True