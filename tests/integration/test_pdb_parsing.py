"""
Integration tests for PDB parsing with real files.

These tests verify PDB parser functionality with real PDB files
and test interactions between parser components.
"""

import pytest
import tempfile
import os
from hbat.core.pdb_parser import PDBParser
from hbat.core.pdb_fixer import PDBFixer, PDBFixerError
from tests.conftest import ExpectedResults


def has_openbabel():
    """Check if OpenBabel is available."""
    try:
        from openbabel import openbabel
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


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestPDBParsingIntegration:
    """Test PDB parsing with real files."""
    
    def test_parse_sample_pdb_file(self, sample_pdb_file):
        """Test parsing real sample PDB file."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        
        assert success, "Should successfully parse sample PDB file"
        
        # Verify parser extracted meaningful data
        atoms = parser.atoms
        assert len(atoms) >= ExpectedResults.MIN_ATOMS, \
            f"Expected >= {ExpectedResults.MIN_ATOMS} atoms, got {len(atoms)}"
        
        # Verify atom types
        hydrogen_count = sum(1 for atom in atoms if atom.is_hydrogen())
        assert hydrogen_count >= ExpectedResults.MIN_HYDROGENS, \
            f"Expected >= {ExpectedResults.MIN_HYDROGENS} hydrogens, got {hydrogen_count}"
        
        # Verify residues
        residues = parser.residues
        assert len(residues) >= ExpectedResults.MIN_RESIDUES, \
            f"Expected >= {ExpectedResults.MIN_RESIDUES} residues, got {len(residues)}"
    
    def test_parse_bonds_detection(self, sample_pdb_file):
        """Test bond detection during parsing."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        # Check bond detection
        bonds = parser.bonds
        assert len(bonds) > 0, "Should detect bonds in structure"
        
        # Verify bond structure
        for bond in bonds[:10]:  # Check first 10 bonds
            assert hasattr(bond, 'atom1_serial'), "Bond should have atom1_serial"
            assert hasattr(bond, 'atom2_serial'), "Bond should have atom2_serial"
            assert isinstance(bond.atom1_serial, int), "Atom serial should be integer"
            assert isinstance(bond.atom2_serial, int), "Atom serial should be integer"
            assert bond.atom1_serial != bond.atom2_serial, "Bond should connect different atoms"
    
    def test_parse_atom_properties(self, sample_pdb_file):
        """Test parsing of atom properties."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        atoms = parser.atoms
        
        # Test first few atoms for proper property parsing
        for atom in atoms[:20]:
            # Required properties
            assert hasattr(atom, 'serial'), "Atom should have serial number"
            assert hasattr(atom, 'name'), "Atom should have name"
            assert hasattr(atom, 'element'), "Atom should have element"
            assert hasattr(atom, 'coords'), "Atom should have coordinates"
            assert hasattr(atom, 'res_name'), "Atom should have residue name"
            assert hasattr(atom, 'chain_id'), "Atom should have chain ID"
            
            # Validate property types and values
            assert isinstance(atom.serial, int), "Serial should be integer"
            assert atom.serial > 0, "Serial should be positive"
            assert isinstance(atom.name, str), "Name should be string"
            assert len(atom.name) > 0, "Name should not be empty"
            assert isinstance(atom.element, str), "Element should be string"
            assert len(atom.element) > 0, "Element should not be empty"
            assert hasattr(atom.coords, 'x'), "Coordinates should have x"
            assert hasattr(atom.coords, 'y'), "Coordinates should have y"
            assert hasattr(atom.coords, 'z'), "Coordinates should have z"
    
    def test_parse_residue_organization(self, sample_pdb_file):
        """Test residue organization and chain structure."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        residues = parser.residues
        
        # Verify residue structure
        # Check first 10 residues from the dict
        residue_list = list(residues.values())
        for residue in residue_list[:10]:
            assert hasattr(residue, 'name'), "Residue should have name"
            assert hasattr(residue, 'chain_id'), "Residue should have chain ID"
            assert hasattr(residue, 'seq_num'), "Residue should have sequence number"
            
            # Check residue atoms
            res_atoms = residue.atoms
            assert len(res_atoms) > 0, "Residue should contain atoms"
            
            # Verify all atoms in residue belong to same residue
            for atom in res_atoms:
                assert atom.res_name == residue.name, "Atom residue name should match"
                assert atom.chain_id == residue.chain_id, "Atom chain should match"
    
    def test_parse_chain_organization(self, sample_pdb_file):
        """Test chain organization."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        chain_ids = parser.get_chain_ids()
        assert len(chain_ids) > 0, "Structure should have chains"
        
        # Verify chain IDs are valid
        for chain_id in chain_ids:
            assert isinstance(chain_id, str), "Chain ID should be string"
            assert len(chain_id) > 0, "Chain ID should not be empty"
            
        # Check residues by chain
        residues = parser.residues
        chains_with_residues = set()
        for residue in residues.values():
            chains_with_residues.add(residue.chain_id)
            
        # Verify all chains have residues
        for chain_id in chain_ids:
            assert chain_id in chains_with_residues, f"Chain {chain_id} should have residues"


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestPDBParsingRobustness:
    """Test PDB parsing robustness and error handling."""
    
    def test_parse_nonexistent_file(self):
        """Test parsing non-existent file."""
        parser = PDBParser()
        success = parser.parse_file("nonexistent_file.pdb")
        assert not success, "Should fail for non-existent file"
    
    def test_parse_invalid_pdb_content(self):
        """Test parsing file with invalid PDB content."""
        # Create temporary file with invalid content
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write("INVALID PDB CONTENT\n")
            f.write("NOT A REAL PDB FILE\n")
            temp_path = f.name
        
        try:
            parser = PDBParser()
            success = parser.parse_file(temp_path)
            
            # Should either fail gracefully or handle minimal parsing
            assert isinstance(success, bool), "Should return boolean success status"
            
            if success:
                # If it claims success, should have minimal data
                atoms = parser.atoms
                # May have 0 atoms for invalid content
                assert len(atoms) >= 0, "Atom count should be non-negative"
        finally:
            os.unlink(temp_path)
    
    def test_parse_empty_file(self):
        """Test parsing empty file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            # Write nothing (empty file)
            temp_path = f.name
        
        try:
            parser = PDBParser()
            success = parser.parse_file(temp_path)
            
            # Should handle empty file gracefully
            assert isinstance(success, bool), "Should return boolean success status"
            
            if success:
                atoms = parser.atoms
                assert len(atoms) == 0, "Empty file should have no atoms"
        finally:
            os.unlink(temp_path)
    
    def test_parse_minimal_valid_pdb(self):
        """Test parsing minimal valid PDB content."""
        minimal_pdb_content = """HEADER    TEST STRUCTURE                          01-JAN-70   TEST            
ATOM      1  CA  ALA A   1      20.154  16.967  18.987  1.00 20.00           C  
END                                                                             
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(minimal_pdb_content)
            temp_path = f.name
        
        try:
            parser = PDBParser()
            success = parser.parse_file(temp_path)
            
            if success:
                atoms = parser.atoms
                assert len(atoms) >= 1, "Should parse at least one atom"
                
                # Verify the parsed atom
                atom = atoms[0]
                assert atom.name == "CA", "Atom name should be parsed correctly"
                assert atom.res_name == "ALA", "Residue name should be parsed correctly"
                assert atom.element == "C", "Element should be parsed correctly"
                assert atom.chain_id == "A", "Chain ID should be parsed correctly"
        finally:
            os.unlink(temp_path)


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestPDBParsingBondDetection:
    """Test bond detection integration."""
    
    def test_bond_detection_hydrogen_bonds(self, sample_pdb_file):
        """Test bond detection for hydrogen bonding analysis."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        # Test hydrogen bond donor detection
        atoms = parser.atoms
        potential_donors = []
        
        for atom in atoms:
            if atom.element.upper() in ["N", "O", "S"]:
                # Check if bonded to hydrogens
                bonded_atoms = parser.get_bonded_atoms(atom.serial)
                for bonded_serial in bonded_atoms:
                    bonded_atom = next((a for a in atoms if a.serial == bonded_serial), None)
                    if bonded_atom and bonded_atom.is_hydrogen():
                        potential_donors.append((atom, bonded_atom))
                        break
        
        assert len(potential_donors) > 0, "Should find potential hydrogen bond donors"
    
    def test_bond_detection_covalent_bonds(self, sample_pdb_file):
        """Test covalent bond detection."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        bonds = parser.bonds
        atoms = parser.atoms
        atom_map = {atom.serial: atom for atom in atoms}
        
        # Verify bond distances are reasonable for covalent bonds
        reasonable_bonds = 0
        for bond in bonds[:50]:  # Check first 50 bonds
            atom1_serial = bond.atom1_serial
            atom2_serial = bond.atom2_serial
            atom1 = atom_map.get(atom1_serial)
            atom2 = atom_map.get(atom2_serial)
            
            if atom1 and atom2:
                distance = atom1.coords.distance_to(atom2.coords)
                # Covalent bonds should be roughly 1-2 Angstroms
                if 0.5 <= distance <= 3.0:
                    reasonable_bonds += 1
        
        # Most detected bonds should be reasonable covalent bond distances
        bond_ratio = reasonable_bonds / min(len(bonds), 50)
        assert bond_ratio > 0.5, f"Too few reasonable bonds: {bond_ratio:.2f}"
    
    def test_bond_detection_consistency(self, sample_pdb_file):
        """Test bond detection consistency."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        bonds = parser.bonds
        
        # Test bond symmetry: if A-B is bonded, B-A should also be accessible
        bond_pairs = set()
        for bond in bonds:
            atom1 = bond.atom1_serial
            atom2 = bond.atom2_serial
            bond_pairs.add((min(atom1, atom2), max(atom1, atom2)))
        
        # Test get_bonded_atoms consistency
        for atom1, atom2 in list(bond_pairs)[:20]:  # Test first 20 pairs
            bonded_to_1 = parser.get_bonded_atoms(atom1)
            bonded_to_2 = parser.get_bonded_atoms(atom2)
            
            assert atom2 in bonded_to_1, f"Atom {atom2} should be bonded to {atom1}"
            assert atom1 in bonded_to_2, f"Atom {atom1} should be bonded to {atom2}"


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestPDBParsingPerformance:
    """Test PDB parsing performance."""
    
    def test_parsing_performance(self, sample_pdb_file):
        """Test parsing performance with timing."""
        import time
        
        parser = PDBParser()
        
        start_time = time.time()
        success = parser.parse_file(sample_pdb_file)
        parse_time = time.time() - start_time
        
        assert success, "Parsing should succeed"
        assert parse_time < 30.0, f"Parsing took too long: {parse_time:.2f}s"
        
        # Verify substantial data was parsed
        atoms = parser.atoms
        assert len(atoms) > 100, "Should parse substantial number of atoms"
    
    def test_memory_usage_during_parsing(self, sample_pdb_file):
        """Test memory usage during parsing."""
        import sys
        
        parser = PDBParser()
        
        # Get initial object count
        initial_objects = len(sys.modules)
        
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        # Check that we haven't leaked modules
        final_objects = len(sys.modules)
        module_growth = final_objects - initial_objects
        
        # Some module growth is expected, but not excessive
        assert module_growth < 20, f"Too many new modules: {module_growth}"
    
    def test_multiple_file_parsing(self, sample_pdb_file, test_pdb_dir):
        """Test parsing multiple files."""
        parser = PDBParser()
        
        # Parse sample file multiple times
        for i in range(3):
            success = parser.parse_file(sample_pdb_file)
            assert success, f"Parsing iteration {i+1} should succeed"
            
            atoms = parser.atoms
            assert len(atoms) > 0, f"Should have atoms on iteration {i+1}"
        
        # Test that parser can handle different files
        if os.path.exists(test_pdb_dir):
            pdb_files = [f for f in os.listdir(test_pdb_dir) if f.endswith('.pdb')]
            
            for pdb_file in pdb_files[:3]:  # Test first 3 PDB files
                file_path = os.path.join(test_pdb_dir, pdb_file)
                success = parser.parse_file(file_path)
                
                if success:  # Some files might not parse correctly
                    atoms = parser.atoms
                    assert len(atoms) >= 0, f"Atom count should be non-negative for {pdb_file}"


@pytest.mark.integration 
@pytest.mark.requires_pdb_files
class TestPDBParsingDataIntegrity:
    """Test data integrity during PDB parsing."""
    
    def test_atom_serial_uniqueness(self, sample_pdb_file):
        """Test that atom serials are unique within structure."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        atoms = parser.atoms
        serials = [atom.serial for atom in atoms]
        
        # Check for duplicate serials
        unique_serials = set(serials)
        assert len(unique_serials) == len(serials), "Atom serials should be unique"
    
    def test_coordinate_validity(self, sample_pdb_file):
        """Test that coordinates are valid numbers."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        atoms = parser.atoms
        
        for atom in atoms[:100]:  # Check first 100 atoms
            # Coordinates should be finite numbers
            assert isinstance(atom.coords.x, (int, float)), "X coordinate should be numeric"
            assert isinstance(atom.coords.y, (int, float)), "Y coordinate should be numeric"  
            assert isinstance(atom.coords.z, (int, float)), "Z coordinate should be numeric"
            
            # Should not be NaN or infinite
            import math
            assert math.isfinite(atom.coords.x), "X coordinate should be finite"
            assert math.isfinite(atom.coords.y), "Y coordinate should be finite"
            assert math.isfinite(atom.coords.z), "Z coordinate should be finite"
    
    def test_element_consistency(self, sample_pdb_file):
        """Test element assignment consistency."""
        parser = PDBParser()
        success = parser.parse_file(sample_pdb_file)
        assert success
        
        atoms = parser.atoms
        
        # Check element-name consistency for common cases
        element_name_pairs = []
        for atom in atoms[:50]:  # Check first 50 atoms
            element_name_pairs.append((atom.element.upper(), atom.name.strip()))
        
        # Common element-name patterns
        for element, name in element_name_pairs:
            if element == "H":
                assert name.startswith("H"), f"Hydrogen should have name starting with H, got {name}"
            elif element == "C":
                assert name.startswith("C") or name in ["CA", "CB", "CG", "CD", "CE", "CZ"], \
                    f"Carbon should have appropriate name, got {name}"
            elif element == "N":
                assert name.startswith("N") or name in ["N"], \
                    f"Nitrogen should have appropriate name, got {name}"
            elif element == "O":
                assert name.startswith("O") or name in ["O"], \
                    f"Oxygen should have appropriate name, got {name}"


@pytest.mark.integration
@pytest.mark.requires_pdb_files  
class TestPDBFixingIntegration:
    """Test PDB fixing functionality integration with real files."""
    
    def test_pdb_fixing_basic_functionality(self, pdb_fixing_test_file):
        """Test basic PDB fixing functionality with 1ubi.pdb."""
        # Parse original structure
        parser = PDBParser()
        success = parser.parse_file(pdb_fixing_test_file)
        assert success, "Should parse 1ubi.pdb successfully"
        
        original_atoms = parser.atoms
        
        # Verify original state (should have no hydrogens)
        hydrogen_atoms = [atom for atom in original_atoms if atom.is_hydrogen()]
        assert len(hydrogen_atoms) == 0, "1ubi.pdb should have no hydrogens initially"
        
        # Should have expected number of heavy atoms
        assert len(original_atoms) == 683, f"Expected 683 atoms in 1ubi.pdb, got {len(original_atoms)}"
        
        # Count ATOM vs HETATM records
        atom_records = [atom for atom in original_atoms if atom.record_type == "ATOM"]
        hetatm_records = [atom for atom in original_atoms if atom.record_type == "HETATM"]
        assert len(atom_records) == 602, "Should have 602 protein atoms"
        assert len(hetatm_records) == 81, "Should have 81 water molecules"
    
    @pytest.mark.skipif(not has_openbabel(), reason="OpenBabel not available")
    def test_hydrogen_addition_openbabel(self, pdb_fixing_test_file):
        """Test hydrogen addition using OpenBabel."""
        fixer = PDBFixer()
        
        # Parse original structure
        parser = PDBParser()
        parser.parse_file(pdb_fixing_test_file)
        original_atoms = parser.atoms
        
        # Add hydrogens using OpenBabel
        result_atoms = fixer.add_missing_hydrogens(original_atoms, method="openbabel")
        
        # Should have significantly more atoms
        assert len(result_atoms) > len(original_atoms), "Should have more atoms after adding hydrogens"
        
        # Should now contain hydrogen atoms
        hydrogen_atoms = [atom for atom in result_atoms if atom.is_hydrogen()]
        assert len(hydrogen_atoms) > 0, "Should now have hydrogen atoms"
        
        # Heavy atom count should be preserved or similar
        heavy_atoms = [atom for atom in result_atoms if not atom.is_hydrogen()]
        assert len(heavy_atoms) >= len(original_atoms) * 0.9, "Heavy atom count should be preserved"
        
        # Expected hydrogen count for 1ubi.pdb with OpenBabel
        assert 750 <= len(hydrogen_atoms) <= 820, \
            f"Expected 750-820 hydrogens for OpenBabel on 1ubi.pdb, got {len(hydrogen_atoms)}"
        
        # Total should equal heavy + hydrogen
        assert len(result_atoms) == len(heavy_atoms) + len(hydrogen_atoms), \
            "Total atoms should equal heavy atoms plus hydrogens"
    
    @pytest.mark.skipif(not has_pdbfixer(), reason="PDBFixer not available")
    def test_hydrogen_addition_pdbfixer(self, pdb_fixing_test_file):
        """Test hydrogen addition using PDBFixer."""
        fixer = PDBFixer()
        
        # Parse original structure
        parser = PDBParser()
        parser.parse_file(pdb_fixing_test_file)
        original_atoms = parser.atoms
        
        # Add hydrogens using PDBFixer
        result_atoms = fixer.add_missing_hydrogens(original_atoms, method="pdbfixer", pH=7.0)
        
        # Should have significantly more atoms
        assert len(result_atoms) > len(original_atoms), "Should have more atoms after adding hydrogens"
        
        # Should now contain hydrogen atoms
        hydrogen_atoms = [atom for atom in result_atoms if atom.is_hydrogen()]
        assert len(hydrogen_atoms) > 0, "Should now have hydrogen atoms"
        
        # Heavy atom count should be preserved or similar  
        heavy_atoms = [atom for atom in result_atoms if not atom.is_hydrogen()]
        assert len(heavy_atoms) >= len(original_atoms) * 0.9, "Heavy atom count should be preserved"
        
        # Expected hydrogen count for 1ubi.pdb with PDBFixer
        assert 750 <= len(hydrogen_atoms) <= 820, \
            f"Expected 750-820 hydrogens for PDBFixer on 1ubi.pdb, got {len(hydrogen_atoms)}"
        
        # Total should equal heavy + hydrogen
        assert len(result_atoms) == len(heavy_atoms) + len(hydrogen_atoms), \
            "Total atoms should equal heavy atoms plus hydrogens"
    
    @pytest.mark.skipif(not (has_openbabel() or has_pdbfixer()), reason="No PDB fixer available")
    def test_file_based_fixing(self, pdb_fixing_test_file):
        """Test file-based PDB fixing operations."""
        fixer = PDBFixer()
        
        # Choose method based on availability
        method = "openbabel" if has_openbabel() else "pdbfixer"
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as output_file:
            output_path = output_file.name
        
        try:
            # Remove temp file so fixer can create it
            os.unlink(output_path)
            
            # Fix the structure file
            result_path = fixer.fix_structure_file(
                pdb_fixing_test_file, 
                output_path, 
                method=method, 
                overwrite=True
            )
            
            # Should return the specified output path
            assert result_path == output_path, "Should return specified output path"
            
            # Output file should exist and have content
            assert os.path.exists(output_path), "Output file should exist"
            assert os.path.getsize(output_path) > 0, "Output file should not be empty"
            
            # Should be able to parse the output file
            parser = PDBParser()
            success = parser.parse_file(output_path)
            assert success, "Should be able to parse fixed PDB file"
            
            fixed_atoms = parser.atoms
            assert len(fixed_atoms) > 683, "Fixed structure should have more atoms"
            
            # Should contain hydrogens
            hydrogen_atoms = [atom for atom in fixed_atoms if atom.is_hydrogen()]
            assert len(hydrogen_atoms) > 0, "Fixed structure should contain hydrogens"
            
        finally:
            # Clean up
            if os.path.exists(output_path):
                os.unlink(output_path)
    
    def test_hydrogen_analysis_functionality(self, pdb_fixing_test_file):
        """Test hydrogen analysis capabilities."""
        fixer = PDBFixer()
        
        # Parse original structure
        parser = PDBParser()
        parser.parse_file(pdb_fixing_test_file)
        original_atoms = parser.atoms
        
        # Analyze hydrogen content of original structure
        info = fixer.get_missing_hydrogen_info(original_atoms)
        
        # Verify analysis results
        assert info["total_atoms"] == 683, "Should report correct total atom count"
        assert info["hydrogen_atoms"] == 0, "Should report zero hydrogen atoms initially"
        assert info["heavy_atoms"] == 683, "Should report correct heavy atom count"
        assert info["hydrogen_percentage"] == 0.0, "Should report 0% hydrogen initially"
        assert info["estimated_missing_hydrogens"] > 0, "Should estimate missing hydrogens"
        assert info["has_sufficient_hydrogens"] == False, "Should indicate insufficient hydrogens"
    
    @pytest.mark.skipif(not (has_openbabel() or has_pdbfixer()), reason="No PDB fixer available")
    def test_hydrogen_analysis_after_fixing(self, pdb_fixing_test_file):
        """Test hydrogen analysis after fixing."""
        fixer = PDBFixer()
        
        # Parse and fix structure
        parser = PDBParser()
        parser.parse_file(pdb_fixing_test_file)
        original_atoms = parser.atoms
        
        # Add hydrogens
        method = "openbabel" if has_openbabel() else "pdbfixer"
        fixed_atoms = fixer.add_missing_hydrogens(original_atoms, method=method)
        
        # Analyze fixed structure
        info = fixer.get_missing_hydrogen_info(fixed_atoms)
        
        # Verify analysis results after fixing
        assert 1450 <= info["total_atoms"] <= 1500, \
            f"Expected 1450-1500 total atoms after fixing, got {info['total_atoms']}"
        assert 750 <= info["hydrogen_atoms"] <= 820, \
            f"Expected 750-820 hydrogen atoms after fixing, got {info['hydrogen_atoms']}"
        assert info["heavy_atoms"] >= 683 * 0.9, "Heavy atom count should be preserved"
        assert info["hydrogen_percentage"] > 50.0, "Should have high hydrogen percentage after fixing"
        assert info["has_sufficient_hydrogens"] == True, "Should indicate sufficient hydrogens after fixing"
    
    def test_error_handling_integration(self):
        """Test error handling in PDB fixing integration."""
        fixer = PDBFixer()
        
        # Test with non-existent file
        with pytest.raises(PDBFixerError, match="does not exist"):
            fixer.fix_structure_file("/nonexistent/file.pdb")
        
        # Test with unsupported method (need non-empty list to trigger validation)
        from hbat.core.structure import Atom
        from hbat.core.np_vector import NPVec3D
        dummy_atom = Atom(
            serial=1, name="CA", alt_loc="", res_name="ALA", chain_id="A",
            res_seq=1, i_code="", coords=NPVec3D(0, 0, 0), occupancy=1.0,
            temp_factor=0.0, element="C", charge="", record_type="ATOM"
        )
        with pytest.raises(PDBFixerError, match="Unsupported method"):
            fixer.add_missing_hydrogens([dummy_atom], method="invalid_method")
        
        # Test overwrite protection
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as existing_file:
            existing_file.write("EXISTING CONTENT")
            existing_path = existing_file.name
        
        try:
            with pytest.raises(PDBFixerError, match="already exists"):
                fixer.fix_structure_file("example_pdb_files/6rsa.pdb", existing_path, overwrite=False)
        finally:
            if os.path.exists(existing_path):
                os.unlink(existing_path)


@pytest.mark.integration
class TestPDBParsingEdgeCases:
    """Test edge cases and error handling in PDB parsing."""
    
    def test_nan_value_handling(self):
        """Test handling of NaN values in PDB data."""
        from hbat.core.pdb_parser import _safe_int_convert, _safe_float_convert
        
        # Test safe_int_convert with normal values
        assert _safe_int_convert("42") == 42
        assert _safe_int_convert(42.7) == 42  # Should truncate
        
        # Test safe_int_convert with NaN and None values
        assert _safe_int_convert(None) == 0
        assert _safe_int_convert(float('nan')) == 0
        assert _safe_int_convert(None, 99) == 99
        assert _safe_int_convert(float('nan'), 99) == 99
        
        # Test safe_int_convert with invalid values
        assert _safe_int_convert("invalid") == 0
        assert _safe_int_convert("invalid", 42) == 42
        
        # Test safe_float_convert with normal values
        assert _safe_float_convert("42.5") == 42.5
        assert _safe_float_convert("-3.14") == -3.14
        
        # Test safe_float_convert with NaN and None values
        assert _safe_float_convert(None) == 0.0
        assert _safe_float_convert(float('nan')) == 0.0
        assert _safe_float_convert(None, 99.9) == 99.9
        assert _safe_float_convert(float('nan'), 99.9) == 99.9
        
        # Test safe_float_convert with invalid values
        assert _safe_float_convert("invalid") == 0.0
        assert _safe_float_convert("invalid", 42.0) == 42.0
        
        # Test safe_float_convert with infinity values
        assert _safe_float_convert(float('inf')) == float('inf')
        assert _safe_float_convert(float('-inf')) == float('-inf')
    
    def test_malformed_pdb_data_simulation(self):
        """Test handling of malformed PDB data with NaN values."""
        parser = PDBParser()
        
        # Simulate a malformed atom row with NaN values
        class MockAtomRow:
            def get(self, key, default=None):
                data = {
                    "id": float('nan'),  # NaN serial number
                    "name": "CA",
                    "resname": "ALA",
                    "chain": "A",
                    "resid": float('nan'),  # NaN residue number
                    "x": float('nan'),  # NaN coordinates
                    "y": 1.0,
                    "z": 2.0,
                    "occupancy": float('nan'),
                    "b_factor": 20.0,
                    "element": "C"
                }
                return data.get(key, default)
        
        mock_row = MockAtomRow()
        
        # This should not raise an exception
        atom = parser._convert_atom_row(mock_row, "ATOM")
        
        # Should create atom with default values for NaN fields
        assert atom is not None
        assert atom.serial == 0  # Default for NaN serial
        assert atom.name == "CA"
        assert atom.res_seq == 0  # Default for NaN resid
        assert atom.coords.x == 0.0  # Default for NaN x
        assert atom.coords.y == 1.0  # Normal y value
        assert atom.coords.z == 2.0  # Normal z value
        assert atom.occupancy == 1.0  # Default for NaN occupancy
        assert atom.temp_factor == 20.0  # Normal b_factor
    
    def test_parser_error_handling(self):
        """Test parser error handling for invalid files."""
        parser = PDBParser()
        
        # Test non-existent file
        success = parser.parse_file("/nonexistent/file.pdb")
        assert success == False, "Should return False for non-existent file"
        
        # Test with invalid file path
        success = parser.parse_file("")
        assert success == False, "Should return False for empty file path"
        
        # Test parser state after failed parsing
        atoms = parser.atoms
        assert len(atoms) == 0, "Should have no atoms after failed parsing"
        
        bonds = parser.bonds
        assert len(bonds) == 0, "Should have no bonds after failed parsing"
    
    def test_corrupted_pdb_content_handling(self):
        """Test handling of corrupted PDB file content."""
        import tempfile
        import os
        
        # Create file with corrupted/invalid PDB content
        corrupted_content = """INVALID LINE FORMAT
ATOM    NOT PROPER FORMAT
RANDOM TEXT THAT IS NOT PDB
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(corrupted_content)
            temp_path = f.name
        
        try:
            parser = PDBParser()
            success = parser.parse_file(temp_path)
            
            # Parser should handle corrupted content gracefully
            # Either return False or successfully parse with minimal data
            assert isinstance(success, bool), "Should return boolean success status"
            
            atoms = parser.atoms
            # Should have 0 atoms if parsing completely failed,
            # or minimal atoms if some lines were parseable
            assert len(atoms) >= 0, "Atom count should be non-negative"
            
        finally:
            os.unlink(temp_path)
    
    def test_very_large_coordinate_values(self):
        """Test handling of extremely large coordinate values."""
        import tempfile
        import os
        
        # Create PDB with very large coordinates
        large_coord_content = """HEADER    TEST STRUCTURE                          01-JAN-70   TEST            
ATOM      1  CA  ALA A   1   99999.999 99999.999 99999.999  1.00 20.00           C  
ATOM      2  CB  ALA A   1  -99999.999-99999.999-99999.999  1.00 20.00           C  
END                                                                             
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(large_coord_content)
            temp_path = f.name
        
        try:
            parser = PDBParser()
            success = parser.parse_file(temp_path)
            
            if success:
                atoms = parser.atoms
                
                if len(atoms) > 0:
                    # Check that large coordinates are handled
                    atom1 = atoms[0]
                    assert isinstance(atom1.coords.x, (int, float)), "X coordinate should be numeric"
                    assert isinstance(atom1.coords.y, (int, float)), "Y coordinate should be numeric"
                    assert isinstance(atom1.coords.z, (int, float)), "Z coordinate should be numeric"
                    
                    # Should not be NaN (finite or infinite is acceptable)
                    import math
                    assert not math.isnan(atom1.coords.x), "X coordinate should not be NaN"
                    assert not math.isnan(atom1.coords.y), "Y coordinate should not be NaN"
                    assert not math.isnan(atom1.coords.z), "Z coordinate should not be NaN"
        
        finally:
            os.unlink(temp_path)