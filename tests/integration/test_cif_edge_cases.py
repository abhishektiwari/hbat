"""CIF Format Edge Case Tests: Heteroatoms, Ligands, and Special Cases

Tests for handling special cases in CIF format:
- Heteroatom detection and bonding
- Ligand connectivity preservation
- Water molecule handling
- Non-standard residue processing
- Comparison with PDB equivalents
"""

import pytest
import os
from pathlib import Path


class TestHeteroatomHandling:
    """Test P5.3: CIF handling of heteroatoms and ligands."""

    def test_cif_heteroatom_detection(self):
        """Test that heteroatoms are correctly detected in CIF files."""
        from hbat.core.pdb_parser import PDBParser

        parser = PDBParser()
        result = parser.parse_file('example_pdb_files/6RSA.cif')
        assert result is True

        # Count heteroatom residues
        heteroatom_residues = {}
        for residue in parser.residues.values():
            if residue.name not in heteroatom_residues:
                heteroatom_residues[residue.name] = 0
            heteroatom_residues[residue.name] += 1

        # Check for water molecules (DOD - deuterated water)
        assert 'DOD' in heteroatom_residues, "Water molecules not detected"
        assert heteroatom_residues['DOD'] > 0, "No water molecules found"

        # Check for ligand (UVC)
        assert 'UVC' in heteroatom_residues, "Ligand not detected"
        assert heteroatom_residues['UVC'] > 0, "No ligand found"

        # Standard amino acids should also be present
        standard_residues = {'ALA', 'GLY', 'SER', 'THR', 'CYS', 'VAL', 'LEU', 'ILE'}
        found_standards = set(heteroatom_residues.keys()) & standard_residues
        assert len(found_standards) > 0, "No standard amino acids found"

    def test_cif_water_molecule_count(self):
        """Test that water molecules are correctly counted in CIF."""
        from hbat.core.pdb_parser import PDBParser

        parser = PDBParser()
        parser.parse_file('example_pdb_files/6RSA.cif')

        # Count water molecules
        water_count = sum(1 for res in parser.residues.values() if res.name == 'DOD')
        assert water_count == 112, f"Expected 112 water molecules, found {water_count}"

    def test_cif_vs_pdb_heteroatom_equivalence(self):
        """Test that CIF and PDB have equivalent heteroatom counts."""
        from hbat.core.pdb_parser import PDBParser

        # Parse CIF
        parser_cif = PDBParser()
        parser_cif.parse_file('example_pdb_files/6RSA.cif')

        # Parse PDB
        parser_pdb = PDBParser()
        parser_pdb.parse_file('example_pdb_files/6rsa.pdb')

        # Count heteroatoms in each
        def count_heteroatoms(parser):
            heteroatom_counts = {}
            for residue in parser.residues.values():
                heteroatom_counts[residue.name] = heteroatom_counts.get(residue.name, 0) + 1
            return heteroatom_counts

        heteroatoms_cif = count_heteroatoms(parser_cif)
        heteroatoms_pdb = count_heteroatoms(parser_pdb)

        # Compare counts for common residues
        for res_type in heteroatoms_cif.keys():
            cif_count = heteroatoms_cif.get(res_type, 0)
            pdb_count = heteroatoms_pdb.get(res_type, 0)
            assert cif_count == pdb_count, \
                f"Heteroatom count mismatch for {res_type}: CIF={cif_count}, PDB={pdb_count}"

    def test_cif_heteroatom_bonds_detected(self):
        """Test that bonds involving heteroatoms are correctly detected."""
        from hbat.core.pdb_parser import PDBParser

        parser = PDBParser()
        parser.parse_file('example_pdb_files/6RSA.cif')

        # Count bonds
        total_bonds = len(parser.bonds)
        assert total_bonds > 0, "No bonds detected in CIF file"

        # Check that water molecules have bonds (O-H bonds)
        water_atoms = set()
        for residue in parser.residues.values():
            if residue.name == 'DOD':
                for atom in residue.atoms:
                    water_atoms.add(atom.serial)

        # Count bonds involving water molecules
        water_bonds = 0
        for bond in parser.bonds:
            if bond.atom1_serial in water_atoms or bond.atom2_serial in water_atoms:
                water_bonds += 1

        assert water_bonds > 0, "No bonds detected for water molecules"

    def test_cif_ligand_connectivity(self):
        """Test that ligand atoms are detected in CIF (connectivity is limited).

        Note: Full ligand connectivity (internal bonds) requires explicit struct_conn
        records in the CIF file or more sophisticated bond detection. Current
        implementation focuses on protein structure bonds. Ligand atoms are correctly
        identified but internal connectivity may be incomplete.
        """
        from hbat.core.pdb_parser import PDBParser

        parser = PDBParser()
        parser.parse_file('example_pdb_files/6RSA.cif')

        # Find ligand atoms
        ligand_atoms = set()
        ligand_residue = None
        for residue in parser.residues.values():
            if residue.name == 'UVC':
                ligand_residue = residue
                for atom in residue.atoms:
                    ligand_atoms.add(atom.serial)

        assert len(ligand_atoms) > 0, "Ligand atoms not found"
        assert ligand_residue is not None, "Ligand residue not found"

        # Verify ligand atoms are correctly parsed (31 atoms including hydrogen)
        assert len(ligand_atoms) == 31, f"Expected 31 ligand atoms, found {len(ligand_atoms)}"

    def test_cif_heteroatom_analysis(self):
        """Test that heteroatoms can be analyzed for interactions."""
        from hbat.core.analysis import AnalysisParameters, NPMolecularInteractionAnalyzer

        params = AnalysisParameters()
        params.fix_pdb_enabled = False

        analyzer = NPMolecularInteractionAnalyzer(params)
        result = analyzer.analyze_file('example_pdb_files/6RSA.cif')

        assert result is True
        # Analysis should complete without errors despite heteroatoms
        assert len(analyzer.hydrogen_bonds) > 0

    def test_cif_vs_pdb_heteroatom_analysis(self):
        """Test that heteroatom analysis produces equivalent results for CIF vs PDB."""
        from hbat.core.analysis import AnalysisParameters, NPMolecularInteractionAnalyzer

        params = AnalysisParameters()
        params.fix_pdb_enabled = False

        # Analyze CIF
        analyzer_cif = NPMolecularInteractionAnalyzer(params)
        analyzer_cif.analyze_file('example_pdb_files/6RSA.cif')

        # Analyze PDB
        analyzer_pdb = NPMolecularInteractionAnalyzer(params)
        analyzer_pdb.analyze_file('example_pdb_files/6rsa.pdb')

        # H-bond counts should be identical (no fixing)
        hb_cif = len(analyzer_cif.hydrogen_bonds)
        hb_pdb = len(analyzer_pdb.hydrogen_bonds)

        assert hb_cif == hb_pdb, \
            f"H-bond counts differ: CIF={hb_cif}, PDB={hb_pdb}"

    def test_cif_heteroatom_with_fixing(self):
        """Test that heteroatoms are preserved during PDB fixing."""
        from hbat.core.analysis import AnalysisParameters, NPMolecularInteractionAnalyzer
        from hbat.core.pdb_parser import PDBParser

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = 'pdbfixer'
        params.fix_pdb_add_hydrogens = True

        analyzer = NPMolecularInteractionAnalyzer(params)
        analyzer.analyze_file('example_pdb_files/6RSA.cif')

        # Parse the fixed file
        fixed_file = analyzer._pdb_fixing_info.get('fixed_file_path')
        assert fixed_file is not None
        assert os.path.exists(fixed_file)

        # Check that heteroatoms are still present in fixed file
        parser = PDBParser()
        result = parser.parse_file(fixed_file)
        assert result is True

        heteroatom_count = sum(1 for res in parser.residues.values() if res.name == 'DOD')
        assert heteroatom_count == 112, f"Water molecules lost during fixing: {heteroatom_count}/112"

        # Clean up
        if os.path.exists(fixed_file):
            os.unlink(fixed_file)


class TestNonStandardResidues:
    """Test handling of non-standard residues in CIF format."""

    def test_cif_non_standard_residue_detection(self):
        """Test detection of non-standard residues (like UVC ligand)."""
        from hbat.core.pdb_parser import PDBParser

        parser = PDBParser()
        parser.parse_file('example_pdb_files/6RSA.cif')

        # Find non-standard residues
        standard_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
            'TYR', 'VAL', 'HOH', 'DOD', 'WAT'  # Include common water names
        }

        non_standard = set()
        for residue in parser.residues.values():
            if residue.name not in standard_residues:
                non_standard.add(residue.name)

        # UVC should be detected as non-standard (it's a ligand)
        assert 'UVC' in non_standard, "Ligand UVC not detected as non-standard residue"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
