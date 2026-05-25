"""Format Equivalence Tests: PDB vs CIF

Comprehensive tests verifying that PDB and CIF formats produce equivalent results
for structure parsing and analysis. These tests ensure that scientists can use
either format without data loss or analysis differences.

Coverage:
- Parsing equivalence: atom counts, coordinates, residue structures, bonds
- Analysis equivalence: hydrogen bonds, π interactions, all interaction types
- Edge cases: heteroatoms, disulfide bonds, non-standard residues
- End-to-end workflows: CLI processing with both formats
"""

import pytest
import os
from pathlib import Path


class TestParsingEquivalence:
    """Verify that CIF and PDB parsing produces equivalent results."""

    def test_cif_vs_pdb_atom_count_equivalence(self):
        """Test that CIF and PDB files have same atom count."""
        from hbat.core.pdb_parser import PDBParser

        # Parse PDB file
        parser_pdb = PDBParser()
        result_pdb = parser_pdb.parse_file("example_pdb_files/6rsa.pdb")
        assert result_pdb is True, "Failed to parse PDB file"

        # Parse CIF file
        parser_cif = PDBParser()
        result_cif = parser_cif.parse_file("example_pdb_files/6RSA.cif")
        assert result_cif is True, "Failed to parse CIF file"

        # Compare atom counts
        atom_count_pdb = len(parser_pdb.atoms)
        atom_count_cif = len(parser_cif.atoms)
        assert atom_count_pdb == atom_count_cif, (
            f"Atom counts differ: PDB={atom_count_pdb}, CIF={atom_count_cif}"
        )

    def test_cif_vs_pdb_coordinate_equivalence(self):
        """Test that coordinates are identical between CIF and PDB."""
        from hbat.core.pdb_parser import PDBParser

        parser_pdb = PDBParser()
        parser_pdb.parse_file("example_pdb_files/6rsa.pdb")

        parser_cif = PDBParser()
        parser_cif.parse_file("example_pdb_files/6RSA.cif")

        # Compare coordinates for first 10 atoms
        for i in range(min(10, len(parser_pdb.atoms))):
            atom_pdb = parser_pdb.atoms[i]
            atom_cif = parser_cif.atoms[i]

            # Check coordinates are very close (within 0.01 Ångströms)
            dx = abs(atom_pdb.coords.x - atom_cif.coords.x)
            dy = abs(atom_pdb.coords.y - atom_cif.coords.y)
            dz = abs(atom_pdb.coords.z - atom_cif.coords.z)

            assert dx < 0.01, f"X-coordinate differs: {dx} Å"
            assert dy < 0.01, f"Y-coordinate differs: {dy} Å"
            assert dz < 0.01, f"Z-coordinate differs: {dz} Å"

    def test_cif_vs_pdb_residue_equivalence(self):
        """Test that residue structures are equivalent where files match.

        Note: PDB and CIF files may differ in how they organize water molecules
        and heteroatoms into chains. This test verifies that where they do match,
        the residue structures are identical.
        """
        from hbat.core.pdb_parser import PDBParser

        parser_pdb = PDBParser()
        parser_pdb.parse_file("example_pdb_files/6rsa.pdb")

        parser_cif = PDBParser()
        parser_cif.parse_file("example_pdb_files/6RSA.cif")

        # Get common chains between files
        chains_pdb = set(atom.chain_id for atom in parser_pdb.atoms)
        chains_cif = set(atom.chain_id for atom in parser_cif.atoms)
        common_chains = chains_pdb & chains_cif

        # Check that at least one common chain exists
        assert len(common_chains) > 0, (
            f"No common chains between PDB and CIF: PDB={chains_pdb}, CIF={chains_cif}"
        )

        # For each common chain, verify protein residues (not water) are equivalent
        for chain in common_chains:
            # Get protein residues (those with typical protein residue names)
            protein_res_names = {
                "ALA",
                "ARG",
                "ASN",
                "ASP",
                "CYS",
                "GLU",
                "GLN",
                "GLY",
                "HIS",
                "ILE",
                "LEU",
                "LYS",
                "MET",
                "PHE",
                "PRO",
                "SER",
                "THR",
                "TRP",
                "TYR",
                "VAL",
            }

            pdb_protein_res = {
                k: v
                for k, v in parser_pdb.residues.items()
                if v.chain_id == chain and v.name in protein_res_names
            }
            cif_protein_res = {
                k: v
                for k, v in parser_cif.residues.items()
                if v.chain_id == chain and v.name in protein_res_names
            }

            # Protein residues should match
            res_count_pdb = len(pdb_protein_res)
            res_count_cif = len(cif_protein_res)
            assert res_count_pdb == res_count_cif, (
                f"Protein residue counts differ on chain {chain}: "
                f"PDB={res_count_pdb}, CIF={res_count_cif}"
            )

            # Check that residue structures match
            for res_key in list(pdb_protein_res.keys())[:5]:  # Check first 5
                if res_key in cif_protein_res:
                    pdb_res = pdb_protein_res[res_key]
                    cif_res = cif_protein_res[res_key]

                    assert pdb_res.name == cif_res.name, (
                        f"Residue name mismatch for {res_key}"
                    )
                    assert len(pdb_res.atoms) == len(cif_res.atoms), (
                        f"Atom count in residue {res_key} differs"
                    )

    def test_cif_vs_pdb_chain_equivalence(self):
        """Test that chain organization is equivalent where files overlap."""
        from hbat.core.pdb_parser import PDBParser

        parser_pdb = PDBParser()
        parser_pdb.parse_file("example_pdb_files/6rsa.pdb")

        parser_cif = PDBParser()
        parser_cif.parse_file("example_pdb_files/6RSA.cif")

        # Extract chain IDs
        chains_pdb = set(atom.chain_id for atom in parser_pdb.atoms)
        chains_cif = set(atom.chain_id for atom in parser_cif.atoms)

        # Check that PDB chains are subset of CIF chains
        # (CIF file might be more complete)
        assert chains_pdb.issubset(chains_cif), (
            f"PDB chains not subset of CIF chains: PDB={chains_pdb}, CIF={chains_cif}"
        )

        # Verify at least one common chain
        assert len(chains_pdb & chains_cif) > 0, f"No common chains between PDB and CIF"

    def test_cif_vs_pdb_bond_count_equivalence(self):
        """Test that detected bonds are equivalent."""
        from hbat.core.pdb_parser import PDBParser

        parser_pdb = PDBParser()
        parser_pdb.parse_file("example_pdb_files/6rsa.pdb")

        parser_cif = PDBParser()
        parser_cif.parse_file("example_pdb_files/6RSA.cif")

        # Compare bond counts
        bond_count_pdb = len(parser_pdb.bonds)
        bond_count_cif = len(parser_cif.bonds)

        # Bonds should be very similar (exact match may vary due to struct_conn)
        # Allow 10% difference as acceptable variance
        max_diff = max(bond_count_pdb, bond_count_cif) * 0.1
        actual_diff = abs(bond_count_pdb - bond_count_cif)

        assert actual_diff <= max_diff, (
            f"Bond counts differ too much: PDB={bond_count_pdb}, CIF={bond_count_cif}, "
            f"difference={actual_diff}, max_allowed={max_diff}"
        )


class TestAnalysisEquivalence:
    """Test P5.2: Verify that analysis results are equivalent between formats per fixing method."""

    def test_cif_vs_pdb_h_bonds_no_fixing(self):
        """Test H-bond equivalence without PDB fixing.

        CIF and PDB should produce identical results without fixing since both
        use the same underlying atom/bond detection on the original structures.
        """
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        # No fixing - original structures only
        params = AnalysisParameters()
        params.fix_pdb_enabled = False

        analyzer_pdb = NPMolecularInteractionAnalyzer(params)
        analyzer_pdb.analyze_file("example_pdb_files/6rsa.pdb")

        analyzer_cif = NPMolecularInteractionAnalyzer(params)
        analyzer_cif.analyze_file("example_pdb_files/6RSA.cif")

        hb_pdb = len(analyzer_pdb.hydrogen_bonds)
        hb_cif = len(analyzer_cif.hydrogen_bonds)

        # Without fixing, results should be identical (0% tolerance)
        assert hb_pdb == hb_cif, (
            f"H-bond counts differ (no fixing): PDB={hb_pdb}, CIF={hb_cif}"
        )

    def test_cif_vs_pdb_h_bonds_pdbfixer(self):
        """Test H-bond equivalence with PDBFixer method.

        PDBFixer produces near-identical results for CIF and PDB formats since it
        applies the same algorithm to both formats. Variance allowed for
        stochastic hydrogen placement and format-specific parsing differences (±10%).
        """
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "pdbfixer"
        params.fix_pdb_add_hydrogens = True

        analyzer_pdb = NPMolecularInteractionAnalyzer(params)
        analyzer_pdb.analyze_file("example_pdb_files/6rsa.pdb")

        analyzer_cif = NPMolecularInteractionAnalyzer(params)
        analyzer_cif.analyze_file("example_pdb_files/6RSA.cif")

        hb_pdb = len(analyzer_pdb.hydrogen_bonds)
        hb_cif = len(analyzer_cif.hydrogen_bonds)

        # PDBFixer: allow ±10% for stochastic hydrogen placement and format-specific parsing
        variance = 0.10
        max_diff = max(hb_pdb, hb_cif) * variance
        actual_diff = abs(hb_pdb - hb_cif)

        assert actual_diff <= max_diff, (
            f"H-bond counts differ (PDBFixer): PDB={hb_pdb}, CIF={hb_cif}, "
            f"difference={actual_diff}, max_allowed={max_diff}"
        )

        # Clean up
        for analyzer in [analyzer_pdb, analyzer_cif]:
            fixed_file = analyzer._pdb_fixing_info.get("fixed_file_path")
            if fixed_file and os.path.exists(fixed_file):
                os.unlink(fixed_file)

    def test_cif_vs_pdb_h_bonds_openbabel(self):
        """Test H-bond compatibility with OpenBabel method.

        OpenBabel uses a different hydrogen placement algorithm than PDBFixer,
        resulting in different (but valid) H-bond counts. Tolerance is higher (±15%)
        to account for algorithmic differences.

        Note: OpenBabel always outputs PDB format (not CIF), even when input is CIF.
        Results are scientifically valid but different from PDBFixer.
        """
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "openbabel"
        params.fix_pdb_add_hydrogens = True

        analyzer_pdb = NPMolecularInteractionAnalyzer(params)
        analyzer_pdb.analyze_file("example_pdb_files/6rsa.pdb")

        analyzer_cif = NPMolecularInteractionAnalyzer(params)
        analyzer_cif.analyze_file("example_pdb_files/6RSA.cif")

        hb_pdb = len(analyzer_pdb.hydrogen_bonds)
        hb_cif = len(analyzer_cif.hydrogen_bonds)

        # OpenBabel: allow ±15% due to different hydrogen placement algorithm
        # Results are valid but different from PDBFixer
        variance = 0.15
        max_diff = max(hb_pdb, hb_cif) * variance
        actual_diff = abs(hb_pdb - hb_cif)

        assert actual_diff <= max_diff, (
            f"H-bond counts differ (OpenBabel): PDB={hb_pdb}, CIF={hb_cif}, "
            f"difference={actual_diff}, max_allowed={max_diff}"
        )

        # Clean up
        for analyzer in [analyzer_pdb, analyzer_cif]:
            fixed_file = analyzer._pdb_fixing_info.get("fixed_file_path")
            if fixed_file and os.path.exists(fixed_file):
                os.unlink(fixed_file)

    def test_cif_vs_pdb_pi_interaction_equivalence(self):
        """Test π interaction equivalence without fixing.

        π interactions depend only on aromatic ring positions, not hydrogens,
        so results should be identical between CIF and PDB without fixing.
        """
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = False

        analyzer_pdb = NPMolecularInteractionAnalyzer(params)
        analyzer_pdb.analyze_file("example_pdb_files/6rsa.pdb")

        analyzer_cif = NPMolecularInteractionAnalyzer(params)
        analyzer_cif.analyze_file("example_pdb_files/6RSA.cif")

        pi_pdb = len(analyzer_pdb.pi_interactions)
        pi_cif = len(analyzer_cif.pi_interactions)

        # π interactions should be identical without fixing
        assert pi_pdb == pi_cif, (
            f"π interaction counts differ: PDB={pi_pdb}, CIF={pi_cif}"
        )

    def test_cif_vs_pdb_carbonyl_interaction_equivalence(self):
        """Test that carbonyl interaction counts are equivalent."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()

        analyzer_pdb = NPMolecularInteractionAnalyzer(params)
        analyzer_pdb.analyze_file("example_pdb_files/6rsa.pdb")

        analyzer_cif = NPMolecularInteractionAnalyzer(params)
        analyzer_cif.analyze_file("example_pdb_files/6RSA.cif")

        # Compare carbonyl interaction counts
        carb_pdb = len(analyzer_pdb.carbonyl_interactions)
        carb_cif = len(analyzer_cif.carbonyl_interactions)

        # Carbonyl interactions should match exactly (not affected by H placement)
        assert carb_pdb == carb_cif, (
            f"Carbonyl interaction counts differ: PDB={carb_pdb}, CIF={carb_cif}"
        )

    def test_cif_vs_pdb_all_interactions_equivalence(self):
        """Test that all interaction types are equivalent."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()

        analyzer_pdb = NPMolecularInteractionAnalyzer(params)
        analyzer_pdb.analyze_file("example_pdb_files/6rsa.pdb")

        analyzer_cif = NPMolecularInteractionAnalyzer(params)
        analyzer_cif.analyze_file("example_pdb_files/6RSA.cif")

        # Get all interaction types
        interaction_types = [
            ("hydrogen_bonds", "Hydrogen Bonds"),
            ("halogen_bonds", "Halogen Bonds"),
            ("pi_interactions", "π Interactions"),
            ("carbonyl_interactions", "Carbonyl Interactions"),
            ("n_pi_interactions", "N-π Interactions"),
        ]

        variance = (
            0.15  # Allow ±15% variance for files with different chain organization
        )

        for attr_name, display_name in interaction_types:
            if hasattr(analyzer_pdb, attr_name) and hasattr(analyzer_cif, attr_name):
                count_pdb = len(getattr(analyzer_pdb, attr_name))
                count_cif = len(getattr(analyzer_cif, attr_name))

                if display_name == "Carbonyl Interactions":
                    # Carbonyl should match exactly
                    assert count_pdb == count_cif, (
                        f"{display_name} counts differ: PDB={count_pdb}, CIF={count_cif}"
                    )
                else:
                    # Others allow variance
                    max_diff = max(count_pdb, count_cif) * variance
                    actual_diff = abs(count_pdb - count_cif)
                    assert actual_diff <= max_diff, (
                        f"{display_name} counts differ: PDB={count_pdb}, CIF={count_cif}"
                    )

    def test_cif_vs_pdb_interaction_details_equivalence(self):
        """Test that both formats detect H-bonds with similar geometric properties.

        Note: Exact residue numbers may differ between files due to different
        chain/water organization, so we test that H-bonds are detected and
        have reasonable geometric properties rather than exact residue matches.
        """
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()

        analyzer_pdb = NPMolecularInteractionAnalyzer(params)
        analyzer_pdb.analyze_file("example_pdb_files/6rsa.pdb")

        analyzer_cif = NPMolecularInteractionAnalyzer(params)
        analyzer_cif.analyze_file("example_pdb_files/6RSA.cif")

        # Verify that both files detect hydrogen bonds
        assert len(analyzer_pdb.hydrogen_bonds) > 0, "No hydrogen bonds detected in PDB"
        assert len(analyzer_cif.hydrogen_bonds) > 0, "No hydrogen bonds detected in CIF"

        # Check that H-bonds in both files have reasonable geometry
        # (distance between 1.5 and 3.5 Å for H...A)
        for hb in analyzer_pdb.hydrogen_bonds[:5]:
            assert 1.5 <= hb.distance <= 3.5, (
                f"PDB H-bond distance out of range: {hb.distance}"
            )

        for hb in analyzer_cif.hydrogen_bonds[:5]:
            assert 1.5 <= hb.distance <= 3.5, (
                f"CIF H-bond distance out of range: {hb.distance}"
            )


class TestWithFixing:
    """Test equivalence when PDB fixing is applied."""

    def test_cif_vs_pdb_equivalence_with_fixing(self):
        """Test that fixed CIF and PDB produce equivalent results."""
        from hbat.core.analysis import (
            AnalysisParameters,
            NPMolecularInteractionAnalyzer,
        )

        params = AnalysisParameters()
        params.fix_pdb_enabled = True
        params.fix_pdb_method = "pdbfixer"
        params.fix_pdb_add_hydrogens = True

        # Analyze fixed PDB
        analyzer_pdb = NPMolecularInteractionAnalyzer(params)
        analyzer_pdb.analyze_file("example_pdb_files/6rsa.pdb")

        # Analyze fixed CIF
        analyzer_cif = NPMolecularInteractionAnalyzer(params)
        analyzer_cif.analyze_file("example_pdb_files/6RSA.cif")

        # After fixing, H-bond counts should be very similar (allow ±10% for format-specific differences)
        hb_pdb = len(analyzer_pdb.hydrogen_bonds)
        hb_cif = len(analyzer_cif.hydrogen_bonds)

        variance = 0.10
        max_diff = max(hb_pdb, hb_cif) * variance
        actual_diff = abs(hb_pdb - hb_cif)

        assert actual_diff <= max_diff, (
            f"Fixed H-bond counts differ: PDB={hb_pdb}, CIF={hb_cif}"
        )

        # Clean up
        for analyzer in [analyzer_pdb, analyzer_cif]:
            fixed_file = analyzer._pdb_fixing_info.get("fixed_file_path")
            if fixed_file and os.path.exists(fixed_file):
                os.unlink(fixed_file)


class TestParsingDetails:
    """Detailed parsing equivalence tests."""

    def test_water_molecule_parsing_equivalence(self):
        """Test that water molecules are parsed identically."""
        from hbat.core.pdb_parser import PDBParser

        parser_pdb = PDBParser()
        parser_pdb.parse_file("example_pdb_files/6rsa.pdb")

        parser_cif = PDBParser()
        parser_cif.parse_file("example_pdb_files/6RSA.cif")

        # Get common chains
        chains_pdb = set(atom.chain_id for atom in parser_pdb.atoms)
        chains_cif = set(atom.chain_id for atom in parser_cif.atoms)
        common_chains = chains_pdb & chains_cif

        # Count water molecules in common chains
        water_pdb = [
            res
            for key, res in parser_pdb.residues.items()
            if res.name == "HOH" and res.chain_id in common_chains
        ]
        water_cif = [
            res
            for key, res in parser_cif.residues.items()
            if res.name == "HOH" and res.chain_id in common_chains
        ]

        assert len(water_pdb) == len(water_cif), (
            f"Water count differs in common chains: PDB={len(water_pdb)}, CIF={len(water_cif)}"
        )

    def test_element_symbol_parsing_equivalence(self):
        """Test that element symbols are parsed correctly."""
        from hbat.core.pdb_parser import PDBParser

        parser_pdb = PDBParser()
        parser_pdb.parse_file("example_pdb_files/6rsa.pdb")

        parser_cif = PDBParser()
        parser_cif.parse_file("example_pdb_files/6RSA.cif")

        # Check element symbols for first 10 atoms
        for i in range(min(10, len(parser_pdb.atoms))):
            atom_pdb = parser_pdb.atoms[i]
            atom_cif = parser_cif.atoms[i]

            assert atom_pdb.element == atom_cif.element, (
                f"Element mismatch at atom {i}: PDB={atom_pdb.element}, "
                f"CIF={atom_cif.element}"
            )

    def test_occupancy_parsing_equivalence(self):
        """Test that occupancy values are parsed correctly."""
        from hbat.core.pdb_parser import PDBParser

        parser_pdb = PDBParser()
        parser_pdb.parse_file("example_pdb_files/6rsa.pdb")

        parser_cif = PDBParser()
        parser_cif.parse_file("example_pdb_files/6RSA.cif")

        # Check occupancy values (may have small floating point differences)
        for i in range(min(10, len(parser_pdb.atoms))):
            occ_pdb = parser_pdb.atoms[i].occupancy
            occ_cif = parser_cif.atoms[i].occupancy

            # Allow 0.01 difference due to rounding
            assert abs(occ_pdb - occ_cif) < 0.01, (
                f"Occupancy mismatch at atom {i}: PDB={occ_pdb}, CIF={occ_cif}"
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
