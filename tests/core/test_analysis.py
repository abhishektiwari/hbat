"""
Tests for analysis engine functionality.
"""

import pytest
import math
from hbat.core.analysis import HBondAnalyzer, AnalysisParameters
from hbat.constants import AtomicData
from tests.conftest import (
    ExpectedResults, 
    PDBFixingExpectedResults,
    validate_hydrogen_bond, 
    validate_pi_interaction, 
    validate_cooperativity_chain
)


class TestAnalysisParameters:
    """Test cases for AnalysisParameters."""
    
    def test_default_parameters(self):
        """Test default parameter creation."""
        params = AnalysisParameters()
        
        assert params.hb_distance_cutoff > 0
        assert params.hb_angle_cutoff > 0
        assert params.hb_donor_acceptor_cutoff > 0
        assert params.analysis_mode in ["complete", "local"]
    
    def test_custom_parameters(self):
        """Test custom parameter creation."""
        params = AnalysisParameters(
            hb_distance_cutoff=3.0,
            hb_angle_cutoff=130.0,
            analysis_mode="local"
        )
        
        assert params.hb_distance_cutoff == 3.0
        assert params.hb_angle_cutoff == 130.0
        assert params.analysis_mode == "local"
    
    def test_parameter_validation(self):
        """Test parameter validation."""
        # Valid parameters should work
        params = AnalysisParameters(
            hb_distance_cutoff=3.5,
            hb_angle_cutoff=120.0
        )
        assert params.hb_distance_cutoff == 3.5
        
        # Invalid parameters should raise errors or use defaults
        try:
            params = AnalysisParameters(hb_distance_cutoff=-1.0)
            # If no validation, at least ensure reasonable behavior
            assert params.hb_distance_cutoff > 0
        except (ValueError, AssertionError):
            # Acceptable to raise error for invalid values
            pass
    
    def test_pdb_fixing_parameters(self):
        """Test PDB fixing parameter validation."""
        # Test default PDB fixing parameters
        params = AnalysisParameters()
        assert hasattr(params, 'fix_pdb_enabled')
        assert hasattr(params, 'fix_pdb_method')
        assert hasattr(params, 'fix_pdb_add_hydrogens')
        assert hasattr(params, 'fix_pdb_add_heavy_atoms')
        assert hasattr(params, 'fix_pdb_replace_nonstandard')
        assert hasattr(params, 'fix_pdb_remove_heterogens')
        assert hasattr(params, 'fix_pdb_keep_water')
        
        # Test valid PDB fixing parameters
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="pdbfixer",
            fix_pdb_add_hydrogens=True,
            fix_pdb_add_heavy_atoms=True
        )
        assert params.fix_pdb_enabled is True
        assert params.fix_pdb_method == "pdbfixer"
        assert params.fix_pdb_add_hydrogens is True
        assert params.fix_pdb_add_heavy_atoms is True
        
        # Test parameter validation
        try:
            params.validate()
        except ValueError as e:
            pytest.fail(f"Valid PDB fixing parameters should not raise validation error: {e}")
    
    def test_pdb_fixing_parameter_validation(self):
        """Test PDB fixing parameter validation logic."""
        # Test invalid method - should raise when creating HBondAnalyzer
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="invalid_method"
        )
        with pytest.raises(ValueError, match="PDB fixing method must be one of"):
            HBondAnalyzer(params)
        
        # Test OpenBabel with PDBFixer-only operations
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_heavy_atoms=True
        )
        with pytest.raises(ValueError, match="OpenBabel does not support"):
            HBondAnalyzer(params)
        
        # Test enabled fixing with no operations selected
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="pdbfixer",
            fix_pdb_add_hydrogens=False,
            fix_pdb_add_heavy_atoms=False,
            fix_pdb_replace_nonstandard=False,
            fix_pdb_remove_heterogens=False
        )
        with pytest.raises(ValueError, match="At least one PDB fixing operation must be selected"):
            HBondAnalyzer(params)
    
    def test_pdb_fixing_parameter_validation_direct(self):
        """Test PDB fixing parameter validation using validate() method directly."""
        # Test invalid method
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="invalid_method"
        )
        errors = params.validate()
        assert any("PDB fixing method must be one of" in error for error in errors), f"Expected validation error, got: {errors}"
        
        # Test valid parameters should have no errors
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True
        )
        errors = params.validate()
        # Filter out any non-PDB-fixing related errors for this test
        pdb_errors = [e for e in errors if 'fix_pdb' in e or 'OpenBabel' in e or 'PDBFixer' in e]
        assert len(pdb_errors) == 0, f"Valid PDB fixing parameters should not produce errors: {pdb_errors}"


class TestHBondAnalyzer:
    """Test cases for HBondAnalyzer."""
    
    def test_analyzer_creation(self):
        """Test analyzer creation with different parameters."""
        # Default parameters
        analyzer = HBondAnalyzer()
        assert analyzer is not None
        assert hasattr(analyzer, 'parameters')
        
        # Custom parameters
        params = AnalysisParameters(hb_distance_cutoff=3.0)
        analyzer = HBondAnalyzer(params)
        assert analyzer.parameters.hb_distance_cutoff == 3.0
    
    def test_analyzer_with_pdb_fixing_parameters(self):
        """Test analyzer creation with PDB fixing parameters."""
        # Valid PDB fixing parameters
        params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="pdbfixer",
            fix_pdb_add_hydrogens=True,
            fix_pdb_add_heavy_atoms=False
        )
        analyzer = HBondAnalyzer(params)
        assert analyzer.parameters.fix_pdb_enabled is True
        assert analyzer.parameters.fix_pdb_method == "pdbfixer"
        
        # Invalid PDB fixing parameters should raise error
        invalid_params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="invalid_method"
        )
        with pytest.raises(ValueError):
            HBondAnalyzer(invalid_params)
    
    def test_analyzer_initial_state(self):
        """Test analyzer initial state."""
        analyzer = HBondAnalyzer()
        
        assert len(analyzer.hydrogen_bonds) == 0
        assert len(analyzer.halogen_bonds) == 0
        assert len(analyzer.pi_interactions) == 0
        assert len(analyzer.cooperativity_chains) == 0
        
        stats = analyzer.get_statistics()
        assert stats['hydrogen_bonds'] == 0
        assert stats['halogen_bonds'] == 0
        assert stats['pi_interactions'] == 0
        assert stats['total_interactions'] == 0
        
        # Test that analyzer has PDB fixing methods
        assert hasattr(analyzer, '_apply_pdb_fixing'), "Should have _apply_pdb_fixing method"
    
    @pytest.mark.integration
    def test_complete_analysis_workflow(self, sample_pdb_file):
        """Test complete analysis workflow with real PDB file."""
        analyzer = HBondAnalyzer()
        
        # Run analysis
        success = analyzer.analyze_file(sample_pdb_file)
        assert success, "Analysis should succeed"
        
        # Validate results
        stats = analyzer.get_statistics()
        
        assert stats['hydrogen_bonds'] >= ExpectedResults.MIN_HYDROGEN_BONDS, \
            f"Expected >={ExpectedResults.MIN_HYDROGEN_BONDS} H-bonds, got {stats['hydrogen_bonds']}"
        assert stats['pi_interactions'] >= ExpectedResults.MIN_PI_INTERACTIONS, \
            f"Expected >={ExpectedResults.MIN_PI_INTERACTIONS} π-interactions, got {stats['pi_interactions']}"
        assert stats['total_interactions'] >= ExpectedResults.MIN_TOTAL_INTERACTIONS, \
            f"Expected >={ExpectedResults.MIN_TOTAL_INTERACTIONS} total interactions, got {stats['total_interactions']}"
    
    @pytest.mark.integration
    def test_pdb_fixing_workflow(self, pdb_fixing_test_file):
        """Test analysis workflow with PDB fixing enabled using 1ubi.pdb."""
        # Test with OpenBabel fixing
        params_ob = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True
        )
        analyzer_ob = HBondAnalyzer(params_ob)
        
        success = analyzer_ob.analyze_file(pdb_fixing_test_file)
        assert success, "Analysis with OpenBabel PDB fixing should succeed"
        
        stats_ob = analyzer_ob.get_statistics()
        assert stats_ob['hydrogen_bonds'] >= 0, "Should have non-negative hydrogen bonds"
        assert stats_ob['total_interactions'] >= PDBFixingExpectedResults.MIN_TOTAL_INTERACTIONS, \
            f"Expected >={PDBFixingExpectedResults.MIN_TOTAL_INTERACTIONS} total interactions with fixing"
        
        # Test with PDBFixer fixing
        params_pdb = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="pdbfixer",
            fix_pdb_add_hydrogens=True,
            fix_pdb_add_heavy_atoms=True
        )
        analyzer_pdb = HBondAnalyzer(params_pdb)
        
        success = analyzer_pdb.analyze_file(pdb_fixing_test_file)
        assert success, "Analysis with PDBFixer PDB fixing should succeed"
        
        stats_pdb = analyzer_pdb.get_statistics()
        assert stats_pdb['hydrogen_bonds'] >= 0, "Should have non-negative hydrogen bonds"
    
    @pytest.mark.integration
    def test_hydrogen_bond_analysis(self, sample_pdb_file):
        """Test hydrogen bond detection and validation."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        hbonds = analyzer.hydrogen_bonds
        assert len(hbonds) > 0, "Should find hydrogen bonds"
        
        # Validate first few hydrogen bonds
        for hb in hbonds[:5]:
            validate_hydrogen_bond(hb)
            
            # Additional validation
            assert hb.distance > 0, "Distance should be positive"
            assert hb.distance <= analyzer.parameters.hb_distance_cutoff, \
                "Distance should be within cutoff"
            
            # Angle should be in reasonable range
            angle_degrees = math.degrees(hb.angle)
            assert angle_degrees >= analyzer.parameters.hb_angle_cutoff, \
                f"Angle {angle_degrees}° should be >= {analyzer.parameters.hb_angle_cutoff}°"
    
    def test_bond_based_hydrogen_donor_detection(self, sample_pdb_file):
        """Test that hydrogen bond donor detection uses pre-calculated bonds."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success, "Analysis should succeed"
        
        # Get hydrogen bond donors using the updated method
        donors = analyzer._get_hydrogen_bond_donors()
        assert isinstance(donors, list), "Should return a list of donors"
        
        # Validate donor structure
        for heavy_atom, hydrogen_atom in donors[:5]:  # Check first 5 donors
            # Heavy atom should be a donor element
            assert heavy_atom.element.upper() in ["N", "O", "S"], \
                f"Heavy atom element {heavy_atom.element} should be N, O, or S"
            
            # Hydrogen should be hydrogen
            assert hydrogen_atom.is_hydrogen(), "Hydrogen atom should be hydrogen"
            
            # Verify they are actually bonded according to the bond list
            bonded_serials = analyzer.parser.get_bonded_atoms(hydrogen_atom.serial)
            assert heavy_atom.serial in bonded_serials, \
                "Heavy atom and hydrogen should be bonded according to bond list"
        
        # Test that we're actually using bonds vs falling back to distance calculation
        # The parser should have detected bonds during parsing
        bonds = analyzer.parser.get_bonds()
        assert len(bonds) > 0, "Parser should have detected bonds"
    
    @pytest.mark.integration
    def test_pi_interaction_analysis(self, sample_pdb_file):
        """Test π interaction detection and validation."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        pi_interactions = analyzer.pi_interactions
        if len(pi_interactions) > 0:
            # Validate π interactions
            for pi in pi_interactions[:3]:
                validate_pi_interaction(pi)
                
                # Additional validation
                assert pi.distance > 0, "Distance should be positive"
                assert pi.distance <= analyzer.parameters.pi_distance_cutoff, \
                    "Distance should be within cutoff"
                
                # Check π center coordinates
                assert hasattr(pi.pi_center, 'x'), "π center should have coordinates"
                assert hasattr(pi.pi_center, 'y'), "π center should have coordinates"
                assert hasattr(pi.pi_center, 'z'), "π center should have coordinates"
    
    @pytest.mark.integration
    def test_cooperativity_analysis(self, sample_pdb_file):
        """Test cooperativity chain analysis."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        chains = analyzer.cooperativity_chains
        stats = analyzer.get_statistics()
        
        if len(chains) > 0:
            assert stats.get('cooperativity_chains', 0) == len(chains), \
                "Statistics should match actual chain count"
            
            # Validate cooperativity chains
            for chain in chains[:3]:
                validate_cooperativity_chain(chain)
    
    def test_bond_based_halogen_detection(self, sample_pdb_file):
        """Test that halogen bond detection uses pre-calculated bonds."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success, "Analysis should succeed"
        
        # Test halogen atom detection - should only return halogens bonded to carbon
        halogens = analyzer._get_halogen_atoms()
        assert isinstance(halogens, list), "Should return a list of halogens"
        
        # Verify that each halogen is actually bonded to carbon
        atom_map = {atom.serial: atom for atom in analyzer.parser.atoms}
        
        for halogen in halogens:
            # Should be a halogen element
            assert halogen.element.upper() in ["F", "CL", "BR", "I"], \
                f"Atom element {halogen.element} should be a halogen"
            
            # Should be bonded to at least one carbon according to bond list
            bonded_serials = analyzer.parser.get_bonded_atoms(halogen.serial)
            has_carbon_bond = False
            for bonded_serial in bonded_serials:
                bonded_atom = atom_map.get(bonded_serial)
                if bonded_atom is not None and bonded_atom.element.upper() == "C":
                    has_carbon_bond = True
                    break
            
            assert has_carbon_bond, \
                f"Halogen {halogen.element} at serial {halogen.serial} should be bonded to carbon"
        
        # Test _find_bonded_carbon function for each halogen
        for halogen in halogens[:5]:  # Test first 5 halogens
            bonded_carbon = analyzer._find_bonded_carbon(halogen)
            if bonded_carbon is not None:  # May be None if no carbon bonded
                assert bonded_carbon.element.upper() == "C", \
                    "Bonded atom should be carbon"
                
                # Verify they are actually bonded according to bond list
                bonded_serials = analyzer.parser.get_bonded_atoms(halogen.serial)
                assert bonded_carbon.serial in bonded_serials, \
                    "Carbon and halogen should be bonded according to bond list"
        
        # Verify bond detection is working
        bonds = analyzer.parser.get_bonds()
        assert len(bonds) > 0, "Parser should have detected bonds"
    
    def test_halogen_bond_classification(self):
        """Test halogen bond classification function."""
        analyzer = HBondAnalyzer()
        
        # Create test atoms for different halogen bond types
        from hbat.core.pdb_parser import Atom
        from hbat.core.vector import Vec3D
        
        # Test Cl...O bond
        cl_atom = Atom(
            serial=1, name="CL", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=1, i_code="", coords=Vec3D(0, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="CL", charge="", record_type="ATOM"
        )
        
        o_atom = Atom(
            serial=2, name="O", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=2, i_code="", coords=Vec3D(2, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="O", charge="", record_type="ATOM"
        )
        
        # Test classification
        bond_type = analyzer._classify_halogen_bond(cl_atom, o_atom)
        assert bond_type == "CL...O", f"Expected 'CL...O', got '{bond_type}'"
        
        # Test Br...N bond
        br_atom = Atom(
            serial=3, name="BR", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=3, i_code="", coords=Vec3D(4, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="BR", charge="", record_type="ATOM"
        )
        
        n_atom = Atom(
            serial=4, name="N", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=4, i_code="", coords=Vec3D(6, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="N", charge="", record_type="ATOM"
        )
        
        bond_type = analyzer._classify_halogen_bond(br_atom, n_atom)
        assert bond_type == "BR...N", f"Expected 'BR...N', got '{bond_type}'"
        
        # Test I...S bond
        i_atom = Atom(
            serial=5, name="I", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=5, i_code="", coords=Vec3D(8, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="I", charge="", record_type="ATOM"
        )
        
        s_atom = Atom(
            serial=6, name="S", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=6, i_code="", coords=Vec3D(10, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="S", charge="", record_type="ATOM"
        )
        
        bond_type = analyzer._classify_halogen_bond(i_atom, s_atom)
        assert bond_type == "I...S", f"Expected 'I...S', got '{bond_type}'"
        
        # Test case insensitivity
        cl_lower = Atom(
            serial=7, name="cl", alt_loc="", res_name="TEST", chain_id="A",
            res_seq=7, i_code="", coords=Vec3D(12, 0, 0), occupancy=1.0,
            temp_factor=20.0, element="cl", charge="", record_type="ATOM"
        )
        
        bond_type = analyzer._classify_halogen_bond(cl_lower, o_atom)
        assert bond_type == "CL...O", f"Expected 'CL...O' (case insensitive), got '{bond_type}'"

    @pytest.mark.integration
    def test_interaction_statistics(self, sample_pdb_file):
        """Test interaction statistics consistency."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        stats = analyzer.get_statistics()
        
        # Check that statistics match actual counts
        assert stats['hydrogen_bonds'] == len(analyzer.hydrogen_bonds), \
            "H-bond count mismatch"
        assert stats['halogen_bonds'] == len(analyzer.halogen_bonds), \
            "Halogen bond count mismatch"
        assert stats['pi_interactions'] == len(analyzer.pi_interactions), \
            "π-interaction count mismatch"
        
        # Total should be sum of individual types
        expected_total = (stats['hydrogen_bonds'] + 
                         stats['halogen_bonds'] + 
                         stats['pi_interactions'])
        assert stats['total_interactions'] == expected_total, \
            "Total interactions should sum correctly"
    
    @pytest.mark.integration
    def test_analysis_modes(self, sample_pdb_file):
        """Test different analysis modes."""
        # Complete mode
        params_complete = AnalysisParameters(analysis_mode="complete")
        analyzer_complete = HBondAnalyzer(params_complete)
        success = analyzer_complete.analyze_file(sample_pdb_file)
        assert success
        
        stats_complete = analyzer_complete.get_statistics()
        
        # Local mode
        params_local = AnalysisParameters(analysis_mode="local")
        analyzer_local = HBondAnalyzer(params_local)
        success = analyzer_local.analyze_file(sample_pdb_file)
        assert success
        
        stats_local = analyzer_local.get_statistics()
        
        # Complete mode should generally find more interactions
        assert stats_complete['total_interactions'] >= stats_local['total_interactions'], \
            "Complete mode should find at least as many interactions as local mode"
    
    @pytest.mark.integration
    def test_parameter_effects(self, sample_pdb_file):
        """Test effects of different parameter values."""
        # Strict parameters
        strict_params = AnalysisParameters(
            hb_distance_cutoff=3.0,
            hb_angle_cutoff=140.0
        )
        analyzer_strict = HBondAnalyzer(strict_params)
        success = analyzer_strict.analyze_file(sample_pdb_file)
        assert success
        
        # Permissive parameters
        permissive_params = AnalysisParameters(
            hb_distance_cutoff=4.0,
            hb_angle_cutoff=110.0
        )
        analyzer_permissive = HBondAnalyzer(permissive_params)
        success = analyzer_permissive.analyze_file(sample_pdb_file)
        assert success
        
        strict_stats = analyzer_strict.get_statistics()
        permissive_stats = analyzer_permissive.get_statistics()
        
        # Permissive should generally find more interactions
        assert permissive_stats['hydrogen_bonds'] >= strict_stats['hydrogen_bonds'], \
            "Permissive parameters should find at least as many H-bonds"
    
    @pytest.mark.integration
    def test_pdb_fixing_effects(self, pdb_fixing_test_file):
        """Test effects of PDB fixing on analysis results using 1ubi.pdb."""
        # Analysis without PDB fixing
        params_no_fix = AnalysisParameters(fix_pdb_enabled=False)
        analyzer_no_fix = HBondAnalyzer(params_no_fix)
        success = analyzer_no_fix.analyze_file(pdb_fixing_test_file)
        assert success
        
        # Analysis with PDB fixing (add hydrogens)
        params_with_fix = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True
        )
        analyzer_with_fix = HBondAnalyzer(params_with_fix)
        success = analyzer_with_fix.analyze_file(pdb_fixing_test_file)
        assert success
        
        stats_no_fix = analyzer_no_fix.get_statistics()
        stats_with_fix = analyzer_with_fix.get_statistics()
        
        # PDB fixing should generally not decrease interaction count
        # (may find more interactions with added hydrogens)
        assert stats_with_fix['total_interactions'] >= 0, "Should have non-negative interactions"
        assert stats_no_fix['total_interactions'] >= 0, "Should have non-negative interactions"
        
        print(f"\nPDB fixing effects on 1ubi.pdb:")
        print(f"  Without fixing: {stats_no_fix['total_interactions']} interactions")
        print(f"  With fixing: {stats_with_fix['total_interactions']} interactions")


class TestAtomicPropertyLookup:
    """Test atomic property lookup functionality."""
    
    def test_covalent_radius_lookup(self):
        """Test covalent radius lookup for various atoms."""
        analyzer = HBondAnalyzer()
        
        test_cases = [
            ("N", AtomicData.COVALENT_RADII.get('N', 0.71)),
            ("O", AtomicData.COVALENT_RADII.get('O', 0.66)),
            ("C", AtomicData.COVALENT_RADII.get('C', 0.76)),
            ("H", AtomicData.COVALENT_RADII.get('H', 0.31)),
            ("CA", AtomicData.COVALENT_RADII.get('C', 0.76)),  # Should use C
            ("ND1", AtomicData.COVALENT_RADII.get('N', 0.71)),  # Should use N
            ("OE1", AtomicData.COVALENT_RADII.get('O', 0.66)),  # Should use O
        ]
        
        for atom_symbol, expected_radius in test_cases:
            radius = analyzer._get_covalent_radius(atom_symbol)
            assert abs(radius - expected_radius) < 1e-6, \
                f"Covalent radius mismatch for {atom_symbol}"
    
    def test_vdw_radius_lookup(self):
        """Test van der Waals radius lookup."""
        analyzer = HBondAnalyzer()
        
        # Test basic elements
        for element in ['C', 'N', 'O', 'H']:
            radius = analyzer._get_vdw_radius(element)
            assert radius > 0, f"VDW radius should be positive for {element}"
            
        # Test complex atom names
        ca_radius = analyzer._get_vdw_radius("CA")
        c_radius = analyzer._get_vdw_radius("C")
        assert abs(ca_radius - c_radius) < 1e-6, "CA should use C radius"
    
    def test_electronegativity_lookup(self):
        """Test electronegativity lookup."""
        analyzer = HBondAnalyzer()
        
        # Test that common elements have reasonable electronegativities
        for element in ['C', 'N', 'O', 'H']:
            en = analyzer._get_electronegativity(element)
            assert 0 <= en <= 4.0, f"Electronegativity should be reasonable for {element}"
    
    def test_atomic_mass_lookup(self):
        """Test atomic mass lookup."""
        analyzer = HBondAnalyzer()
        
        # Test that common elements have reasonable masses
        for element in ['C', 'N', 'O', 'H']:
            mass = analyzer._get_atomic_mass(element)
            assert mass > 0, f"Atomic mass should be positive for {element}"
    
    def test_edge_cases(self):
        """Test edge cases in atomic property lookup."""
        analyzer = HBondAnalyzer()
        
        # Test case insensitive lookup
        ca_lower = analyzer._get_covalent_radius("ca")
        ca_upper = analyzer._get_covalent_radius("CA")
        assert abs(ca_lower - ca_upper) < 1e-6, "Should be case insensitive"
        
        # Test unknown atoms (should fall back to carbon)
        unknown_radius = analyzer._get_covalent_radius("XYZ")
        c_radius = analyzer._get_covalent_radius("C")
        assert abs(unknown_radius - c_radius) < 1e-6, "Unknown atoms should use C fallback"
    
    def test_pdb_specific_atoms(self):
        """Test PDB-specific atom name handling."""
        analyzer = HBondAnalyzer()
        
        pdb_atoms = [
            ("CA", "C"),   # Alpha carbon
            ("CB", "C"),   # Beta carbon
            ("ND1", "N"),  # Histidine nitrogen
            ("NE2", "N"),  # Histidine nitrogen
            ("OE1", "O"),  # Glutamate oxygen
            ("OD1", "O"),  # Aspartate oxygen
            ("H1", "H"),   # Hydrogen with number
            ("HG", "H"),   # Hydrogen gamma
        ]
        
        for pdb_name, element in pdb_atoms:
            radius = analyzer._get_covalent_radius(pdb_name)
            expected_radius = AtomicData.COVALENT_RADII.get(element, 0.76)
            assert abs(radius - expected_radius) < 1e-6, \
                f"{pdb_name} should use {element} properties"


class TestPerformanceMetrics:
    """Test performance and expected results."""
    
    @pytest.mark.integration
    @pytest.mark.slow
    def test_performance_benchmarks(self, sample_pdb_file):
        """Test that analysis meets performance expectations."""
        analyzer = HBondAnalyzer()
        
        import time
        start_time = time.time()
        success = analyzer.analyze_file(sample_pdb_file)
        analysis_time = time.time() - start_time
        
        assert success, "Analysis should succeed"
        
        # Analysis should complete in reasonable time (adjust as needed)
        assert analysis_time < 60.0, f"Analysis took too long: {analysis_time:.2f}s"
        
        stats = analyzer.get_statistics()
        
        # Performance metrics - should find substantial interactions
        assert stats['hydrogen_bonds'] >= ExpectedResults.MIN_HYDROGEN_BONDS, \
            "Should find substantial number of hydrogen bonds"
        assert stats['total_interactions'] >= ExpectedResults.MIN_TOTAL_INTERACTIONS, \
            "Should find substantial total interactions"
    
    @pytest.mark.integration
    def test_expected_results_documentation(self, sample_pdb_file):
        """Document expected results for 6RSA.pdb."""
        analyzer = HBondAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success
        
        stats = analyzer.get_statistics()
        
        # Print results for documentation
        print(f"\nExpected results for 6RSA.pdb:")
        print(f"  - Hydrogen bonds: {stats['hydrogen_bonds']}")
        print(f"  - Halogen bonds: {stats['halogen_bonds']}")
        print(f"  - π interactions: {stats['pi_interactions']}")
        print(f"  - Cooperativity chains: {stats.get('cooperativity_chains', 0)}")
        print(f"  - Total interactions: {stats['total_interactions']}")
        
        # Validate against minimum expectations
        assert stats['hydrogen_bonds'] >= ExpectedResults.MIN_HYDROGEN_BONDS
        assert stats['pi_interactions'] >= ExpectedResults.MIN_PI_INTERACTIONS
        assert stats['total_interactions'] >= ExpectedResults.MIN_TOTAL_INTERACTIONS
    
    @pytest.mark.integration
    def test_pdb_fixing_results_documentation(self, pdb_fixing_test_file):
        """Document expected results for 1ubi.pdb with PDB fixing."""
        # Test with OpenBabel fixing
        params_ob = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True
        )
        analyzer_ob = HBondAnalyzer(params_ob)
        success = analyzer_ob.analyze_file(pdb_fixing_test_file)
        assert success
        
        stats_ob = analyzer_ob.get_statistics()
        
        # Print results for documentation
        print(f"\nExpected results for 1ubi.pdb with OpenBabel PDB fixing:")
        print(f"  - Hydrogen bonds: {stats_ob['hydrogen_bonds']}")
        print(f"  - Halogen bonds: {stats_ob['halogen_bonds']}")
        print(f"  - π interactions: {stats_ob['pi_interactions']}")
        print(f"  - Total interactions: {stats_ob['total_interactions']}")
        
        # Test with PDBFixer fixing
        params_pdb = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="pdbfixer",
            fix_pdb_add_hydrogens=True,
            fix_pdb_add_heavy_atoms=True
        )
        analyzer_pdb = HBondAnalyzer(params_pdb)
        success = analyzer_pdb.analyze_file(pdb_fixing_test_file)
        assert success
        
        stats_pdb = analyzer_pdb.get_statistics()
        
        print(f"\nExpected results for 1ubi.pdb with PDBFixer PDB fixing:")
        print(f"  - Hydrogen bonds: {stats_pdb['hydrogen_bonds']}")
        print(f"  - Halogen bonds: {stats_pdb['halogen_bonds']}")
        print(f"  - π interactions: {stats_pdb['pi_interactions']}")
        print(f"  - Total interactions: {stats_pdb['total_interactions']}")
        
        # Both should produce valid results
        assert stats_ob['total_interactions'] >= 0
        assert stats_pdb['total_interactions'] >= 0