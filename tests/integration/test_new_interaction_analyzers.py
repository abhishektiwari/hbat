"""
Integration tests for new molecular interaction analyzers.

These tests verify that π-π stacking, carbonyl-carbonyl, and n→π* interaction
analyzers work correctly with real PDB structures and integrate properly
with the overall HBAT analysis workflow.
"""

import pytest
from hbat.core.analysis import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters
from hbat.core.interactions import PiPiInteraction, CarbonylInteraction, NPiInteraction


@pytest.mark.integration
class TestPiPiStackingAnalyzer:
    """Test π-π stacking interaction analyzer integration."""
    
    def test_pi_pi_detection_7nwd(self):
        """Test π-π stacking detection with 7NWD.pdb structure."""
        # 7NWD contains aromatic residues suitable for π-π stacking analysis
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success, "Failed to analyze 7NWD.pdb"
        
        # Check that π-π interactions were detected
        assert hasattr(analyzer, 'pi_pi_interactions'), "π-π interactions not found"
        
        pi_pi_interactions = analyzer.pi_pi_interactions
        assert len(pi_pi_interactions) > 0, "No π-π interactions detected in 7NWD"
        
        # Validate first π-π interaction
        first_interaction = pi_pi_interactions[0]
        assert isinstance(first_interaction, PiPiInteraction)
        
        # Check required properties
        assert hasattr(first_interaction, 'ring1_atoms')
        assert hasattr(first_interaction, 'ring2_atoms')
        assert hasattr(first_interaction, 'distance')
        assert hasattr(first_interaction, 'plane_angle')
        assert hasattr(first_interaction, 'stacking_type')
        
        # Validate distance range (typical π-π stacking: 3.3-6.0 Å)
        assert 3.0 <= first_interaction.distance <= 6.5
        
        # Validate angle range (0-180°, angles can be obtuse)
        assert 0.0 <= first_interaction.plane_angle <= 180.0
        
        # Validate stacking type
        valid_types = ['parallel', 'T-shaped', 'offset']
        assert first_interaction.stacking_type in valid_types
        
        print(f"✓ Detected {len(pi_pi_interactions)} π-π interactions in 7NWD")
        print(f"  First interaction: {first_interaction.ring1_residue} - {first_interaction.ring2_residue}")
        print(f"  Distance: {first_interaction.distance:.2f}Å, Type: {first_interaction.stacking_type}")
    
    def test_pi_pi_different_stacking_types(self):
        """Test detection of different π-π stacking types."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        if hasattr(analyzer, 'pi_pi_interactions') and analyzer.pi_pi_interactions:
            stacking_types = {interaction.stacking_type for interaction in analyzer.pi_pi_interactions}
            
            # Should detect at least one type
            assert len(stacking_types) > 0
            
            # All types should be valid
            valid_types = {'parallel', 'T-shaped', 'offset'}
            assert stacking_types.issubset(valid_types)
            
            print(f"✓ π-π stacking types found: {sorted(stacking_types)}")
    
    def test_pi_pi_aromatic_residue_types(self):
        """Test π-π interactions involve aromatic residues."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        if hasattr(analyzer, 'pi_pi_interactions') and analyzer.pi_pi_interactions:
            # Include nucleotides and aromatic amino acids
            aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS', 'DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'T', 'U'}
            
            for interaction in analyzer.pi_pi_interactions:
                # Extract residue types from residue strings
                ring1_type = interaction.ring1_type
                ring2_type = interaction.ring2_type
                
                # Should involve aromatic residues or nucleotides
                assert ring1_type in aromatic_residues or ring2_type in aromatic_residues
            
            print(f"✓ All π-π interactions involve aromatic residues")


@pytest.mark.integration 
class TestCarbonylCarbonylAnalyzer:
    """Test carbonyl-carbonyl interaction analyzer integration."""
    
    def test_carbonyl_detection_with_structure(self):
        """Test carbonyl-carbonyl detection with protein structure."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        # Use 6RSA which has known carbonyl interactions
        success = analyzer.analyze_file('example_pdb_files/6rsa.pdb')
        assert success, "Failed to analyze structure for carbonyl interactions"
        
        # Check that carbonyl interactions were detected
        assert hasattr(analyzer, 'carbonyl_interactions'), "Carbonyl interactions not found"
        
        if analyzer.carbonyl_interactions:
            carbonyl_interactions = analyzer.carbonyl_interactions
            
            # Validate first carbonyl interaction
            first_interaction = carbonyl_interactions[0]
            assert isinstance(first_interaction, CarbonylInteraction)
            
            # Check required properties
            assert hasattr(first_interaction, 'donor_carbon')
            assert hasattr(first_interaction, 'donor_oxygen')
            assert hasattr(first_interaction, 'acceptor_carbon')
            assert hasattr(first_interaction, 'acceptor_oxygen')
            assert hasattr(first_interaction, 'distance')
            assert hasattr(first_interaction, 'burgi_dunitz_angle')
            
            # Validate distance range (typical n→π*: 2.8-4.0 Å)
            assert 2.5 <= first_interaction.distance <= 4.5
            
            # Validate Bürgi-Dunitz angle range (95-125°)
            assert 90.0 <= first_interaction.burgi_dunitz_angle <= 130.0
            
            print(f"✓ Detected {len(carbonyl_interactions)} carbonyl interactions")
            print(f"  First interaction: {first_interaction.donor_residue} - {first_interaction.acceptor_residue}")
            print(f"  Distance: {first_interaction.distance:.2f}Å, Angle: {first_interaction.burgi_dunitz_angle:.1f}°")
        else:
            print("ℹ No carbonyl interactions detected in this structure")
    
    def test_carbonyl_angle_calculations(self):
        """Test that Bürgi-Dunitz angles are calculated correctly."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/6rsa.pdb')
        assert success
        
        if hasattr(analyzer, 'carbonyl_interactions') and analyzer.carbonyl_interactions:
            for interaction in analyzer.carbonyl_interactions:
                # Bürgi-Dunitz angle should be in the expected range
                angle = interaction.burgi_dunitz_angle
                assert 90.0 <= angle <= 130.0, f"Invalid Bürgi-Dunitz angle: {angle:.1f}°"
                
                # Angle property should return radians
                angle_rad = interaction.angle
                expected_rad = angle * 3.14159 / 180.0
                assert abs(angle_rad - expected_rad) < 0.1
            
            print(f"✓ All Bürgi-Dunitz angles in valid range (90-130°)")
    
    def test_carbonyl_backbone_sidechain_classification(self):
        """Test carbonyl interaction backbone/sidechain classification."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/6rsa.pdb')
        assert success
        
        if hasattr(analyzer, 'carbonyl_interactions') and analyzer.carbonyl_interactions:
            valid_classifications = {
                'backbone-backbone', 'backbone-sidechain', 
                'sidechain-backbone', 'sidechain-sidechain'
            }
            
            for interaction in analyzer.carbonyl_interactions:
                assert interaction.interaction_classification in valid_classifications
            
            print(f"✓ All carbonyl interactions have valid classifications")


@pytest.mark.integration
class TestNPiAnalyzer:
    """Test n→π* interaction analyzer integration."""
    
    def test_n_pi_detection_with_structure(self):
        """Test n→π* interaction detection with protein structure."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success, "Failed to analyze structure for n→π* interactions"
        
        # Check that n→π* interactions were detected
        assert hasattr(analyzer, 'n_pi_interactions'), "n→π* interactions not found"
        
        if analyzer.n_pi_interactions:
            n_pi_interactions = analyzer.n_pi_interactions
            
            # Validate first n→π* interaction
            first_interaction = n_pi_interactions[0]
            assert isinstance(first_interaction, NPiInteraction)
            
            # Check required properties
            assert hasattr(first_interaction, 'lone_pair_atom')
            assert hasattr(first_interaction, 'pi_center')
            assert hasattr(first_interaction, 'pi_atoms')
            assert hasattr(first_interaction, 'distance')
            assert hasattr(first_interaction, 'angle_to_plane')
            assert hasattr(first_interaction, 'subtype')
            
            # Validate distance range (typical n→π*: 3.0-5.0 Å)
            assert 2.5 <= first_interaction.distance <= 5.5
            
            # Validate angle to plane range (0-90°)
            assert 0.0 <= first_interaction.angle_to_plane <= 90.0
            
            # Validate donor element
            valid_elements = {'O', 'N', 'S'}
            assert first_interaction.donor_element in valid_elements
            
            print(f"✓ Detected {len(n_pi_interactions)} n→π* interactions")
            print(f"  First interaction: {first_interaction.donor_residue} - {first_interaction.acceptor_residue}")
            print(f"  Distance: {first_interaction.distance:.2f}Å, Angle: {first_interaction.angle_to_plane:.1f}°")
            print(f"  Subtype: {first_interaction.subtype}")
        else:
            print("ℹ No n→π* interactions detected in this structure")
    
    def test_n_pi_donor_element_types(self):
        """Test n→π* interactions with different donor elements."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        if hasattr(analyzer, 'n_pi_interactions') and analyzer.n_pi_interactions:
            donor_elements = {interaction.donor_element for interaction in analyzer.n_pi_interactions}
            
            # Should only have valid donor elements
            valid_elements = {'O', 'N', 'S'}
            assert donor_elements.issubset(valid_elements)
            
            print(f"✓ n→π* donor elements found: {sorted(donor_elements)}")
    
    def test_n_pi_subtype_classification(self):
        """Test n→π* interaction subtype classification."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        if hasattr(analyzer, 'n_pi_interactions') and analyzer.n_pi_interactions:
            subtypes = {interaction.subtype for interaction in analyzer.n_pi_interactions}
            
            # Should have meaningful subtypes
            assert len(subtypes) > 0
            for subtype in subtypes:
                assert len(subtype) > 0  # Non-empty subtype
                assert isinstance(subtype, str)
            
            print(f"✓ n→π* interaction subtypes: {sorted(subtypes)}")


@pytest.mark.integration
class TestNewInteractionsIntegration:
    """Test integration of all new interaction types together."""
    
    def test_all_new_interactions_detected(self):
        """Test that all new interaction types can be detected together."""
        params = AnalysisParameters()
        
        # Test with 7nwd for π-π interactions
        analyzer_pi = MolecularInteractionAnalyzer(parameters=params)
        success_pi = analyzer_pi.analyze_file('example_pdb_files/7nwd.pdb')
        assert success_pi, "Failed to analyze 7nwd.pdb"
        
        # Test with 6rsa for carbonyl interactions  
        analyzer_carbonyl = MolecularInteractionAnalyzer(parameters=params)
        success_carbonyl = analyzer_carbonyl.analyze_file('example_pdb_files/6rsa.pdb')
        assert success_carbonyl, "Failed to analyze 6rsa.pdb"
        
        # Use 7nwd as primary analyzer for summary
        analyzer = analyzer_pi
        
        # Check that analyzer has all new interaction attributes
        assert hasattr(analyzer, 'pi_pi_interactions'), "π-π interactions not available"
        assert hasattr(analyzer, 'carbonyl_interactions'), "Carbonyl interactions not available"
        assert hasattr(analyzer, 'n_pi_interactions'), "n→π* interactions not available"
        
        # Count total new interactions
        pi_pi_count = len(analyzer.pi_pi_interactions) if analyzer.pi_pi_interactions else 0
        carbonyl_count = len(analyzer.carbonyl_interactions) if analyzer.carbonyl_interactions else 0
        n_pi_count = len(analyzer.n_pi_interactions) if analyzer.n_pi_interactions else 0
        
        total_new_interactions = pi_pi_count + carbonyl_count + n_pi_count
        
        print(f"✓ New interaction detection summary for 7NWD:")
        print(f"  π-π stacking: {pi_pi_count}")
        print(f"  Carbonyl-carbonyl: {carbonyl_count}")
        print(f"  n→π* interactions: {n_pi_count}")
        print(f"  Total new interactions: {total_new_interactions}")
        
        # Should detect at least some new interactions in a protein structure
        assert total_new_interactions > 0, "No new interactions detected"
    
    def test_new_interactions_do_not_interfere_with_existing(self):
        """Test that new interaction detection doesn't break existing functionality."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        # Traditional interactions should still work
        assert hasattr(analyzer, 'hydrogen_bonds'), "Hydrogen bonds not available"
        assert hasattr(analyzer, 'halogen_bonds'), "Halogen bonds not available"
        
        # Should detect traditional interactions
        h_bonds = analyzer.hydrogen_bonds if analyzer.hydrogen_bonds else []
        x_bonds = analyzer.halogen_bonds if analyzer.halogen_bonds else []
        
        print(f"✓ Existing interaction detection still works:")
        print(f"  Hydrogen bonds: {len(h_bonds)}")
        print(f"  Halogen bonds: {len(x_bonds)}")
        
        # Should have some traditional interactions in a protein
        assert len(h_bonds) > 0, "No hydrogen bonds detected - existing functionality broken"
    
    def test_interaction_residue_consistency(self):
        """Test that all new interactions have consistent residue identifiers."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        all_interactions = []
        
        # Collect all new interactions
        if analyzer.pi_pi_interactions:
            all_interactions.extend(analyzer.pi_pi_interactions)
        if analyzer.carbonyl_interactions:
            all_interactions.extend(analyzer.carbonyl_interactions)
        if analyzer.n_pi_interactions:
            all_interactions.extend(analyzer.n_pi_interactions)
        
        if all_interactions:
            for interaction in all_interactions:
                # All should have donor and acceptor residues
                assert hasattr(interaction, 'get_donor_residue')
                assert hasattr(interaction, 'get_acceptor_residue')
                
                donor_res = interaction.get_donor_residue()
                acceptor_res = interaction.get_acceptor_residue()
                
                # Residue identifiers should be non-empty strings
                assert isinstance(donor_res, str) and len(donor_res) > 0
                assert isinstance(acceptor_res, str) and len(acceptor_res) > 0

                # Should have proper format (e.g., "A123ALA" or "DG2" for nucleotides)
                # Minimum 3 characters to allow for single-digit residue numbers
                assert len(donor_res) >= 3  # Minimum: "A1X" or "DG2"
                assert len(acceptor_res) >= 3
            
            print(f"✓ All {len(all_interactions)} new interactions have valid residue identifiers")


@pytest.mark.integration
class TestAnalyzerPerformance:
    """Test performance and scalability of new interaction analyzers."""
    
    def test_analysis_completes_reasonably_fast(self):
        """Test that analysis with new interactions completes in reasonable time."""
        import time
        
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        start_time = time.time()
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        end_time = time.time()
        
        analysis_time = end_time - start_time
        
        assert success, "Analysis failed"
        assert analysis_time < 60.0, f"Analysis took too long: {analysis_time:.2f}s"
        
        print(f"✓ Analysis completed in {analysis_time:.2f} seconds")
    
    def test_memory_usage_reasonable(self):
        """Test that new interactions don't cause excessive memory usage."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)
        
        success = analyzer.analyze_file('example_pdb_files/7nwd.pdb')
        assert success
        
        # Count total objects created
        total_interactions = 0
        if analyzer.pi_pi_interactions:
            total_interactions += len(analyzer.pi_pi_interactions)
        if analyzer.carbonyl_interactions:
            total_interactions += len(analyzer.carbonyl_interactions)
        if analyzer.n_pi_interactions:
            total_interactions += len(analyzer.n_pi_interactions)
        
        # Should be reasonable number of interactions
        assert total_interactions < 10000, f"Excessive interactions created: {total_interactions}"
        
        print(f"✓ Created {total_interactions} new interaction objects (reasonable)")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])