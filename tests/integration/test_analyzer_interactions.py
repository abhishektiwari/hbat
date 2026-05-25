"""
Interaction detection integration tests for analyzer.

Tests verify detection and validation of all interaction types:
cooperativity chains, π-π stacking, carbonyl-carbonyl, and n→π* interactions.
"""

import pytest
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters
from hbat.core.interactions import PiPiInteraction, CarbonylInteraction, NPiInteraction


@pytest.mark.integration
@pytest.mark.requires_pdb_files
class TestAnalyzerCooperativityIntegration:
    """Test cooperativity analysis integration."""

    def test_cooperativity_chain_detection(self, sample_pdb_file):
        """Test cooperativity chain detection integration."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success

        chains = analyzer.cooperativity_chains
        assert isinstance(chains, list), "Should return list of chains"

        # If chains are found, validate them
        for chain in chains[:5]:  # Check first 5 chains
            assert hasattr(chain, "interactions"), "Chain should have interactions"
            assert hasattr(chain, "chain_length"), "Chain should have length"
            assert hasattr(chain, "chain_type"), "Chain should have type"

            # Validate chain structure
            assert len(chain.interactions) == chain.chain_length, (
                "Length should match interaction count"
            )
            assert chain.chain_length >= 0, "Chain length should be non-negative"

            # Validate interactions in chain
            for interaction in chain.interactions:
                assert hasattr(interaction, "interaction_type"), (
                    "Chain interaction should have type"
                )
                assert interaction.interaction_type in [
                    "H-Bond",
                    "X-Bond",
                    "π–Inter",
                ], f"Unknown interaction type: {interaction.interaction_type}"

    def test_cooperativity_statistics_integration(self, sample_pdb_file):
        """Test cooperativity statistics integration."""
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(sample_pdb_file)
        assert success

        chains = analyzer.cooperativity_chains
        summary = analyzer.get_summary()

        # Verify statistics consistency
        if "cooperativity_chains" in summary:
            assert summary["cooperativity_chains"]["count"] == len(chains), (
                "Summary should match actual chain count"
            )


@pytest.mark.integration
class TestPiPiStackingAnalyzer:
    """Test π-π stacking interaction analyzer integration."""

    def test_pi_pi_detection_7nwd(self):
        """Test π-π stacking detection with 7NWD.pdb structure."""
        # 7NWD contains aromatic residues suitable for π-π stacking analysis
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success, "Failed to analyze 7NWD.pdb"

        # Check that π-π interactions were detected
        assert hasattr(analyzer, "pi_pi_interactions"), "π-π interactions not found"

        pi_pi_interactions = analyzer.pi_pi_interactions
        assert len(pi_pi_interactions) > 0, "No π-π interactions detected in 7NWD"

        # Validate first π-π interaction
        first_interaction = pi_pi_interactions[0]
        assert isinstance(first_interaction, PiPiInteraction)

        # Check required properties
        assert hasattr(first_interaction, "ring1_atoms")
        assert hasattr(first_interaction, "ring2_atoms")
        assert hasattr(first_interaction, "distance")
        assert hasattr(first_interaction, "plane_angle")
        assert hasattr(first_interaction, "stacking_type")

        # Validate distance range (typical π-π stacking: 3.3-6.0 Å)
        assert 3.0 <= first_interaction.distance <= 6.5

        # Validate angle range (0-180°, angles can be obtuse)
        assert 0.0 <= first_interaction.plane_angle <= 180.0

        # Validate stacking type
        valid_types = ["parallel", "T-shaped", "offset"]
        assert first_interaction.stacking_type in valid_types

        print(f"✓ Detected {len(pi_pi_interactions)} π-π interactions in 7NWD")
        print(
            f"  First interaction: {first_interaction.ring1_residue} - {first_interaction.ring2_residue}"
        )
        print(
            f"  Distance: {first_interaction.distance:.2f}Å, Type: {first_interaction.stacking_type}"
        )

    def test_pi_pi_different_stacking_types(self):
        """Test detection of different π-π stacking types."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        if hasattr(analyzer, "pi_pi_interactions") and analyzer.pi_pi_interactions:
            stacking_types = {
                interaction.stacking_type for interaction in analyzer.pi_pi_interactions
            }

            # Should detect at least one type
            assert len(stacking_types) > 0

            # All types should be valid
            valid_types = {"parallel", "T-shaped", "offset"}
            assert stacking_types.issubset(valid_types)

            print(f"✓ π-π stacking types found: {sorted(stacking_types)}")

    def test_pi_pi_aromatic_residue_types(self):
        """Test π-π interactions involve aromatic residues."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        if hasattr(analyzer, "pi_pi_interactions") and analyzer.pi_pi_interactions:
            # Include nucleotides and aromatic amino acids
            aromatic_residues = {
                "PHE",
                "TYR",
                "TRP",
                "HIS",
                "DA",
                "DC",
                "DG",
                "DT",
                "A",
                "C",
                "G",
                "T",
                "U",
            }

            for interaction in analyzer.pi_pi_interactions:
                # Extract residue types from residue strings
                ring1_type = interaction.ring1_type
                ring2_type = interaction.ring2_type

                # Should involve aromatic residues or nucleotides
                assert (
                    ring1_type in aromatic_residues or ring2_type in aromatic_residues
                )

            print("✓ All π-π interactions involve aromatic residues")


@pytest.mark.integration
class TestCarbonylCarbonylAnalyzer:
    """Test carbonyl-carbonyl interaction analyzer integration."""

    def test_carbonyl_detection_with_structure(self):
        """Test carbonyl-carbonyl detection with protein structure."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        # Use 6RSA which has known carbonyl interactions
        success = analyzer.analyze_file("example_pdb_files/6rsa.pdb")
        assert success, "Failed to analyze structure for carbonyl interactions"

        # Check that carbonyl interactions were detected
        assert hasattr(analyzer, "carbonyl_interactions"), (
            "Carbonyl interactions not found"
        )

        if analyzer.carbonyl_interactions:
            carbonyl_interactions = analyzer.carbonyl_interactions

            # Validate first carbonyl interaction
            first_interaction = carbonyl_interactions[0]
            assert isinstance(first_interaction, CarbonylInteraction)

            # Check required properties
            assert hasattr(first_interaction, "donor_carbon")
            assert hasattr(first_interaction, "donor_oxygen")
            assert hasattr(first_interaction, "acceptor_carbon")
            assert hasattr(first_interaction, "acceptor_oxygen")
            assert hasattr(first_interaction, "distance")
            assert hasattr(first_interaction, "burgi_dunitz_angle")

            # Validate distance range (typical n→π*: 2.8-4.0 Å)
            assert 2.5 <= first_interaction.distance <= 4.5

            # Validate Bürgi-Dunitz angle range (95-125°)
            assert 90.0 <= first_interaction.burgi_dunitz_angle <= 130.0

            print(f"✓ Detected {len(carbonyl_interactions)} carbonyl interactions")
            print(
                f"  First interaction: {first_interaction.donor_residue} - {first_interaction.acceptor_residue}"
            )
            print(
                f"  Distance: {first_interaction.distance:.2f}Å, Angle: {first_interaction.burgi_dunitz_angle:.1f}°"
            )
        else:
            print("ℹ No carbonyl interactions detected in this structure")

    def test_carbonyl_angle_calculations(self):
        """Test that Bürgi-Dunitz angles are calculated correctly."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/6rsa.pdb")
        assert success

        if (
            hasattr(analyzer, "carbonyl_interactions")
            and analyzer.carbonyl_interactions
        ):
            for interaction in analyzer.carbonyl_interactions:
                # Bürgi-Dunitz angle should be in the expected range
                angle = interaction.burgi_dunitz_angle
                assert 90.0 <= angle <= 130.0, (
                    f"Invalid Bürgi-Dunitz angle: {angle:.1f}°"
                )

                # Angle property should return radians
                angle_rad = interaction.angle
                expected_rad = angle * 3.14159 / 180.0
                assert abs(angle_rad - expected_rad) < 0.1

            print("✓ All Bürgi-Dunitz angles in valid range (90-130°)")

    def test_carbonyl_backbone_sidechain_classification(self):
        """Test carbonyl interaction backbone/sidechain classification."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/6rsa.pdb")
        assert success

        if (
            hasattr(analyzer, "carbonyl_interactions")
            and analyzer.carbonyl_interactions
        ):
            valid_classifications = {
                "backbone-backbone",
                "backbone-sidechain",
                "sidechain-backbone",
                "sidechain-sidechain",
            }

            for interaction in analyzer.carbonyl_interactions:
                assert interaction.interaction_classification in valid_classifications

            print("✓ All carbonyl interactions have valid classifications")


@pytest.mark.integration
class TestNPiAnalyzer:
    """Test n→π* interaction analyzer integration."""

    def test_n_pi_detection_with_structure(self):
        """Test n→π* interaction detection with protein structure."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success, "Failed to analyze structure for n→π* interactions"

        # Check that n→π* interactions were detected
        assert hasattr(analyzer, "n_pi_interactions"), "n→π* interactions not found"

        if analyzer.n_pi_interactions:
            n_pi_interactions = analyzer.n_pi_interactions

            # Validate first n→π* interaction
            first_interaction = n_pi_interactions[0]
            assert isinstance(first_interaction, NPiInteraction)

            # Check required properties
            assert hasattr(first_interaction, "lone_pair_atom")
            assert hasattr(first_interaction, "pi_center")
            assert hasattr(first_interaction, "pi_atoms")
            assert hasattr(first_interaction, "distance")
            assert hasattr(first_interaction, "angle_to_plane")
            assert hasattr(first_interaction, "subtype")

            # Validate distance range (typical n→π*: 3.0-5.0 Å)
            assert 2.5 <= first_interaction.distance <= 5.5

            # Validate angle to plane range (0-90°)
            assert 0.0 <= first_interaction.angle_to_plane <= 90.0

            # Validate donor element
            valid_elements = {"O", "N", "S"}
            assert first_interaction.donor_element in valid_elements

            print(f"✓ Detected {len(n_pi_interactions)} n→π* interactions")
            print(
                f"  First interaction: {first_interaction.donor_residue} - {first_interaction.acceptor_residue}"
            )
            print(
                f"  Distance: {first_interaction.distance:.2f}Å, Angle: {first_interaction.angle_to_plane:.1f}°"
            )
            print(f"  Subtype: {first_interaction.subtype}")
        else:
            print("ℹ No n→π* interactions detected in this structure")

    def test_n_pi_donor_element_types(self):
        """Test n→π* interactions with different donor elements."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        if hasattr(analyzer, "n_pi_interactions") and analyzer.n_pi_interactions:
            donor_elements = {
                interaction.donor_element for interaction in analyzer.n_pi_interactions
            }

            # Should only have valid donor elements
            valid_elements = {"O", "N", "S"}
            assert donor_elements.issubset(valid_elements)

            print(f"✓ n→π* donor elements found: {sorted(donor_elements)}")

    def test_n_pi_subtype_classification(self):
        """Test n→π* interaction subtype classification."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        if hasattr(analyzer, "n_pi_interactions") and analyzer.n_pi_interactions:
            subtypes = {
                interaction.subtype for interaction in analyzer.n_pi_interactions
            }

            # Should have meaningful subtypes
            assert len(subtypes) > 0
            for subtype in subtypes:
                assert len(subtype) > 0  # Non-empty subtype
                assert isinstance(subtype, str)

            print(f"✓ n→π* interaction subtypes: {sorted(subtypes)}")


@pytest.mark.integration
class TestWaterBridgeDetection:
    """Test water bridge interaction detection and validation."""

    def test_water_bridge_detection_with_structure(self):
        """Test water bridge detection with protein structure."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success, "Failed to analyze structure for water bridges"

        # Check that water bridges were detected
        assert hasattr(analyzer, "water_bridges"), "Water bridges not found"

        water_bridges = analyzer.water_bridges
        assert isinstance(water_bridges, list), "Should return list of water bridges"

        # If water bridges are found, validate them
        for wb in water_bridges[:5]:  # Check first 5
            assert hasattr(wb, "donor_atom"), "Water bridge should have donor atom"
            assert hasattr(wb, "acceptor_atom"), "Water bridge should have acceptor atom"
            assert hasattr(wb, "bridge_path"), "Water bridge should have bridge path"
            assert hasattr(wb, "total_distance"), "Water bridge should have total distance"

            # Validate distances
            assert wb.total_distance > 0, "Total distance should be positive"

            # Validate bridge path (list of HydrogenBond objects)
            assert isinstance(wb.bridge_path, list), "Bridge path should be a list"
            assert len(wb.bridge_path) > 0, "Bridge path should have at least one bond"

        print(f"✓ Detected {len(water_bridges)} water bridges in structure")

    def test_water_bridge_properties_validation(self):
        """Test that water bridges have valid properties."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        water_bridges = analyzer.water_bridges

        if water_bridges:
            for wb in water_bridges:
                # Test donor and acceptor residue format
                donor_res = wb.get_donor_residue()
                acceptor_res = wb.get_acceptor_residue()

                assert isinstance(donor_res, str) and len(donor_res) > 0
                assert isinstance(acceptor_res, str) and len(acceptor_res) > 0

                # Test interaction type
                assert wb.get_interaction_type() == "water_bridge"

                # Test geometry (should have valid distances)
                da_distance = wb.get_donor_acceptor_distance()
                assert da_distance > 0, "Donor-acceptor distance should be positive"

            print(f"✓ All {len(water_bridges)} water bridges have valid properties")

    def test_water_bridge_geometry_constraints(self):
        """Test water bridge geometry follows expected constraints."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        water_bridges = analyzer.water_bridges

        if water_bridges:
            for wb in water_bridges:
                # Total distance should be reasonable (sum of bond distances)
                total = wb.total_distance
                bridge_length = len(wb.bridge_path)

                # Each hop should have reasonable distance
                # Typical hydrogen bond: 2.5-3.5 Å
                # Water bridges usually 2-4 hops: 5-14 Å total
                assert 4.0 <= total <= 15.0, (
                    f"Water bridge total distance {total:.2f}Å should be reasonable"
                )

                # Should have at least 1 bond in path
                assert bridge_length >= 1, "Bridge should have at least one bond"

            print(f"✓ All water bridges follow geometry constraints")

    def test_water_bridge_residue_consistency(self):
        """Test water bridge residue identifiers are consistent."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        water_bridges = analyzer.water_bridges

        if water_bridges:
            for wb in water_bridges:
                # Both donor and acceptor should be valid residues
                donor_res = wb.get_donor_residue()
                acceptor_res = wb.get_acceptor_residue()

                # Format check: should have chain:residue:resnum format or similar
                assert ":" in donor_res or len(donor_res) >= 3
                assert ":" in acceptor_res or len(acceptor_res) >= 3

                # Should not be the same
                assert donor_res != acceptor_res, (
                    "Donor and acceptor residues should be different"
                )

            print(f"✓ All water bridge residues are consistent")


@pytest.mark.integration
class TestLigandInteractionDetection:
    """Test ligand interaction detection and validation."""

    def test_ligand_interaction_detection(self):
        """Test ligand interaction detection with structure."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        # Use a structure likely to have ligands
        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success, "Failed to analyze structure"

        # Check that ligand interactions object exists
        assert hasattr(analyzer, "ligand_interactions"), (
            "Ligand interactions not found"
        )

        ligand_interactions = analyzer.ligand_interactions
        assert ligand_interactions is not None

        # Check interactions list
        assert hasattr(ligand_interactions, "interactions")
        interactions_list = ligand_interactions.interactions
        assert isinstance(interactions_list, list)

        # If ligands are found, validate them
        if interactions_list:
            for lig_int in interactions_list[:5]:  # Check first 5
                assert hasattr(lig_int, "donor") or hasattr(lig_int, "acceptor")

            print(f"✓ Detected {len(interactions_list)} ligand interactions")

    def test_ligand_info_structure(self):
        """Test that ligand_info has proper structure."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        ligand_interactions = analyzer.ligand_interactions
        assert hasattr(ligand_interactions, "ligand_info")

        ligand_info = ligand_interactions.ligand_info
        assert isinstance(ligand_info, dict)

        # Check structure of ligand info entries
        for ligand_key, lig_data in ligand_info.items():
            # Key should be in format like "A:GTP:1" or similar
            assert isinstance(ligand_key, str)
            assert len(ligand_key) > 0

            # Value should be a dictionary with expected fields
            assert isinstance(lig_data, dict)
            assert "count" in lig_data, "Ligand info should have count"
            assert "chain" in lig_data, "Ligand info should have chain"
            assert "name" in lig_data, "Ligand info should have name"
            assert "seq" in lig_data, "Ligand info should have seq"

            # Values should have correct types
            assert isinstance(lig_data["count"], int)
            assert lig_data["count"] > 0

        if ligand_info:
            print(f"✓ Ligand info structure valid for {len(ligand_info)} ligands")

    def test_ligand_classification(self):
        """Test ligand classification and identification."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        ligand_interactions = analyzer.ligand_interactions
        ligand_info = ligand_interactions.ligand_info

        if ligand_info:
            for ligand_key, lig_data in ligand_info.items():
                # Ligand name should be non-empty
                name = lig_data.get("name", "")
                assert len(name) > 0, "Ligand should have a name"

                # Chain should be valid (single character or empty)
                chain = lig_data.get("chain", "")
                assert isinstance(chain, str)

                # Sequence number should be valid
                seq = lig_data.get("seq", 0)
                assert isinstance(seq, int)
                assert seq >= 0 or seq == 0, "Sequence number should be non-negative"

            print(f"✓ All ligands properly classified")

    def test_ligand_residue_identification(self):
        """Test ligand residue identification."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        ligand_interactions = analyzer.ligand_interactions

        if ligand_interactions.interactions:
            for lig_int in ligand_interactions.interactions:
                # Should be able to get donor and acceptor residues
                if hasattr(lig_int, "get_donor_residue"):
                    donor_res = lig_int.get_donor_residue()
                    assert isinstance(donor_res, str)
                    assert len(donor_res) > 0

                if hasattr(lig_int, "get_acceptor_residue"):
                    acceptor_res = lig_int.get_acceptor_residue()
                    assert isinstance(acceptor_res, str)
                    assert len(acceptor_res) > 0

            print(f"✓ All ligand interactions have valid residue IDs")


@pytest.mark.integration
class TestWaterBridgeLigandIntegration:
    """Test integration of water bridges and ligand interactions."""

    def test_water_bridges_with_ligand_interactions(self):
        """Test that water bridges work correctly with ligand interactions."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        # Both should be available
        assert hasattr(analyzer, "water_bridges")
        assert hasattr(analyzer, "ligand_interactions")

        water_bridges = analyzer.water_bridges
        ligand_interactions = analyzer.ligand_interactions

        # Get counts
        wb_count = len(water_bridges)
        lig_int_count = (
            len(ligand_interactions.interactions)
            if ligand_interactions
            else 0
        )

        print(f"✓ Water bridges: {wb_count}, Ligand interactions: {lig_int_count}")

        # Both should be valid lists (even if empty)
        assert isinstance(water_bridges, list)
        assert ligand_interactions is not None

    def test_all_interaction_types_together(self):
        """Test that all interaction types can coexist without interference."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        # Count all interaction types
        all_counts = {
            "hydrogen_bonds": len(analyzer.hydrogen_bonds),
            "halogen_bonds": len(analyzer.halogen_bonds),
            "pi_interactions": len(analyzer.pi_interactions),
            "pi_pi_interactions": len(analyzer.pi_pi_interactions),
            "carbonyl_interactions": len(analyzer.carbonyl_interactions),
            "n_pi_interactions": len(analyzer.n_pi_interactions),
            "water_bridges": len(analyzer.water_bridges),
            "ligand_interactions": (
                len(analyzer.ligand_interactions.interactions)
                if analyzer.ligand_interactions
                else 0
            ),
            "cooperativity_chains": len(analyzer.cooperativity_chains),
        }

        # All should be non-negative
        for interaction_type, count in all_counts.items():
            assert count >= 0, f"{interaction_type} count should be non-negative"

        # At least some should be detected
        total_interactions = sum(all_counts.values())
        assert total_interactions > 0, "Should detect some interactions"

        print("✓ All interaction types coexist without interference:")
        for interaction_type, count in all_counts.items():
            print(f"  {interaction_type}: {count}")


@pytest.mark.integration
class TestNewInteractionsIntegration:
    """Test integration of all new interaction types together."""

    def test_all_new_interactions_detected(self):
        """Test that all new interaction types can be detected together."""
        params = AnalysisParameters()

        # Test with 7nwd for π-π interactions
        analyzer_pi = MolecularInteractionAnalyzer(parameters=params)
        success_pi = analyzer_pi.analyze_file("example_pdb_files/7nwd.pdb")
        assert success_pi, "Failed to analyze 7nwd.pdb"

        # Test with 6rsa for carbonyl interactions
        analyzer_carbonyl = MolecularInteractionAnalyzer(parameters=params)
        success_carbonyl = analyzer_carbonyl.analyze_file("example_pdb_files/6rsa.pdb")
        assert success_carbonyl, "Failed to analyze 6rsa.pdb"

        # Use 7nwd as primary analyzer for summary
        analyzer = analyzer_pi

        # Check that analyzer has all new interaction attributes
        assert hasattr(analyzer, "pi_pi_interactions"), "π-π interactions not available"
        assert hasattr(analyzer, "carbonyl_interactions"), (
            "Carbonyl interactions not available"
        )
        assert hasattr(analyzer, "n_pi_interactions"), "n→π* interactions not available"

        # Count total new interactions
        pi_pi_count = (
            len(analyzer.pi_pi_interactions) if analyzer.pi_pi_interactions else 0
        )
        carbonyl_count = (
            len(analyzer.carbonyl_interactions) if analyzer.carbonyl_interactions else 0
        )
        n_pi_count = (
            len(analyzer.n_pi_interactions) if analyzer.n_pi_interactions else 0
        )

        total_new_interactions = pi_pi_count + carbonyl_count + n_pi_count

        print("✓ New interaction detection summary for 7NWD:")
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

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
        assert success

        # Traditional interactions should still work
        assert hasattr(analyzer, "hydrogen_bonds"), "Hydrogen bonds not available"
        assert hasattr(analyzer, "halogen_bonds"), "Halogen bonds not available"

        # Should detect traditional interactions
        h_bonds = analyzer.hydrogen_bonds if analyzer.hydrogen_bonds else []
        x_bonds = analyzer.halogen_bonds if analyzer.halogen_bonds else []

        print("✓ Existing interaction detection still works:")
        print(f"  Hydrogen bonds: {len(h_bonds)}")
        print(f"  Halogen bonds: {len(x_bonds)}")

        # Should have some traditional interactions in a protein
        assert len(h_bonds) > 0, (
            "No hydrogen bonds detected - existing functionality broken"
        )

    def test_interaction_residue_consistency(self):
        """Test that all new interactions have consistent residue identifiers."""
        params = AnalysisParameters()
        analyzer = MolecularInteractionAnalyzer(parameters=params)

        success = analyzer.analyze_file("example_pdb_files/7nwd.pdb")
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
                assert hasattr(interaction, "get_donor_residue")
                assert hasattr(interaction, "get_acceptor_residue")

                donor_res = interaction.get_donor_residue()
                acceptor_res = interaction.get_acceptor_residue()

                # Residue identifiers should be non-empty strings
                assert isinstance(donor_res, str) and len(donor_res) > 0
                assert isinstance(acceptor_res, str) and len(acceptor_res) > 0

                # Should have proper format (e.g., "A123ALA" or "DG2" for nucleotides)
                # Minimum 3 characters to allow for single-digit residue numbers
                assert len(donor_res) >= 3  # Minimum: "A1X" or "DG2"
                assert len(acceptor_res) >= 3

            print(
                f"✓ All {len(all_interactions)} new interactions have valid residue identifiers"
            )
