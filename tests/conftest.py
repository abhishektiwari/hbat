"""
Pytest configuration for HBAT tests.
"""

import pytest
import sys
import os

# Add the parent directory to Python path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "gui: marks tests that require GUI components"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests that require sample files"
    )

def find_sample_pdb_file():
    """Find the sample PDB file regardless of working directory."""
    sample_paths = [
        os.path.join(os.path.dirname(__file__), "..", "example_pdb_files", "6rsa.pdb"),
        "../example_pdb_files/6rsa.pdb",  # When running from tests/
        "example_pdb_files/6rsa.pdb",     # When running from project root
    ]
    
    for path in sample_paths:
        if os.path.exists(path):
            return path
    return None

@pytest.fixture
def sample_pdb_file():
    """Provide path to sample PDB file."""
    file_path = find_sample_pdb_file()
    if not file_path:
        pytest.skip("Sample PDB file (6rsa.pdb) not found")
    return file_path

@pytest.fixture
def analyzer():
    """Provide a configured HBondAnalyzer instance."""
    from hbat.core.analysis import HBondAnalyzer
    return HBondAnalyzer()

@pytest.fixture
def analysis_parameters():
    """Fixture providing standard analysis parameters."""
    from hbat.core.analysis import AnalysisParameters
    return AnalysisParameters(
        hb_distance_cutoff=3.5,
        hb_angle_cutoff=120.0,
        hb_donor_acceptor_cutoff=4.0,
        analysis_mode="complete"
    )

@pytest.fixture
def strict_analysis_parameters():
    """Fixture providing strict analysis parameters."""
    from hbat.core.analysis import AnalysisParameters
    return AnalysisParameters(
        hb_distance_cutoff=3.0,
        hb_angle_cutoff=130.0,
        hb_donor_acceptor_cutoff=3.7,
        analysis_mode="complete"
    )

# Constants for expected results with 6RSA.pdb
class ExpectedResults:
    """Expected results for 6RSA.pdb analysis."""
    MIN_HYDROGEN_BONDS = 100
    MIN_PI_INTERACTIONS = 5
    MIN_TOTAL_INTERACTIONS = 50
    MIN_COOPERATIVITY_CHAINS = 5
    MIN_ATOMS = 2000
    MIN_HYDROGENS = 1000
    MIN_RESIDUES = 100

# Test data validation utilities
def validate_interaction_attributes(interaction, interaction_type):
    """Validate that an interaction has required attributes."""
    required_attrs = ['distance', 'angle', 'interaction_type']
    
    for attr in required_attrs:
        assert hasattr(interaction, attr), f"Interaction missing {attr} attribute"
    
    assert interaction.interaction_type == interaction_type, \
        f"Expected {interaction_type}, got {interaction.interaction_type}"
    
    assert interaction.distance > 0, "Distance should be positive"
    assert 0 <= interaction.angle <= 3.14159, "Angle should be in radians [0, π]"

def validate_hydrogen_bond(hbond):
    """Validate hydrogen bond specific attributes."""
    validate_interaction_attributes(hbond, "hydrogen_bond")
    
    required_attrs = ['donor', 'hydrogen', 'acceptor', 'donor_residue', 'acceptor_residue']
    for attr in required_attrs:
        assert hasattr(hbond, attr), f"H-bond missing {attr} attribute"

def validate_pi_interaction(pi_interaction):
    """Validate π interaction specific attributes."""
    validate_interaction_attributes(pi_interaction, "pi_interaction")
    
    required_attrs = ['donor', 'hydrogen', 'pi_center', 'donor_residue', 'pi_residue']
    for attr in required_attrs:
        assert hasattr(pi_interaction, attr), f"π-interaction missing {attr} attribute"

def validate_cooperativity_chain(chain):
    """Validate cooperativity chain attributes."""
    required_attrs = ['interactions', 'chain_length', 'chain_type']
    for attr in required_attrs:
        assert hasattr(chain, attr), f"Chain missing {attr} attribute"
    
    assert len(chain.interactions) > 0, "Chain should have at least one interaction"
    assert chain.chain_length == len(chain.interactions), "Length should match interactions count"
    
    # Validate chain_type format - should be like "H-Bond -> X-Bond" or "H-Bond -> H-Bond"
    valid_components = {"H-Bond", "X-Bond", "π-Int", "Unknown", "Empty"}
    
    if chain.chain_type == "Empty":
        assert len(chain.interactions) == 0, "Empty chain should have no interactions"
    else:
        # Split by " -> " and validate each component
        components = chain.chain_type.split(" -> ")
        assert all(comp in valid_components for comp in components), \
            f"Invalid chain type components in '{chain.chain_type}'. Valid components: {valid_components}"
        assert len(components) == len(chain.interactions), \
            f"Chain type components ({len(components)}) should match number of interactions ({len(chain.interactions)})"
    
    # Validate each interaction in the chain
    for interaction in chain.interactions:
        assert hasattr(interaction, 'interaction_type'), "Chain interaction missing type"
        assert interaction.interaction_type in ["hydrogen_bond", "halogen_bond", "pi_interaction"], \
            f"Unknown interaction type: {interaction.interaction_type}"