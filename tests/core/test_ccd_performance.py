"""
Performance test for CCD bond data loading and lookup.

This test measures the performance of loading CCD bond data into memory
and performing lookups for specific residues like PHE (phenylalanine).
"""

import time
import pytest
from typing import Dict, List

from hbat.ccd.ccd_analyzer import CCDDataManager


def format_time(seconds: float) -> str:
    """Format time in a human-readable way."""
    if seconds < 0.001:
        return f"{seconds * 1000000:.1f} μs"
    elif seconds < 1:
        return f"{seconds * 1000:.1f} ms"
    else:
        return f"{seconds:.2f} s"


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.ccd
@pytest.mark.performance
def test_ccd_data_manager_initialization():
    """Test CCDDataManager initialization performance."""
    start_time = time.perf_counter()
    ccd_manager = CCDDataManager(ccd_folder="~/.hbat/ccd-data")
    init_time = time.perf_counter() - start_time
    
    # Initialization should be very fast (< 1ms)
    assert init_time < 0.001, f"Initialization took too long: {format_time(init_time)}"
    assert ccd_manager is not None
    print(f"✅ Initialization time: {format_time(init_time)}")


@pytest.mark.integration 
@pytest.mark.slow
@pytest.mark.ccd
def test_ccd_files_availability():
    """Test CCD file availability check performance."""
    ccd_manager = CCDDataManager(ccd_folder="~/.hbat/ccd-data")
    
    start_time = time.perf_counter()
    files_ready = ccd_manager.ensure_files_exist()
    files_time = time.perf_counter() - start_time
    
    assert files_ready, "CCD files should be available"
    print(f"✅ File check time: {format_time(files_time)}")


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.ccd
@pytest.mark.performance
def test_ccd_bond_data_loading():
    """Test bond data loading performance."""
    ccd_manager = CCDDataManager(ccd_folder="~/.hbat/ccd-data")
    
    # Ensure files exist first
    assert ccd_manager.ensure_files_exist(), "CCD files must be available"
    
    start_time = time.perf_counter()
    success = ccd_manager.load_bonds_data()
    load_time = time.perf_counter() - start_time
    
    assert success, "Bond data loading should succeed"
    # Allow up to 30 seconds for loading large dataset
    assert load_time < 30.0, f"Bond data loading took too long: {format_time(load_time)}"
    print(f"✅ Bond data loading time: {format_time(load_time)}")


@pytest.mark.integration
@pytest.mark.slow  
@pytest.mark.ccd
@pytest.mark.performance
def test_ccd_atom_data_loading():
    """Test atom data loading performance."""
    ccd_manager = CCDDataManager(ccd_folder="~/.hbat/ccd-data")
    
    # Ensure files exist first
    assert ccd_manager.ensure_files_exist(), "CCD files must be available"
    
    start_time = time.perf_counter()
    success = ccd_manager.load_atoms_data()
    load_time = time.perf_counter() - start_time
    
    assert success, "Atom data loading should succeed"
    # Allow up to 20 seconds for loading large dataset
    assert load_time < 20.0, f"Atom data loading took too long: {format_time(load_time)}"
    print(f"✅ Atom data loading time: {format_time(load_time)}")


@pytest.mark.integration
@pytest.mark.ccd
@pytest.mark.performance
def test_phe_residue_bond_lookup():
    """Test PHE residue bond lookup performance and correctness."""
    ccd_manager = CCDDataManager(ccd_folder="~/.hbat/ccd-data")
    
    # Setup: ensure data is loaded
    assert ccd_manager.ensure_files_exist(), "CCD files must be available"
    assert ccd_manager.load_bonds_data(), "Bond data must be loaded"
    
    # Test single PHE lookup
    start_time = time.perf_counter()
    phe_bonds = ccd_manager.get_component_bonds("PHE")
    lookup_time = time.perf_counter() - start_time
    
    # Verify results
    assert len(phe_bonds) > 0, "PHE should have bonds"
    assert len(phe_bonds) >= 20, f"PHE should have at least 20 bonds, found {len(phe_bonds)}"
    
    # Performance assertion: lookup should be very fast (< 1ms)
    assert lookup_time < 0.001, f"PHE lookup took too long: {format_time(lookup_time)}"
    
    # Verify bond structure
    sample_bond = phe_bonds[0]
    required_keys = ['atom_id_1', 'atom_id_2', 'comp_id', 'value_order']
    for key in required_keys:
        assert key in sample_bond, f"Bond should contain {key}"
    
    assert sample_bond['comp_id'] == 'PHE', "Bond should belong to PHE"
    
    print(f"✅ PHE lookup time: {format_time(lookup_time)}")
    print(f"✅ PHE bonds found: {len(phe_bonds)}")


@pytest.mark.integration
@pytest.mark.ccd
@pytest.mark.performance
def test_multiple_residue_lookup_performance():
    """Test performance of multiple residue lookups."""
    ccd_manager = CCDDataManager(ccd_folder="~/.hbat/ccd-data")
    
    # Setup
    assert ccd_manager.ensure_files_exist(), "CCD files must be available"
    assert ccd_manager.load_bonds_data(), "Bond data must be loaded"
    
    test_residues = ["PHE", "TYR", "TRP", "ALA", "GLY"]
    
    start_time = time.perf_counter()
    lookup_results = {}
    for residue in test_residues:
        bonds = ccd_manager.get_component_bonds(residue)
        lookup_results[residue] = len(bonds)
    
    total_time = time.perf_counter() - start_time
    avg_time = total_time / len(test_residues)
    
    # All residues should have bonds
    for residue, bond_count in lookup_results.items():
        assert bond_count > 0, f"{residue} should have bonds"
    
    # Average lookup should be very fast
    assert avg_time < 0.001, f"Average lookup too slow: {format_time(avg_time)}"
    
    print(f"✅ Multiple lookup avg time: {format_time(avg_time)}")
    print(f"✅ Lookup results: {lookup_results}")


@pytest.mark.integration
@pytest.mark.ccd
@pytest.mark.performance
def test_phe_atom_specific_lookups():
    """Test PHE aromatic ring atom lookups."""
    ccd_manager = CCDDataManager(ccd_folder="~/.hbat/ccd-data")
    
    # Setup
    assert ccd_manager.ensure_files_exist(), "CCD files must be available"
    assert ccd_manager.load_atoms_data(), "Atom data must be loaded"
    assert ccd_manager.load_bonds_data(), "Bond data must be loaded"
    
    # PHE aromatic ring atoms
    phe_aromatic_atoms = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
    
    start_time = time.perf_counter()
    atom_results = {}
    for atom_id in phe_aromatic_atoms:
        atom_data = ccd_manager.get_atom_by_id("PHE", atom_id)
        bonds_involving = ccd_manager.get_bonds_involving_atom("PHE", atom_id)
        atom_results[atom_id] = {
            "found": atom_data is not None,
            "bonds_count": len(bonds_involving)
        }
    
    total_time = time.perf_counter() - start_time
    avg_time = total_time / len(phe_aromatic_atoms)
    
    # All aromatic atoms should be found
    for atom_id, result in atom_results.items():
        assert result["found"], f"PHE atom {atom_id} should be found"
        assert result["bonds_count"] > 0, f"PHE atom {atom_id} should have bonds"
    
    # Lookups should be fast
    assert avg_time < 0.01, f"Atom lookup too slow: {format_time(avg_time)}"
    
    print(f"✅ Atom lookup avg time: {format_time(avg_time)}")
    print(f"✅ Atom results: {atom_results}")


@pytest.mark.integration
@pytest.mark.ccd
def test_phe_component_summary():
    """Test PHE component summary generation."""
    ccd_manager = CCDDataManager(ccd_folder="~/.hbat/ccd-data")
    
    # Setup
    assert ccd_manager.ensure_files_exist(), "CCD files must be available"
    assert ccd_manager.load_atoms_data(), "Atom data must be loaded"
    assert ccd_manager.load_bonds_data(), "Bond data must be loaded"
    
    start_time = time.perf_counter()
    phe_summary = ccd_manager.get_component_summary("PHE")
    summary_time = time.perf_counter() - start_time
    
    # Verify summary content
    assert phe_summary["component_id"] == "PHE"
    assert phe_summary["available"], "PHE should be available"
    assert phe_summary["atom_count"] > 0, "PHE should have atoms"
    assert phe_summary["bond_count"] > 0, "PHE should have bonds"
    assert "bond_orders" in phe_summary
    assert "atoms" in phe_summary
    
    # PHE should have expected structure
    assert phe_summary["atom_count"] >= 20, f"PHE should have at least 20 atoms, found {phe_summary['atom_count']}"
    assert phe_summary["bond_count"] >= 20, f"PHE should have at least 20 bonds, found {phe_summary['bond_count']}"
    
    # Summary generation should be fast
    assert summary_time < 0.01, f"Summary generation too slow: {format_time(summary_time)}"
    
    print(f"✅ Summary generation time: {format_time(summary_time)}")
    print(f"✅ PHE summary: {phe_summary['atom_count']} atoms, {phe_summary['bond_count']} bonds")


def test_ccd_performance_comprehensive():
    """
    Comprehensive CCD performance test (for manual execution).
    
    This test measures:
    1. Time to initialize CCDDataManager
    2. Time to load bond data into memory
    3. Time to lookup bond information for PHE residue
    4. Time for multiple lookups to test caching efficiency
    """
    print("🧪 CCD Performance Test")
    print("=" * 50)
    
    # Test 1: Initialize CCDDataManager
    print("\n📋 Test 1: Initializing CCDDataManager")
    start_time = time.perf_counter()
    
    ccd_manager = CCDDataManager(ccd_folder="~/.hbat/ccd-data")
    
    init_time = time.perf_counter() - start_time
    print(f"   ⏱️  Initialization time: {format_time(init_time)}")
    
    # Test 2: Ensure files exist (download if needed)
    print("\n📥 Test 2: Ensuring CCD files exist")
    start_time = time.perf_counter()
    
    files_ready = ccd_manager.ensure_files_exist()
    
    files_time = time.perf_counter() - start_time
    print(f"   ⏱️  File check/download time: {format_time(files_time)}")
    print(f"   ✅ Files ready: {files_ready}")
    
    if not files_ready:
        print("❌ CCD files not available. Cannot continue test.")
        return
    
    # Test 3: Load bond data into memory
    print("\n💾 Test 3: Loading bond data into memory")
    start_time = time.perf_counter()
    
    bond_load_success = ccd_manager.load_bonds_data()
    
    bond_load_time = time.perf_counter() - start_time
    print(f"   ⏱️  Bond data loading time: {format_time(bond_load_time)}")
    print(f"   ✅ Bond data loaded: {bond_load_success}")
    
    if not bond_load_success:
        print("❌ Failed to load bond data. Cannot continue test.")
        return
    
    # Test 4: Load atom data into memory
    print("\n💾 Test 4: Loading atom data into memory")
    start_time = time.perf_counter()
    
    atom_load_success = ccd_manager.load_atoms_data()
    
    atom_load_time = time.perf_counter() - start_time
    print(f"   ⏱️  Atom data loading time: {format_time(atom_load_time)}")
    print(f"   ✅ Atom data loaded: {atom_load_success}")
    
    # Test 5: Single PHE lookup
    print("\n🔍 Test 5: Single PHE residue lookup")
    start_time = time.perf_counter()
    
    phe_bonds = ccd_manager.get_component_bonds("PHE")
    
    single_lookup_time = time.perf_counter() - start_time
    print(f"   ⏱️  Single lookup time: {format_time(single_lookup_time)}")
    print(f"   📊 PHE bonds found: {len(phe_bonds)}")
    
    # Display some PHE bond details
    if phe_bonds:
        print("   📝 Sample PHE bonds:")
        for i, bond in enumerate(phe_bonds[:3]):  # Show first 3 bonds
            atom1 = bond.get('atom_id_1', 'N/A')
            atom2 = bond.get('atom_id_2', 'N/A')
            order = bond.get('value_order', 'N/A')
            aromatic = bond.get('pdbx_aromatic_flag', 'N')
            print(f"      {i+1}. {atom1} - {atom2} ({order}, aromatic: {aromatic})")
    
    # Test 6: Multiple lookups for caching efficiency
    print("\n🔄 Test 6: Multiple lookups (caching test)")
    test_residues = ["PHE", "TYR", "TRP", "ALA", "GLY"]
    
    start_time = time.perf_counter()
    
    lookup_results = {}
    for residue in test_residues:
        bonds = ccd_manager.get_component_bonds(residue)
        lookup_results[residue] = len(bonds)
    
    multiple_lookup_time = time.perf_counter() - start_time
    avg_lookup_time = multiple_lookup_time / len(test_residues)
    
    print(f"   ⏱️  Total time for {len(test_residues)} lookups: {format_time(multiple_lookup_time)}")
    print(f"   ⏱️  Average lookup time: {format_time(avg_lookup_time)}")
    print("   📊 Lookup results:")
    for residue, bond_count in lookup_results.items():
        print(f"      {residue}: {bond_count} bonds")
    
    # Test 7: Atom-specific lookups
    print("\n🎯 Test 7: Atom-specific lookups for PHE")
    phe_atoms = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]  # Aromatic ring atoms
    
    start_time = time.perf_counter()
    
    atom_lookup_results = {}
    for atom_id in phe_atoms:
        atom_data = ccd_manager.get_atom_by_id("PHE", atom_id)
        bonds_involving_atom = ccd_manager.get_bonds_involving_atom("PHE", atom_id)
        atom_lookup_results[atom_id] = {
            "found": atom_data is not None,
            "bonds_count": len(bonds_involving_atom)
        }
    
    atom_lookup_time = time.perf_counter() - start_time
    avg_atom_lookup_time = atom_lookup_time / len(phe_atoms)
    
    print(f"   ⏱️  Total time for {len(phe_atoms)} atom lookups: {format_time(atom_lookup_time)}")
    print(f"   ⏱️  Average atom lookup time: {format_time(avg_atom_lookup_time)}")
    print("   📊 Atom lookup results:")
    for atom_id, result in atom_lookup_results.items():
        status = "✅" if result["found"] else "❌"
        print(f"      {atom_id}: {status} found, {result['bonds_count']} bonds")
    
    # Test 8: Component summary
    print("\n📈 Test 8: Component summary generation")
    start_time = time.perf_counter()
    
    phe_summary = ccd_manager.get_component_summary("PHE")
    
    summary_time = time.perf_counter() - start_time
    print(f"   ⏱️  Summary generation time: {format_time(summary_time)}")
    print("   📊 PHE Summary:")
    print(f"      Atoms: {phe_summary['atom_count']}")
    print(f"      Bonds: {phe_summary['bond_count']}")
    print(f"      Aromatic bonds: {phe_summary['aromatic_bonds']}")
    print(f"      Bond orders: {phe_summary['bond_orders']}")
    print(f"      Available atoms: {', '.join(phe_summary['atoms'][:10])}...")  # Show first 10
    
    # Performance Summary
    print("\n📊 Performance Summary")
    print("=" * 50)
    total_time = (init_time + files_time + bond_load_time + atom_load_time + 
                  single_lookup_time + multiple_lookup_time + atom_lookup_time + summary_time)
    
    print(f"   🚀 Total test time: {format_time(total_time)}")
    print(f"   💾 Data loading time: {format_time(bond_load_time + atom_load_time)}")
    print(f"   🔍 Lookup performance: {format_time(avg_lookup_time)} avg per residue")
    print(f"   🎯 Atom lookup performance: {format_time(avg_atom_lookup_time)} avg per atom")
    
    # Performance analysis
    print("\n🎯 Performance Analysis")
    print("=" * 50)
    
    if bond_load_time < 5.0:
        print("   ✅ Excellent: Bond data loading < 5 seconds")
    elif bond_load_time < 15.0:
        print("   ⚠️  Good: Bond data loading < 15 seconds")
    else:
        print("   ❌ Slow: Bond data loading > 15 seconds")
    
    if avg_lookup_time < 0.001:
        print("   ✅ Excellent: Average lookup < 1ms (sub-millisecond)")
    elif avg_lookup_time < 0.01:
        print("   ⚠️  Good: Average lookup < 10ms")
    else:
        print("   ❌ Slow: Average lookup > 10ms")
    
    if avg_atom_lookup_time < 0.001:
        print("   ✅ Excellent: Atom lookup < 1ms")
    elif avg_atom_lookup_time < 0.005:
        print("   ⚠️  Good: Atom lookup < 5ms")
    else:
        print("   ❌ Slow: Atom lookup > 5ms")
    
    print(f"\n🎉 Test completed successfully!")
    return phe_summary


def main():
    """Run the CCD performance test."""
    try:
        result = test_ccd_performance()
        return 0
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())