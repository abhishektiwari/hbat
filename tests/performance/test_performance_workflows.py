"""
Performance tests using large structures.

These tests verify performance characteristics and scalability
of molecular interaction analysis using large PDB files like 4jsv.pdb.
"""

import pytest
import time
import tempfile
import os
import json
import psutil
from hbat.core.analyzer import MolecularInteractionAnalyzer
from hbat.constants.parameters import AnalysisParameters


@pytest.mark.performance
@pytest.mark.slow
@pytest.mark.requires_pdb_files
@pytest.mark.integration
class TestLargeStructurePerformance:
    """Test performance with large structures using 4jsv.pdb."""
    
    @pytest.fixture
    def large_pdb_file(self, test_pdb_dir):
        """Provide path to large PDB file (4jsv.pdb)."""
        large_file_path = os.path.join(test_pdb_dir, "4jsv.pdb")
        if not os.path.exists(large_file_path):
            pytest.skip("Large PDB file (4jsv.pdb) not found")
        return large_file_path
    
    def test_large_structure_analysis_performance(self, large_pdb_file):
        """Test analysis performance with large structure (4jsv.pdb)."""
        analyzer = MolecularInteractionAnalyzer()
        
        # Measure analysis time
        start_time = time.time()
        success = analyzer.analyze_file(large_pdb_file)
        analysis_time = time.time() - start_time
        
        assert success, "Large structure analysis should succeed"
        
        # Performance requirements (adjust based on expectations)
        assert analysis_time < 300.0, f"Analysis took too long: {analysis_time:.2f}s (should be < 5 minutes)"
        
        # Verify substantial results for large structure
        summary = analyzer.get_summary()
        assert summary['hydrogen_bonds']['count'] > 100, "Large structure should have many H-bonds"
        assert summary['total_interactions'] > 200, "Large structure should have many interactions"
        
        # Log performance metrics
        print(f"\n4jsv.pdb Performance Metrics:")
        print(f"  Analysis time: {analysis_time:.2f}s")
        print(f"  Hydrogen bonds: {summary['hydrogen_bonds']['count']}")
        print(f"  Halogen bonds: {summary['halogen_bonds']}")
        print(f"  π interactions: {summary['pi_interactions']}")
        print(f"  Total interactions: {summary['total_interactions']}")
        print(f"  Rate: {summary['total_interactions']/analysis_time:.1f} interactions/second")
    
    def test_large_structure_memory_usage(self, large_pdb_file):
        """Test memory usage with large structure."""
        process = psutil.Process(os.getpid())
        
        # Measure initial memory
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        analyzer = MolecularInteractionAnalyzer()
        success = analyzer.analyze_file(large_pdb_file)
        assert success
        
        # Measure peak memory
        peak_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_increase = peak_memory - initial_memory
        
        # Memory usage should be reasonable (adjust threshold as needed)
        assert memory_increase < 1000, f"Memory usage too high: {memory_increase:.1f}MB"
        
        print(f"\n4jsv.pdb Memory Usage:")
        print(f"  Initial memory: {initial_memory:.1f}MB")
        print(f"  Peak memory: {peak_memory:.1f}MB")
        print(f"  Memory increase: {memory_increase:.1f}MB")
        
        # Clean up to test memory release
        del analyzer
        
        # Force garbage collection
        import gc
        gc.collect()
        
        # Check memory after cleanup
        final_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_retained = final_memory - initial_memory
        
        # Memory release check (informational only due to Python GC behavior and system variations)
        memory_release_ratio = memory_retained / memory_increase
        print(f"  Memory release ratio: {memory_release_ratio:.2f} (lower is better)")
        
        # Informational check - memory behavior varies significantly across systems
        if memory_retained > memory_increase * 0.9:
            print(f"  Warning: Memory usage grew beyond peak (ratio: {memory_release_ratio:.2f})")
        elif memory_retained < memory_increase * 0.5:
            print(f"  Good: Significant memory released (ratio: {memory_release_ratio:.2f})")
        else:
            print(f"  Normal: Moderate memory retention (ratio: {memory_release_ratio:.2f})")
        
        print(f"  Final memory: {final_memory:.1f}MB")
        print(f"  Memory retained: {memory_retained:.1f}MB")
    
    def test_large_structure_scalability_comparison(self, large_pdb_file, sample_pdb_file):
        """Test scalability by comparing large vs small structure performance."""
        # Analyze small structure
        small_analyzer = MolecularInteractionAnalyzer()
        small_start = time.time()
        small_success = small_analyzer.analyze_file(sample_pdb_file)
        small_time = time.time() - small_start
        small_stats = small_analyzer.get_summary()
        
        assert small_success
        
        # Analyze large structure
        large_analyzer = MolecularInteractionAnalyzer()
        large_start = time.time()
        large_success = large_analyzer.analyze_file(large_pdb_file)
        large_time = time.time() - large_start
        large_stats = large_analyzer.get_summary()
        
        assert large_success
        
        # Calculate scalability metrics
        size_ratio = large_stats['total_interactions'] / max(small_stats['total_interactions'], 1)
        time_ratio = large_time / max(small_time, 0.01)
        
        # Performance should scale reasonably (not exponentially)
        efficiency_ratio = time_ratio / size_ratio
        assert efficiency_ratio < 15.0, f"Performance scaling poor: {efficiency_ratio:.2f}x"
        
        print(f"\nScalability Analysis:")
        print(f"  Small structure (6rsa.pdb):")
        print(f"    Time: {small_time:.2f}s")
        print(f"    Interactions: {small_stats['total_interactions']}")
        print(f"  Large structure (4jsv.pdb):")
        print(f"    Time: {large_time:.2f}s")
        print(f"    Interactions: {large_stats['total_interactions']}")
        print(f"  Scaling factors:")
        print(f"    Size ratio: {size_ratio:.1f}x")
        print(f"    Time ratio: {time_ratio:.1f}x")
        print(f"    Efficiency ratio: {efficiency_ratio:.2f}x")


@pytest.mark.performance
@pytest.mark.slow
@pytest.mark.requires_pdb_files
@pytest.mark.integration
class TestLargeStructureParameterPerformance:
    """Test performance impact of different parameters on large structures."""
    
    @pytest.fixture
    def large_pdb_file(self, test_pdb_dir):
        """Provide path to large PDB file (4jsv.pdb)."""
        large_file_path = os.path.join(test_pdb_dir, "4jsv.pdb")
        if not os.path.exists(large_file_path):
            pytest.skip("Large PDB file (4jsv.pdb) not found")
        return large_file_path
    
    def test_analysis_mode_performance_impact(self, large_pdb_file):
        """Test performance impact of analysis modes on large structures."""
        # Test complete mode
        complete_params = AnalysisParameters(analysis_mode="complete")
        complete_analyzer = MolecularInteractionAnalyzer(complete_params)
        
        complete_start = time.time()
        complete_success = complete_analyzer.analyze_file(large_pdb_file)
        complete_time = time.time() - complete_start
        complete_stats = complete_analyzer.get_summary()
        
        assert complete_success
        
        # Test local mode
        local_params = AnalysisParameters(analysis_mode="local")
        local_analyzer = MolecularInteractionAnalyzer(local_params)
        
        local_start = time.time()
        local_success = local_analyzer.analyze_file(large_pdb_file)
        local_time = time.time() - local_start
        local_stats = local_analyzer.get_summary()
        
        assert local_success
        
        # Local mode should generally be faster
        time_speedup = complete_time / local_time
        
        print(f"\nAnalysis Mode Performance (4jsv.pdb):")
        print(f"  Complete mode:")
        print(f"    Time: {complete_time:.2f}s")
        print(f"    Interactions: {complete_stats['total_interactions']}")
        print(f"  Local mode:")
        print(f"    Time: {local_time:.2f}s")
        print(f"    Interactions: {local_stats['total_interactions']}")
        print(f"  Speedup: {time_speedup:.2f}x")
        
        # Both should complete in reasonable time
        assert complete_time < 300.0, "Complete mode should finish within 5 minutes"
        assert local_time < 200.0, "Local mode should finish within 3.5 minutes"
    
    def test_parameter_strictness_performance_impact(self, large_pdb_file):
        """Test performance impact of parameter strictness."""
        # Strict parameters (fewer interactions to compute)
        strict_params = AnalysisParameters(
            hb_distance_cutoff=2.5,
            hb_angle_cutoff=150.0,
            analysis_mode="local"
        )
        strict_analyzer = MolecularInteractionAnalyzer(strict_params)
        
        strict_start = time.time()
        strict_success = strict_analyzer.analyze_file(large_pdb_file)
        strict_time = time.time() - strict_start
        strict_stats = strict_analyzer.get_summary()
        
        assert strict_success
        
        # Permissive parameters (more interactions to compute)
        permissive_params = AnalysisParameters(
            hb_distance_cutoff=4.0,
            hb_angle_cutoff=100.0,
            analysis_mode="complete"
        )
        permissive_analyzer = MolecularInteractionAnalyzer(permissive_params)
        
        permissive_start = time.time()
        permissive_success = permissive_analyzer.analyze_file(large_pdb_file)
        permissive_time = time.time() - permissive_start
        permissive_stats = permissive_analyzer.get_summary()
        
        assert permissive_success
        
        # Verify expected behavior
        assert permissive_stats['hydrogen_bonds']['count'] >= strict_stats['hydrogen_bonds']['count'], \
            "Permissive parameters should find more interactions"
        
        print(f"\nParameter Strictness Performance (4jsv.pdb):")
        print(f"  Strict parameters:")
        print(f"    Time: {strict_time:.2f}s")
        print(f"    H-bonds: {strict_stats['hydrogen_bonds']['count']}")
        print(f"    Total: {strict_stats['total_interactions']}")
        print(f"  Permissive parameters:")
        print(f"    Time: {permissive_time:.2f}s")
        print(f"    H-bonds: {permissive_stats['hydrogen_bonds']['count']}")
        print(f"    Total: {permissive_stats['total_interactions']}")
        
        # Both should complete reasonably
        assert strict_time < 180.0, "Strict analysis should be fast"
        assert permissive_time < 400.0, "Permissive analysis should complete"


@pytest.mark.performance
@pytest.mark.slow
@pytest.mark.requires_pdb_files
@pytest.mark.integration
class TestLargeStructureOutputPerformance:
    """Test output generation performance with large structures."""
    
    @pytest.fixture
    def large_pdb_file(self, test_pdb_dir):
        """Provide path to large PDB file (4jsv.pdb)."""
        large_file_path = os.path.join(test_pdb_dir, "4jsv.pdb")
        if not os.path.exists(large_file_path):
            pytest.skip("Large PDB file (4jsv.pdb) not found")
        return large_file_path
    
    def test_large_structure_json_export_performance(self, large_pdb_file):
        """Test JSON export performance with large results."""
        analyzer = MolecularInteractionAnalyzer()
        
        # Run analysis
        analysis_start = time.time()
        success = analyzer.analyze_file(large_pdb_file)
        analysis_time = time.time() - analysis_start
        
        assert success
        
        # Test JSON export performance
        export_start = time.time()
        
        # Create comprehensive results for export
        results = {
            'metadata': {
                'file': large_pdb_file,
                'analysis_time': analysis_time,
                'timestamp': time.time()
            },
            'parameters': {
                'hb_distance_cutoff': analyzer.parameters.hb_distance_cutoff,
                'hb_angle_cutoff': analyzer.parameters.hb_angle_cutoff,
                'analysis_mode': analyzer.parameters.analysis_mode
            },
            'statistics': analyzer.get_summary(),
            'summary': analyzer.get_summary(),
            'interactions': {
                'hydrogen_bonds': [
                    {
                        'donor_residue': hb.donor_residue,
                        'acceptor_residue': hb.acceptor_residue,
                        'distance': hb.distance,
                        'angle': hb.angle,
                        'bond_type': hb.bond_type
                    }
                    for hb in analyzer.hydrogen_bonds[:1000]  # Limit for performance
                ],
                'pi_interactions': [
                    {
                        'donor_residue': pi.donor_residue,
                        'pi_residue': pi.acceptor_residue,
                        'distance': pi.distance,
                        'angle': pi.angle
                    }
                    for pi in analyzer.pi_interactions[:100]  # Limit for performance
                ]
            }
        }
        
        # Export to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(results, f, indent=2)
            export_path = f.name
        
        export_time = time.time() - export_start
        
        try:
            # Verify export file
            file_size = os.path.getsize(export_path) / 1024 / 1024  # MB
            
            # Export should be reasonably fast
            assert export_time < 30.0, f"JSON export too slow: {export_time:.2f}s"
            
            # File should be reasonable size
            assert file_size < 100.0, f"Export file too large: {file_size:.1f}MB"
            
            # Verify file can be read back
            read_start = time.time()
            with open(export_path, 'r') as f:
                loaded_results = json.load(f)
            read_time = time.time() - read_start
            
            assert read_time < 10.0, f"JSON read too slow: {read_time:.2f}s"
            assert loaded_results['statistics']['total_interactions'] > 0
            
            print(f"\nJSON Export Performance (4jsv.pdb):")
            print(f"  Export time: {export_time:.2f}s")
            print(f"  File size: {file_size:.1f}MB")
            print(f"  Read time: {read_time:.2f}s")
            print(f"  H-bonds exported: {len(results['interactions']['hydrogen_bonds'])}")
            
        finally:
            os.unlink(export_path)
    
    def test_large_structure_summary_generation_performance(self, large_pdb_file):
        """Test summary generation performance."""
        analyzer = MolecularInteractionAnalyzer()
        
        # Run analysis
        success = analyzer.analyze_file(large_pdb_file)
        assert success
        
        # Test summary generation performance
        summary_times = []
        stats_times = []
        
        # Run multiple iterations to get stable timing
        for _ in range(5):
            # Test summary generation
            summary_start = time.time()
            summary = analyzer.get_summary()
            summary_times.append(time.time() - summary_start)
            
            # Test statistics generation
            stats_start = time.time()
            stats = analyzer.get_summary()
            stats_times.append(time.time() - stats_start)
            
            # Verify consistency
            assert summary['total_interactions'] == stats['total_interactions']
        
        avg_summary_time = sum(summary_times) / len(summary_times)
        avg_stats_time = sum(stats_times) / len(stats_times)
        
        # Summary generation should be very fast
        assert avg_summary_time < 0.1, f"Summary generation too slow: {avg_summary_time:.3f}s"
        assert avg_stats_time < 0.1, f"Statistics generation too slow: {avg_stats_time:.3f}s"
        
        print(f"\nSummary Generation Performance (4jsv.pdb):")
        print(f"  Average summary time: {avg_summary_time*1000:.1f}ms")
        print(f"  Average stats time: {avg_stats_time*1000:.1f}ms")
        print(f"  Total interactions: {stats['total_interactions']}")


@pytest.mark.e2e
@pytest.mark.performance
@pytest.mark.slow  
@pytest.mark.requires_pdb_files
class TestLargeStructureStressTest:
    """Stress tests using large structures to test system limits."""
    
    @pytest.fixture
    def large_pdb_file(self, test_pdb_dir):
        """Provide path to large PDB file (4jsv.pdb)."""
        large_file_path = os.path.join(test_pdb_dir, "4jsv.pdb")
        if not os.path.exists(large_file_path):
            pytest.skip("Large PDB file (4jsv.pdb) not found")
        return large_file_path
    
    def test_multiple_large_structure_analyses(self, large_pdb_file):
        """Test multiple consecutive analyses of large structure."""
        num_runs = 2
        times = []
        
        for run in range(num_runs):
            analyzer = MolecularInteractionAnalyzer()
            
            start_time = time.time()
            success = analyzer.analyze_file(large_pdb_file)
            run_time = time.time() - start_time
            
            assert success, f"Run {run+1} should succeed"
            times.append(run_time)
            
            # Verify consistent results
            stats = analyzer.get_summary()
            assert stats['total_interactions'] > 100, f"Run {run+1} should find interactions"
            
            print(f"Run {run+1}: {run_time:.2f}s, {stats['total_interactions']} interactions")
            
            # Clean up
            del analyzer
        
        # Performance should be consistent
        avg_time = sum(times) / len(times)
        max_time = max(times)
        min_time = min(times)
        time_variance = max_time - min_time
        
        # Variance should be reasonable (< 100% of average for large structures)
        assert time_variance < avg_time * 1.0, f"Performance too variable: {time_variance:.2f}s"
        
        print(f"\nMultiple Run Performance (4jsv.pdb):")
        print(f"  Average time: {avg_time:.2f}s")
        print(f"  Min time: {min_time:.2f}s")
        print(f"  Max time: {max_time:.2f}s")
        print(f"  Variance: {time_variance:.2f}s")
    
    def test_large_structure_with_pdb_fixing_performance(self, test_pdb_dir):
        """Test PDB fixing performance with large structure."""
        # Use 1ubi.pdb for fixing test (larger than 6rsa but fixable)
        fixing_file = os.path.join(test_pdb_dir, "1ubi.pdb")
        if not os.path.exists(fixing_file):
            pytest.skip("PDB fixing test file (1ubi.pdb) not found")
        
        # Test with PDB fixing enabled
        fixing_params = AnalysisParameters(
            fix_pdb_enabled=True,
            fix_pdb_method="openbabel",
            fix_pdb_add_hydrogens=True
        )
        fixing_analyzer = MolecularInteractionAnalyzer(fixing_params)
        
        fixing_start = time.time()
        fixing_success = fixing_analyzer.analyze_file(fixing_file)
        fixing_time = time.time() - fixing_start
        
        assert fixing_success, "PDB fixing analysis should succeed"
        
        # Test without fixing for comparison
        no_fixing_params = AnalysisParameters(fix_pdb_enabled=False)
        no_fixing_analyzer = MolecularInteractionAnalyzer(no_fixing_params)
        
        no_fixing_start = time.time()
        no_fixing_success = no_fixing_analyzer.analyze_file(fixing_file)
        no_fixing_time = time.time() - no_fixing_start
        
        assert no_fixing_success, "No fixing analysis should succeed"
        
        # PDB fixing adds overhead but should still be reasonable
        fixing_overhead = fixing_time - no_fixing_time
        
        print(f"\nPDB Fixing Performance (1ubi.pdb):")
        print(f"  Without fixing: {no_fixing_time:.2f}s")
        print(f"  With fixing: {fixing_time:.2f}s")
        print(f"  Overhead: {fixing_overhead:.2f}s")
        
        # Both should complete in reasonable time
        assert no_fixing_time < 60.0, "Analysis without fixing should be fast"
        assert fixing_time < 180.0, "Analysis with fixing should complete"
        assert fixing_overhead < 120.0, "PDB fixing overhead should be reasonable"


@pytest.mark.e2e
@pytest.mark.performance
@pytest.mark.slow
class TestPerformanceBenchmarks:
    """Performance benchmarks and regression testing."""
    
    def test_performance_baseline_documentation(self, test_pdb_dir):
        """Document performance baselines for different structure sizes."""
        test_files = [
            ("6rsa.pdb", "Small structure"),
            ("4jsv.pdb", "Large structure"),
            ("1ubi.pdb", "Medium structure")
        ]
        
        results = {}
        
        for filename, description in test_files:
            file_path = os.path.join(test_pdb_dir, filename)
            if not os.path.exists(file_path):
                print(f"Skipping {filename} - file not found")
                continue
            
            analyzer = MolecularInteractionAnalyzer()
            
            # Measure performance
            start_time = time.time()
            success = analyzer.analyze_file(file_path)
            analysis_time = time.time() - start_time
            
            if success:
                stats = analyzer.get_summary()
                
                results[filename] = {
                    'description': description,
                    'time': analysis_time,
                    'interactions': stats['total_interactions'],
                    'hydrogen_bonds': stats['hydrogen_bonds']['count'],
                    'rate': stats['total_interactions'] / analysis_time
                }
        
        # Document results
        print(f"\nPerformance Baseline Documentation:")
        print("=" * 80)
        print(f"{'File':<15} {'Description':<15} {'Time(s)':<10} {'H-bonds':<10} {'Total':<10} {'Rate':<10}")
        print("-" * 80)
        
        for filename, data in results.items():
            print(f"{filename:<15} {data['description']:<15} {data['time']:<10.2f} "
                  f"{data['hydrogen_bonds']:<10} {data['interactions']:<10} {data['rate']:<10.1f}")
        
        print("=" * 80)
        
        # Basic performance assertions
        for filename, data in results.items():
            if "Small" in data['description']:
                assert data['time'] < 30.0, f"{filename} should analyze quickly"
            elif "Large" in data['description']:
                assert data['time'] < 300.0, f"{filename} should complete within 5 minutes"
                assert data['interactions'] > 100, f"{filename} should find many interactions"
            
            # All should have reasonable interaction rates
            assert data['rate'] > 1.0, f"{filename} should process at least 1 interaction/second"