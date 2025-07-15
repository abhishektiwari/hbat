#!/usr/bin/env python3
"""
Test runner for HBAT test suite.

This script can run tests in multiple ways:
1. Import and run the main test_hbat module directly
2. Use pytest if available for more advanced test reporting
3. Run specific test categories

Usage:
    python run_tests.py              # Run all tests
    python run_tests.py --fast       # Skip slow tests
    python run_tests.py --pytest     # Force pytest usage
    python run_tests.py --legacy     # Use legacy test runner
"""

import sys
import os
import argparse

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def run_with_pytest(args):
    """Run tests using pytest."""
    try:
        import pytest
        
        # Determine test directory
        test_dir = os.path.dirname(os.path.abspath(__file__))
        
        pytest_args = ["-v"]
        
        # Add specific test categories if requested (new structure)
        if hasattr(args, 'unit') and args.unit:
            pytest_args.extend(["-m", "unit", os.path.join(test_dir, "unit")])
        elif hasattr(args, 'integration') and args.integration:
            pytest_args.extend(["-m", "integration", os.path.join(test_dir, "integration")])
        elif hasattr(args, 'e2e') and args.e2e:
            pytest_args.extend(["-m", "e2e", os.path.join(test_dir, "e2e")])
        elif hasattr(args, 'gui') and args.gui:
            pytest_args.append(os.path.join(test_dir, "gui"))
        elif hasattr(args, 'cli') and args.cli:
            pytest_args.append(os.path.join(test_dir, "cli"))
        else:
            # Run all tests from all directories
            pytest_args.extend([
                os.path.join(test_dir, "unit"),
                os.path.join(test_dir, "integration"),
                os.path.join(test_dir, "e2e"),
                os.path.join(test_dir, "performance"),
                os.path.join(test_dir, "cli")
            ])
            # Only add GUI tests if not skipping them
            if not (hasattr(args, 'no_gui') and args.no_gui):
                pytest_args.append(os.path.join(test_dir, "gui"))
        
        # Add markers for different test types
        if args.fast:
            pytest_args.extend(["-m", "unit and not slow"])
        elif hasattr(args, 'performance') and args.performance:
            pytest_args.extend(["-m", "performance", os.path.join(test_dir, "performance")])
        
        if args.verbose:
            pytest_args.append("-s")
            
        if args.coverage:
            pytest_args.extend(["--cov=hbat", "--cov-report=html", "--cov-report=term"])
        
        # Skip GUI tests if requested or if display not available
        if hasattr(args, 'no_gui') and args.no_gui:
            pytest_args.extend(["-m", "not gui"])
        
        # Skip tests requiring PDB files if requested
        if hasattr(args, 'no_pdb') and args.no_pdb:
            pytest_args.extend(["-m", "not requires_pdb_files"])
        
        print("Running tests with pytest...")
        print(f"Test command: pytest {' '.join(pytest_args)}")
        exit_code = pytest.main(pytest_args)
        return exit_code
        
    except ImportError:
        print("pytest not available, falling back to legacy test runner...")
        return run_legacy_tests(args)

def run_legacy_tests(args=None):
    """Run tests using the legacy test_hbat module."""
    try:
        # For now, legacy test runner runs all tests regardless of selection
        # since the legacy test doesn't support modular selection
        if args and (args.gui or args.cli):
            print(f"Note: Legacy test runner does not support modular test selection.")
            print(f"Running all available tests...")
        
        from tests.test_hbat import main as test_main
        return test_main()
    except ImportError as e:
        print(f"Error importing test module: {e}")
        return 1

def main():
    """Main test runner."""
    parser = argparse.ArgumentParser(
        description="HBAT Test Runner - Refactored test suite with modular organization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                    # Run all tests with pytest
  %(prog)s --unit             # Run only unit tests (fast)
  %(prog)s --integration      # Run only integration tests  
  %(prog)s --e2e              # Run only end-to-end tests
  %(prog)s --performance      # Run only performance tests
  %(prog)s --fast             # Run only fast unit tests
  %(prog)s --cli              # Run only CLI tests
  %(prog)s --gui              # Run only GUI tests
  %(prog)s --no-gui           # Skip GUI tests
  %(prog)s --no-pdb           # Skip tests requiring PDB files
  %(prog)s --coverage         # Generate coverage report
  %(prog)s --legacy           # Use legacy test runner
        """
    )
    
    # Test runner options
    parser.add_argument("--pytest", action="store_true", help="Force pytest usage")
    parser.add_argument("--legacy", action="store_true", help="Use legacy test runner")
    
    # Test selection options (new structure)
    parser.add_argument("--unit", action="store_true", help="Run only unit tests (fast, isolated)")
    parser.add_argument("--integration", action="store_true", help="Run only integration tests")
    parser.add_argument("--e2e", action="store_true", help="Run only end-to-end tests")
    parser.add_argument("--performance", action="store_true", help="Run only performance tests")
    
    parser.add_argument("--gui", action="store_true", help="Run only GUI tests")  
    parser.add_argument("--cli", action="store_true", help="Run only CLI tests")
    
    # Test filtering options
    parser.add_argument("--no-gui", action="store_true", help="Skip GUI tests")
    parser.add_argument("--no-pdb", action="store_true", help="Skip tests requiring PDB files")
    
    # Test execution options
    parser.add_argument("--fast", action="store_true", help="Skip slow and integration tests")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--coverage", action="store_true", help="Generate coverage report")
    
    args = parser.parse_args()
    
    print("HBAT Test Suite Runner - Refactored")
    print("=" * 50)
    print("Test organization:")
    print("  - tests/unit/         : Unit tests (fast, isolated components)")
    print("  - tests/integration/  : Integration tests (component interactions)")  
    print("  - tests/e2e/          : End-to-end tests (complete workflows)")
    print("  - tests/performance/  : Performance and benchmark tests")
    print("  - tests/cli/          : Command-line interface tests")
    print("  - tests/gui/          : Graphical interface tests")
    print("  - tests/conftest.py   : Shared fixtures and utilities")
    print()
    
    # Validate argument combinations
    test_selectors = [args.unit, args.integration, args.e2e, args.performance, 
                     args.gui, args.cli]
    if sum(test_selectors) > 1:
        print("Error: Cannot specify multiple test selection options")
        return 1
    
    if args.legacy:
        print("Using legacy test runner...")
        return run_legacy_tests(args)
    elif args.pytest:
        return run_with_pytest(args)
    else:
        # Try pytest first, fall back to legacy
        try:
            import pytest
            return run_with_pytest(args)
        except ImportError:
            print("pytest not available, using legacy test runner...")
            return run_legacy_tests(args)

if __name__ == "__main__":
    sys.exit(main())