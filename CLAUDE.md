# Development Guidelines

## Core Development Rules

1. Testing Requirements

    - When updating existing functions or classes must update the existing test cases
    - When adding a new feature must create a new test case

2.  Documentation Requirements

    - Every time you change modules make sure corresponding.rst files in docs/source/modules are updated.

3. Code Quality

   - Type hints required for all code in `core` module
   - All public APIs in `core`, `cli`, and `gui` module must have Sphinx docstring format
   - Run formatters before type checks
   - Format code using `make format`
   - Perform type check using `make type-check`
   - Fix type checking errors by adding type hint
   - For typing support install required package type stubs if available
   - When refactoring modules/packages, make sure tests are refactored and matched with similar structure

4. Requirements management

    - After installing new packages add them to requirements files

5. For testing use following pdb files in example_pdb_files folder
    - 6RSA.pdb and 2IZF.pdb for general hydrogen bonds. 6RSA and 2IZF come with Hydrogen atoms so no PDB fixing required.
    - 4X21.pdb, 4LAZ.pdb and 4UB7.pdb for halogen bonds - they provide <5 Halogen bonds.
    - 1UBI.pdb and 1BHL.pdb lack hydrogen atoms. Add missing hydrogen atoms (i.e. PDB fixing) and then perform hydrogen bond analysis. Test both OpenBabel and PDB fixer for fixing missing hydrogen atoms.
