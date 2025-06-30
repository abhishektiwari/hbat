# Development Guidelines

## Core Development Rules

1. Testing Requirements

    - When updating existing functions or classes must update the existing test cases
    - When adding a new feature must create a new test case

2.  Documentation Requirements

    - When changing in `core` module must update `logic.rst`
    - When adding new file to `hbat` or `tests` must update `development.rst`

3. Code Quality

   - Type hints required for all code in `core` module
   - All public APIs in `core`, `cli`, and `gui` module must have Sphinx docstring format
   - Run formatters before type checks
   - Format code using `make format`
   - Perform type check using `make type-check`
   - Fix type checking errors by adding type hint
   - For typing support install required package type stubs if available

4. Requirements management

    - After installing new packages add them to requirements files
