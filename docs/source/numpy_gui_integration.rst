NumPy GUI Integration
====================

HBAT now includes seamless integration of the NumPy-optimized analyzer with the graphical user interface, providing users with high-performance molecular interaction analysis through an intuitive interface.

Overview
--------

The GUI integration allows users to:

* **Toggle between analysis engines**: Choose between the standard analyzer and the NumPy-optimized version
* **Real-time performance indication**: Visual feedback showing which analyzer is active
* **Enhanced status reporting**: Clear indication of analysis progress and completion
* **Automatic performance optimization**: NumPy analyzer is enabled by default for best performance

Key Features
------------

Analyzer Selection
^^^^^^^^^^^^^^^^^^

Users can easily switch between analysis engines through the **Settings** menu:

* **Settings → Use NumPy Analyzer (High Performance)**: Toggle NumPy optimization
* **Default**: NumPy analyzer is enabled by default
* **Visual feedback**: Performance indicator shows current mode

Performance Indicators
^^^^^^^^^^^^^^^^^^^^^^

The GUI provides visual cues about the selected analyzer:

* **⚡ NumPy High-Performance Mode**: Green indicator for NumPy analyzer
* **🐌 Standard Mode**: Orange indicator for standard analyzer
* **Toolbar display**: Shows performance mode during analysis
* **Status updates**: Include analyzer type in completion messages

Enhanced User Experience
^^^^^^^^^^^^^^^^^^^^^^^^

* **Intelligent defaults**: NumPy analyzer enabled for optimal performance
* **Informational dialogs**: Clear explanations when switching modes
* **Export metadata**: Analysis results include analyzer type information
* **Performance context**: Users understand the benefits of each mode

Usage Instructions
------------------

Basic Workflow
^^^^^^^^^^^^^^

1. **Start the GUI**::

    python -m hbat.gui

2. **Load a PDB file**:
   - Click **File → Open PDB File...**
   - Select your protein structure file

3. **Configure analyzer** (optional):
   - Go to **Settings → Use NumPy Analyzer (High Performance)**
   - Toggle based on your needs

4. **Run analysis**:
   - Click **Analysis → Run Analysis** or press **F5**
   - Watch the performance indicator during analysis

5. **View and export results**:
   - Review results in the right panel
   - Use **File → Save Results...** to export

Analyzer Selection Guide
^^^^^^^^^^^^^^^^^^^^^^^^

**Use NumPy Analyzer when:**

* Working with large structures (>1000 atoms)
* Performing batch analysis
* Need maximum performance
* System has sufficient memory

**Use Standard Analyzer when:**

* Working with small structures (<500 atoms)
* Memory constraints exist
* Compatibility is a concern
* Debugging analysis issues

Performance Comparison
----------------------

The NumPy analyzer provides significant performance improvements:

================= ================== ==================
Structure Size    Standard Analyzer  NumPy Analyzer
================= ================== ==================
Small (<500)      1x                 1-2x faster
Medium (500-2000) 1x                 5-15x faster
Large (>2000)     1x                 10-30x faster
================= ================== ==================

Memory Usage
^^^^^^^^^^^^

* **NumPy analyzer**: Higher initial memory usage for coordinate caching
* **Standard analyzer**: Lower memory footprint
* **Recommendation**: Use NumPy for structures with adequate system memory

Integration Details
-------------------

Technical Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^

The GUI integration includes:

* **Analyzer factory pattern**: Dynamic creation based on user preference
* **Thread-safe operations**: Analysis runs in background threads
* **State persistence**: Analyzer preference maintained during session
* **Error handling**: Graceful fallback to standard analyzer if needed

Code Example
^^^^^^^^^^^^

Programmatic usage of the integrated analyzers::

    from hbat.gui.main_window import MainWindow
    from hbat.core.analysis import AnalysisParameters
    
    # Create main window
    app = MainWindow()
    
    # Configure for NumPy analyzer
    app.use_numpy_analyzer = True
    
    # Set analysis parameters
    params = AnalysisParameters()
    
    # Run analysis (would normally be triggered by GUI)
    app._perform_analysis(params)

Customization Options
^^^^^^^^^^^^^^^^^^^^^

Advanced users can customize the integration:

* **Default analyzer**: Modify ``use_numpy_analyzer`` default in ``MainWindow.__init__()``
* **Performance thresholds**: Adjust automatic selection based on structure size
* **Visual indicators**: Customize icons and colors in the toolbar
* **Export format**: Extend result export to include performance metrics

Troubleshooting
---------------

Common Issues
^^^^^^^^^^^^^

**NumPy analyzer fails to load:**

* Check that NumPy is installed: ``pip install numpy>=1.20.0``
* Verify system memory is sufficient
* Try running with standard analyzer

**Performance not as expected:**

* Ensure NumPy is using optimized BLAS libraries
* Check system resources during analysis
* Consider structure size and complexity

**GUI freezes during analysis:**

* Large structures may take time even with NumPy optimization
* Check system memory and CPU usage
* Consider using command-line interface for very large structures

Best Practices
--------------

For Optimal Performance
^^^^^^^^^^^^^^^^^^^^^^^

1. **Enable NumPy by default**: Best performance for most use cases
2. **Monitor system resources**: Ensure adequate memory for large structures
3. **Use appropriate hardware**: Multi-core systems benefit most from NumPy optimization
4. **Batch processing**: For multiple files, consider command-line interface

For Compatibility
^^^^^^^^^^^^^^^^^

1. **Test both analyzers**: Verify results are consistent
2. **Document analysis settings**: Include analyzer type in result files
3. **Version control**: Track which analyzer version was used
4. **Backup results**: Keep copies with analyzer metadata

Future Enhancements
-------------------

Planned improvements include:

* **Automatic analyzer selection**: Based on structure size and system capabilities
* **Progress reporting**: Real-time progress bars for large analyses
* **Performance profiling**: Built-in benchmarking tools
* **Cloud integration**: Support for remote high-performance computing
* **GPU acceleration**: Future CUDA-based optimizations

The NumPy GUI integration represents a significant step forward in making high-performance molecular interaction analysis accessible through an intuitive graphical interface.