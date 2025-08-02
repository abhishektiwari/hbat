"""
Unit tests for GraphViz detection utilities.

Tests the GraphVizDetector class and related functionality for detecting
GraphViz installation, version information, and available engines.
"""

import os
import subprocess
import sys
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pytest

from hbat.utilities.graphviz_utils import (
    GRAPHVIZ_EXECUTABLES,
    PLATFORM_PATHS,
    GraphVizDetector,
    get_graphviz_info,
)


class TestGraphVizDetector:
    """Test cases for GraphVizDetector class."""

    def setup_method(self):
        """Set up test fixtures before each test method."""
        # Clear cache before each test
        GraphVizDetector.clear_cache()

    def teardown_method(self):
        """Clean up after each test method."""
        # Clear cache after each test
        GraphVizDetector.clear_cache()

    def test_is_graphviz_available_found_in_path(self):
        """Test GraphViz detection when dot is found in PATH."""
        with patch.object(GraphVizDetector, '_which', return_value='/usr/bin/dot'):
            result = GraphVizDetector.is_graphviz_available()
            assert result is True

    def test_is_graphviz_available_not_found(self):
        """Test GraphViz detection when dot is not found anywhere."""
        with patch.object(GraphVizDetector, '_which', return_value=None), \
             patch.object(GraphVizDetector, '_check_executable_in_directory', return_value=False):
            result = GraphVizDetector.is_graphviz_available()
            assert result is False

    def test_is_graphviz_available_found_in_platform_path(self):
        """Test GraphViz detection when dot is found in platform-specific path."""
        with patch.object(GraphVizDetector, '_which', return_value=None), \
             patch.object(GraphVizDetector, '_check_executable_in_directory', return_value=True), \
             patch('sys.platform', 'darwin'):
            result = GraphVizDetector.is_graphviz_available()
            assert result is True

    def test_is_graphviz_available_caching(self):
        """Test that detection results are cached."""
        with patch.object(GraphVizDetector, '_which', return_value='/usr/bin/dot') as mock_which:
            # First call
            result1 = GraphVizDetector.is_graphviz_available()
            assert result1 is True
            
            # Second call should use cache
            result2 = GraphVizDetector.is_graphviz_available()
            assert result2 is True
            
            # _which should only be called once due to caching
            assert mock_which.call_count == 1

    def test_get_graphviz_version_success(self):
        """Test successful GraphViz version detection."""
        mock_result = Mock()
        mock_result.returncode = 0
        mock_result.stdout = ""
        mock_result.stderr = "dot - graphviz version 2.44.1 (20200629.0800)"
        
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=True), \
             patch('subprocess.run', return_value=mock_result):
            version = GraphVizDetector.get_graphviz_version()
            assert version == "2.44.1"

    def test_get_graphviz_version_not_available(self):
        """Test version detection when GraphViz is not available."""
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=False):
            version = GraphVizDetector.get_graphviz_version()
            assert version is None

    def test_get_graphviz_version_subprocess_error(self):
        """Test version detection when subprocess fails."""
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=True), \
             patch('subprocess.run', side_effect=subprocess.SubprocessError("Command failed")):
            version = GraphVizDetector.get_graphviz_version()
            assert version is None

    def test_get_graphviz_version_timeout(self):
        """Test version detection when subprocess times out."""
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=True), \
             patch('subprocess.run', side_effect=subprocess.TimeoutExpired("dot", 5)):
            version = GraphVizDetector.get_graphviz_version()
            assert version is None

    def test_get_graphviz_version_caching(self):
        """Test that version results are cached."""
        mock_result = Mock()
        mock_result.returncode = 0
        mock_result.stdout = ""
        mock_result.stderr = "dot - graphviz version 2.44.1 (20200629.0800)"
        
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=True), \
             patch('subprocess.run', return_value=mock_result) as mock_subprocess:
            # First call
            version1 = GraphVizDetector.get_graphviz_version()
            assert version1 == "2.44.1"
            
            # Second call should use cache
            version2 = GraphVizDetector.get_graphviz_version()
            assert version2 == "2.44.1"
            
            # subprocess.run should only be called once due to caching
            assert mock_subprocess.call_count == 1

    def test_get_available_engines_success(self):
        """Test successful engine detection."""
        def mock_check_path(engine):
            return engine in ['dot', 'neato', 'fdp']
        
        def mock_check_directory(engine, path):
            return False  # Don't find anything in platform-specific paths
        
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=True), \
             patch.object(GraphVizDetector, '_check_executable_in_path', side_effect=mock_check_path), \
             patch.object(GraphVizDetector, '_check_executable_in_directory', side_effect=mock_check_directory):
            engines = GraphVizDetector.get_available_engines()
            assert 'dot' in engines
            assert 'neato' in engines
            assert 'fdp' in engines
            assert 'circo' not in engines

    def test_get_available_engines_not_available(self):
        """Test engine detection when GraphViz is not available."""
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=False):
            engines = GraphVizDetector.get_available_engines()
            assert engines == []

    def test_get_available_engines_platform_paths(self):
        """Test engine detection in platform-specific paths."""
        def mock_check_path(engine):
            return False  # Not in PATH
        
        def mock_check_directory(engine, directory):
            return engine == 'dot' and '/usr/local/bin' in directory
        
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=True), \
             patch.object(GraphVizDetector, '_check_executable_in_path', side_effect=mock_check_path), \
             patch.object(GraphVizDetector, '_check_executable_in_directory', side_effect=mock_check_directory), \
             patch('sys.platform', 'darwin'):
            engines = GraphVizDetector.get_available_engines()
            assert 'dot' in engines

    def test_get_available_engines_caching(self):
        """Test that engine detection results are cached."""
        def mock_check_path(engine):
            return engine == 'dot'
        
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=True), \
             patch.object(GraphVizDetector, '_check_executable_in_path', side_effect=mock_check_path) as mock_check:
            # First call
            engines1 = GraphVizDetector.get_available_engines()
            assert 'dot' in engines1
            
            # Second call should use cache
            engines2 = GraphVizDetector.get_available_engines()
            assert engines1 == engines2
            
            # Should only check once due to caching
            assert mock_check.call_count == len(GRAPHVIZ_EXECUTABLES)

    def test_validate_engine_valid(self):
        """Test engine validation for valid engine."""
        with patch.object(GraphVizDetector, 'get_available_engines', return_value=['dot', 'neato']):
            assert GraphVizDetector.validate_engine('dot') is True
            assert GraphVizDetector.validate_engine('neato') is True

    def test_validate_engine_invalid(self):
        """Test engine validation for invalid engine."""
        with patch.object(GraphVizDetector, 'get_available_engines', return_value=['dot', 'neato']):
            assert GraphVizDetector.validate_engine('invalid_engine') is False

    def test_get_engine_path_in_path(self):
        """Test getting engine path when in PATH."""
        with patch.object(GraphVizDetector, '_which', return_value='/usr/bin/dot'):
            path = GraphVizDetector.get_engine_path('dot')
            assert path == '/usr/bin/dot'

    def test_get_engine_path_not_found(self):
        """Test getting engine path when not found."""
        with patch.object(GraphVizDetector, '_which', return_value=None), \
             patch('pathlib.Path.exists', return_value=False):
            path = GraphVizDetector.get_engine_path('dot')
            assert path is None

    @pytest.mark.skip("Complex mocking - covered by integration tests")
    def test_get_engine_path_platform_specific(self):
        """Test getting engine path from platform-specific location."""
        # This test is covered by integration tests and real system behavior
        pass

    def test_clear_cache(self):
        """Test cache clearing functionality."""
        # Test that clear_cache runs without error and clears all cache variables
        GraphVizDetector.clear_cache()
        
        # Verify cache variables are cleared
        assert GraphVizDetector._detection_cache is None
        assert GraphVizDetector._version_cache is None  
        assert GraphVizDetector._engines_cache is None

    def test_check_executable_in_path_found(self):
        """Test private method _check_executable_in_path when found."""
        with patch.object(GraphVizDetector, '_which', return_value='/usr/bin/dot'):
            result = GraphVizDetector._check_executable_in_path('dot')
            assert result is True

    def test_check_executable_in_path_not_found(self):
        """Test private method _check_executable_in_path when not found."""
        with patch.object(GraphVizDetector, '_which', return_value=None):
            result = GraphVizDetector._check_executable_in_path('dot')
            assert result is False

    def test_check_executable_in_directory_found(self):
        """Test private method _check_executable_in_directory when found."""
        def mock_path_constructor(*args, **kwargs):
            mock_path = Mock()
            mock_path.exists.return_value = True
            mock_path.is_file.return_value = True
            mock_path.__truediv__ = lambda self, other: mock_path
            return mock_path
        
        with patch('hbat.utilities.graphviz_utils.Path', side_effect=mock_path_constructor):
            result = GraphVizDetector._check_executable_in_directory('dot', '/usr/bin')
            assert result is True

    def test_check_executable_in_directory_not_found(self):
        """Test private method _check_executable_in_directory when not found."""
        mock_path_instance = Mock()
        mock_path_instance.exists.return_value = False
        mock_path_instance.is_file.return_value = False
        
        mock_path_constructor = Mock()
        mock_path_constructor.__truediv__ = Mock(return_value=mock_path_instance)
        
        with patch('hbat.utilities.graphviz_utils.Path', return_value=mock_path_constructor):
            result = GraphVizDetector._check_executable_in_directory('dot', '/usr/bin')
            assert result is False

    def test_check_executable_in_directory_windows(self):
        """Test private method _check_executable_in_directory on Windows."""
        mock_path = Mock()
        mock_path.with_suffix.return_value = mock_path
        mock_path.exists.return_value = True
        mock_path.is_file.return_value = True
        mock_path.__truediv__ = lambda self, other: mock_path
        
        def mock_path_constructor(*args, **kwargs):
            return mock_path
        
        with patch('hbat.utilities.graphviz_utils.Path', side_effect=mock_path_constructor), \
             patch('sys.platform', 'win32'):
            result = GraphVizDetector._check_executable_in_directory('dot', r'C:\Program Files\Graphviz\bin')
            assert result is True
            mock_path.with_suffix.assert_called_with('.exe')

    def test_which_found(self):
        """Test private method _which when executable is found."""
        mock_path = Mock()
        mock_path.exists.return_value = True
        mock_path.is_file.return_value = True
        mock_path.__str__ = Mock(return_value='/usr/bin/dot')
        mock_path.__truediv__ = lambda self, other: mock_path
        
        def mock_path_constructor(path):
            return mock_path
        
        with patch('hbat.utilities.graphviz_utils.Path', side_effect=mock_path_constructor), \
             patch.dict(os.environ, {'PATH': '/usr/bin:/usr/local/bin'}):
            result = GraphVizDetector._which('dot')
            assert str(result) == '/usr/bin/dot'

    def test_which_not_found(self):
        """Test private method _which when executable is not found."""
        mock_path_instance = Mock()
        mock_path_instance.exists.return_value = False
        mock_path_instance.is_file.return_value = False
        
        mock_path_constructor = Mock()
        mock_path_constructor.__truediv__ = Mock(return_value=mock_path_instance)
        
        with patch('hbat.utilities.graphviz_utils.Path', return_value=mock_path_constructor), \
             patch.dict(os.environ, {'PATH': '/usr/bin:/usr/local/bin'}):
            result = GraphVizDetector._which('nonexistent')
            assert result is None

    def test_which_windows(self):
        """Test private method _which on Windows platform."""
        mock_path = Mock()
        mock_path.exists.return_value = True
        mock_path.is_file.return_value = True
        mock_path.__str__ = Mock(return_value=r'C:\Program Files\Graphviz\bin\dot.exe')
        mock_path.__truediv__ = lambda self, other: mock_path
        
        def mock_path_constructor(path):
            return mock_path
        
        with patch('hbat.utilities.graphviz_utils.Path', side_effect=mock_path_constructor), \
             patch('sys.platform', 'win32'), \
             patch.dict(os.environ, {'PATH': r'C:\Program Files\Graphviz\bin'}):
            result = GraphVizDetector._which('dot')
            assert '.exe' in str(result)


class TestGetGraphVizInfo:
    """Test cases for get_graphviz_info function."""

    def test_get_graphviz_info_complete(self):
        """Test get_graphviz_info with complete information."""
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=True), \
             patch.object(GraphVizDetector, 'get_graphviz_version', return_value='2.44.1'), \
             patch.object(GraphVizDetector, 'get_available_engines', return_value=['dot', 'neato']):
            info = get_graphviz_info()
            
            assert info['available'] is True
            assert info['version'] == '2.44.1'
            assert info['engines'] == ['dot', 'neato']
            assert 'platform' in info

    def test_get_graphviz_info_not_available(self):
        """Test get_graphviz_info when GraphViz is not available."""
        with patch.object(GraphVizDetector, 'is_graphviz_available', return_value=False), \
             patch.object(GraphVizDetector, 'get_graphviz_version', return_value=None), \
             patch.object(GraphVizDetector, 'get_available_engines', return_value=[]):
            info = get_graphviz_info()
            
            assert info['available'] is False
            assert info['version'] is None
            assert info['engines'] == []
            assert 'platform' in info


class TestConstants:
    """Test cases for module constants."""

    def test_graphviz_executables_defined(self):
        """Test that GRAPHVIZ_EXECUTABLES is properly defined."""
        assert isinstance(GRAPHVIZ_EXECUTABLES, list)
        assert len(GRAPHVIZ_EXECUTABLES) > 0
        assert 'dot' in GRAPHVIZ_EXECUTABLES
        assert 'neato' in GRAPHVIZ_EXECUTABLES

    def test_platform_paths_defined(self):
        """Test that PLATFORM_PATHS is properly defined."""
        assert isinstance(PLATFORM_PATHS, dict)
        assert 'win32' in PLATFORM_PATHS
        assert 'darwin' in PLATFORM_PATHS
        assert 'linux' in PLATFORM_PATHS
        
        # Each platform should have a list of paths
        for platform, paths in PLATFORM_PATHS.items():
            assert isinstance(paths, list)
            assert len(paths) > 0