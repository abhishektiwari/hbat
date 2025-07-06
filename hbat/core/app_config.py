"""
HBAT Application Configuration and Data Management

This module handles application-specific configuration, including creating and
managing the user's .hbat directory for storing CCD data, application state,
and other persistent data.
"""

import os
import json
from pathlib import Path
from typing import Dict, Any, Optional
from datetime import datetime


class HBATConfig:
    """
    Manages HBAT application configuration and user data directory.
    
    This class handles:
    - Creating and managing the user's ~/.hbat directory
    - Downloading and storing CCD data files
    - Managing application state and preferences
    - Providing paths for various data storage needs
    """
    
    def __init__(self):
        """Initialize HBAT configuration manager."""
        self.user_home = Path.home()
        self.hbat_dir = self.user_home / ".hbat"
        self.ccd_dir = self.hbat_dir / "ccd-data"
        self.config_file = self.hbat_dir / "config.json"
        self.state_file = self.hbat_dir / "app_state.json"
        self.logs_dir = self.hbat_dir / "logs"
        self.cache_dir = self.hbat_dir / "cache"
        self.presets_dir = self.hbat_dir / "presets"
        
        # Default configuration
        self.default_config = {
            "version": "1.0",
            "created": datetime.now().isoformat(),
            "last_updated": datetime.now().isoformat(),
            "ccd_data": {
                "auto_update": False,
                "last_download": None,
                "files_present": False
            },
            "preferences": {
                "theme": "default",
                "verbose_output": False,
                "auto_save_results": True,
                "default_output_format": "json"
            },
            "paths": {
                "last_pdb_directory": None,
                "last_output_directory": None
            }
        }
    
    def ensure_hbat_directory(self) -> bool:
        """
        Ensure the .hbat directory structure exists.
        
        Creates the main .hbat directory and all subdirectories if they don't exist.
        
        Returns:
            True if directory structure is ready, False if creation failed
        """
        try:
            # Create main .hbat directory
            self.hbat_dir.mkdir(exist_ok=True)
            print(f"📁 HBAT directory: {self.hbat_dir}")
            
            # Create subdirectories
            subdirs = [
                self.ccd_dir,
                self.logs_dir, 
                self.cache_dir,
                self.presets_dir
            ]
            
            for subdir in subdirs:
                subdir.mkdir(exist_ok=True)
                print(f"   📂 Created: {subdir.name}/")
            
            # Initialize configuration if it doesn't exist
            if not self.config_file.exists():
                self._create_initial_config()
                print(f"   ⚙️  Created initial configuration")
            
            # Initialize state file if it doesn't exist
            if not self.state_file.exists():
                self._create_initial_state()
                print(f"   📊 Created application state file")
            
            return True
            
        except Exception as e:
            print(f"❌ Error creating HBAT directory: {e}")
            return False
    
    def _create_initial_config(self):
        """Create initial configuration file."""
        with open(self.config_file, 'w') as f:
            json.dump(self.default_config, f, indent=2)
    
    def _create_initial_state(self):
        """Create initial application state file."""
        initial_state = {
            "first_run": datetime.now().isoformat(),
            "last_used": datetime.now().isoformat(),
            "run_count": 0,
            "recent_files": [],
            "recent_presets": [],
            "window_geometry": None,
            "last_analysis_parameters": None
        }
        
        with open(self.state_file, 'w') as f:
            json.dump(initial_state, f, indent=2)
    
    def load_config(self) -> Dict[str, Any]:
        """
        Load configuration from config file.
        
        Returns:
            Configuration dictionary
        """
        try:
            if self.config_file.exists():
                with open(self.config_file, 'r') as f:
                    return json.load(f)
            else:
                return self.default_config.copy()
        except Exception as e:
            print(f"⚠️  Error loading config, using defaults: {e}")
            return self.default_config.copy()
    
    def save_config(self, config: Dict[str, Any]) -> bool:
        """
        Save configuration to config file.
        
        Args:
            config: Configuration dictionary to save
            
        Returns:
            True if successful, False otherwise
        """
        try:
            config["last_updated"] = datetime.now().isoformat()
            with open(self.config_file, 'w') as f:
                json.dump(config, f, indent=2)
            return True
        except Exception as e:
            print(f"❌ Error saving config: {e}")
            return False
    
    def load_state(self) -> Dict[str, Any]:
        """
        Load application state from state file.
        
        Returns:
            Application state dictionary
        """
        try:
            if self.state_file.exists():
                with open(self.state_file, 'r') as f:
                    return json.load(f)
            else:
                return {}
        except Exception as e:
            print(f"⚠️  Error loading state: {e}")
            return {}
    
    def save_state(self, state: Dict[str, Any]) -> bool:
        """
        Save application state to state file.
        
        Args:
            state: Application state dictionary to save
            
        Returns:
            True if successful, False otherwise
        """
        try:
            state["last_used"] = datetime.now().isoformat()
            with open(self.state_file, 'w') as f:
                json.dump(state, f, indent=2)
            return True
        except Exception as e:
            print(f"❌ Error saving state: {e}")
            return False
    
    def update_run_count(self):
        """Update the application run count."""
        state = self.load_state()
        state["run_count"] = state.get("run_count", 0) + 1
        self.save_state(state)
    
    def add_recent_file(self, file_path: str, max_recent: int = 10):
        """
        Add a file to the recent files list.
        
        Args:
            file_path: Path to the file to add
            max_recent: Maximum number of recent files to keep
        """
        state = self.load_state()
        recent_files = state.get("recent_files", [])
        
        # Remove if already in list
        if file_path in recent_files:
            recent_files.remove(file_path)
        
        # Add to beginning
        recent_files.insert(0, file_path)
        
        # Trim to max_recent
        recent_files = recent_files[:max_recent]
        
        state["recent_files"] = recent_files
        self.save_state(state)
    
    def get_recent_files(self) -> list:
        """Get list of recent files."""
        state = self.load_state()
        return state.get("recent_files", [])
    
    def get_ccd_data_path(self) -> str:
        """
        Get the path for CCD data storage.
        
        Returns:
            Path to the CCD data directory
        """
        return str(self.ccd_dir)
    
    def is_first_run(self) -> bool:
        """
        Check if this is the first time HBAT is being run.
        
        Returns:
            True if this is the first run, False otherwise
        """
        return not self.config_file.exists()
    
    def update_ccd_status(self, files_present: bool, last_download: Optional[str] = None):
        """
        Update CCD data status in configuration.
        
        Args:
            files_present: Whether CCD files are present
            last_download: ISO timestamp of last download (optional)
        """
        config = self.load_config()
        config["ccd_data"]["files_present"] = files_present
        if last_download:
            config["ccd_data"]["last_download"] = last_download
        self.save_config(config)
    
    def get_preference(self, key: str, default: Any = None) -> Any:
        """
        Get a user preference value.
        
        Args:
            key: Preference key
            default: Default value if preference not found
            
        Returns:
            Preference value or default
        """
        config = self.load_config()
        return config.get("preferences", {}).get(key, default)
    
    def set_preference(self, key: str, value: Any):
        """
        Set a user preference value.
        
        Args:
            key: Preference key
            value: Preference value
        """
        config = self.load_config()
        if "preferences" not in config:
            config["preferences"] = {}
        config["preferences"][key] = value
        self.save_config(config)
    
    def cleanup_old_files(self, days_old: int = 30):
        """
        Clean up old cache and log files.
        
        Args:
            days_old: Files older than this many days will be removed
        """
        import time
        
        cutoff_time = time.time() - (days_old * 24 * 60 * 60)
        
        for directory in [self.cache_dir, self.logs_dir]:
            if directory.exists():
                for file_path in directory.iterdir():
                    if file_path.is_file() and file_path.stat().st_mtime < cutoff_time:
                        try:
                            file_path.unlink()
                            print(f"🗑️  Cleaned up old file: {file_path.name}")
                        except Exception as e:
                            print(f"⚠️  Could not remove {file_path.name}: {e}")
    
    def get_info(self) -> Dict[str, Any]:
        """
        Get comprehensive information about the HBAT installation.
        
        Returns:
            Dictionary with installation and configuration information
        """
        config = self.load_config()
        state = self.load_state()
        
        return {
            "hbat_directory": str(self.hbat_dir),
            "ccd_data_directory": str(self.ccd_dir),
            "config_exists": self.config_file.exists(),
            "state_exists": self.state_file.exists(),
            "is_first_run": self.is_first_run(),
            "run_count": state.get("run_count", 0),
            "last_used": state.get("last_used"),
            "ccd_files_present": config.get("ccd_data", {}).get("files_present", False),
            "recent_files_count": len(state.get("recent_files", [])),
            "preferences": config.get("preferences", {}),
            "directory_sizes": self._get_directory_sizes()
        }
    
    def _get_directory_sizes(self) -> Dict[str, str]:
        """Get sizes of various directories in human-readable format."""
        def get_dir_size(path: Path) -> str:
            if not path.exists():
                return "0 B"
            
            total_size = sum(f.stat().st_size for f in path.rglob('*') if f.is_file())
            
            # Convert to human readable
            for unit in ['B', 'KB', 'MB', 'GB']:
                if total_size < 1024.0:
                    return f"{total_size:.1f} {unit}"
                total_size /= 1024.0
            return f"{total_size:.1f} TB"
        
        return {
            "total": get_dir_size(self.hbat_dir),
            "ccd_data": get_dir_size(self.ccd_dir),
            "cache": get_dir_size(self.cache_dir),
            "logs": get_dir_size(self.logs_dir)
        }


# Global configuration instance
_hbat_config = None

def get_hbat_config() -> HBATConfig:
    """
    Get the global HBAT configuration instance.
    
    Returns:
        HBATConfig instance
    """
    global _hbat_config
    if _hbat_config is None:
        _hbat_config = HBATConfig()
    return _hbat_config


def initialize_hbat_environment(verbose: bool = True) -> bool:
    """
    Initialize the HBAT environment on first run.
    
    This function should be called when HBAT CLI or GUI starts up.
    
    Args:
        verbose: Whether to print initialization messages
        
    Returns:
        True if initialization successful, False otherwise
    """
    config = get_hbat_config()
    
    if config.is_first_run() and verbose:
        print("🚀 Welcome to HBAT!")
        print("   Setting up your personal HBAT environment...")
    
    # Ensure directory structure exists
    success = config.ensure_hbat_directory()
    
    if success:
        # Update run count
        config.update_run_count()
        
        if verbose:
            info = config.get_info()
            print(f"✅ HBAT environment ready")
            print(f"   📁 Data directory: {info['hbat_directory']}")
            print(f"   🗃️  CCD data path: {info['ccd_data_directory']}")
            
            if info['is_first_run']:
                print("   🆕 This is your first time running HBAT!")
            else:
                print(f"   📊 Run count: {info['run_count']}")
    
    return success