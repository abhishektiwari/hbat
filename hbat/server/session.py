"""
Session management for HBAT web server.

This module provides session-based file isolation for multi-user scenarios,
ensuring that each user's uploaded files and analysis results are kept separate.
"""

import shutil
import uuid
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional


class SessionManager:
    """Manages user sessions with isolated file storage."""

    def __init__(self, base_dir: Path, session_timeout_hours: int = 24):
        """Initialize session manager.

        :param base_dir: Base directory for all session data
        :param session_timeout_hours: Hours before a session is considered expired
        """
        self.base_dir = base_dir
        self.session_timeout_hours = session_timeout_hours
        self.base_dir.mkdir(exist_ok=True)

    def create_session(self) -> str:
        """Create a new session with a unique UUID.

        :returns: Session ID (UUID)
        :rtype: str
        """
        session_id = str(uuid.uuid4())
        session_dir = self.get_session_dir(session_id)
        session_dir.mkdir(parents=True, exist_ok=True)

        # Create a timestamp file to track session creation
        timestamp_file = session_dir / ".session_timestamp"
        timestamp_file.write_text(datetime.now().isoformat())

        return session_id

    def get_session_dir(self, session_id: str) -> Path:
        """Get the directory path for a session.

        :param session_id: Session ID (UUID)
        :returns: Path to session directory
        :rtype: Path
        """
        return self.base_dir / session_id

    def session_exists(self, session_id: str) -> bool:
        """Check if a session exists and is valid.

        :param session_id: Session ID (UUID)
        :returns: True if session exists and is valid
        :rtype: bool
        """
        session_dir = self.get_session_dir(session_id)
        if not session_dir.exists():
            return False

        timestamp_file = session_dir / ".session_timestamp"
        if not timestamp_file.exists():
            return False

        # Check if session has expired
        try:
            timestamp = datetime.fromisoformat(timestamp_file.read_text())
            expiry_time = timestamp + timedelta(hours=self.session_timeout_hours)
            return datetime.now() < expiry_time
        except (ValueError, OSError):
            return False

    def cleanup_expired_sessions(self) -> int:
        """Clean up expired sessions.

        :returns: Number of sessions cleaned up
        :rtype: int
        """
        if not self.base_dir.exists():
            return 0

        cleaned_count = 0
        for session_dir in self.base_dir.iterdir():
            if not session_dir.is_dir():
                continue

            timestamp_file = session_dir / ".session_timestamp"
            if not timestamp_file.exists():
                # No timestamp file, remove the directory
                try:
                    shutil.rmtree(session_dir)
                    cleaned_count += 1
                except OSError:
                    pass
                continue

            try:
                timestamp = datetime.fromisoformat(timestamp_file.read_text())
                expiry_time = timestamp + timedelta(hours=self.session_timeout_hours)

                if datetime.now() >= expiry_time:
                    # Session has expired, remove it
                    shutil.rmtree(session_dir)
                    cleaned_count += 1
            except (ValueError, OSError):
                # Invalid timestamp or error reading, remove the directory
                try:
                    shutil.rmtree(session_dir)
                    cleaned_count += 1
                except OSError:
                    pass

        return cleaned_count

    def delete_session(self, session_id: str) -> bool:
        """Delete a specific session and all its files.

        :param session_id: Session ID (UUID)
        :returns: True if session was deleted, False if it didn't exist
        :rtype: bool
        """
        session_dir = self.get_session_dir(session_id)
        if not session_dir.exists():
            return False

        try:
            shutil.rmtree(session_dir)
            return True
        except OSError:
            return False

    def get_session_file_path(
        self, session_id: str, filename: str, create_dirs: bool = True
    ) -> Path:
        """Get the full path for a file within a session.

        :param session_id: Session ID (UUID)
        :param filename: Name of the file
        :param create_dirs: Whether to create session directory if it doesn't exist
        :returns: Full path to the file
        :rtype: Path
        """
        session_dir = self.get_session_dir(session_id)
        if create_dirs and not session_dir.exists():
            session_dir.mkdir(parents=True, exist_ok=True)

        return session_dir / filename
