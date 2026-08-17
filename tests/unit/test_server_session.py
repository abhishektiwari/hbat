"""Unit tests for server session file isolation."""

import pytest

from hbat.server.session import SessionManager


@pytest.mark.unit
@pytest.mark.parametrize(
    "filename",
    ["../outside.pdb", "../../outside.pdb", "/tmp/outside.pdb", r"..\outside.pdb"],
)
def test_session_file_path_rejects_path_traversal(tmp_path, filename):
    """Client-controlled filenames cannot escape the session directory."""
    manager = SessionManager(tmp_path / "sessions")
    session_id = manager.create_session()

    with pytest.raises(ValueError, match="simple name"):
        manager.get_session_file_path(session_id, filename)


def test_session_file_path_accepts_simple_filename(tmp_path):
    """A normal uploaded filename remains inside its session directory."""
    manager = SessionManager(tmp_path / "sessions")
    session_id = manager.create_session()

    path = manager.get_session_file_path(session_id, "structure.pdb")

    assert path == manager.get_session_dir(session_id) / "structure.pdb"
