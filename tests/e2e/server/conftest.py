"""Shared NiceGUI fixtures and isolated storage for server E2E tests."""

# ruff: noqa: I001

import atexit
import os
import shutil
import tempfile
from pathlib import Path

import pytest


pytest.importorskip("nicegui")

SERVER_TEST_UPLOADS_DIR = Path(tempfile.mkdtemp(prefix="hbat-server-e2e-"))
os.environ["HBAT_UPLOADS_DIR"] = str(SERVER_TEST_UPLOADS_DIR)
os.environ["HBAT_ANALYTICS_ENABLED"] = "false"
os.environ["HBAT_SESSION_CLEANUP_ENABLED"] = "false"
atexit.register(shutil.rmtree, SERVER_TEST_UPLOADS_DIR, ignore_errors=True)

# Import fixtures directly so the user and screen helpers can coexist without
# registering their duplicate pytest command-line hooks as global plugins.
from nicegui.testing.general_fixtures import (  # noqa: F401
    nicegui_reset_globals,
    pytest_addoption,
    pytest_configure,
)
from nicegui.testing.user_plugin import create_user, user  # noqa: F401
