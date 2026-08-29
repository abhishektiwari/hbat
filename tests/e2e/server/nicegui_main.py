"""Isolated NiceGUI entry point shared by user and screen tests."""

import importlib
import os

if screen_port := os.getenv("NICEGUI_SCREEN_TEST_PORT"):
    os.environ["HBAT_PORT"] = screen_port

os.environ["HBAT_ENV"] = "production"
os.environ["HBAT_ANALYTICS_ENABLED"] = "false"
os.environ["HBAT_SESSION_CLEANUP_ENABLED"] = "false"

importlib.import_module("hbat.server.app").create_app()
