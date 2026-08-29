"""Isolated NiceGUI entry point shared by user and screen tests."""

import importlib
import os

from nicegui import Client

if screen_port := os.getenv("NICEGUI_SCREEN_TEST_PORT"):
    os.environ["HBAT_PORT"] = screen_port

os.environ["HBAT_ENV"] = "production"
os.environ["HBAT_ANALYTICS_ENABLED"] = "false"
os.environ["HBAT_SESSION_CLEANUP_ENABLED"] = "false"

importlib.import_module("hbat.server.app").create_app()

# NiceGUI's user-test cleanup removes the module (and all parent modules) that
# owns each registered page. The production page is defined in
# ``hbat.server.app``, so leaving its original owner here would remove ``hbat``
# from ``sys.modules`` after every simulated-user test and corrupt later tests'
# package attribute resolution. Mark the registered page as owned by this test
# entry point so cleanup remains confined to test state.
for page in Client.page_routes:
    if page.__module__.startswith("hbat.server"):
        page.__module__ = "tests.e2e.server.nicegui_main"
