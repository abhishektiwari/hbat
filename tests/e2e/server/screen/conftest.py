"""NiceGUI real-browser fixtures, scoped to screen tests only."""

# ruff: noqa: I001

import importlib.util
import os
import shutil
from pathlib import Path

import pytest


MACOS_CHROME = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
if not os.getenv("CHROME_BINARY_LOCATION") and MACOS_CHROME.is_file():
    os.environ["CHROME_BINARY_LOCATION"] = str(MACOS_CHROME)

HAS_CHROME = bool(
    os.getenv("CHROME_BINARY_LOCATION")
    or shutil.which("google-chrome")
    or shutil.which("google-chrome-stable")
    or shutil.which("chromium")
    or shutil.which("chromium-browser")
    or shutil.which("chrome")
)

if importlib.util.find_spec("selenium") and HAS_CHROME:
    from nicegui.testing.screen import Screen
    from nicegui.testing.screen_plugin import (  # noqa: F401
        nicegui_chrome_options,
        nicegui_driver,
        nicegui_remove_all_screenshots,
        pytest_configure,
        screen,
    )

    @pytest.hookimpl(tryfirst=True, hookwrapper=True)
    def pytest_runtest_makereport(item, call):
        """Keep NiceGUI's root screenshot only when a test does not fully pass."""
        outcome = yield
        report = outcome.get_result()
        setattr(item, f"rep_{report.when}", report)
        if report.when == "teardown" and report.passed:
            redundant_screenshot = Screen.SCREENSHOT_DIR / f"{item.name}.png"
            redundant_screenshot.unlink(missing_ok=True)
else:

    @pytest.fixture
    def screen():
        """Explain the optional dependency when browser tests are unavailable."""
        pytest.skip("Selenium and Chrome are required for NiceGUI screen tests")


@pytest.fixture(autouse=True)
def large_browser_viewport(screen):
    """Use an exact desktop viewport for screen tests and screenshots."""
    screen.selenium.set_window_rect(width=1440, height=1000)
    screen.selenium.execute_cdp_cmd(
        "Emulation.setDeviceMetricsOverride",
        {
            "width": 1440,
            "height": 1000,
            "deviceScaleFactor": 1,
            "mobile": False,
        },
    )
