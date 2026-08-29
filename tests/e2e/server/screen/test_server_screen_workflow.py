"""Real-browser smoke test for the HBAT NiceGUI server."""

import importlib.util
import os
import shutil
from pathlib import Path

import pytest

from tests.fixtures.fixed_pdb_cases import (
    FIXED_PDB_CASES,
    PRIMARY_INTERACTION_ATTRIBUTES,
)


def _has_chrome() -> bool:
    """Return whether a local Chrome/Chromium executable is available."""
    return importlib.util.find_spec("selenium") is not None and bool(
        os.getenv("CHROME_BINARY_LOCATION")
        or shutil.which("google-chrome")
        or shutil.which("google-chrome-stable")
        or shutil.which("chromium")
        or shutil.which("chromium-browser")
        or shutil.which("chrome")
        or Path(
            "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
        ).is_file()
    )


pytestmark = [
    pytest.mark.e2e,
    pytest.mark.server,
    pytest.mark.browser,
    pytest.mark.slow,
    pytest.mark.requires_pdb_files,
    pytest.mark.nicegui_main_file("tests/e2e/server/nicegui_main.py"),
    pytest.mark.skipif(not _has_chrome(), reason="Chrome or Chromium is not installed"),
]


@pytest.mark.parametrize(
    "fixed_pdb_case",
    FIXED_PDB_CASES,
    ids=[case["name"] for case in FIXED_PDB_CASES],
)
def test_fixed_pdb_workflow_in_real_browser(screen, fixed_pdb_case):
    """Validate each fixed PDB from upload through exact browser-rendered results."""
    case = fixed_pdb_case
    screen.open("/", timeout=10)
    screen.allowed_js_errors.append("socket.dev/api/badge")

    file_input = screen.find_by_css(
        '[data-testid="structure-upload"] input[type="file"]'
    )
    file_input.send_keys(case["file"])
    screen.wait_for_js(
        "document.querySelector('[data-testid=\"edit-pdb-fixing\"]') !== null",
        True,
        timeout=20,
    )

    screen.find_by_css('[data-testid="edit-pdb-fixing"]').click()
    fixing_switch = screen.find_by_css('[data-testid="fix-pdb-enabled"]')
    assert fixing_switch.get_attribute("aria-checked") == "true"
    fixing_switch.click()
    screen.wait_for_js(
        "document.querySelector('[data-testid=\"fix-pdb-enabled\"]')"
        '.getAttribute("aria-checked")',
        "false",
        timeout=5,
    )
    screen.find_by_css('[data-testid="save-pdb-fixing"]').click()
    analyze_button = screen.find_by_css('[data-testid="analyze"]')
    screen.selenium.execute_script(
        "arguments[0].scrollIntoView({block: 'center'});", analyze_button
    )
    screen.wait(0.5)  # allow the parameter drawer's close animation to finish
    analyze_button.click()

    screen.wait_for_js(
        "document.querySelector('[data-testid=\"results-ready\"]')?.innerText.includes("
        f"'{Path(case['file']).name}')",
        True,
        timeout=120,
    )

    for interaction_type, expected_count in case["expected_counts"].items():
        count = screen.find_by_css(f'[data-testid="summary-{interaction_type}-count"]')
        assert count.text == str(expected_count)

    expected_total = sum(
        case["expected_counts"][interaction_type]
        for interaction_type in PRIMARY_INTERACTION_ATTRIBUTES
    )
    total = screen.find_by_css('[data-testid="summary-total_interactions-count"]')
    assert total.text == str(expected_total)

    expected_interactions = set(case["expected_interactions"])
    for interaction_type in PRIMARY_INTERACTION_ATTRIBUTES:
        selector = f'[data-testid="results-tab-{interaction_type}"]'
        tab_exists = screen.selenium.execute_script(
            "return document.querySelector(arguments[0]) !== null;", selector
        )
        assert tab_exists is (interaction_type in expected_interactions)

    # Leave the page at the top of the results so NiceGUI's teardown screenshot
    # captures the interaction tabs below the fixed application header.
    summary_tab = screen.find_by_css('[data-testid="results-tab-summary"]')
    screen.selenium.execute_script(
        "arguments[0].scrollIntoView({block: 'start'}); window.scrollBy(0, -88);",
        summary_tab,
    )
    screen.wait(0.25)
