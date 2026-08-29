"""End-to-end HBAT server tests using NiceGUI's simulated user."""

import asyncio
import importlib
import json
import os
from pathlib import Path

import pytest

from tests.e2e.server.helpers import (
    analyze_fixed_structure,
    assert_expected_ui_results,
)
from tests.fixtures.fixed_pdb_cases import (
    FIXED_PDB_CASES,
    PRIMARY_INTERACTION_ATTRIBUTES,
)

pytestmark = [
    pytest.mark.e2e,
    pytest.mark.server,
    pytest.mark.requires_pdb_files,
    pytest.mark.nicegui_main_file("tests/e2e/server/nicegui_main.py"),
]


@pytest.mark.parametrize(
    "fixed_pdb_case",
    FIXED_PDB_CASES,
    ids=[case["name"] for case in FIXED_PDB_CASES],
)
async def test_fixed_pdb_analysis_workflow(user, fixed_pdb_case):
    """Upload, configure, analyze, and render exact results for each fixed PDB."""
    await user.open("/")
    await analyze_fixed_structure(user, fixed_pdb_case)
    await assert_expected_ui_results(user, fixed_pdb_case)


async def test_json_export_matches_rendered_analysis(user):
    """The downloaded JSON contains the same counts and interaction content as UI."""
    case = next(case for case in FIXED_PDB_CASES if case["name"] == "7nwd.pdb")
    await user.open("/")
    await analyze_fixed_structure(user, case)
    await assert_expected_ui_results(user, case)

    user.find(marker="nav-export").click()
    await user.should_see(marker="export-json", retries=20)
    download = asyncio.create_task(user.download.next(timeout=10))
    await asyncio.sleep(0)
    user.find(marker="export-json").click()
    response = await download
    assert response.status_code == 200

    exported = json.loads(response.content)
    assert exported["metadata"]["input_file"] == case["name"]
    assert {
        name: exported["summary"][name]["count"] for name in case["expected_counts"]
    } == case["expected_counts"]

    export_keys = {"pi_pi_interactions": "pi_pi_stacking"}
    for interaction_type in PRIMARY_INTERACTION_ATTRIBUTES:
        export_key = export_keys.get(interaction_type, interaction_type)
        expected_count = case["expected_counts"][interaction_type]
        if expected_count:
            assert len(exported[export_key]) == expected_count
        else:
            assert export_key not in exported
    assert (
        len(exported["cooperativity_chains"])
        == case["expected_counts"]["cooperativity_chains"]
    )


async def test_two_users_have_isolated_files_and_results(user, create_user):
    """Concurrent browser sessions cannot overwrite or display each other's PDB."""
    first_case = next(case for case in FIXED_PDB_CASES if case["name"] == "7nwd.pdb")
    second_case = next(case for case in FIXED_PDB_CASES if case["name"] == "1ubi.pdb")
    sessions_dir = Path(os.environ["HBAT_UPLOADS_DIR"]) / "sessions"
    before = set(sessions_dir.iterdir()) if sessions_dir.exists() else set()

    second_user = create_user()
    await user.open("/")
    await second_user.open("/")
    await analyze_fixed_structure(user, first_case)
    await analyze_fixed_structure(second_user, second_case)

    await assert_expected_ui_results(user, first_case)
    await assert_expected_ui_results(second_user, second_case)
    await user.should_not_see(content=second_case["name"])
    await second_user.should_not_see(content=first_case["name"])

    after = set(sessions_dir.iterdir())
    new_sessions = after - before
    assert len(new_sessions) == 2
    assert {
        path.name for directory in new_sessions for path in directory.glob("*.pdb")
    } == {first_case["name"], second_case["name"]}


async def test_user_simulation_cleanup_preserves_hbat_package():
    """NiceGUI cleanup must not detach HBAT subpackages used by later tests."""
    from nicegui.testing.user_simulation import user_simulation

    subpackages = {
        name: importlib.import_module(f"hbat.{name}")
        for name in ("core", "gui", "utilities", "visualization")
    }
    main_file = Path(__file__).parents[1] / "nicegui_main.py"

    async with user_simulation(main_file=main_file):
        pass

    hbat_package = importlib.import_module("hbat")
    for name, subpackage in subpackages.items():
        assert getattr(hbat_package, name) is subpackage
