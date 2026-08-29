"""Helpers shared by NiceGUI server workflow tests."""

from pathlib import Path

from nicegui import ui
from nicegui.elements.upload_files import SmallFileUpload

from tests.fixtures.fixed_pdb_cases import PRIMARY_INTERACTION_ATTRIBUTES


async def upload_fixed_structure(user, case) -> None:
    """Upload one pre-fixed structure through the real upload component."""
    upload_match = user.find(kind=ui.upload, marker="structure-upload")
    assert len(upload_match.elements) == 1
    upload = next(iter(upload_match.elements))
    pdb_path = Path(case["file"])
    await upload.handle_uploads(
        [
            SmallFileUpload(
                name=case["name"],
                content_type="chemical/x-pdb",
                _data=pdb_path.read_bytes(),
            )
        ]
    )
    await user.should_see(marker="upload-status", content=case["name"], retries=100)


async def disable_pdb_fixing(user) -> None:
    """Disable fixing because the regression structures are already fixed."""
    user.find(marker="edit-pdb-fixing").click()
    await user.should_see(marker="fix-pdb-enabled", retries=20)

    switch_match = user.find(kind=ui.switch, marker="fix-pdb-enabled")
    assert len(switch_match.elements) == 1
    switch = next(iter(switch_match.elements))
    assert switch.value is True
    switch_match.click()
    assert switch.value is False
    user.find(marker="save-pdb-fixing").click()


async def analyze_fixed_structure(user, case) -> None:
    """Run the complete upload/configure/analyze workflow for a fixed PDB."""
    await upload_fixed_structure(user, case)
    await disable_pdb_fixing(user)
    user.find(marker="analyze").click()
    await user.should_see(marker="results-ready", content=case["name"], retries=1200)


async def assert_expected_ui_results(user, case) -> None:
    """Assert exact summary counts and every expected non-empty result tab."""
    expected_counts = case["expected_counts"]
    for interaction_type, expected_count in expected_counts.items():
        await user.should_see(
            marker=f"summary-{interaction_type}-count",
            content=str(expected_count),
            retries=20,
        )

    expected_total = sum(
        expected_counts[interaction_type]
        for interaction_type in PRIMARY_INTERACTION_ATTRIBUTES
    )
    await user.should_see(
        marker="summary-total_interactions-count",
        content=str(expected_total),
        retries=20,
    )

    for interaction_type in case["expected_interactions"]:
        await user.should_see(marker=f"results-tab-{interaction_type}", retries=20)
