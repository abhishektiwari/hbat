"""Tests for optional Google Analytics behavior in the web server."""

import asyncio
from unittest.mock import AsyncMock, patch

import pytest

from hbat.server.app import (
    HBATWebApp,
    get_google_analytics_head_html,
    is_analytics_enabled,
)


@pytest.mark.unit
def test_analytics_is_disabled_by_default(monkeypatch):
    """Analytics is opt-in when the environment variable is absent."""
    monkeypatch.delenv("HBAT_ANALYTICS_ENABLED", raising=False)

    assert is_analytics_enabled() is False
    assert get_google_analytics_head_html() is None


@pytest.mark.unit
@pytest.mark.parametrize("value", ["0", "false", "no", "off", "FALSE"])
def test_analytics_disabled_values(monkeypatch, value):
    """False-like environment values disable analytics."""
    monkeypatch.setenv("HBAT_ANALYTICS_ENABLED", value)

    assert is_analytics_enabled() is False
    assert get_google_analytics_head_html() is None


@pytest.mark.unit
@pytest.mark.parametrize("value", ["1", "true", "yes", "on", "TRUE"])
def test_analytics_enabled_values(monkeypatch, value):
    """True-like environment values enable analytics."""
    monkeypatch.setenv("HBAT_ANALYTICS_ENABLED", value)

    assert is_analytics_enabled() is True
    head_html = get_google_analytics_head_html()
    assert head_html is not None
    assert "googletagmanager.com/gtag/js" in head_html
    assert "G-Y4J82QZJ50" in head_html


@pytest.mark.unit
def test_tracking_methods_are_noops_when_analytics_is_disabled(monkeypatch):
    """Disabled analytics does not execute browser JavaScript."""
    monkeypatch.delenv("HBAT_ANALYTICS_ENABLED", raising=False)
    web_app = HBATWebApp()

    with patch("hbat.server.app.ui.run_javascript", new_callable=AsyncMock) as run_js:
        asyncio.run(web_app._track_analysis_completion(1.25))
        asyncio.run(web_app._track_export("json"))

    run_js.assert_not_awaited()


@pytest.mark.unit
def test_tracking_methods_send_events_when_analytics_is_enabled(monkeypatch):
    """Enabled analytics sends completion and export events."""
    monkeypatch.setenv("HBAT_ANALYTICS_ENABLED", "true")
    web_app = HBATWebApp()
    web_app.current_file = "1abc.pdb"

    with patch("hbat.server.app.ui.run_javascript", new_callable=AsyncMock) as run_js:
        asyncio.run(web_app._track_analysis_completion(1.25))
        asyncio.run(web_app._track_export("json"))

    assert run_js.await_count == 2
    javascript = "\n".join(call.args[0] for call in run_js.await_args_list)
    assert "analysis_completed" in javascript
    assert "analysis_export" in javascript
