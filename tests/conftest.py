"""Pytest configuration and fixtures for microDNA tests."""

import pytest
from pathlib import Path

# Root of the project
PROJECT_ROOT = Path(__file__).parent.parent


@pytest.fixture
def sample_data_dir() -> Path:
    """Return the path to the sample data directory."""
    return PROJECT_ROOT / "data"


@pytest.fixture
def src_scripts_dir() -> Path:
    """Return the path to the Python scripts directory."""
    return PROJECT_ROOT / "src" / "microdna"
