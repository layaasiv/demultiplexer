from pathlib import Path

import pytest


def pytest_addoption(parser):
    parser.addoption("--input", action="store", default=None)


@pytest.fixture
def input_path(request):
    value = request.config.getoption("--input")
    if value is not None:
        return Path(value)
    return Path(__file__).parent / "data"
