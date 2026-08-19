import pytest

def pytest_addoption(parser):
    parser.addoption("--input", action="store", default=None)
    parser.addoption("--output", action="store", default=None)

@pytest.fixture
def input_path(request):
    value = request.config.getoption("--input")
    if value is None:
        pytest.fail("--input is required for this test.")
    return request.config.getoption("--input")

@pytest.fixture
def output_path(request):
    value = request.config.getoption("--output")
    if value is None:
        pytest.fail("--output is required for this test.")
    return value
