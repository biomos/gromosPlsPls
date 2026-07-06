import os, pytest
from pathlib import Path


@pytest.fixture(scope="session")
def bin_dir():
    b = os.environ.get("BIN_PATH")
    if b:
        p = Path(b)
        if p.is_dir():
            return p
    for cand in [
        "build/bin",
        "BUILD/bin",
        "BUILD_deb12/bin",
        "BUILD_deb13/bin",
        "BUILD_ubu22/bin",
        "BUILD_ubu24/bin",
    ]:
        p = Path(cand)
        if p.is_dir():
            return p
    raise RuntimeError(
        "Set BIN_PATH to the directory with gromos++ programs, "
        "e.g. export BIN_PATH=$PWD/BUILD_*/bin"
    )


@pytest.fixture(scope="session")
def test_data_dir():
    t = os.environ.get("TEST_REPO")
    if t:
        return Path(t) / "gromosPlsPls_tests"
    for cand in ["../../gromos_test_files"]:
        p = Path(cand)
        if p.is_dir():
            return p / "gromosPlsPls_tests"
    raise RuntimeError(
        "Set TEST_REPO to the gromos_test_files directory, "
        "e.g. export TEST_REPO=$PWD/../../gromos_test_files"
    )
