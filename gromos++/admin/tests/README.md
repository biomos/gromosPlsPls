# GROMOS++ simulation tests

Numerical regression tests for gromos++ programs (rmsd, frameout, ...).

## Running locally

```bash
export BIN_PATH=$PWD/BUILD_*/bin
export TEST_REPO=/path/to/gromos_test_files
cd admin/tests && pytest -v
```

Without env vars, conftest.py tries common BUILD_*/bin directories
and ../../gromos_test_files as fallbacks.

## Adding a new test

Create `test_<name>.py` in this directory. Use the `bin_dir` and
`test_data_dir` fixtures from conftest.py. Call the program with
`subprocess.run()`, parse output, and `assert` expected values.

Input files go in `gromos_test_files/gromosPlsPls_tests/<name>/`.

## Tolerance

Numerical assertions use 0.005 (5-decimal precision of the `rmsd` program).
