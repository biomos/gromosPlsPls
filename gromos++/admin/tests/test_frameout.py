import subprocess

_FRAME_ATOMS = 578  # ATOM lines expected from dna.top


def test_cnf_to_pdb(bin_dir, test_data_dir, tmp_path):
    td = test_data_dir / "frameout"
    r = subprocess.run(
        [
            str(bin_dir / "frameout"),
            "@topo", str(td / "dna.top"),
            "@pbc", "v",
            "@notimeblock",
            "@outformat", "pdb",
            "@single",
            "@traj", str(td / "1naj.cnf"),
        ],
        capture_output=True, text=True,
        cwd=tmp_path,
    )
    assert r.returncode == 0, \
        f"frameout exited with code {r.returncode}\nstderr:\n{r.stderr}"

    pdb = tmp_path / "FRAME_00001.pdb"
    assert pdb.is_file(), f"Output file {pdb} was not created"

    atoms = 0
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM"):
            atoms += 1
            parts = line.split()
            assert len(parts) >= 8, f"Malformed ATOM line: {line!r}"
            x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
            assert abs(x) > 1e-6 or abs(y) > 1e-6 or abs(z) > 1e-6, \
                f"All-zero coordinates on line:\n{line}"

    assert atoms == _FRAME_ATOMS, \
        f"Expected {_FRAME_ATOMS} ATOM lines, got {atoms}"


def test_cnf_to_cif(bin_dir, test_data_dir, tmp_path):
    td = test_data_dir / "frameout"
    r = subprocess.run(
        [
            str(bin_dir / "frameout"),
            "@topo", str(td / "dna.top"),
            "@pbc", "v",
            "@notimeblock",
            "@outformat", "cif",
            "@single",
            "@traj", str(td / "1naj.cnf"),
        ],
        capture_output=True, text=True,
        cwd=tmp_path,
    )
    assert r.returncode == 0, \
        f"frameout exited with code {r.returncode}\nstderr:\n{r.stderr}"

    cif = tmp_path / "FRAME_00001.cif"
    assert cif.is_file(), f"Output file {cif} was not created"

    atoms = 0
    for line in cif.read_text().splitlines():
        if line.startswith("ATOM"):
            atoms += 1
            cols = line.split()
            assert len(cols) >= 10, f"Malformed atom_site line: {line!r}"
            x, y, z = float(cols[7]), float(cols[8]), float(cols[9])
            assert abs(x) > 1e-6 or abs(y) > 1e-6 or abs(z) > 1e-6, \
                f"All-zero coordinates on line:\n{line}"

    assert atoms == _FRAME_ATOMS, \
        f"Expected {_FRAME_ATOMS} atom_site rows, got {atoms}"
