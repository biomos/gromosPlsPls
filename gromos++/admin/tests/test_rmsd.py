import subprocess

_TRAJ_CNT = 5

_EXPECTED = [
    ("md_1.cnf", 0.39255),
    ("md_2.cnf", 0.36502),
    ("md_3.cnf", 0.34109),
    ("md_4.cnf", 0.51512),
    ("md_5.cnf", 0.39492),
]

_TOL = 0.005


def test_rmsd(bin_dir, test_data_dir):
    td = test_data_dir / "rmsd"

    traj = sorted(td.glob("md_*.cnf"))
    assert len(traj) == _TRAJ_CNT, \
        f"Expected {_TRAJ_CNT} md_*.cnf files in {td}, found {len(traj)}"

    cmd = [
        str(bin_dir / "rmsd"),
        "@topo", str(td / "dna_ion.top"),
        "@pbc", "r", "cog",
        "@reftopo", str(td / "dna.top"),
        "@ref", str(td / "ref.cnf"),
        "@atomsrmsd", "1-2:res(1-12:a)",
        "@traj",
    ] + [str(f) for f in traj]

    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, \
        f"rmsd exited with code {r.returncode}\nstderr:\n{r.stderr}"

    lines = [ln for ln in r.stdout.splitlines()
             if ln.strip() and not ln.strip().startswith("#")]
    assert len(lines) == len(_EXPECTED), \
        f"Expected {len(_EXPECTED)} data lines, got {len(lines)}:\n{r.stdout}"

    for (fname, exp_val), line in zip(_EXPECTED, lines):
        parts = line.strip().split()
        assert len(parts) == 2, \
            f"Unexpected line format for {fname}: {line!r}"
        val = float(parts[1])
    assert abs(val - exp_val) < _TOL, \
        f"{fname}: RMSD {val:.5f} differs from expected {exp_val:.5f} (tolerance {_TOL})"
