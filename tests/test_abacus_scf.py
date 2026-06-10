import os
import tempfile
from pathlib import Path
import unittest

from valkyrie.entrypoints import abacus_scf


MINIMAL_POSCAR = """In2
1.0
3.25947054 0.0 0.0
0.0 3.25947054 0.0
0.0 0.0 4.98579870
In
2
direct
0.0 0.0 0.0
0.5 0.5 0.5
"""

MINIMAL_INPUT = """INPUT_PARAMETERS
pseudo_dir /pseudo
orbital_dir /orbital
calculation scf
ecutwfc 120
basis_type lcao
noncolin 1
lspinorb 1
nspin 4
"""

MINIMAL_KPT = """K_POINTS
0
Gamma
12 15 1 0 0 0
"""

MINIMAL_STRU = """ATOMIC_SPECIES
In 114.818 In.upf

NUMERICAL_ORBITAL
In.orb

LATTICE_CONSTANT
1.8897261258369282

LATTICE_VECTORS
1 0 0
0 1 0
0 0 1

ATOMIC_POSITIONS
Direct

In
0.0
1
0 0 0 1 1 1 mag 0 0 0 sc 1 1 1
"""


class AbacusScfTest(unittest.TestCase):
    def test_parse_key_value_map(self):
        self.assertEqual(
            abacus_scf.parse_key_value_map(["V=V.upf", "Se=Se.upf"]),
            {"V": "V.upf", "Se": "Se.upf"},
        )

    def test_read_abacus_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "INPUT"
            path.write_text(MINIMAL_INPUT)
            params = abacus_scf.read_abacus_input(path)
        self.assertEqual(params["calculation"], "scf")
        self.assertEqual(params["ecutwfc"], 120)
        self.assertEqual(params["nspin"], 4)

    def test_read_abacus_kpt(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "KPT"
            path.write_text(MINIMAL_KPT)
            params = abacus_scf.read_abacus_kpt(path)
        self.assertEqual(params["kmode"], "Gamma")
        self.assertEqual(params["kpts"], [12, 15, 1])
        self.assertEqual(params["koffset"], [0, 0, 0])

    def test_read_stru_pp_basis(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "STRU"
            path.write_text(MINIMAL_STRU)
            pp, basis = abacus_scf.read_stru_pp_basis(path)
        self.assertEqual(pp, {"In": "In.upf"})
        self.assertEqual(basis, {"In": "In.orb"})

    def test_main_generates_abacus_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = Path.cwd()
            tmp = Path(tmpdir)
            (tmp / "POSCAR").write_text(MINIMAL_POSCAR)
            ref = tmp / "1"
            ref.mkdir()
            (ref / "INPUT").write_text(MINIMAL_INPUT)
            (ref / "KPT").write_text(MINIMAL_KPT)
            (ref / "STRU").write_text(MINIMAL_STRU)
            try:
                os.chdir(tmp)
                abacus_scf.main(notSub=True)
            finally:
                os.chdir(cwd)

            for name in ["INPUT", "KPT", "STRU", "job"]:
                self.assertTrue((tmp / name).exists(), name)
            self.assertIn("calculation", (tmp / "INPUT").read_text())
            self.assertIn("scf", (tmp / "INPUT").read_text())
            self.assertIn("K_POINTS", (tmp / "KPT").read_text())
            self.assertIn("12 15 1 0 0 0", (tmp / "KPT").read_text())
            self.assertIn("ATOMIC_SPECIES", (tmp / "STRU").read_text())


if __name__ == "__main__":
    unittest.main()
