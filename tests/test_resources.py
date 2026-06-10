from pathlib import Path
import unittest

from valkyrie.resources import python_script, shell_template


class ResourcesTest(unittest.TestCase):
    def test_shell_template_exists(self):
        self.assertTrue(Path(shell_template("job_vasp_relax")).is_file())
        self.assertTrue(Path(shell_template("job_abacus_scf")).is_file())

    def test_python_script_exists(self):
        self.assertTrue(Path(python_script("vasp_dos_sum_spd.py")).is_file())


if __name__ == "__main__":
    unittest.main()
