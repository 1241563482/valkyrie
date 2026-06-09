import tempfile
import unittest
from pathlib import Path

from valkyrie.entrypoints.status import collect_status, parse_qe_output, parse_vasp_outcar


class StatusTest(unittest.TestCase):
    def test_parse_vasp_outcar_negative_energy(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outcar = Path(tmpdir) / "OUTCAR"
            outcar.write_text(
                "free  energy   TOTEN  =       -123.456789 eV\n"
                "reached required accuracy\n"
                "General timing and accounting informations for this job:\n"
            )
            result = parse_vasp_outcar(outcar)
        self.assertTrue(result["completed"])
        self.assertTrue(result["converged"])
        self.assertEqual(result["final_energy"], {"value": -123.456789, "unit": "eV"})
        self.assertEqual(result["steps"], 1)

    def test_parse_qe_output_negative_energy(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "scf.out"
            output.write_text(
                "iteration #  1\n"
                "!    total energy              =    -50.123456 Ry\n"
                "convergence has been achieved\n"
                "JOB DONE.\n"
            )
            result = parse_qe_output(output)
        self.assertTrue(result["completed"])
        self.assertTrue(result["converged"])
        self.assertEqual(result["final_energy"], {"value": -50.123456, "unit": "Ry"})
        self.assertEqual(result["steps"], 1)

    def test_collect_status_no_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            status = collect_status(tmpdir)
        self.assertEqual(status["state"], "no_output")


if __name__ == "__main__":
    unittest.main()
