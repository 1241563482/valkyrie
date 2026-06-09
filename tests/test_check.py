import tempfile
import unittest
from pathlib import Path
from unittest import mock

from valkyrie.entrypoints.check import (
    check_executable,
    check_local_inputs,
    command_executable,
    run_checks,
    summarize,
)


class CheckTest(unittest.TestCase):
    def test_command_executable_skips_mpi_launcher(self):
        self.assertEqual(command_executable("mpiexec.hydra vasp_std"), "vasp_std")
        self.assertEqual(command_executable("mpiexec.hydra /path/to/pw.x -npool 4"), "/path/to/pw.x")

    def test_check_executable_with_mocked_path(self):
        with mock.patch("shutil.which", return_value="/usr/bin/tool"):
            item = check_executable("tool test", "tool")
        self.assertEqual(item["status"], "ok")

    def test_missing_poscar_does_not_crash(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            items = check_local_inputs(tmpdir)
            names = [item["name"] for item in items]
            self.assertIn("POSCAR", names)

    def test_summarize_strict_counts_warnings_as_failure(self):
        summary = summarize([{"status": "warn"}], strict=True)
        self.assertTrue(summary["failed"])

    def test_run_checks_returns_items(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "POTCAR").write_text("ENMAX = 400;\n")
            checks = run_checks(tmpdir)
        self.assertTrue(checks)


if __name__ == "__main__":
    unittest.main()
