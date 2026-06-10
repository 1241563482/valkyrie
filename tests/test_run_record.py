import argparse
import json
import tempfile
import unittest
from pathlib import Path

from valkyrie.run_record import (
    build_run_record,
    load_history,
    load_latest_run,
    record_invocation,
)


class RunRecordTest(unittest.TestCase):
    def test_build_run_record_is_json_serializable(self):
        args = argparse.Namespace(command="vasp_scf", no_record=False, nodes=24, pot=["Li_sv", "O"])
        record = build_run_record(args, argv=["vasp_scf", "-ns"], cwd=".")
        self.assertEqual(record["command"], "vasp_scf")
        json.dumps(record)

    def test_record_invocation_writes_latest_and_history(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="vasp_scf", no_record=False, nodes=24)
            record = record_invocation(args, argv=["vasp_scf"], cwd=tmpdir)
            self.assertIsNotNone(record)
            metadata = Path(tmpdir) / ".valkyrie"
            self.assertTrue((metadata / "run.json").exists())
            self.assertTrue((metadata / "run_history.jsonl").exists())
            self.assertEqual(load_latest_run(tmpdir)["command"], "vasp_scf")
            self.assertEqual(len(load_history(tmpdir)), 1)

    def test_abacus_scf_is_recordable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="abacus_scf", no_record=False, nodes=24)
            record = record_invocation(args, argv=["abacus_scf", "-ns"], cwd=tmpdir)
            self.assertIsNotNone(record)
            self.assertEqual(load_latest_run(tmpdir)["command"], "abacus_scf")

    def test_utility_commands_are_skipped(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="status", no_record=False)
            self.assertIsNone(record_invocation(args, cwd=tmpdir))
            self.assertFalse((Path(tmpdir) / ".valkyrie").exists())

    def test_load_missing_metadata_gracefully(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self.assertIsNone(load_latest_run(tmpdir))
            self.assertEqual(load_history(tmpdir), [])


if __name__ == "__main__":
    unittest.main()
