import tempfile
from pathlib import Path
import unittest

from valkyrie.entrypoints.keep import plan_removals


class KeepTest(unittest.TestCase):
    def test_plan_removals_keeps_patterns_and_skips_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            keep = tmp_path / "INCAR"
            delete = tmp_path / "WAVECAR"
            directory = tmp_path / "tmp"
            keep.write_text("keep")
            delete.write_text("delete")
            directory.mkdir()

            removals = plan_removals(tmp_path.iterdir(), ["INCAR"])

            self.assertIn(delete, removals)
            self.assertNotIn(keep, removals)
            self.assertNotIn(directory, removals)


if __name__ == "__main__":
    unittest.main()
