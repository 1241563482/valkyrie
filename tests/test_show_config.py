import json
import unittest

from valkyrie.entrypoints.show_config import collect_config


class ShowConfigTest(unittest.TestCase):
    def test_collect_config_has_expected_sections(self):
        data = collect_config()
        self.assertIn("version", data)
        self.assertIn("scheduler", data)
        self.assertIn("job", data)
        self.assertIn("executables", data)
        self.assertIn("paths", data)
        self.assertIn("run_vasp", data["executables"])
        self.assertIn("run_pw", data["executables"])
        self.assertIn("run_abacus", data["executables"])
        self.assertIn("abacus_pseudo_dir", data["paths"])
        self.assertIn("abacus_orbital_dir", data["paths"])

    def test_collect_config_is_json_serializable(self):
        json.dumps(collect_config())


if __name__ == "__main__":
    unittest.main()
