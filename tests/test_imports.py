import importlib
import unittest


class ImportsTest(unittest.TestCase):
    def test_core_modules_import(self):
        modules = [
            "valkyrie",
            "valkyrie.config",
            "valkyrie.resources",
            "valkyrie.entrypoints.main",
            "valkyrie.abacus_ase_io",
            "valkyrie.entrypoints.abacus_scf",
            "valkyrie.entrypoints.qe_scf",
            "valkyrie.entrypoints.dmft",
            "valkyrie.shell_scripts.job",
        ]
        for module in modules:
            with self.subTest(module=module):
                importlib.import_module(module)


if __name__ == "__main__":
    unittest.main()
