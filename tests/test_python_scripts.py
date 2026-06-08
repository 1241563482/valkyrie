import importlib
import unittest


class PythonScriptsTest(unittest.TestCase):
    def test_helper_scripts_import_without_running_workflows(self):
        modules = [
            "dmft_get_command",
            "get_free_energy",
            "qe_elph_alpha2F",
            "qe_elph_sumlambda",
            "vasp_cohp",
            "vasp_cohp_plot",
            "vasp_dos_sum_spd",
            "vasp_elf_electride",
            "vasp_md_MDdata",
            "vasp_md_reprocess",
            "vasp_ph_loop_in",
        ]
        for module in modules:
            with self.subTest(module=module):
                importlib.import_module(f"valkyrie.python_scripts.{module}")


if __name__ == "__main__":
    unittest.main()
