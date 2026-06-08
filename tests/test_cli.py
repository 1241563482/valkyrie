import unittest

from valkyrie.entrypoints.main import load_command, parse_args


class CliTest(unittest.TestCase):
    def test_qescf_alias_dispatches_to_qe_scf(self):
        self.assertEqual(load_command("qescf").__module__, "valkyrie.entrypoints.qe_scf")
        self.assertEqual(load_command("qe_scf").__module__, "valkyrie.entrypoints.qe_scf")

    def test_qe_scf_help_exists(self):
        with self.assertRaises(SystemExit) as ctx:
            parse_args(["qe_scf", "--help"])
        self.assertEqual(ctx.exception.code, 0)

    def test_nodes_default_is_int(self):
        args = parse_args(["vasp_scf"])
        self.assertEqual(args.nodes, 24)


if __name__ == "__main__":
    unittest.main()
