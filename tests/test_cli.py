import unittest

from valkyrie.entrypoints.main import load_command, parse_args


class CliTest(unittest.TestCase):
    def test_qe_scf_dispatches_to_qe_scf(self):
        self.assertEqual(load_command("qe_scf").__module__, "valkyrie.entrypoints.qe_scf")

    def test_abacus_scf_dispatches_to_abacus_scf(self):
        self.assertEqual(load_command("abacus_scf").__module__, "valkyrie.entrypoints.abacus_scf")

    def test_qe_scf_help_exists(self):
        with self.assertRaises(SystemExit) as ctx:
            parse_args(["qe_scf", "--help"])
        self.assertEqual(ctx.exception.code, 0)

    def test_abacus_scf_help_exists(self):
        with self.assertRaises(SystemExit) as ctx:
            parse_args(["abacus_scf", "--help"])
        self.assertEqual(ctx.exception.code, 0)

    def test_nodes_default_is_int(self):
        args = parse_args(["vasp_scf"])
        self.assertEqual(args.nodes, 24)

    def test_utility_commands_dispatch(self):
        self.assertEqual(load_command("check").__module__, "valkyrie.entrypoints.check")
        self.assertEqual(load_command("show-config").__module__, "valkyrie.entrypoints.show_config")
        self.assertEqual(load_command("status").__module__, "valkyrie.entrypoints.status")

    def test_utility_command_args(self):
        check_args = parse_args(["check", "--json"])
        self.assertTrue(check_args.json_output)
        self.assertFalse(hasattr(check_args, "nodes"))

        config_args = parse_args(["show-config", "--json"])
        self.assertTrue(config_args.json_output)

        status_args = parse_args(["status", "--engine", "vasp", "--history", "2"])
        self.assertEqual(status_args.engine, "vasp")
        self.assertEqual(status_args.history, 2)

    def test_no_record_global_arg(self):
        args = parse_args(["--no-record", "vasp_scf"])
        self.assertTrue(args.no_record)


if __name__ == "__main__":
    unittest.main()
