"""Show Valkyrie configuration."""

import json

from .. import __version__
from .. import config
from ..resources import resource_dir


def collect_config():
    """Collect active Valkyrie configuration in a JSON-serializable dict."""
    return {
        "version": __version__,
        "scheduler": {
            "comment_startswith": config.comment_startswith,
            "core_number_startswith": config.core_number_startswith,
            "node_name_startswith": config.node_name_startswith,
            "sub_command": config.sub_command,
        },
        "job": {
            "job_head": config.job_head,
            "environment": config.environment,
        },
        "executables": {
            "run_vasp": config.run_vasp,
            "run_vasp_opt": config.run_vasp_opt,
            "run_vasp_gam": config.run_vasp_gam,
            "run_pw": config.run_pw,
            "run_abacus": config.run_abacus,
        },
        "paths": {
            "qe_pot_path": config.qe_pot_path,
            "abacus_pseudo_dir": config.abacus_pseudo_dir,
            "abacus_orbital_dir": config.abacus_orbital_dir,
            "shell_scripts": resource_dir("shell_scripts"),
            "python_scripts": resource_dir("python_scripts"),
        },
    }


def print_text(data):
    """Print config in a readable grouped format."""
    print(f"Valkyrie version: {data['version']}")
    print("\nScheduler:")
    for key, value in data["scheduler"].items():
        print(f"  {key}: {value}")
    print("\nJob:")
    for key, value in data["job"].items():
        print(f"  {key}:")
        for line in value.strip("\n").splitlines():
            print(f"    {line}")
    print("\nExecutables:")
    for key, value in data["executables"].items():
        print(f"  {key}: {value}")
    print("\nPaths:")
    for key, value in data["paths"].items():
        print(f"  {key}: {value}")


def main(*args, json_output=False, **kwargs):
    data = collect_config()
    if json_output:
        print(json.dumps(data, indent=2, sort_keys=True))
    else:
        print_text(data)
    return 0


if __name__ == "__main__":
    main()
