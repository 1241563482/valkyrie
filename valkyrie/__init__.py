from .config import (
    environment,
    job_head,
    abacus_orbital_dir,
    abacus_pseudo_dir,
    qe_pot_path,
    run_abacus,
    run_pw,
    run_vasp,
    run_vasp_gam,
    run_vasp_opt,
)
from .resources import resource_dir


__version__ = "0.1"
__picture__ = rf"""
                 _  _                  _
  /\   /\  __ _ | || | __ _   _  _ __ (_)  ___
  \ \ / / / _` || || |/ /| | | || '__|| | / _ \
   \ V / | (_| || ||   < | |_| || |   | ||  __/
    \_/   \__,_||_||_|\_\ \__, ||_|   |_| \___|
            Version {__version__.ljust(5)}  |___/
                                by Yijie Zhu
""".strip("\n")

# Backwards-compatible resource directory variables used by existing entrypoints.
__shell__ = resource_dir("shell_scripts")
__python__ = resource_dir("python_scripts")
