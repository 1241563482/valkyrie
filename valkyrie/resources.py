"""Helpers for locating and copying package resources."""

from importlib.resources import files
from pathlib import Path
import shutil


_PACKAGE = "valkyrie"


def resource_path(*parts):
    """Return a filesystem path string for a resource inside ``valkyrie``."""
    return str(files(_PACKAGE).joinpath(*parts))


def resource_dir(name):
    """Return a package resource directory path."""
    return resource_path(name)


def shell_template(name):
    """Return the path to a shell/job template bundled with Valkyrie."""
    return resource_path("shell_scripts", name)


def python_script(name):
    """Return the path to a bundled helper script."""
    return resource_path("python_scripts", name)


def copy_python_script(name, destination=None):
    """Copy a bundled helper script into the working directory.

    Parameters
    ----------
    name
        File name under ``valkyrie/python_scripts``.
    destination
        Output path. If omitted, the original file name is used.
    """
    src = Path(python_script(name))
    dst = Path(destination or name)
    shutil.copy2(src, dst)
    return dst


def copy_shell_script(name, destination=None):
    """Copy a bundled shell helper into the working directory."""
    src = Path(shell_template(name))
    dst = Path(destination or name)
    shutil.copy2(src, dst)
    return dst
