"""Run metadata recording helpers for Valkyrie CLI commands."""

import getpass
import json
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import __version__


RECORDABLE_COMMANDS = {
    "vasp_relax",
    "vasp_scf",
    "vasp_band",
    "vasp_dos",
    "vasp_elf",
    "vasp_ph",
    "vasp_md",
    "vasp_cohp",
    "qe_relax",
    "qescf",
    "qe_scf",
    "dmft",
}

METADATA_DIR = ".valkyrie"
LATEST_FILE = "run.json"
HISTORY_FILE = "run_history.jsonl"


def is_recordable(command):
    """Return whether a command should create run metadata."""
    return command in RECORDABLE_COMMANDS


def _json_safe(value):
    """Convert argparse values to JSON-serializable objects."""
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    return str(value)


def build_run_record(args, argv=None, cwd=None):
    """Build a serializable run metadata record from parsed CLI args."""
    command = getattr(args, "command", None)
    parameters = {
        key: _json_safe(value)
        for key, value in vars(args).items()
        if key not in {"command", "no_record"}
    }
    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    cwd_path = Path(cwd or Path.cwd()).resolve()
    return {
        "schema_version": 1,
        "recorded_at": now,
        "cwd": str(cwd_path),
        "command": command,
        "argv": list(argv) if argv is not None else sys.argv[1:],
        "parameters": parameters,
        "valkyrie_version": __version__,
        "python": sys.version.split()[0],
        "user": getpass.getuser(),
        "host": socket.gethostname(),
    }


def metadata_dir(path=None):
    """Return the metadata directory for a calculation path."""
    return Path(path or Path.cwd()) / METADATA_DIR


def record_invocation(args, argv=None, cwd=None):
    """Record the latest and historical CLI invocation.

    Returns the record dictionary when a record was written, otherwise ``None``.
    """
    command = getattr(args, "command", None)
    if not is_recordable(command) or getattr(args, "no_record", False):
        return None

    record = build_run_record(args, argv=argv, cwd=cwd)
    directory = metadata_dir(cwd)
    directory.mkdir(exist_ok=True)

    latest_path = directory / LATEST_FILE
    history_path = directory / HISTORY_FILE
    latest_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    with history_path.open("a") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")
    return record


def load_latest_run(path=None):
    """Load latest run metadata, or return ``None`` when it is unavailable."""
    latest_path = metadata_dir(path) / LATEST_FILE
    if not latest_path.exists():
        return None
    try:
        return json.loads(latest_path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def load_history(path=None, limit=None):
    """Load run history records, newest last unless ``limit`` is provided."""
    history_path = metadata_dir(path) / HISTORY_FILE
    if not history_path.exists():
        return []
    records = []
    try:
        for line in history_path.read_text().splitlines():
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    except OSError:
        return []
    if limit is not None and limit > 0:
        return records[-limit:]
    return records
