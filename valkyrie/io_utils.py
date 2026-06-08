"""Small filesystem and process helpers used by entrypoints."""

from pathlib import Path
import subprocess


def run_command(command, *, input=None, stdout=None, stderr=None):
    """Run an external command and fail if it exits non-zero."""
    return subprocess.run(command, input=input, text=True, check=True, stdout=stdout, stderr=stderr)


def comment_out_key(path, key):
    """Comment out lines containing *key* in a text file."""
    file_path = Path(path)
    lines = file_path.read_text().splitlines(keepends=True)
    updated = []
    for line in lines:
        if key in line and not line.lstrip().startswith("#"):
            updated.append("#" + line)
        else:
            updated.append(line)
    file_path.write_text("".join(updated))


def remove_files(*patterns):
    """Remove files matching glob patterns in the current directory."""
    for pattern in patterns:
        for path in Path(".").glob(pattern):
            if path.is_file() or path.is_symlink():
                path.unlink(missing_ok=True)
