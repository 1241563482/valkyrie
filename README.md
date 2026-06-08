# Valkyrie

[中文文档](README-cn.md)

Valkyrie is a command-line toolkit for preparing and managing first-principles calculation workflows. It focuses on generating input files, job scripts, and helper post-processing scripts for VASP, Quantum ESPRESSO, phonopy, lobster/COHP, and related workflows on an LSF-style cluster.

> Valkyrie automates the repetitive setup around calculations; it does not replace the underlying simulation codes. You still need working installations of VASP, Quantum ESPRESSO, phonopy, vaspkit, lobster, WIEN2k/EDMFT, and a compatible scheduler where relevant.

## Features

- Generate VASP inputs for relaxation, SCF, band structure, DOS, ELF, phonon, AIMD, and COHP workflows.
- Generate Quantum ESPRESSO relaxation and SCF inputs from `POSCAR`.
- Generate LSF job scripts from bundled templates.
- Copy task-specific helper scripts into the calculation directory.
- Provide safe cleanup helpers for calculation directories.
- Package resources through `pyproject.toml`, so shell templates and helper scripts are available after installation.

## Installation

Valkyrie requires Python 3.8 or newer.

```bash
cd /path/to/valkyrie
python -m pip install -e .
```

For this project, the package has been tested in the `magus` conda environment:

```bash
conda activate magus
python -m pip install -e .
```

If your shell does not activate the expected environment, use the environment's Python explicitly:

```bash
/fsb/home/jiansun/js_zhuyj/.conda/envs/magus/bin/python -m pip install -e .
```

## External Requirements

Depending on the command, Valkyrie expects the following external tools to be available on `PATH` or configured in `valkyrie/config.py`:

- VASP executables, such as `vasp_std`, `vasp_gam`, or `vasprelax`
- Quantum ESPRESSO `pw.x`
- `vaspkit`
- `phonopy`
- `lobster`
- LSF `bsub`
- WIEN2k/EDMFT tools for `dmft`

Cluster-specific defaults, executable commands, queue settings, and pseudopotential paths are currently defined in:

```text
valkyrie/config.py
```

You should review this file before using Valkyrie on a new machine or cluster.

## Quick Start

Most calculation commands assume that the current directory contains a valid `POSCAR` file. VASP workflows also need a usable `POTCAR` or a working `vaspkit` setup for automatic POTCAR generation.

Always start with `-ns`/`--notSub` when testing a command. This generates input files without submitting the job.

```bash
val vasp_relax -ns
val vasp_scf -ns
val vasp_dos -ns
val qe_relax -ns
```

After checking the generated files, run the same command without `-ns` to submit the job:

```bash
val vasp_relax
```

## Command Overview

```bash
val --help
```

Available subcommands:

| Command | Description |
| --- | --- |
| `keep` | Keep selected files and delete other regular files in the current directory. |
| `vasp_relax` | Generate VASP relaxation inputs and job script. |
| `vasp_scf` | Generate VASP SCF inputs and job script. |
| `vasp_band` | Generate VASP band structure workflow inputs. |
| `vasp_dos` | Generate VASP DOS workflow inputs and helper script. |
| `vasp_elf` | Generate VASP ELF calculation inputs and helper script. |
| `vasp_ph` | Generate VASP + phonopy finite-displacement workflow inputs. |
| `vasp_md` | Generate VASP AIMD supercell inputs and helper scripts. |
| `vasp_cohp` | Generate VASP + lobster COHP workflow inputs. |
| `qe_relax` | Generate Quantum ESPRESSO variable-cell relaxation input. |
| `qescf` | Legacy alias for QE SCF input generation. |
| `qe_scf` | Generate Quantum ESPRESSO SCF input. |
| `dmft` | Prepare a WIEN2k/EDMFT-style workflow. |

## Common Job Options

Most calculation commands accept:

```bash
-ns, --notSub       Generate inputs only; do not submit the job
-q,  --queue        LSF queue name
-n,  --nodes        Number of cores
-c,  --comment      Job name/comment
```

Example:

```bash
val vasp_relax -ns -q 9242opa! -n 24 -c LiIn-relax
```

## VASP Workflows

### Relaxation

```bash
val vasp_relax -ns
val vasp_relax -ns --encut 520 --pot Li_sv In_d
val vasp_relax -ns -p 10
val vasp_relax -ns --optcell
```

Useful options:

- `--encut`: set `ENCUT` manually. If omitted or `0`, Valkyrie reads `ENMAX` from `POTCAR` and uses `1.5 * max(ENMAX)`.
- `--pot`: specify POTCAR labels for `vaspkit`; default is `auto`.
- `-p/--pressure`: pressure in GPa.
- `--optcell`: generate constrained 2D-style relaxation settings.
- `--spin`: enable spin polarization.
- `--fun gga|ggau|hse`: choose functional mode.
- `-u ATOM U`: enable GGA+U for a specific atom.
- `-fd/--fermiDirac T`: use Fermi-Dirac smearing at temperature `T`.

### SCF, Band, DOS, ELF

```bash
val vasp_scf -ns
val vasp_band -ns --fun gga
val vasp_band -ns --fun hse
val vasp_dos -ns
val vasp_elf -ns
```

Band and DOS workflows generate multiple INCAR files and may call `vaspkit` to generate K-path or DOS-related files.

### Phonon

```bash
val vasp_ph -ns -d 2 2 2 -k 3 3 3
```

- `-d/--dim`: supercell size for phonopy.
- `-k/--kpoints`: VASP KPOINTS mesh for displaced structures.

### AIMD

```bash
val vasp_md -ns -t 300 1000 -d 3 3 3 -k 1 1 1
```

- `-t/--temperature`: `TEBEG` and `TEEND`.
- `-d/--dim`: supercell size.
- `-k/--kpoints`: KPOINTS mesh. `1 1 1` uses the gamma executable configured in `valkyrie/config.py`.

### COHP

```bash
val vasp_cohp -ns
```

This prepares VASP/lobster inputs and copies helper scripts for COHP analysis.

## Quantum ESPRESSO Workflows

### Relaxation

```bash
val qe_relax -ns --encut 100
val qe_relax -ns --encut 120 --spin
val qe_relax -ns --optcell
```

### SCF

Both command names are supported:

```bash
val qe_scf -ns --encut 100
val qescf -ns --encut 100
```

`qescf` is kept as a backward-compatible alias. Prefer `qe_scf` for new workflows.

## Cleanup Helper

The `keep` command deletes regular files in the current directory that do not match the keep patterns. Directories are skipped.

Preview before deleting:

```bash
val keep --dry-run --keepFiles POSCAR POTCAR INCAR job "*.py"
```

Use a predefined task profile:

```bash
val keep --dry-run --task vasp_relax
val keep --task vasp_relax
```

## Project Layout

```text
valkyrie/
├── entrypoints/       # CLI subcommands
├── input/             # Input-file generation helpers
├── python_scripts/    # Helper scripts copied into calculation directories
├── shell_scripts/     # Job templates and scheduler helper module
├── config.py          # Cluster/site defaults
├── resources.py       # Package resource access helpers
└── io_utils.py        # Small IO/process helpers

examples/              # Example VASP input files
tests/                 # Unit tests and legacy fixtures
```

## Development

Run tests with the standard library test runner:

```bash
python -m unittest discover -s tests -v
```

Build a wheel and check packaged resources:

```bash
python -m pip wheel . --no-deps --wheel-dir /tmp/valkyrie-wheel
```

Install in editable mode:

```bash
python -m pip install -e .
```

## Notes and Limitations

- Many defaults are site-specific and reflect the original cluster environment. Review `valkyrie/config.py` before production use.
- Some workflows invoke external tools such as `vaspkit`, `phonopy`, `bsub`, or `init_lapw`; these must be installed separately.
- `-ns/--notSub` is recommended for checking generated inputs before submitting jobs.
- The project currently targets practical calculation setup rather than being a fully general workflow manager.

## License

MIT
