# Valkyrie

[English README](README.md)

Valkyrie 是一个面向第一性原理计算工作流的命令行工具，主要用于生成输入文件、作业脚本和后处理辅助脚本。目前重点支持 VASP、Quantum ESPRESSO、phonopy、lobster/COHP、WIEN2k/EDMFT 等常见工作流，尤其适合 LSF 集群环境。

> Valkyrie 负责自动化计算前后的重复准备工作；它并不替代 VASP、Quantum ESPRESSO、phonopy、vaspkit、lobster、WIEN2k/EDMFT 等实际计算程序。相关外部软件仍需用户自行安装和配置。

## 功能特点

- 生成 VASP relax、SCF、band、DOS、ELF、phonon、AIMD、COHP 等计算输入。
- 从 `POSCAR` 生成 Quantum ESPRESSO relax / SCF 输入文件。
- 根据内置模板生成 LSF 作业脚本。
- 自动复制任务相关的后处理脚本到计算目录。
- 提供计算目录清理工具。
- 使用 `pyproject.toml` 进行现代 Python 打包，job 模板和辅助脚本会随包一起安装。

## 安装

Valkyrie 需要 Python 3.8 或更新版本。

```bash
cd /path/to/valkyrie
python -m pip install -e .
```

本项目已在 `magus` conda 环境中测试：

```bash
conda activate magus
python -m pip install -e .
```

如果 shell 中的 `conda activate` 没有切到预期环境，可以直接使用环境里的 Python：

```bash
/fsb/home/jiansun/js_zhuyj/.conda/envs/magus/bin/python -m pip install -e .
```

## 外部依赖

不同子命令会调用不同外部程序。根据实际任务，可能需要：

- VASP 可执行文件，例如 `vasp_std`、`vasp_gam`、`vasprelax`
- Quantum ESPRESSO `pw.x`
- `vaspkit`
- `phonopy`
- `lobster`
- LSF `bsub`
- WIEN2k/EDMFT 相关工具

集群相关默认配置、程序路径、提交命令、赝势路径等目前写在：

```text
valkyrie/config.py
```

如果要在新机器或新集群上使用，建议先检查并修改这个文件。

## 快速开始

大多数计算命令默认当前目录存在有效的 `POSCAR`。VASP 工作流还需要可用的 `POTCAR`，或者需要 `vaspkit` 能自动生成 POTCAR。

调试时建议先加 `-ns` / `--notSub`，只生成输入文件，不提交作业：

```bash
val vasp_relax -ns
val vasp_scf -ns
val vasp_dos -ns
val qe_relax -ns
```

检查生成的 `INCAR`、`KPOINTS`、`job` 等文件没有问题后，再去掉 `-ns` 提交：

```bash
val vasp_relax
```

## 命令概览

查看帮助：

```bash
val --help
```

当前子命令包括：

| 命令 | 说明 |
| --- | --- |
| `keep` | 保留指定文件，删除当前目录中其他普通文件。 |
| `vasp_relax` | 生成 VASP 结构优化输入和作业脚本。 |
| `vasp_scf` | 生成 VASP SCF 输入和作业脚本。 |
| `vasp_band` | 生成 VASP 能带计算工作流输入。 |
| `vasp_dos` | 生成 VASP DOS 工作流输入和后处理脚本。 |
| `vasp_elf` | 生成 VASP ELF 计算输入和分析脚本。 |
| `vasp_ph` | 生成 VASP + phonopy 有限位移声子工作流输入。 |
| `vasp_md` | 生成 VASP AIMD 超胞输入和后处理脚本。 |
| `vasp_cohp` | 生成 VASP + lobster COHP 工作流输入。 |
| `qe_relax` | 生成 Quantum ESPRESSO 变胞优化输入。 |
| `qescf` | QE SCF 输入生成的旧命令别名。 |
| `qe_scf` | 生成 Quantum ESPRESSO SCF 输入。 |
| `dmft` | 准备 WIEN2k/EDMFT 风格的工作流。 |

## 通用作业参数

大多数计算命令都支持：

```bash
-ns, --notSub       只生成输入，不提交作业
-q,  --queue        LSF 队列名
-n,  --nodes        核数
-c,  --comment      作业名/注释
```

例子：

```bash
val vasp_relax -ns -q 9242opa! -n 24 -c LiIn-relax
```

## VASP 工作流

### 结构优化

```bash
val vasp_relax -ns
val vasp_relax -ns --encut 520 --pot Li_sv In_d
val vasp_relax -ns -p 10
val vasp_relax -ns --optcell
```

常用选项：

- `--encut`：手动指定 `ENCUT`。如果为 `0` 或不指定，会从 `POTCAR` 读取 `ENMAX`，并使用 `1.5 * max(ENMAX)`。
- `--pot`：指定 `vaspkit` 使用的 POTCAR 类型，默认 `auto`。
- `-p/--pressure`：压力，单位 GPa。
- `--optcell`：生成 2D/受限晶胞优化相关设置。
- `--spin`：开启自旋极化。
- `--fun gga|ggau|hse`：选择泛函模式。
- `-u ATOM U`：为指定元素启用 GGA+U。
- `-fd/--fermiDirac T`：使用温度为 `T` 的 Fermi-Dirac smearing。

### SCF、能带、DOS、ELF

```bash
val vasp_scf -ns
val vasp_band -ns --fun gga
val vasp_band -ns --fun hse
val vasp_dos -ns
val vasp_elf -ns
```

其中 band 和 DOS 工作流会生成多个 INCAR 文件，并可能调用 `vaspkit` 生成 K-path 或 DOS 相关文件。

### 声子

```bash
val vasp_ph -ns -d 2 2 2 -k 3 3 3
```

- `-d/--dim`：phonopy 超胞大小。
- `-k/--kpoints`：位移结构使用的 VASP KPOINTS 网格。

### AIMD

```bash
val vasp_md -ns -t 300 1000 -d 3 3 3 -k 1 1 1
```

- `-t/--temperature`：对应 VASP 的 `TEBEG` 和 `TEEND`。
- `-d/--dim`：超胞大小。
- `-k/--kpoints`：KPOINTS 网格。`1 1 1` 时会使用 [valkyrie/config.py](valkyrie/config.py) 里配置的 gamma 版本 VASP 命令。

### COHP

```bash
val vasp_cohp -ns
```

该命令会准备 VASP/lobster 输入，并复制 COHP 分析脚本。

## Quantum ESPRESSO 工作流

### 结构优化

```bash
val qe_relax -ns --encut 100
val qe_relax -ns --encut 120 --spin
val qe_relax -ns --optcell
```

### SCF

新旧两个命令都可用：

```bash
val qe_scf -ns --encut 100
val qescf -ns --encut 100
```

`qescf` 是为了兼容旧用法保留的别名。新的工作流建议使用 `qe_scf`。

## 清理计算目录

`keep` 命令会删除当前目录中不匹配保留规则的普通文件。目录会被跳过，不会递归删除。

删除前建议先 dry-run：

```bash
val keep --dry-run --keepFiles POSCAR POTCAR INCAR job "*.py"
```

使用预定义任务规则：

```bash
val keep --dry-run --task vasp_relax
val keep --task vasp_relax
```

## 项目结构

```text
valkyrie/
├── entrypoints/       # CLI 子命令入口
├── input/             # 输入文件生成辅助模块
├── python_scripts/    # 会复制到计算目录的后处理脚本
├── shell_scripts/     # job 模板和调度器辅助模块
├── config.py          # 集群/站点默认配置
├── resources.py       # 包内资源访问工具
└── io_utils.py        # 小型 IO / 外部命令工具

examples/              # 示例 VASP 输入文件
tests/                 # 单元测试和旧测试 fixtures
```

## 开发与测试

运行测试：

```bash
python -m unittest discover -s tests -v
```

构建 wheel：

```bash
python -m pip wheel . --no-deps --wheel-dir /tmp/valkyrie-wheel
```

editable 安装：

```bash
python -m pip install -e .
```

## 注意事项和限制

- 很多默认值是针对原集群环境设置的，正式使用前请检查 [valkyrie/config.py](valkyrie/config.py)。
- 部分命令会调用外部程序，例如 `vaspkit`、`phonopy`、`bsub`、`init_lapw`，这些程序需要用户自行安装。
- 建议每次新任务先使用 `-ns/--notSub` 检查生成文件，再提交作业。
- Valkyrie 当前定位是实用型输入/作业准备工具，不是完整通用的工作流管理系统。

## 许可证

MIT
