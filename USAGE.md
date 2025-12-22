# Cheq CLI: User Manual

Welcome to the Cheq command-line interface! This manual explains how to install, configure, and run the `cheq` CLI to perform charge equilibration (QEq) calculations on molecular geometries.

## Table of Contents

- [Cheq CLI: User Manual](#cheq-cli-user-manual)
  - [Table of Contents](#table-of-contents)
  - [Preparing Input Data](#preparing-input-data)
    - [XYZ File Requirements](#xyz-file-requirements)
    - [Using Standard Input](#using-standard-input)
  - [Core Charge Equilibration (`cheq`)](#core-charge-equilibration-cheq)
    - [Basic Usage](#basic-usage)
    - [Argument Reference](#argument-reference)
      - [Positional Argument](#positional-argument)
      - [Output Options](#output-options)
      - [Calculation Options](#calculation-options)
      - [Solver Options](#solver-options)
      - [General Arguments](#general-arguments)
  - [Custom Parameter Files (`params.toml`)](#custom-parameter-files-paramstoml)
    - [File Structure](#file-structure)
    - [Example Parameter File](#example-parameter-file)
    - [Detailed Parameter Descriptions](#detailed-parameter-descriptions)
  - [Practical Examples (Use Cases)](#practical-examples-use-cases)
    - [Example 1: Quick Check with Defaults](#example-1-quick-check-with-defaults)
    - [Example 2: Structured Output to Disk](#example-2-structured-output-to-disk)
    - [Example 3: Custom Total Charge and Parameters](#example-3-custom-total-charge-and-parameters)
    - [Example 4: Machine-Readable JSON for Pipelines](#example-4-machine-readable-json-for-pipelines)
  - [Argument Reference Table](#argument-reference-table)

---

## Preparing Input Data

Cheq consumes structures stored in the standard XYZ text format. Each calculation requires coordinates and element identity for every atom in the system.

### XYZ File Requirements

- Line 1: integer count of atoms.
- Line 2: free-form comment (reported in outputs; can hold metadata such as a system name).
- Remaining lines: one atom per line with either an element symbol or atomic number followed by `x y z` coordinates in ångströms, separated by whitespace.
- Additional columns after `z` are ignored; charges from previous calculations are not reused.
- The parser stops after reading the declared number of atoms. Extra lines are ignored; missing lines trigger a descriptive parse error.

### Using Standard Input

- Pass `-` as the input path to stream XYZ data from standard input.
- Example: `cat water.xyz | cheq - --format json`.
- Progress messages and errors are still written to standard error, so command pipelines remain clean.

---

## Core Charge Equilibration (`cheq`)

### Basic Usage

```sh
cheq INPUT [options]
```

- `INPUT` is the XYZ file to process (or `-` for stdin).
- By default, results print as a table to standard output while a spinner indicates progress.
- Calculations rely on built-in parameter data (`resources/qeq.data.toml`). Use `--params` to override with your own TOML file.

### Argument Reference

Cheq exposes a single command with one positional argument and several option groups. Flags can be combined in any order; long names use `--kebab-case` variants.

#### Positional Argument

- `INPUT` (**required**): Path to an XYZ file or `-` for stdin.

#### Output Options

- `-o, --output <FILE>`: Write results to `FILE` instead of stdout. Directories must exist beforehand.
- `-f, --format <pretty|xyz|csv|json>`: Select output representation. `pretty` renders a Unicode table, `xyz` appends charges as a fifth column, `csv` emits comma-separated rows, and `json` returns structured data. Default: `pretty`.
- `-p, --precision <INT>`: Number of decimal places for floating-point output. Applies to coordinates, charges, and metadata. Default: `6`.

#### Calculation Options

- `-P, --params <FILE>`: Load atomic parameters from a TOML file. When omitted, the bundled QEq set is used.
- `-q, --total-charge <FLOAT>`: Constrain the total molecular charge (in electrons). Default: `0.0`.

#### Solver Options

- `--tolerance <FLOAT>`: Convergence threshold for the RMS change in charges between iterations. Default: `1e-6`.
- `--max-iterations <INT>`: Maximum number of solver iterations before aborting. Default: `100`.
- `--lambda-scale <FLOAT>`: Multiplier applied to the screening length in the Coulomb operator. Default: `0.5`.
- `--hydrogen-scf <BOOL>`: Enable (true) or disable (false) the hydrogen hardness self-consistency update each iteration. Default: `true`.
- `--basis <sto|gto>`: Basis functions to use for Coulomb integrals. `sto` uses exact Slater-Type Orbitals, `gto` uses approximate Gaussian-Type Orbitals. Default: `sto`.
- `--damping <auto|fixed|none>`: Damping strategy for the SCF iteration. `auto` adjusts damping based on convergence, `fixed` uses a constant factor, `none` disables damping. Default: `auto`.
- `--damping-factor <FLOAT>`: Damping factor (0.0-1.0). Used as the fixed value for `fixed` strategy or the initial value for `auto` strategy. Default: `0.4`.

#### General Arguments

- `-h, --help`: Print help information, including all options grouped as shown here.
- `-V, --version`: Display the CLI version.

---

## Custom Parameter Files (`params.toml`)

Cheq parameters map atomic numbers (or symbols) to electronegativity, hardness, covalent radius, and principal quantum number. Provide overrides when experimenting with alternative datasets or extending the element coverage.

### File Structure

- Single top-level table: `[elements]`.
- Keys may be a quoted atomic number (`"8"`) or an element symbol (`"O"`). Mixed styles are allowed.
- Each entry contains four required fields:
  - `chi`: electronegativity (electron volts).
  - `j`: chemical hardness (electron volts).
  - `radius`: covalent radius (ångströms).
  - `n`: principal quantum number (unsigned integer).

### Example Parameter File

```toml
[elements]
"H" = { chi = 2.20, j = 13.60, radius = 0.37, n = 1 }
"O" = { chi = 3.44, j = 13.62, radius = 0.60, n = 2 }
"Na" = { chi = 0.93, j = 5.14, radius = 1.86, n = 3 }
```

Save this as `params.toml` and reference it via `cheq water.xyz --params params.toml`.

### Detailed Parameter Descriptions

- **Electronegativity (`chi`)**: Drives charge transfer direction; higher values attract electron density.
- **Hardness (`j`)**: Resists charge flow; large values dampen partial charges.
- **Radius (`radius`)**: Sets the screening distance for Coulomb interactions.
- **Principal Quantum Number (`n`)**: Distinguishes valence shells with similar radii but different electronic structure.

The solver requires parameters for every element present in `INPUT`. Missing entries raise a descriptive `Calculation error` pointing to the absent atomic number.

---

## Practical Examples (Use Cases)

Assume `examples/water.xyz` contains a three-atom water geometry.

### Example 1: Quick Check with Defaults

```sh
cheq examples/water.xyz
```

> Emits a pretty table to stdout using bundled parameters, zero total charge, and six decimal places.

### Example 2: Structured Output to Disk

```sh
cheq examples/water.xyz \
    --format csv \
    --output results/water.csv \
    --precision 4
```

> Produces a CSV file ready for spreadsheets or plotting tools. Ensure `results/` exists before running the command.

### Example 3: Custom Total Charge and Parameters

```sh
cheq models/cation.xyz \
    --total-charge 1.0 \
    --params data/transition-metals.toml \
    --tolerance 1e-8 \
    --max-iterations 150
```

> Forces a +1 charge state, uses custom parameters, and tightens convergence criteria to improve accuracy for challenging systems.

### Example 4: Machine-Readable JSON for Pipelines

```sh
cat trajectories/frame_050.xyz | \
    cheq - --format json --precision 8 > frame_050.json
```

> Streams an XYZ frame into Cheq, collects JSON output, and can be chained into downstream analysis scripts or databases.

---

## Argument Reference Table

| CLI Argument (Short) | CLI Argument (Long) | Value Type | Default      | Description                                                   |
| :------------------- | :------------------ | :--------- | :----------- | :------------------------------------------------------------ |
| _(positional)_       | `INPUT`             | File Path  | **Required** | XYZ file to process, or `-` for stdin.                        |
| `-o`                 | `--output`          | File Path  | stdout       | Destination for formatted results.                            |
| `-f`                 | `--format`          | Enum       | `pretty`     | Output style: `pretty`, `xyz`, `csv`, or `json`.              |
| `-p`                 | `--precision`       | Integer    | `6`          | Decimal digits for floating-point values.                     |
| `-P`                 | `--params`          | File Path  | bundled set  | Optional TOML parameter file.                                 |
| `-q`                 | `--total-charge`    | Float      | `0.0`        | Target total charge applied during equilibration.             |
| _(none)_             | `--tolerance`       | Float      | `1e-6`       | RMS change threshold for solver convergence.                  |
| _(none)_             | `--max-iterations`  | Integer    | `100`        | Hard cap on solver iterations.                                |
| _(none)_             | `--lambda-scale`    | Float      | `0.5`        | Screening length multiplier in Coulomb term.                  |
| _(none)_             | `--hydrogen-scf`    | Bool       | `true`       | Toggle hydrogen hardness self-consistency.                    |
| _(none)_             | `--basis`           | Enum       | `sto`        | Coulomb integral basis: `sto` (Slater) or `gto` (Gaussian).   |
| _(none)_             | `--damping`         | Enum       | `auto`       | SCF damping strategy: `auto`, `fixed`, or `none`.             |
| _(none)_             | `--damping-factor`  | Float      | `0.4`        | Damping factor (0.0–1.0) for `fixed` or initial `auto` value. |
| `-h`                 | `--help`            | Flag       | (N/A)        | Display contextual help and exit.                             |
| `-V`                 | `--version`         | Flag       | (N/A)        | Show version information and exit.                            |
