# Spectral Analysis of Numerical Schemes

This repository contains the developments and numerical experiments associated with a study of the **spectral properties of numerical schemes for the linear convection equation**.

The objective is to investigate how spatial and temporal discretization affect the numerical solution, with particular emphasis on **numerical dispersion, dissipation, phase velocity, and group velocity**. The repository is intended as a research and exploration environment, combining analytical developments with Python-based numerical experiments and visualizations. The main focus is the **excellent Lax-Wendroff scheme**.

## Repository structure

```text
.
├── src/              # Python source code and numerical experiments
├── notebooks/        # Exploratory notebooks
├── doc/              # LaTeX report and documentation
│   └── Images/       # Figures used in the report
├── pyproject.toml    # Project dependencies and configuration
└── uv.lock           # Locked Python environment
```

## Installation

The project uses [uv](https://docs.astral.sh/uv/) for Python environment and dependency management.

Clone the repository and synchronize the environment:

```bash
git clone git@github.com:nicolamarque/<repository-name>.git
cd <repository-name>
uv sync
```

The project environment is created automatically in `.venv/`.

## Usage

Python scripts can be run with:

```bash
uv run python src/<script>.py
```

Jupyter notebooks can be launched from the project environment with:

```bash
uv run jupyter lab
```

The generated figures are used in the LaTeX report located in `doc/`.

## Report

The report documents the analytical and numerical developments carried out in the project.

The current version is **v1.1**.

The report source files are available in `doc/`, together with the figures used to produce the document.

## Status

This is an ongoing personal research project. The repository currently focuses on the spectral analysis of classical numerical schemes and will progressively include additional schemes, numerical experiments, and extensions of the analysis... if I find some time to do this!

## License

This project is licensed under the MIT License. See the LICENSE file for details.