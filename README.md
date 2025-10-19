# CMOD - Chromatic Surface Brightness MODulation Analysis Pipeline

[![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/victoralonsorodriguez/CMOD.svg)](https://github.com/victoralonsorodriguez/CMOD/tags)[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Analysis pipeline to study the CMOD effect in nearby galaxies, based on surface brightness profile decomposition using [Galfit](https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html)/[Imfit](https://www.mpe.mpg.de/~erwin/code/imfit/).

**Status:** Code refactored (v1.0.0). Ready for implementing new features (e.g., resampling).

---

## Installation

1. **Clone the repository:**

```bash
    git clone [https://github.com/victoralonsorodriguez/CMOD.git](https://github.com/victoralonsorodriguez/CMOD.git)
    cd CMOD
```

2. **Create the Conda environment:** Requires [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

```bash
    conda env create -f environment.yml
    conda activate cmod_env
```

This will install all Python dependencies and the `cmod` package in editable mode.

3. **External Dependencies:** This pipeline requires **Galfit** and/or **Imfit** to be installed on your system and available in the system's `PATH`. Follow the installation instructions from their respective websites.

---

## Usage

The main pipeline is executed via the `scripts/run_cmod.py` script.

```bash
conda activate cmod_env
python scripts/run_cmod.py [OPTIONS]
```

## Configuration

The pipeline's behavior (which galaxies to analyze, which program to use, initial parameters, etc.) is primarily controlled through a configuration file. Pass the path using the `-conf` flag:

```bash
python scripts/run_cmod.py -conf config/your_config.txt
```

(We need to define and document the format of this configuration file).

Specific parameters can be overridden via command-line flags. Run `python scripts/run_cmod.py --help` to see available options (to be implemented in `io.py`).

## Project Structure

* `src/cmod/`: Main source code for the Python package.
    * `io.py`: Data and configuration reading/writing.
    * `processing.py`: Image processing (PSF, isophotes, initial conditions).
    * `fitting.py`: Script generation and execution for Galfit/Imfit.
    * `photometry.py`: Photometric conversions.
    * `cosmology.py`: Cosmological conversions (scales, z-lambda).
    * `plotting.py`: Functions for generating plots.
    * `resampling.py`: (Future) Functions for observation simulation.
    * `pipeline.py`: Functions orchestrating the analysis pipeline.
    * `utils.py`: Helper functions.
* `scripts/`: Executable scripts.
    * `run_cmod.py`: Main entry point.
* `data/`: Input data (original FITS, PSFs). *(Note: Add to .gitignore if large)*.
* `config/`: Example/template configuration files.
* `notebooks/`: (Optional) Jupyter notebooks for exploration or visualization.
* `legacy_scripts/`: Archived old scripts (.sh, .py).

---

## Contact

* Name: Victor Alonso Rodriguez
* email: victoralonsorodriguez61@gmail.com


