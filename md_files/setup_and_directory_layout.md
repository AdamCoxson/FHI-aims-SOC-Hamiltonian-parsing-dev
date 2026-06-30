# FHI-aims SOC Extraction and DeepH Workflow Development

This directory (`FHIaims_SOC_extract_dev`) covers the development, testing, and documentation for extracting noncollinear Spin-Orbit Coupling (SOC) Hamiltonians
from FHI-aims and processing them into Hamiltonians compatible with the DeepH pipeline. 

**This instruction file will cover**:
- Building a modified FHI-aims with SOC Hamiltonian file dumping.
- How to run a simulation workflow to reproduce MoS2 results.
- The Directory Layout.

## Making a modified FHI-aims build with SOC Hamiltonian dumping.
Go to `.../FHIaims_SOC_extract_dev/modified_source_files` and you will see a modified version of `calculate_second_variational_soc.f90` and a `ReadMe.txt` showing the
relevant code block and how to modify the source code.

- Make a copy of you current FHI-aims directory for development, and substitute this file into src/soc.
- Use cmake and make to form a new FHI-aims build with the modified source code.

## Simulation to reproduce MoS2 results.
- Create a relevant environment
Environment dependences: ASE, pyfhiaims, deephdock (python 3.12 or higher), pytest (for unit tests in another file checking hermicity)

```bash
conda create -n fhi-soc-dock
conda activate fhi-soc-dock
conda install -c conda-forge python=3.14.6
conda install -c conda-forge ase
conda install -c conda-forge pyfhiaims
pip install deepx-dock
```

- Go to `.../FHIaims_SOC_extract_dev/2D_materials_data/MoS2_reproduce/run_workflow_unitcell.py` and update the number of cores and relevant paths and directories
at the beginning of main. Run this workflow, it may take some time depending on resources (only about 2 mins for my 32 core AMD-thread ripper)

```bash
python 2D_materials_data/MoS2_reproduce/run_workflow_unitcell.py
```

The soc_intermediate is a single-point scf calculation of the relaxed structures. Here we can see that the DUMP SOC MATRIX block was activated and we now have two
extra files "realspace_soc_matrix.out" and "realspace_soc_unfmttd.out". The Human readable realspace matrix is ~ 4x the size as rs_hamiltonian.out. However
rs_hamiltonian.out has 1 column while realspace_soc_matrix.out has 3 columns of the 3 different Pi matrices.

- For band structure plotting, I installed and used clims
```bash
clims-aimsplot --band --species_dos --total_dos --emin -5 --emax 5 --save_only
clims-aimsplot --band  --emin -5 --emax 5 --save_only
```

- Go to `/.../FHIaims_SOC_extract_dev/DeepH_soc_process/run_aims_soc_to_deeph.py`, change base dir and run.
uncomment `base_dir = current_dir.parent / '2D_materials_data/MoS2_reproduce'` to change base dir

run using:

```bash
python DeepH_soc_process/run_aims_soc_to_deeph.py
```

- This will create a new directory called `.../FHIaims_SOC_extract_dev/2D_materials_data/MoS2_reproduce/deepH_files`


## Directory Layout

The workspace is organized into several key subdirectories:

### 1. `md_files/` (Documentation)
Contains in-depth markdown documentation explaining the physics, methodology, and modifications made throughout this project.
- **modified_deeph_code_for_SOC.md**: Details modifications to DeepH's parser for SOC.
- **compare_parsers_methodology.md**: Explains how the DeepH-based parsing route was mathematically verified against an explicit custom Python parser.
- **extract_realspace_soc.md**: A guide on how the real-space SOC Hamiltonian ($\Pi$ matrix) is identified, dumped from FHI-aims, and processed in Python.
- **full_workflow_summary.md**: A comprehensive overview of the full pipeline (from Fortran source code to HDF5 tensors).

### 2. `modified_source_files/` (FHI-aims Source Modifications)
Contains the specific Fortran files modified within FHI-aims to support SOC extraction.
- `calculate_second_variational_soc.f90`: The patched source file that natively dumps the three spatial components of the SOC operator ($\pi_x, \pi_y, \pi_z$) to `realspace_soc_matrix.out` alongside standard scalar Hamiltonians.
- `ReadMe.txt`: Notes pertaining to the FHI-aims source patches.

### 3. `DeepH_soc_process/` (DeepH Conversion Scripts)
Contains the production-ready python parsers designed to ingest FHI-aims outputs and convert them into the block-sparse format utilized by DeepH.
- `aims_soc_to_deeph.py`: The modified DeepH extraction tool capable of handling the noncollinear $2N \times 2N$ complex SOC Hamiltonian.
- `run_aims_soc_to_deeph.py`: An execution wrapper to run the extraction locally.

### 4. `script_testing/` (Validation and Unit Tests)
A dedicated testing suite verifying that matrix parsing and mathematical reconstruction (including Hermiticity and Parity constraints) are perfectly maintained.
- `compare_parsers.py`: Compares the output of `aims_soc_to_deeph.py` against a custom mathematical script (`build_soc_hamiltonian.py`) to ensure 100% numerical consistency.
- `unit_tests/`: Pytest suite to rigorously validate the parser logic and matrix parity constraints (`symm_signs`).

### 5. `2D_materials_data/` (Automated Test Workflows)
End-to-end Python automation scripts (`run_workflow_unitcell.py`) and associated data evaluating SOC properties on canonical 2D materials. These scripts:
- Perform structural relaxations (e.g., `relaxation_light`).
- Evaluate Spin-Orbit Coupling single points (e.g., `soc_intermediate`).
- Trigger the custom Fortran dump to produce SOC matrices for the DeepH parser.

See the two directories for graphene and MoS2. I used clims to plot the band structures and obviously MoS2 has non-negligible SOC effects. These tests are only a unit cell so it runs quite quickly.

## Summary

This package provides everything needed to reproduce, validate, and utilize the FHI-aims to DeepH SOC data pipeline. For an in-depth understanding of the mathematical operations involved, refer to the documentation inside `md_files/`.
