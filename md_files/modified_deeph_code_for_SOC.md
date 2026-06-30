# DeepH Hamiltonian Processing Overview

This document provides a high-level overview of how DeepH processes Hamiltonians from FHI-aims outputs, and details the specific modifications made to the original parser to support Noncollinear Spin-Orbit Coupling (SOC).

## How DeepH Processes Hamiltonians

DeepH uses a data processing pipeline to map real-space data from density functional theory (DFT) codes into its sparse HDF5 format. The workflow for parsing Hamiltonians follows these general steps:

1. **High-Level Abstraction**: The entry point is a translator class (e.g., `PeriodicAimsDataTranslator`) which controls the data conversion via a master function (e.g., `transfer_one_aims_to_deeph`).
2. **Raw Matrix Ingestion**: FHI-aims outputs scalar matrices ($H_{scalar}$), and in the SOC case, the Pauli matrices ($\pi_x, \pi_y, \pi_z$) as raw 1D array lists.
3. **Assembly and Mapping**: DeepH relies on mapping utilities (such as `_trans_mxs_to_R_dict` and `_trans_R_dict_to_entries`). These functions loop over the FHI-aims sparsity index file (`rs_indices.out`), apply spherical harmonic parity factors, perform unit conversions (`HARTREE_TO_EV`), and assert spatial Hermiticity constraints.
4. **Tensor Recombination**: The reconstructed blocks are passed to concatenation routines to build the final block-sparse tensors. For collinear spin, they are placed diagonally. For noncollinear SOC, they are combined into complex $2 \times 2$ spin-orbital sub-blocks.
5. **Data Export**: The final tensors, alongside auxiliary structures like `info.json`, `POSCAR`, and overlaps, are written to sparse `.h5` formats suitable for machine learning workflows.

Our customized noncollinear SOC parser (`aims_soc_to_deeph.py`) builds upon this pipeline by introducing explicit Pauli tensor algebra while preserving the internal mapping logic of DeepH.

## Change Log: `aims_soc_to_deeph.py` vs `aims_to_deeph.py`

The following modifications were made to the original DeepH dock parser (`aims_to_deeph.py`) to create the SOC variant (`aims_soc_to_deeph.py`), ensuring the $2N \times 2N$ complex SOC Hamiltonian is extracted and reconstructed.

### 1. Control Parameter Parsing (`_read_ctrl_in`)
- **Original (`aims_to_deeph.py`)**: Only detected collinear spin settings via `spin collinear`.
- **Modified (`aims_soc_to_deeph.py`)**: Added detection for `include_spin_orbit` in `control.in`. When present, it correctly flags the data as SOC-enabled (`spin_orbit = True`).

### 2. Reading SOC Pauli Matrices (`_read_soc_matrices`)
- **Original**: Only read the scalar `rs_hamiltonian.out` (and its spin up/down counterparts).
- **Modified**: Introduced the `_read_soc_matrices` function to parse `realspace_soc_matrix.out`. FHI-aims dumps the three spatial $\pi_x, \pi_y, \pi_z$ matrices inside this file. The script parses the array and returns the individual matrix components, correctly applying the unit conversion (`HARTREE_TO_EV`).

### 3. Data Structure Packing (`analysis_data`)
- **Original**: Packed only the `Ovlp` and `H_scalar` matrices into the translation processing list (`self.mx_lst`).
- **Modified**: When `spin_orbit` is True, it packs `Ovlp`, `H_scalar`, and all three Pauli matrices (`pi_x, pi_y, pi_z`) into `self.mx_lst`.

### 4. Spin-Orbit Tensor Recombination (`_soc_concatenate`)
- **Original**: Only featured `_spin_concatenate` which placed $H_{up}$ and $H_{dn}$ on the diagonal of a block matrix for collinear calculations.
- **Modified**: Added the complex concatenation routine `_soc_concatenate`. After the spatial matrices are mapped by the DeepH parsing algorithms, this routine combines them into the full $2N \times 2N$ complex block-sparse array using the following mathematical formulation:
  - $H_{\uparrow\uparrow} = H_{scalar} - i \pi_z$
  - $H_{\downarrow\downarrow} = H_{scalar} + i \pi_z$
  - $H_{\uparrow\downarrow} = -i \pi_x - \pi_y$
  - $H_{\downarrow\uparrow} = -i \pi_x + \pi_y$

### 5. Ensuring Spatial Hermiticity (`symm_signs`)
- **Original**: DeepH executes a spatial Hermiticity check inside `_trans_mxs_to_R_dict` that mirrors matrices across $\mathbf{R}$ and $-\mathbf{R}$. For collinear scalar outputs, this works via symmetric pairing.
- **Modified**: To construct a valid complex $2N \times 2N$ Hamiltonian, the spatial components must obey explicit parity constraints: the scalar matrices are symmetric, but the SOC Pauli matrices ($\pi_x, \pi_y, \pi_z$) are anti-symmetric. We enforce this by setting `symm_signs = [1, 1, -1, -1, -1]`. Furthermore, because FHI-aims computes the SOC real-space integrations numerically on a grid, a rigid Hermiticity threshold can fail due to grid noise. We exposed an adjustable internal DeepH validation `tolerance` (default `1e-3`) to map the SOC matrices while maintaining Hermiticity validations and accounting for grid noise across different materials.

### 6. Dual Export Naming Convention
- **Original**: For collinear models, `export_H0` looked for separate files like `rs_hamiltonian0.out` to isolate the base Hamiltonian.
- **Modified**: To support various training workflows, the parsing pipeline generates two distinct files. It outputs the complex SOC tensor as `soc-hamiltonian.h5` and the purely scalar base Hamiltonian (mapped diagonally into a $2N \times 2N$ real block) as `hamiltonian.h5`. This allows extraction of both the SOC and non-SOC representations in a single run.

### 7. Package Portability Adjustments
- **Modified**: To allow `aims_soc_to_deeph.py` to act as an independent, locally executed wrapper script overriding the installed package, the `transfer_all_aims_to_deeph` function and the import for `parallel_map` were removed. This prevents conflicts with the installed `deepx_dock` environment, ensuring that the local file pulls all its other dependencies (like constants and math utilities) directly from the active Conda package.
