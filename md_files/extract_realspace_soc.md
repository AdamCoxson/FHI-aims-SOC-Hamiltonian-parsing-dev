# Extraction of the Real-Space SOC Hamiltonian ($\Pi$ Matrix)

The real-space SOC Hamiltonian is calculated before being transformed into the state-space basis. This matrix is referred to as the $\Pi$ matrix in the literature, and in the FHI-aims source code, it is represented by the variable `soc_matrix`. This document details how this matrix was extracted and dumped to a file for later post-processing (e.g., using Python) as part of our ML Hamiltonian workflow.

## 1. Extraction Location

The $\Pi$ matrix is allocated and constructed in the main SOC driver routine in the FHI-aims codebase:
**File:** `soc_source/calculate_second_variational_soc.f90`

The extraction logic was added directly after the following subroutine call (around line 664) where the `soc_matrix` is populated:
```fortran
  ! Fill in the matrix for the real-space SOC matrix elements (\Pi from the
  ! Phys Rev Mater. paper)
  call get_timestamps(time_matrix, clock_time_matrix)
  ! Step 3 in Section III.3 from Huhn and Blum, Phys. Rev. Mater. (2017)

  call integrate_soc_matrix (rho,hartree_potential,partition_tab,soc_matrix)
  *** ADDED HERE ***
```

Directly after this call, the `soc_matrix` variable contains the compressed real-space SOC Hamiltonian. The dump was placed here before any subsequent rotation (e.g., by `rotate_Pi_columns` if `calculate_mae` is true).

## 2. Format of `soc_matrix`

`soc_matrix` is a 2D array allocated as:
`real*8, dimension(ld_soc_matrix, 3) :: soc_matrix`

- **First dimension (`ld_soc_matrix`)**: Represents the packed pairs of basis functions. Because the Hamiltonian matrix is highly sparse in real-space, FHI-aims packs the non-zero matrix elements into a contiguous 1D array per operator. `ld_soc_matrix` is either `n_hamiltonian_matrix_size` (full size) or `batch_perm(n_bp_integ)%n_local_matrix_size` (if load balancing is used).
- **Second dimension (`3`)**: Represents the three spatial coordinates ($x, y, z$) of the spin-orbit coupling operator.

## 3. Implementation of the Matrix Dump

To dump this matrix, an unformatted (binary) write statement was injected to handle the large matrix efficiently. 

The write statement was restricted to a single MPI task (`myid == 0`), assuming runs without `use_local_index` (where every process has the full array). The following code block was added to the FHI-aims source right after `call integrate_soc_matrix(...)`:

```fortran
  ! --- DUMP SOC MATRIX ---
  ! --- MODIFICATION BY Adam Coxson
  if (myid == 0) then
    ! Unformatted dump (smaller file, faster read in Python)
    open(unit=999, file="realspace_soc_matrix_unfmttd.out", form='unformatted', status='replace')
    ! Write the dimensions first so the reader script knows the size
    write(999) ld_soc_matrix
    write(999) 3
    ! Write the actual matrix block
    write(999) soc_matrix(1:ld_soc_matrix, 1:3)
    close(999)

    ! Alternative: Formatted dump (easier to debug, but larger file)
    open(unit=998, file="realspace_soc_matrix.out", status='replace')
    write(998, *) ld_soc_matrix, 3
    do i = 1, ld_soc_matrix
      write(998, '(3E24.15)') soc_matrix(i, 1:3)
      end do
    close(998)
  end if
  ! -----------------------
```

## 4. Reading the Matrix in Python

The unformatted (binary) dump is loaded into Python using `scipy.io.FortranFile` during post-processing. Because it is dumped as a compressed vector of overlapping basis pairs, the sparsity pattern array is required to map it back to basis pairs. The formatted dump is also performed and this file is utilised in the down stream tasks in the example workflows.

The following Python snippet demonstrates how the binary data is read in to python:

```python
import numpy as np
from scipy.io import FortranFile

# Read the dumped binary file
f = FortranFile('realspace_soc_matrix_unfmttd.out', 'r')

# FHI-aims wrote unformatted records. The dimensions are read first:
ld_soc_matrix = f.read_ints(dtype=np.int32)[0] # read_ints() advances pointer to next Fortran record
num_coords = f.read_ints(dtype=np.int32)[0]

# Then the actual array is read:
soc_matrix_packed = f.read_reals(dtype=np.float64)
# Reshape taking Fortran column-major order into account
soc_matrix_packed = soc_matrix_packed.reshape((num_coords, ld_soc_matrix)).T
f.close()

print(f"Loaded packed SOC matrix of shape: {soc_matrix_packed.shape}")
```

Alternatively, the formatted dump (`realspace_soc_matrix.out`) can be read directly into Python using `numpy`:

```python
import numpy as np

# Read the SOC matrix, skipping the first line (dimensions)
soc_data = np.loadtxt('realspace_soc_matrix.out', skiprows=1)
Pi_x = soc_data[:, 0]
Pi_y = soc_data[:, 1]
Pi_z = soc_data[:, 2]

print(f"Loaded formatted SOC matrix of shape: {soc_data.shape}")
```

## 5. Integration and Hamiltonian Construction

The parsing of the dumped `soc_matrix` can be performed using two methods within this workflow: via a modified **DeepH Dock** or using a standalone Python script (`build_soc_hamiltonian.py`).

### Method A: DeepH Dock Integration

The parsing of the dumped `soc_matrix` into **DeepH Dock** utilises the following formats:

- **Sparsity Pattern:** The `soc_matrix` shares the exact same sparsity pattern and indexing format as the scalar-relativistic Hamiltonian.
- **Index File:** DeepH Dock natively parses `rs_indices.out` to map the 1D compressed vector to the basis pairs. Because the indexing is identical, the standard `rs_indices.out` generated by FHI-aims (via `output_rs_matrices plain/hdf5`) was used to unpack the `soc_matrix_packed` array.
- **Channel Count:** The extracted `soc_matrix` contains **3 channels** representing the $x, y, z$ spatial components of the SOC operator $\vec{\Pi}$.

To integrate this with the DeepH conversion scripts, the 3-channel numpy array (`soc_matrix_packed`) was loaded, and the same CSR/COO matrix construction routine used for `rs_hamiltonian.out` was applied, using the indices directly from `rs_indices.out`.

This relies on modifying and overwriting functions from DeepH, which can be found in `aims_soc_to_deeph.py`.

### Method B: Standalone Construction Script

Alternatively, the formatted output can be processed using the standalone script `script_testing/build_soc_hamiltonian.py`. This script constructs the full Spin-Orbit Coupled Hamiltonian in the Localized Atomic Orbital (LAO) basis directly from the real-space output matrices.

The script operates by:
1. Parsing the `rs_indices.out` file to retrieve the sparsity pattern and mapping logic.
2. Loading the scalar-relativistic Hamiltonian (`rs_hamiltonian.out`) and the formatted SOC matrix (`realspace_soc_matrix.out`).
3. Constructing the final complex spin-dependent Hamiltonian matrix in a spin-blocked format by mapping the spatial operators and verifying its Hermiticity.

You can run this script directly to yield a finalized `H_soc_rs.npy` Hamiltonian for further analysis:
```bash
python script_testing/build_soc_hamiltonian.py \
  --dir path/to/aims/output \
  --indices rs_indices.out \
  --hamiltonian rs_hamiltonian.out \
  --soc realspace_soc_matrix.out \
  --output H_soc_rs.npy
```

**Validation:**
Unit testing using `pytest` has confirmed that when the 3-channel matrix is packed alongside the scalar field with properly mapped anti-symmetric parity (`symm_signs=[1, 1, -1, -1, -1]`), the resulting complex $2N \times 2N$ Hamiltonian satisfies both Spatial Hermiticity ($H(\mathbf{R}) = H(-\mathbf{R})^\dagger$) and Time-Reversal Symmetry. The finalized outputs are dumped as `soc-hamiltonian.h5` and `hamiltonian.h5` for immediate use.
