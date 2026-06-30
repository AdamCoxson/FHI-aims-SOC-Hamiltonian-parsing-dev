# Full Workflow: Extracting and Processing FHI-aims SOC Hamiltonians for DeepH

This document serves as a comprehensive summary of the entire pipeline developed to extract noncollinear Spin-Orbit Coupling (SOC) Hamiltonians from FHI-aims and seamlessly process them into the DeepH machine learning ecosystem.

---

## 1. Modification of the FHI-aims Fortran Source Code

In standard periodic Density Functional Theory (DFT) calculations using FHI-aims, Spin-Orbit Coupling is typically applied perturbatively or self-consistently onto the localized atomic orbital (LAO) basis. To utilize these matrices for machine learning, we needed raw access to the underlying spatial real-space components of the SOC operator before the generalized eigenvalue solver diagonalizes the system. Consequently, we modified the FHI-aims Fortran source code to actively dump these operator matrices during the real-space integration phase. 

Specifically, alongside the standard real-space scalar Hamiltonian ($H_{scalar}$) written to `rs_hamiltonian.out`, the modified code now dumps the three spatial components of the SOC operator:
- $\pi_x$ (the $x$-component)
- $\pi_y$ (the $y$-component)
- $\pi_z$ (the $z$-component)

which are written to a new output file called `realspace_soc_matrix.out`.

## 2. Automating the Test Case (`run_workflow_unitcell.py`)

**Note, you must edit the paths and directories in `run_workflow_unitcell.py` to suit your system.**

To validate the modified Fortran code, a Python automation script (`2D_materials_data/MoS2_unitcell/run_workflow_unitcell.py`) orchestrates the end-to-end DFT calculation.

1. **Structure Generation**: It utilizes ASE (Atomic Simulation Environment) and PyFHIaims to generate a 2D layer object of Molybdenum Disulfide (MoS2).
2. **Slab Relaxation**: It executes a geometric structural relaxation using `light` species defaults and Van der Waals corrections to better capture inter-layer interactions.
3. **Spin-Orbit Single Point Evaluation**: It then initiates a follow-up single-point calculation using `tight` species defaults with the Spin-Orbit Coupling (`include_spin_orbit`) keyword explicitly enabled.
4. **Data Extraction**: During this final noncollinear subroutine, the hardcoded Fortran block we introduced is actively triggered, prompting FHI-aims to physically dump the relevant SOC and scalar matrices to the working directory.


## 3. Dumping the Real-Space Hamiltonians

Once a calculation concludes, the FHI-aims working directory contains the necessary raw outputs:
- **`rs_indices.out`**: Defines the sparse non-zero elements, mapping them to specific basis functions (orbitals) and spatial displacement vectors ($\mathbf{R}$).
- **`rs_hamiltonian.out`**: Contains the $N \times N$ scalar Hamiltonian ($H_{scalar}$) evaluating the core electrostatic, kinetic, and non-relativistic interaction terms.
- **`realspace_soc_matrix.out`**: Contains the decoupled $\pi_x, \pi_y, \pi_z$ matrices corresponding to the SOC angular momentum operators, 3 columns of x, y, z pauli pi matrices.
- **`rs_overlap.out`**: Contains the $N \times N$ spatial overlap matrix ($S$).

Because these are evaluated on a strictly localized, non-orthogonal atomic orbital basis, the matrices reflect the raw generalized eigenvalue problem $HC = SCE$. 

## 4. DeepH Based Processing (`aims_soc_to_deeph.py`)

DeepH requires the Hamiltonian to be formatted as a complex, $2N \times 2N$ block-sparse tensor stored in chunked `.h5` files. To make this compatible we copied and adapted the DeepH data translator, `aims_to_deeph.py`, to form `aims_soc_to_deeph.py`, which now acts as a standalone parser for FHI-aims SOC Hamiltonians.

### Mathematical Reconstruction of the SOC Tensor
The parser reads the four scalar fields ($H_{scalar}, \pi_x, \pi_y, \pi_z$) mapping them across the sparsity indices defined in `rs_indices.out`. 

It then combines them into the full Noncollinear $2N \times 2N$ complex block matrix using the standard Pauli matrix formulations. For any orbital pair $(i, j)$ at displacement $\mathbf{R}$, the $2 \times 2$ spin-subblock is constructed as:

$$
H_{i,j}(\mathbf{R}) = 
\begin{pmatrix}
H_{\uparrow\uparrow} & H_{\uparrow\downarrow} \\
H_{\downarrow\uparrow} & H_{\downarrow\downarrow}
\end{pmatrix}
$$

Where the individual tensor elements are reconstructed as:
$$ H_{\uparrow\uparrow} = H_{scalar} - i \pi_z $$
$$ H_{\downarrow\downarrow} = H_{scalar} + i \pi_z $$
$$ H_{\uparrow\downarrow} = -i \pi_x - \pi_y $$
$$ H_{\downarrow\uparrow} = -i \pi_x + \pi_y $$

### Resolving Spatial Hermiticity
A major challenge arose regarding DeepH's native spatial Hermiticity assumptions: $H(\mathbf{R}) = H(-\mathbf{R})^\dagger$. 

When equating the real and imaginary parts of the $2 \times 2$ block tensor, it dictates strict physical rules for the raw outputs: the scalar Hamiltonian must be geometrically symmetric, while the spatial SOC Pauli matrices ($\pi_x, \pi_y, \pi_z$) must be geometrically anti-symmetric.

We resolved Hermiticity conflicts by configuring the parser with `symm_signs=[1, 1, -1, -1, -1]`, forcing DeepH to accurately interpret the spatial SOC matrices with antisymmetric parity. Additionally, because the SOC matrices are integrated numerically over a real-space grid in FHI-aims, discrete grid artifacts lead to some noise ($\sim 10^{-6}$). By including an adjustable DeepH Hermiticity `tolerance` argument (default `1e-3`), the pipeline can account for the material-dependent grid noise and validate the parsed data.

### Dual output of SOC and non-SOC Hamiltonians
Finally, to support diverse model training objectives, the parser was configured to simultaneously output two formats:
1. **`soc-hamiltonian.h5`**: The full $2N \times 2N$ complex SOC Hamiltonian representing the fully reconstructed equations above.
2. **`hamiltonian.h5`**: The purely non-SOC baseline, generated by feeding $H_{scalar}$ directly into the diagonal spin blocks as a real matrix.
Note, these are named `hamiltonian.h5` and `hamiltonian0.h5` in the original DeepH data translator.


