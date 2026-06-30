# Methodology: Two Routes to DeepH Sparse HDF5

The `compare_parsers.py` script validates our noncollinear spin-orbit coupling (SOC) extraction process by constructing the identical $2N \times 2N$ complex block-sparse Hamiltonian using two independent scripts: the modified DeepH Dock functions and standalone `build_soc_hamiltonian.py` script.

Despite their different architectures, both routes yield $100\%$ identical numerical outputs, proving that our matrix formulations are accurate and compatible with the `deepx_dock` environment.

---

## Route 1: "The DeepH Way" (`run_deeph_way`)

This route utilizes the native, object-oriented ecosystem of the DeepH data translation module.

### How it works:
1. **High-Level Abstraction**: It instantiates the `PeriodicAimsDataTranslator` from `aims_soc_to_deeph.py` and triggers `transfer_one_aims_to_deeph`.
2. **Internal DeepH Mapping**: It relies heavily on deeply nested `deepx_dock` utilities (like `_trans_mxs_to_R_dict` and `_trans_R_dict_to_entries`).
3. **Array Packaging**: The FHI-aims outputs ($H_{scalar}, \pi_x, \pi_y, \pi_z$) are passed as a raw list of 1D arrays into the mapping utility. 
4. **Automated Assembly**: The deep-level functions automatically loop over the FHI-aims sparsity indexing (`rs_indices.out`), apply spherical harmonic parity factors (`phase_factor`), apply the unit conversions (`HARTREE_TO_EV`), and historically execute spatial Hermiticity assertions.
5. **Tensor Recombination**: The recombined blocks are eventually passed to `_soc_concatenate` and `_spin_concatenate` to build the complex tensors, which are then natively dumped to `.h5` files.

### Pros:
- Deeply integrated into the standard machine-learning pipeline.
- Automatically handles auxiliary file creation (`info.json`, `POSCAR`, `overlap.h5`).

---

## Route 2: "The Alternative Way" (`custom_alternative_way`)

This route bypasses DeepH's native Hamiltonian mapping ecosystem entirely, replacing it with explicit, transparent, custom Python loop logic to verify the mathematical formulations independently.

### How it works:
1. **Hybrid Foundation**: It uses DeepH strictly to construct auxiliary non-Hamiltonian files (like `POSCAR` and `info.json`), then handles the Hamiltonian autonomously.
2. **Raw Matrix Ingestion**: It reads `rs_hamiltonian.out` and `realspace_soc_matrix.out` independently and applies the `HARTREE_TO_EV` scalar directly to the flat arrays.
3. **Explicit Spatial Block Construction**: 
   - Instead of using DeepH's recursive abstractions, it executes a manual `for` loop over `n_cells` ($\mathbf{R}$) and `n_basis` (orbitals).
   - It builds a dictionary `mx_R` where keys are the spatial vectors and atom pairs `(R_1, R_2, R_3, i_atom, j_atom)`.
   - The matrix elements are actively sorted and placed into multi-dimensional NumPy arrays representing specific atomic sub-blocks.
4. **Manual Pauli Tensor Application**:
   - Once the spatial blocks are assembled, it explicitly reconstructs the $2 \times 2$ spin-orbital sub-blocks mathematically:
     - $H_{\uparrow\uparrow} = H_{scalar} - 1j \cdot \pi_z$
     - $H_{\downarrow\downarrow} = H_{scalar} + 1j \cdot \pi_z$
     - $H_{\uparrow\downarrow} = -1j \cdot \pi_x - \pi_y$
     - $H_{\downarrow\uparrow} = -1j \cdot \pi_x + \pi_y$
5. **Sparse Packaging**: It manually computes the `chunk_boundaries`, `chunk_shapes`, and `atom_pairs` required by DeepH and uses the standard `h5py` library to dump the arrays to `soc-hamiltonian.h5` and `hamiltonian.h5`.

### Pros:
- Mathematically transparent and highly readable.
- Easy to debug and verify because the matrix tensor math ($H_{\uparrow\downarrow}$, etc.) is applied to fully assembled $N \times N$ block matrices rather than scattered 1D arrays.
- Proves conclusively that the native DeepH route is interpreting the SOC data correctly.
