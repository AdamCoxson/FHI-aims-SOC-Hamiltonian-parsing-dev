#!/usr/bin/env python3
"""
Build Spin-Orbit Coupled Hamiltonian in the Localized Atomic Orbital (LAO) basis
from FHI-aims real-space output matrices.
"""

import os
import argparse
import numpy as np

def parse_indices(indices_path):
    print(f"Parsing index file: {indices_path}")
    
    n_ham = None
    n_cells = None
    n_basis = None
    cell_index = []
    idx_ham_1 = []
    idx_ham_2 = []
    column_idx = []
    
    with open(indices_path, 'r') as f:
        lines = f.readlines()
        
    state = "HEADER"
    for line in lines:
        line_str = line.strip()
        if not line_str:
            continue
            
        if state == "HEADER":
            if line_str.startswith('n_hamiltonian_matrix_size:'):
                n_ham = int(line_str.split()[-1])
            elif line_str.startswith('n_cells_in_hamiltonian:'):
                n_cells = int(line_str.split()[-1])
            elif line_str.startswith('n_basis:'):
                n_basis = int(line_str.split()[-1])
            elif line_str == 'cell_index':
                state = "CELL_INDEX"
        elif state == "CELL_INDEX":
            if line_str == 'index_hamiltonian(1,:,:)':
                state = "HAM_1"
            else:
                cell_index.append([int(x) for x in line_str.split()])
        elif state == "WAIT_HAM_1":
            if line_str.startswith('index_hamiltonian(1,'):
                state = "HAM_1"
        elif state == "HAM_1":
            if line_str.startswith('index_hamiltonian(2,'):
                state = "HAM_2"
            else:
                idx_ham_1.extend([int(x) for x in line_str.split()])
        elif state == "HAM_2":
            if line_str == 'column_index_hamiltonian':
                state = "COLUMN_INDEX"
            else:
                idx_ham_2.extend([int(x) for x in line_str.split()])
        elif state == "COLUMN_INDEX":
            column_idx.extend([int(x) for x in line_str.split()])
            
    cell_index = np.array(cell_index).reshape(-1, 3)
    idx_ham_1 = np.array(idx_ham_1).reshape(n_cells, n_basis)
    idx_ham_2 = np.array(idx_ham_2).reshape(n_cells, n_basis)
    column_idx = np.array(column_idx)
    
    assert len(cell_index) == n_cells, f"Expected {n_cells} cells, got {len(cell_index)}"
    assert len(column_idx) == n_ham, f"Expected {n_ham} column indices, got {len(column_idx)}"
    
    print(f"Successfully parsed indices:")
    print(f"  n_basis: {n_basis}")
    print(f"  n_cells: {n_cells}")
    print(f"  n_hamiltonian_matrix_size: {n_ham}")
    
    return n_ham, n_cells, n_basis, cell_index, idx_ham_1, idx_ham_2, column_idx

def main():
    parser = argparse.ArgumentParser(description="Construct real-space SOC Hamiltonian from FHI-aims output.")
    current_dir = os.path.dirname(os.path.abspath(__file__))
    default_dir = os.path.join(current_dir, "..", "2D_materials_data", "MoS2_unitcell", "soc_intermediate")
    parser.add_argument("--dir", type=str, default=default_dir, help="Directory containing aims real-space output files")
    parser.add_argument("--indices", type=str, default="rs_indices.out", help="Path to rs_indices.out")
    parser.add_argument("--hamiltonian", type=str, default="rs_hamiltonian.out", help="Path to rs_hamiltonian.out")
    parser.add_argument("--soc", type=str, default="realspace_soc_matrix.out", help="Path to realspace_soc_matrix.out")
    parser.add_argument("--output", type=str, default="H_soc_rs.npy", help="Path to save the constructed Hamiltonian (.npy)")
    args = parser.parse_args()
    
    if args.dir:
        indices_path = os.path.join(args.dir, args.indices)
        ham_path = os.path.join(args.dir, args.hamiltonian)
        soc_path = os.path.join(args.dir, args.soc)
    else:
        indices_path = args.indices
        ham_path = args.hamiltonian
        soc_path = args.soc
        
    # Check that paths exist
    for p in [indices_path, ham_path, soc_path]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Required file not found: {p}")
            
    # Parse indices
    n_ham, n_cells, n_basis, cell_index, idx_ham_1, idx_ham_2, column_idx = parse_indices(indices_path)
    
    # Read scalar Hamiltonian
    print(f"Reading scalar Hamiltonian: {ham_path}")
    H_scalar = np.loadtxt(ham_path)
    assert len(H_scalar) == n_ham, f"Expected H_scalar size {n_ham}, got {len(H_scalar)}"
    
    # Read SOC matrix
    print(f"Reading SOC matrix: {soc_path}")
    with open(soc_path, 'r') as f:
        first_line = f.readline().strip().split()
        ld_soc_matrix = int(first_line[0])
        n_cols = int(first_line[1])
        assert ld_soc_matrix == n_ham, f"SOC matrix size {ld_soc_matrix} does not match n_hamiltonian_matrix_size {n_ham}"
        assert n_cols == 3, f"SOC matrix should have 3 columns, got {n_cols}"
        
    soc_data = np.loadtxt(soc_path, skiprows=1)
    assert soc_data.shape == (n_ham, 3), f"Expected SOC data shape {(n_ham, 3)}, got {soc_data.shape}"
    Pi_x = soc_data[:, 0]
    Pi_y = soc_data[:, 1]
    Pi_z = soc_data[:, 2]
    
    # Identify the central cell index (0, 0, 0)
    cell_idx_0 = None
    for idx, c in enumerate(cell_index):
        if np.all(c == 0):
            cell_idx_0 = idx
            break
            
    if cell_idx_0 is None:
        raise ValueError("Could not find the central cell (0, 0, 0) in cell_index.")
    print(f"Central cell index (0,0,0) is at position: {cell_idx_0} (1-based index: {cell_idx_0 + 1})")
    
    # Construct the full Hamiltonian matrix
    # Spin Blocked format (Option A)
    # Indices 0 to n_basis-1: spin-up
    # Indices n_basis to 2*n_basis-1: spin-down
    H_final = np.zeros((2 * n_basis, 2 * n_basis), dtype=complex)
    
    print("Constructing spin-dependent Hamiltonian matrix...")
    for i_basis in range(n_basis):
        start_idx = idx_ham_1[cell_idx_0, i_basis] - 1
        end_idx = idx_ham_2[cell_idx_0, i_basis] - 1
        
        # If start_idx + 1 == 0, it means no elements for this cell/basis function
        if start_idx + 1 == 0:
            continue
            
        for k in range(start_idx, end_idx + 1):
            j_basis = column_idx[k] - 1
            
            h_scal = H_scalar[k]
            pi_x = Pi_x[k]
            pi_y = Pi_y[k]
            pi_z = Pi_z[k]
            
            # H_block = [ [h - i*pi_z,  -i*pi_x - pi_y],
            #             [-i*pi_x + pi_y, h + i*pi_z] ]
            H_up_up = h_scal - 1j * pi_z
            H_dn_dn = h_scal + 1j * pi_z
            H_up_dn = -1j * pi_x - pi_y
            H_dn_up = -1j * pi_x + pi_y
            
            # Set elements in Blocked format
            # Row i_basis, Col j_basis
            H_final[i_basis, j_basis] = H_up_up
            H_final[i_basis + n_basis, j_basis + n_basis] = H_dn_dn
            H_final[i_basis, j_basis + n_basis] = H_up_dn
            H_final[i_basis + n_basis, j_basis] = H_dn_up
            
            # If off-diagonal in spatial indices, mirror using Hermitian conjugate
            if j_basis < i_basis:
                H_final[j_basis, i_basis] = H_up_up.conjugate()
                H_final[j_basis + n_basis, i_basis + n_basis] = H_dn_dn.conjugate()
                H_final[j_basis, i_basis + n_basis] = H_dn_up.conjugate()
                H_final[j_basis + n_basis, i_basis] = H_up_dn.conjugate()
            elif j_basis > i_basis:
                # This shouldn't normally happen for i_cell=1 (only lower triangular stored)
                print(f"Warning: encountered j_basis ({j_basis}) > i_basis ({i_basis}) in sparse packing.")
                
    # Verify Hermiticity
    print("Verifying Hermiticity of the constructed Hamiltonian...")
    hermitian_diff = np.max(np.abs(H_final - H_final.conj().T))
    print(f"Maximum absolute difference between H and H^dagger: {hermitian_diff:.2e}")
    if hermitian_diff < 1e-12:
        print("Verification SUCCESS: The matrix is perfectly Hermitian!")
    else:
        print("Warning: The matrix is NOT perfectly Hermitian within numerical precision.")
        
    # Save output
    np.save(args.output, H_final)
    print(f"Successfully saved dense SOC Hamiltonian of shape {H_final.shape} to {args.output}")

if __name__ == "__main__":
    main()
