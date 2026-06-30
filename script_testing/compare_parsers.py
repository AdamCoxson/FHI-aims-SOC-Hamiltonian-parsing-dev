import sys
import numpy as np
import h5py
import shutil
from pathlib import Path

from aims_soc_to_deeph import PeriodicAimsDataTranslator, _read_ctrl_in, _parse_struct, _parse_basis, _read_mx_indices, _check_and_fix_basis_idx, _fix_loadtxt
import time

HARTREE_TO_EV = 27.2113845

def run_deeph_way(aims_dir, deeph_dir):
    print("--- Running DeepH Way ---")
    out_path = Path(deeph_dir)
    if out_path.exists():
        shutil.rmtree(out_path)
    out_path.mkdir(parents=True)
    
    t0 = time.time()
    translator = PeriodicAimsDataTranslator(aims_data_dir=aims_dir, deeph_data_dir=out_path, export_H=True, export_H0=True)
    translator.transfer_one_aims_to_deeph('soc_intermediate', Path(aims_dir), Path(deeph_dir), export_H=True, export_H0=True)
    print(f"DeepH Way finished in {time.time()-t0:.2f} s")

def custom_alternative_way(aims_dir, alt_dir):
    print("--- Running Alternative Way ---")
    aims_path = Path(aims_dir) / 'soc_intermediate'
    out_path = Path(alt_dir) / 'soc_intermediate'
    
    if out_path.exists():
        shutil.rmtree(out_path)
    out_path.mkdir(parents=True)
    
    # 1. Use DeepH dock functions to get basis mapping and dump non-Hamiltonian files
    # (POSCAR, info.json, overlap.h5)
    translator = PeriodicAimsDataTranslator(aims_data_dir=aims_dir, deeph_data_dir=alt_dir, export_H=False, export_H0=False)
    translator.transfer_one_aims_to_deeph('soc_intermediate', Path(aims_dir), Path(alt_dir), export_H=False, export_H0=False)
    
    # Extract basis and structural mapping exactly like DeepH does
    ctrl_params, spinful, spin_orbit, aims_data_type = _read_ctrl_in(aims_path)
    is_periodic, lat, site_positions, element, species, sort_idxs, total_occ_num = _parse_struct(aims_path)
    _check_and_fix_basis_idx(aims_path)
    phase_factor, orbit_quantity_list, atom_elem_dict, elem_orb_map, basis_trans_index, N_atom, N_orb, sub_idx = _parse_basis(
        aims_path, len(element), species, sort_idxs
    )
    n_ham_size, n_cells, n_basis, cell_indices, start_idx_matrix, end_idx_matrix, col_idx = _read_mx_indices(aims_path)
    
    # 2. Build our own Hamiltonian parser (The Alternative Way)
    # Read files
    h_scal = _fix_loadtxt(aims_path / "rs_hamiltonian.out", dtype=np.float64) * HARTREE_TO_EV
    with open(aims_path / "realspace_soc_matrix.out", 'r') as f:
        soc_lines = f.readlines()
    if len(soc_lines[0].split()) == 2:
        soc_lines = soc_lines[1:]
    soc_data = np.array([[float(x) for x in line.split()] for line in soc_lines if line.strip()], dtype=np.float64) * HARTREE_TO_EV
    if h_scal.shape[0] == n_ham_size + 1: h_scal = h_scal[:-1]
    if soc_data.shape[0] == n_ham_size + 1: soc_data = soc_data[:-1, :]
    p_x = soc_data[:, 0]
    p_y = soc_data[:, 1]
    p_z = soc_data[:, 2]
    
    # Output structure
    mx_R = {}
    symm_signs = [1, -1, -1, -1] # [h_scal, p_x, p_y, p_z] (signs for spatial Hermiticity logic)
    
    print("Constructing block matrices...")
    t0 = time.time()
    for idx_R in range(n_cells):
        nR = tuple(cell_indices[idx_R])
        for idx_i in range(n_basis):
            i_atom = sub_idx[idx_i, 0]
            idx_row_sub = sub_idx[idx_i, 2]
            
            start_idx = start_idx_matrix[idx_R, idx_i]
            end_idx = end_idx_matrix[idx_R, idx_i]
            if start_idx < 0 or end_idx < 0 or end_idx < start_idx - 0.1:
                continue
                
            for idx_val in range(start_idx, end_idx):
                idx_j = col_idx[idx_val]
                j_atom = sub_idx[idx_j, 0]
                idx_col_sub = sub_idx[idx_j, 2]
                
                # Setup keys
                R_key = (nR[0], nR[1], nR[2], i_atom, j_atom)
                R_hermi_key = (-nR[0], -nR[1], -nR[2], j_atom, i_atom)
                
                if R_key not in mx_R:
                    mx_R[R_key] = [np.zeros((orbit_quantity_list[i_atom], orbit_quantity_list[j_atom]), dtype=np.float64) for _ in range(4)]
                if R_hermi_key not in mx_R:
                    mx_R[R_hermi_key] = [np.zeros((orbit_quantity_list[j_atom], orbit_quantity_list[i_atom]), dtype=np.float64) for _ in range(4)]
                    
                # Values with DeepH phase factors
                phase = phase_factor[idx_i] * phase_factor[idx_j]
                vals = [h_scal[idx_val]*phase, p_x[idx_val]*phase, p_y[idx_val]*phase, p_z[idx_val]*phase]
                
                for k in range(4):
                    mx_R[R_key][k][idx_row_sub, idx_col_sub] = vals[k]
                    # Hermiticity
                    expected_val = vals[k] * symm_signs[k]
                    existing_val = mx_R[R_hermi_key][k][idx_col_sub, idx_row_sub]
                    if abs(existing_val) <= 1e-10:
                        mx_R[R_hermi_key][k][idx_col_sub, idx_row_sub] = expected_val
                    else:
                        if abs(existing_val - expected_val) >= 1e-3:
                            print(f"Hermitian check failed for k={k} at R={nR}, atom pair=({i_atom},{j_atom})")
                            assert False, f"Hermitian check failed at R={nR}"
                    
    print(f"Matrix parsing took {time.time()-t0:.2f} s")
    
    # 3. Concatenate and Spinor mapping
    print("Concatenating into chunked DeepH arrays...")
    num_pairs = len(mx_R)
    atom_pairs = np.zeros((num_pairs, 5), dtype=int)
    chunk_shapes = np.zeros((num_pairs, 2), dtype=int)
    chunk_boundaries = np.zeros((num_pairs + 1,), dtype=int)
    chunk_boundaries_spin = np.zeros((num_pairs + 1,), dtype=int)
    chunk_shapes_spin = np.zeros((num_pairs, 2), dtype=int)
    
    entries_h0 = []
    entries_soc = []
    
    for idx_pair, (R_key, mats) in enumerate(mx_R.items()):
        atom_pairs[idx_pair] = R_key
        n_rows, n_cols = mats[0].shape
        
        chunk_shapes[idx_pair] = [n_rows, n_cols]
        chunk_boundaries[idx_pair+1] = chunk_boundaries[idx_pair] + (n_rows * n_cols)
        
        chunk_shapes_spin[idx_pair] = [2*n_rows, 2*n_cols]
        chunk_boundaries_spin[idx_pair+1] = chunk_boundaries_spin[idx_pair] + (4 * n_rows * n_cols)
        
        # Build H0 spin mapping (Block diagonal)
        b_h = mats[0]
        block_h0 = np.zeros((2*n_rows, 2*n_cols), dtype=np.complex128)
        block_h0[0:n_rows, 0:n_cols] = b_h
        block_h0[n_rows:, n_cols:] = b_h
        entries_h0.append(block_h0.flatten())
        
        # Build SOC mapping
        b_px, b_py, b_pz = mats[1], mats[2], mats[3]
        H_up_up = b_h - 1j * b_pz
        H_dn_dn = b_h + 1j * b_pz
        H_up_dn = -1j * b_px - b_py
        H_dn_up = -1j * b_px + b_py
        
        block_soc = np.zeros((2*n_rows, 2*n_cols), dtype=np.complex128)
        block_soc[0:n_rows, 0:n_cols] = H_up_up
        block_soc[n_rows:, n_cols:] = H_dn_dn
        block_soc[0:n_rows, n_cols:] = H_up_dn
        block_soc[n_rows:, 0:n_cols] = H_dn_up
        entries_soc.append(block_soc.flatten())
        
    entries_h0_flat = np.concatenate(entries_h0)
    entries_soc_flat = np.concatenate(entries_soc)
    
    # 4. Dump to HDF5
    print("Dumping matrices to HDF5...")
    with h5py.File(out_path / 'hamiltonian.h5', 'w') as f:
        f.create_dataset('atom_pairs', data=atom_pairs, dtype='i4')
        f.create_dataset('chunk_boundaries', data=chunk_boundaries_spin, dtype='i4')
        f.create_dataset('chunk_shapes', data=chunk_shapes_spin, dtype='i4')
        f.create_dataset('entries', data=entries_h0_flat)
        
    with h5py.File(out_path / 'soc-hamiltonian.h5', 'w') as f:
        f.create_dataset('atom_pairs', data=atom_pairs, dtype='i4')
        f.create_dataset('chunk_boundaries', data=chunk_boundaries_spin, dtype='i4')
        f.create_dataset('chunk_shapes', data=chunk_shapes_spin, dtype='i4')
        f.create_dataset('entries', data=entries_soc_flat)

    print("Alternative Way finished!")

def compare_results(deeph_dir, alt_dir):
    import json
    from scipy.stats import pearsonr

    print("\n--- Comparing Results ---")
    deeph_path = Path(deeph_dir) / 'soc_intermediate'
    alt_path = Path(alt_dir) / 'soc_intermediate'
    
    poscar_file = alt_path / 'POSCAR'
    with open(poscar_file, 'r') as f:
        lines = f.readlines()
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    elements = []
    for s, c in zip(species, counts):
        elements.extend([s] * c)
    
    for filename in ['hamiltonian.h5', 'soc-hamiltonian.h5']:
        print(f"Comparing {filename}...")
        try:
            with h5py.File(deeph_path / filename, 'r') as f1, h5py.File(alt_path / filename, 'r') as f2:
                # Check atom pairs
                np.testing.assert_array_equal(f1['atom_pairs'][:], f2['atom_pairs'][:], err_msg="atom_pairs mismatch")
                # Check shapes and boundaries
                np.testing.assert_array_equal(f1['chunk_boundaries'][:], f2['chunk_boundaries'][:], err_msg="chunk_boundaries mismatch")
                np.testing.assert_array_equal(f1['chunk_shapes'][:], f2['chunk_shapes'][:], err_msg="chunk_shapes mismatch")
                # Check entries (with small tolerance for floating point ops order)
                np.testing.assert_allclose(f1['entries'][:], f2['entries'][:], rtol=1e-10, atol=1e-10, err_msg="entries mismatch")
            print(f"[OK] {filename} matches perfectly between DeepH and Alternative methods!")
        except Exception as e:
            print(f"[ERROR] {filename} mismatch: {e}")

    print("\n--- Extracting Test Blocks and Computing Statistics ---")
    test_blocks_dict = {}
    summary_md = ["# DeepH vs Alternative Parser: Statistical Summary\n",
                  "Output from compare_parsers.py\nThis document summarizes the block-level validation between the two independent parsers.\n"]
                  
    def safe_pearsonr(x, y):
        if np.std(x) == 0 or np.std(y) == 0:
            if np.allclose(x, y, atol=1e-12):
                return 1.0
            return np.nan
        return pearsonr(x, y)[0]

    for filename in ['hamiltonian.h5', 'soc-hamiltonian.h5']:
        is_soc = (filename == 'soc-hamiltonian.h5')
        summary_md.append(f"## {filename}\n")
        
        with h5py.File(deeph_path / filename, 'r') as f1, h5py.File(alt_path / filename, 'r') as f2:
            atom_pairs = f1['atom_pairs'][:]
            chunk_boundaries = f1['chunk_boundaries'][:]
            chunk_shapes = f1['chunk_shapes'][:]
            entries1 = f1['entries'][:]
            entries2 = f2['entries'][:]
            
            target_counts = {'S-S': 10, 'Mo-S': 10, 'Mo-Mo': 10}
            collected_blocks = {'S-S': [], 'Mo-S': [], 'Mo-Mo': []}
            
            for idx in range(len(atom_pairs)):
                pair = atom_pairs[idx]
                elem1 = elements[pair[3]]
                elem2 = elements[pair[4]]
                
                pair_name = f"{elem1}-{elem2}"
                
                if pair_name in target_counts and len(collected_blocks[pair_name]) < target_counts[pair_name]:
                    shape = tuple(chunk_shapes[idx])
                    start, end = chunk_boundaries[idx], chunk_boundaries[idx+1]
                    block1 = entries1[start:end].reshape(shape)
                    block2 = entries2[start:end].reshape(shape)
                    collected_blocks[pair_name].append((block1, block2))
            
            # """
            # # Dummy list aggregation demonstrating flattened list handling
            # # This approach doesn't require shapes to be identical:
            # mixed_blocks_1 = []
            # mixed_blocks_2 = []
            # for b1, b2 in collected_blocks['Mo-S'] + collected_blocks['Mo-Mo']:
            #     mixed_blocks_1.append(b1.flatten())
            #     mixed_blocks_2.append(b2.flatten())
            # all_flat_1 = np.concatenate(mixed_blocks_1)
            # all_flat_2 = np.concatenate(mixed_blocks_2)
            # mae_flat_mixed = np.mean(np.abs(all_flat_1 - all_flat_2))
            # """
            
            for pair_name, blocks in collected_blocks.items():
                if len(blocks) == 0:
                    continue
                b1_array = np.array([b[0] for b in blocks])
                b2_array = np.array([b[1] for b in blocks])
                test_blocks_dict[f"{filename}_{pair_name}_deeph"] = b1_array
                test_blocks_dict[f"{filename}_{pair_name}_alt"] = b2_array
                
                summary_md.append(f"### {pair_name} Pairwise Blocks\n")
                if not is_soc:
                    flat1 = b1_array.real.flatten()
                    flat2 = b2_array.real.flatten()
                    mae_flat = np.mean(np.abs(flat1 - flat2))
                    r_flat = safe_pearsonr(flat1, flat2)
                    
                    mean1 = np.mean(b1_array.real, axis=(1, 2))
                    mean2 = np.mean(b2_array.real, axis=(1, 2))
                    mae_mean = np.mean(np.abs(mean1 - mean2))
                    r_mean = safe_pearsonr(mean1, mean2)
                    
                    summary_md.append("| Statistic | Method | Real ($\\Re(H)$) |\n")
                    summary_md.append("| :--- | :--- | :--- |\n")
                    summary_md.append(f"| MAE | Flattened | {mae_flat:.3e} |\n")
                    summary_md.append(f"| Pearson $r$ | Flattened | {r_flat:.6f} |\n")
                    summary_md.append(f"| MAE | Block Average | {mae_mean:.3e} |\n")
                    summary_md.append(f"| Pearson $r$ | Block Average | {r_mean:.6f} |\n\n")
                else:
                    metrics = []
                    for comp_name, func in [('Magnitude ($|H|$)', np.abs), ('Real ($\\Re(H)$)', np.real), ('Imaginary ($\\Im(H)$)', np.imag)]:
                        comp1 = func(b1_array)
                        comp2 = func(b2_array)
                        
                        flat1 = comp1.flatten()
                        flat2 = comp2.flatten()
                        mae_flat = np.mean(np.abs(flat1 - flat2))
                        r_flat = safe_pearsonr(flat1, flat2)
                        
                        mean1 = np.mean(comp1, axis=(1, 2))
                        mean2 = np.mean(comp2, axis=(1, 2))
                        mae_mean = np.mean(np.abs(mean1 - mean2))
                        r_mean = safe_pearsonr(mean1, mean2)
                        
                        metrics.append((comp_name, mae_flat, r_flat, mae_mean, r_mean))
                    
                    summary_md.append("| Statistic | Method | Magnitude ($|H|$) | Real ($\\Re(H)$) | Imaginary ($\\Im(H)$) |\n")
                    summary_md.append("| :--- | :--- | :--- | :--- | :--- |\n")
                    summary_md.append(f"| MAE | Flattened | {metrics[0][1]:.3e} | {metrics[1][1]:.3e} | {metrics[2][1]:.3e} |\n")
                    summary_md.append(f"| Pearson $r$ | Flattened | {metrics[0][2]:.6f} | {metrics[1][2]:.6f} | {metrics[2][2]:.6f} |\n")
                    summary_md.append(f"| MAE | Block Average | {metrics[0][3]:.3e} | {metrics[1][3]:.3e} | {metrics[2][3]:.3e} |\n")
                    summary_md.append(f"| Pearson $r$ | Block Average | {metrics[0][4]:.6f} | {metrics[1][4]:.6f} | {metrics[2][4]:.6f} |\n\n")

    npz_path = Path(alt_dir) / 'test_blocks.npz'
    np.savez(npz_path, **test_blocks_dict)
    print(f"Saved extracted blocks to {npz_path}")
    
    current_dir = Path(__file__).resolve().parent
    md_path = current_dir.parent / 'md_files' / 'statistics_summary.md'
    with open(md_path, 'w') as f:
        f.write("\n".join(summary_md))
    print(f"Saved markdown summary to {md_path}")

if __name__ == '__main__':
    current_dir = Path(__file__).resolve().parent
    base_dir = str(current_dir.parent / '2D_materials_data' / 'MoS2_unitcell')
    deeph_dir = f'{base_dir}/LAO_processed_data/deeph_way'
    alt_dir = f'{base_dir}/LAO_processed_data/alternative_way'
    
    run_deeph_way(base_dir, deeph_dir)
    custom_alternative_way(base_dir, alt_dir)
    compare_results(deeph_dir, alt_dir)
