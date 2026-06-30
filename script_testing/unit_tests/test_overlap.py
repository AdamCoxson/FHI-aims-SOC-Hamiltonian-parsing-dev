import pytest
import h5py
import numpy as np
from pathlib import Path


def load_overlap_blocks(filepath):
    """Utility to load overlap h5 file into a dictionary of block matrices."""
    blocks = {}
    with h5py.File(filepath, 'r') as f:
        atom_pairs = f['atom_pairs'][:]
        chunk_boundaries = f['chunk_boundaries'][:]
        chunk_shapes = f['chunk_shapes'][:]
        entries = f['entries'][:]
        
        for idx, pair in enumerate(atom_pairs):
            key = tuple(pair)
            shape = tuple(chunk_shapes[idx])
            start, end = chunk_boundaries[idx], chunk_boundaries[idx+1]
            block = entries[start:end].reshape(shape)
            blocks[key] = block
            
    return blocks

@pytest.fixture
def overlap_matrix(data_dir):
    h5_path = data_dir / "overlap.h5"
    if not h5_path.exists():
        pytest.skip(f"Test data not found at {h5_path}. Run parser first.")
    return load_overlap_blocks(h5_path)

def test_overlap_positive_definite(overlap_matrix):
    """
    Test 3: Overlap Matrix (S) Positive Definiteness
    We construct the Gamma-point (k=0) overlap matrix by summing S(R) over all R,
    and assert it is symmetric positive definite.
    """
    assert len(overlap_matrix) > 0, "No blocks loaded from the .h5 file."
    
    # 1. Determine matrix dimensions (number of atoms and orbitals per atom)
    max_atom = -1
    for key in overlap_matrix.keys():
        max_atom = max(max_atom, key[3], key[4])
        
    num_atoms = max_atom + 1
    
    orb_counts = {}
    for key, block in overlap_matrix.items():
        i_atom, j_atom = key[3], key[4]
        n_rows, n_cols = block.shape
        if i_atom not in orb_counts:
            orb_counts[i_atom] = n_rows
        if j_atom not in orb_counts:
            orb_counts[j_atom] = n_cols
            
    # Calculate starting indices for each atom's block in the full matrix
    total_orbs = sum(orb_counts[i] for i in range(num_atoms))
    offsets = {}
    current = 0
    for i in range(num_atoms):
        offsets[i] = current
        current += orb_counts[i]
        
    # 2. Construct S(k=0) = sum_R S(R)
    S_gamma = np.zeros((total_orbs, total_orbs), dtype=np.float64)
    
    for key, block in overlap_matrix.items():
        _, _, _, i_atom, j_atom = key
        i_start = offsets[i_atom]
        i_end = i_start + orb_counts[i_atom]
        j_start = offsets[j_atom]
        j_end = j_start + orb_counts[j_atom]
        
        S_gamma[i_start:i_end, j_start:j_end] += block
        
    # 3. Check if S_gamma is symmetric
    np.testing.assert_allclose(
        S_gamma, 
        S_gamma.T, 
        atol=1e-8, 
        err_msg="Gamma-point Overlap matrix is not symmetric"
    )
    
    # 4. Check if S_gamma is positive definite via Cholesky decomposition
    try:
        np.linalg.cholesky(S_gamma)
    except np.linalg.LinAlgError:
        pytest.fail("Cholesky decomposition failed: Overlap matrix is not positive definite.")
        
    # 5. Check eigenvalues directly as well
    eigenvalues = np.linalg.eigvalsh(S_gamma)
    assert np.all(eigenvalues > 0), f"Found non-positive eigenvalues in Overlap matrix: min(eig)={np.min(eigenvalues)}"
