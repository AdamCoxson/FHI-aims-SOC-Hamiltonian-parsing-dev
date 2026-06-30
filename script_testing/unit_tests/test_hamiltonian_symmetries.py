import pytest
import h5py
import numpy as np
from pathlib import Path



def load_h5_blocks(filepath):
    """Utility to load h5 files into a dictionary of block matrices."""
    blocks = {}
    with h5py.File(filepath, 'r') as f:
        atom_pairs = f['atom_pairs'][:]
        chunk_boundaries = f['chunk_boundaries'][:]
        chunk_shapes = f['chunk_shapes'][:]
        entries = f['entries'][:]
        
        for idx, pair in enumerate(atom_pairs):
            key = tuple(pair) # (R1, R2, R3, i_atom, j_atom)
            shape = tuple(chunk_shapes[idx])
            start, end = chunk_boundaries[idx], chunk_boundaries[idx+1]
            block = entries[start:end].reshape(shape)
            blocks[key] = block
            
    return blocks

@pytest.fixture
def soc_hamiltonian(data_dir):
    h5_path = data_dir / "soc-hamiltonian.h5"
    if not h5_path.exists():
        pytest.skip(f"Test data not found at {h5_path}. Run parser first.")
    return load_h5_blocks(h5_path)

def test_hermiticity(soc_hamiltonian):
    """
    Test 1: Final Tensor Hermiticity Check
    Verifies that H_{ij}(R) = H_{ji}(-R)^dagger
    """
    assert len(soc_hamiltonian) > 0, "No blocks loaded from the .h5 file."
    
    missing_mirrors = 0
    tested_mirrors = 0
    
    for key, block in soc_hamiltonian.items():
        R1, R2, R3, i_atom, j_atom = key
        mirror_key = (-R1, -R2, -R3, j_atom, i_atom)
        
        if mirror_key not in soc_hamiltonian:
            missing_mirrors += 1
            continue
            
        mirror_block = soc_hamiltonian[mirror_key]
        tested_mirrors += 1
        
        # Check H(R) == H(-R)^dagger
        np.testing.assert_allclose(
            block, 
            mirror_block.conj().T, 
            atol=1e-4, 
            rtol=1e-4,
            err_msg=f"Hermiticity failed for block {key} vs {mirror_key}"
        )
        
    print(f"\nTested {tested_mirrors} explicit mirror pairs. Missing mirrors: {missing_mirrors}.")
    assert tested_mirrors > 0, "No mirror pairs were found to test Hermiticity."
    
def test_time_reversal_symmetry(soc_hamiltonian):
    """
    Test 2: Time-Reversal Symmetry (TRS) Validation
    T H(R) T^-1 = H(R) where T = -i sigma_y K
    This implies the following block constraints:
      H_up_up = H_dn_dn^*
      H_up_dn = -H_dn_up^*
    """
    assert len(soc_hamiltonian) > 0, "No blocks loaded."
    
    for key, block in soc_hamiltonian.items():
        # block is 2N x 2M
        n_rows, n_cols = block.shape
        N = n_rows // 2
        M = n_cols // 2
        
        H_up_up = block[:N, :M]
        H_dn_dn = block[N:, M:]
        H_up_dn = block[:N, M:]
        H_dn_up = block[N:, :M]
        
        # Test diagonal blocks
        np.testing.assert_allclose(
            H_up_up, 
            H_dn_dn.conj(), 
            atol=1e-10, 
            err_msg=f"TRS failed for diagonal spin blocks at {key}"
        )
        
        # Test off-diagonal blocks
        np.testing.assert_allclose(
            H_up_dn, 
            -H_dn_up.conj(), 
            atol=1e-10, 
            err_msg=f"TRS failed for off-diagonal spin blocks at {key}"
        )
