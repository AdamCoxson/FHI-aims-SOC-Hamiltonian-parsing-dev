import sys
import shutil
from pathlib import Path
import time

# Import the custom DeepH SOC parser translator
from aims_soc_to_deeph import PeriodicAimsDataTranslator

def parse_soc_hamiltonian(aims_dir_str, out_dir_str, tolerance=1e-3):
    """
    Parses FHI-aims outputs to DeepH format using the noncollinear SOC logic.
    
    Args:
        aims_dir_str: The full path to the FHI-aims directory containing outputs.
        out_dir_str: The exact directory where the DeepH output files will be saved.
        tolerance: Float threshold for absorbing Hermiticity grid noise (default 1e-3).
    """
    print(f"--- Running DeepH SOC Parser ---")
    
    aims_path_obj = Path(aims_dir_str)
    out_path_obj = Path(out_dir_str)
    
    print(f"Reading from: {aims_path_obj}")
    print(f"Writing to:   {out_path_obj}")
    
    if out_path_obj.exists():
        print("Cleaning existing output directory...")
        shutil.rmtree(out_path_obj)
    out_path_obj.mkdir(parents=True, exist_ok=True)
    
    t0 = time.time()
    
    translator = PeriodicAimsDataTranslator(
        aims_data_dir=str(aims_path_obj), 
        deeph_data_dir=out_path_obj, 
        export_H=True, 
        export_H0=True,
        tolerance=tolerance
    )
    
    # By passing an empty string as dir_name, we force it to output exactly to the provided paths
    translator.transfer_one_aims_to_deeph(
        "", 
        aims_path_obj, 
        out_path_obj, 
        export_H=True, 
        export_H0=True,
        tolerance=tolerance
    )
    
    print(f"DeepH SOC Parsing finished in {time.time()-t0:.2f} seconds!")


if __name__ == '__main__':
    # ==========================================
    # USER CONFIGURABLE DIRECTORIES
    # ==========================================
    
    # Resolve paths dynamically relative to the script location
    current_dir = Path(__file__).resolve().parent
    base_dir = current_dir.parent / '2D_materials_data/MoS2_unitcell'
    #base_dir = current_dir.parent / '2D_materials_data/graphene_monolayer_unitcell'
    #base_dir = current_dir.parent / '2D_materials_data/MoS2_reproduce'
    
    # Input directory containing FHI-aims files (control.in, rs_hamiltonian.out, realspace_soc_matrix.out, etc.)
    AIMS_DIR = str(base_dir / 'soc_intermediate')
    
    # Exact output directory. 
    OUTPUT_DIR = str(base_dir / 'deepH_files')
    
    parse_soc_hamiltonian(AIMS_DIR, OUTPUT_DIR)
