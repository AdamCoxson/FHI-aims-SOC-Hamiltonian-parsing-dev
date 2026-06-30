import os
import subprocess
import shutil
import time

from pyfhiaims.geometry import AimsGeometry
from pyfhiaims.control import AimsControl
from ase.build import graphene, mx2


def build_bilayer_graphene(supercell_xyz=(3,3,1)):
    # Build a primitive unit cell of graphene
    prim_mono = graphene(size=(1, 1, 1))
    
    # Create the bilayer structure
    prim_bi = prim_mono.copy()
    layer2 = prim_mono.copy()
    
    # For AA stacking, we don't shift the second layer fractionally
    # Just shift in z direction
    layer2.positions[:, 2] += 3.35
    prim_bi.extend(layer2)
    prim_bi.center(vacuum=10, axis=2)
    
    # Create the 3x3x1 supercell
    supercell = prim_bi.repeat(supercell_xyz)
    return supercell

def build_monolayer_graphene(supercell_xyz=(3,3,1)):
    # Build a primitive unit cell of graphene
    prim_mono = graphene(size=(1, 1, 1), vacuum=10)

    # Create the 3x3x1 supercell
    supercell = prim_mono.repeat(supercell_xyz)
    return supercell

def build_monolayer_mos2(supercell_xyz=(3,3,1)):
    # Build a primitive unit cell of MoS2
    prim_mono = mx2(formula='MoS2', kind='2H', a=3.18, thickness=3.19, size=(1, 1, 1), vacuum=10)

    # Create the supercell
    supercell = prim_mono.repeat(supercell_xyz)
    return supercell


if __name__ == "__main__":
    total_start_time = time.time()
    timing_data = {}
    num_cores = 8 # I used 28
    
    print(f"Starting FHI-aims workflow using {num_cores} cores.\n")

    # --- BLOCK 1: Initialization and Setup ---
    t_start = time.time()
    
    # Paths based on relevant_filepaths.txt
    fhipath=["/home/name/.../FHIaims_modded_with_SOC",""][0]
    light_basis = fhipath+"/species_defaults/defaults_2020/light"
    intermediate_basis = fhipath+"/species_defaults/defaults_2020/intermediate"
    tight_basis = fhipath+"/species_defaults/defaults_2020/tight"

    maindir=["/home/name/.../FHIaims_SOC_extract_dev/2D_materials_data",""][0]
    # Run or skip the relaxation calculation if you have already done it
    skip_relaxation = False
    run_command = f"mpirun -n {num_cores} "+fhipath+"/build/your_fhi_aims_exe > aims.out 2>&1"
    # For me it is
    #run_command = f"mpirun -n {num_cores} "+fhipath+"/build/aims.260530.scalapack.mpi.x > aims.out 2>&1"


    relax_dir = maindir+"/MoS2_reproduce/relaxation_light"
    soc_dir = maindir+"/MoS2_reproduce/soc_intermediate"
    nosoc_dir = maindir+"/MoS2_reproduce/nosoc_intermediate"
    os.makedirs(relax_dir, exist_ok=True)
    os.makedirs(soc_dir, exist_ok=True)
    os.makedirs(nosoc_dir, exist_ok=True)
    dummy=''
    
    # write geometry for a primitive unit cell (1x1x1)
    atoms = build_monolayer_mos2((1,1,1))
    aims_geometry = AimsGeometry.from_atoms(atoms)
    
    aims_geometry.write_file(relax_dir+"/geometry.in")
    
    # Step 1: Relaxation in graphene_nosoc_relaxation
    print("Setting up relaxation...")
    
    # Prepare the control parameters directly as dictionary to pass into Aims
    relax_params = dict(
        xc='pbe',
        vdw_ts=dummy,
        relativistic='atomic_zora scalar',
        k_grid=(8, 8, 1),
        sc_iter_limit=100,
        relax_geometry='bfgs 5e-3',
        relax_unit_cell='slab',
        species_dir=light_basis,
    )
    control_in = AimsControl(parameters=relax_params)
    control_in.write_file(aims_geometry, writer=relax_dir, overwrite=True)
    
    t_end = time.time()
    timing_data["Initialization & Setup"] = t_end - t_start
    print(f"-> Initialization & Setup finished in {timing_data['Initialization & Setup']:.2f} seconds.\n")

    # --- BLOCK 2: Relaxation Execution ---
    t_start = time.time()
    

    
    os.chdir(relax_dir)
    if skip_relaxation is True:
        print("Skipping relaxation.")
    else:
        print("Running relaxation...")
        subprocess.run(run_command, shell=True, check=True)
        print("Relaxation finished.")
        
    t_end = time.time()
    timing_data["Relaxation Run"] = t_end - t_start
    print(f"-> Relaxation Run finished in {timing_data['Relaxation Run']:.2f} seconds.\n")

    # --- BLOCK 3: Parsing Geometry and Preparing SOC/No-SOC Inputs ---
    t_start = time.time()
    
    # Step 2: Read the relaxed geometry
    if os.path.exists("geometry.in.next_step"):
        geom_relaxed = AimsGeometry.from_file(relax_dir+"/geometry.in.next_step")
        geom_relaxed.write_file(soc_dir+"/geometry.in")
        geom_relaxed.write_file(nosoc_dir+"/geometry.in")
        
        if os.path.exists(os.path.join(relax_dir, "hessian.aims")):
            shutil.copy(os.path.join(relax_dir, "hessian.aims"), os.path.join(soc_dir, "hessian.aims"))
            shutil.copy(os.path.join(relax_dir, "hessian.aims"), os.path.join(nosoc_dir, "hessian.aims"))
            print("Copied hessian.aims to soc and nosoc directories.")
    else:
        print("Relaxation failed. No relaxed geometry found.")
        exit(1)


    # Output band parameters
    bp = atoms.cell.bandpath('GMKG', npoints=50)
    band_lines = []
    path_string = bp.path
    special_points = bp.special_points

    for i in range(len(path_string) - 1):
        start_pt = path_string[i]
        end_pt = path_string[i+1]
        
        if start_pt == ',' or end_pt == ',':
            continue
            
        start_name = 'Gamma' if start_pt == 'G' else start_pt
        end_name = 'Gamma' if end_pt == 'G' else end_pt
        
        start_coord = special_points[start_pt]
        end_coord = special_points[end_pt]
        
        npoints = 50
        
        line = f"band {start_coord[0]:.8f} {start_coord[1]:.8f} {start_coord[2]:.8f} {end_coord[0]:.8f} {end_coord[1]:.8f} {end_coord[2]:.8f} {npoints} {start_name} {end_name}"
        band_lines.append(line)
    
    # NO SOC ###############################################################
    nosoc_params = dict(
        xc='pbe',
        relativistic='atomic_zora scalar',
        k_grid=(8, 8, 1),
        output_rs_matrices='plain',
        species_dir=tight_basis,
    )

    control_in_nosoc = AimsControl(parameters=nosoc_params)
    control_in_nosoc.outputs.append("dos -20 10 15001 0.05")
    for line in band_lines:
        control_in_nosoc.outputs.append(line)
    control_in_nosoc.write_file(geom_relaxed, writer=nosoc_dir, overwrite=True)

    # SOC ###############################################################
    soc_params = dict(
        xc='pbe',
        relativistic='atomic_zora scalar',
        k_grid=(8, 8, 1),
        output_rs_matrices='plain',
        include_spin_orbit=dummy,
        species_dir=tight_basis
    )

    control_in_soc = AimsControl(parameters=soc_params)
    control_in_soc.outputs.append("dos -20 10 15001 0.05")
    for line in band_lines:
        control_in_soc.outputs.append(line)
    control_in_soc.write_file(geom_relaxed, writer=soc_dir, overwrite=True)

    t_end = time.time()
    timing_data["SOC/No-SOC Prep"] = t_end - t_start
    print(f"-> SOC/No-SOC Input Prep finished in {timing_data['SOC/No-SOC Prep']:.2f} seconds.\n")

    # --- BLOCK 4: No-SOC Execution ---
    t_start = time.time()
    
    print("Running nosoc...")
    os.makedirs(nosoc_dir, exist_ok=True)
    os.chdir(nosoc_dir)
    subprocess.run(run_command, shell=True, check=True)
    print("nosoc finished.")
    
    t_end = time.time()
    timing_data["No-SOC Run"] = t_end - t_start
    print(f"-> No-SOC Run finished in {timing_data['No-SOC Run']:.2f} seconds.\n")

    # --- BLOCK 5: SOC Execution ---
    t_start = time.time()
    
    print("Running soc...")
    os.makedirs(soc_dir, exist_ok=True)
    os.chdir(soc_dir)
    subprocess.run(run_command, shell=True, check=True)
    print("soc finished.")
    
    t_end = time.time()
    timing_data["SOC Run"] = t_end - t_start
    print(f"-> SOC Run finished in {timing_data['SOC Run']:.2f} seconds.\n")
    
    # --- SUMMARY ---
    total_end_time = time.time()
    total_duration = total_end_time - total_start_time
    
    print("=" * 50)
    print(" TIMING SUMMARY")
    print("=" * 50)
    print(f" Cores Used : {num_cores}")
    print("-" * 50)
    for section, duration in timing_data.items():
        print(f" {section:<25}: {duration:>10.2f} s")
    print("-" * 50)
    print(f" {'Total Execution Time':<25}: {total_duration:>10.2f} s")
    print("=" * 50)
    print("Finished calculations.")