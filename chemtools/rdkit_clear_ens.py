#!/usr/bin/env python3

import argparse
import os

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from scipy import sparse
from scipy.io import mmwrite
from tqdm import tqdm

# ==========================================
# 1. DATA LOADING (XYZ -> Graph)
# ==========================================

def load_xyz_ensemble(file_path, charge=0):
    """
    Reads a multi-frame XYZ file and infers bonds automatically.
    
    It uses the FIRST frame to determine the molecular graph (bonds/topology).
    This topology is then enforced on all subsequent frames to ensure 
    RMSD compatibility.
    """
    structures = []
    
    # Read the whole file into memory
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    i = 0
    total_lines = len(lines)
    template_mol = None
    
    print(f"Loading XYZ file: {file_path}...")
    
    # We use a manual loop to handle the variable-length blocks
    pbar = tqdm(total=total_lines, desc="Parsing XYZ")
    
    while i < total_lines:
        try:
            # Parse header
            line = lines[i].strip()
            if not line: break # End of file
            
            n_atoms = int(line)
            
            # Extract the XYZ block for this frame
            block_end = i + n_atoms + 2
            xyz_block = "".join(lines[i : block_end])
            
            # Create a raw molecule from coordinates (no bonds yet)
            raw_mol = Chem.MolFromXYZBlock(xyz_block)
            
            if template_mol is None:
                # --- FIRST FRAME: DETERMINE TOPOLOGY ---
                try:
                    rdDetermineBonds.DetermineBonds(raw_mol, charge=charge)
                    template_mol = raw_mol
                    # print(f"Topology determined. Atoms: {template_mol.GetNumAtoms()}, Bonds: {template_mol.GetNumBonds()}")
                except Exception as e:
                    print(f"Critical Error: Failed to determine bonds on first frame: {e}")
                    return []
                
                structures.append({'mol': template_mol})
            
            else:
                # --- SUBSEQUENT FRAMES: ENFORCE TOPOLOGY ---
                # Create a fresh copy of the template (with bonds)
                new_mol = Chem.Mol(template_mol)
                new_mol.RemoveAllConformers()
                
                # Steal the conformer (coordinates) from the raw XYZ parse
                conf = raw_mol.GetConformer()
                new_mol.AddConformer(conf)
                
                structures.append({'mol': new_mol})
            
            # Update progress
            pbar.update(block_end - i)
            i = block_end
            
        except (ValueError, IndexError) as e:
            print(f"Error parsing XYZ at line {i}: {e}")
            break
            
    pbar.close()
    print(f"Successfully loaded {len(structures)} conformers.")
    return structures

# ==========================================
# 2. HELPER FUNCTIONS
# ==========================================

def prepare_rmsd_structures(structures, include_hydrogens=True):
    """
    Prepares the molecule list for RMSD calculation.
    If include_hydrogens is False, returns a new list of molecules 
    with Hydrogens removed.
    """
    if include_hydrogens:
        return structures
    
    # print("Stripping hydrogens for Heavy-Atom RMSD...")
    stripped_structures = []
    # No tqdm here as it's usually fast, but could be added for massive sets
    for s in structures:
        no_h_mol = Chem.RemoveHs(s['mol'])
        stripped_structures.append({'mol': no_h_mol})
        
    return stripped_structures

def build_multi_conformer_molecule(structures):
    """
    Combines individual molecules into one object, STRICTLY enforcing
    that Conformer ID matches the list index (0 to N-1).
    """
    if not structures: return None
    
    # 1. Create a container with NO conformers initially
    # We take the first molecule as a template for atoms/bonds, 
    # but to verify we strip any existing conformers to start clean.
    multi_conf_mol = Chem.Mol(structures[0]['mol'])
    multi_conf_mol.RemoveAllConformers()
    
    # 2. Add conformers one by one, forcing their ID to match the loop index
    for i, s in enumerate(tqdm(structures, desc="Building Ensemble")):
        # Get the conformer from the individual molecule
        original_conf = s['mol'].GetConformer(0)
        
        # Create a new conformer object to ensure we don't mess up the original
        new_conf = Chem.Conformer(original_conf)
        
        # FORCE the ID to be 'i' (0, 1, 2...)
        new_conf.SetId(i)
        
        # Add to the multi-molecule
        # addConformer(conf, assignId=False) respects the ID we just set
        multi_conf_mol.AddConformer(new_conf, assignId=False)
        
    return multi_conf_mol

def get_symmetry_complexity(mol):
    return len(mol.GetSubstructMatches(mol, uniquify=False))

def export_clusters_to_xyz(structures, clusters, output_dir="cluster_inspection"):
    """
    Exports each cluster to a separate multi-frame XYZ file.
    
    Parameters:
    -----------
    structures : list
        The original list of all conformers (dictionaries with 'mol').
    clusters : dict
        The dictionary returned by filter_conformers.
        Format: {leader_idx: [leader_idx, member_idx1, member_idx2...]}
    output_dir : str
        Directory where files will be saved.
    """
    # 1. Create directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    
    print(f"Exporting {len(clusters)} clusters to '{output_dir}/'...")
    
    # 2. Iterate over each cluster
    for leader_idx, members in clusters.items():
        filename = os.path.join(output_dir, f"cluster_{leader_idx}_size_{len(members)}.xyz")
        
        with open(filename, "w") as f:
            # Write every member of this cluster into the same file
            for member_idx in members:
                mol = structures[member_idx]['mol']
                
                # Get the XYZ block
                # Note: We use the conformer ID -1 (default) because these mols 
                # only have 1 conformer each after our loader.
                xyz_block = Chem.MolToXYZBlock(mol)
                
                # RDKit's XYZ block includes the atom count and a blank comment line.
                # We can inject a custom comment to help you identify the frame.
                lines = xyz_block.split('\n')
                if len(lines) > 1:
                    # Line 0 is atom count, Line 1 is comment
                    is_leader = "(Leader/Kept)" if member_idx == leader_idx else "(Discarded)"
                    lines[1] = f"Conformer_ID: {member_idx} {is_leader}"
                    
                f.write("\n".join(lines))
                # MolToXYZBlock sometimes misses the final newline or adds too many,
                # but standard XYZ parsers usually handle it. 
                # Just ensuring separation:
                if not lines[-1].strip() == "":
                    f.write("\n")
                    
    print("Export complete.")

def export_matrix_to_mtx(matrix, filename):
    """
    Exports the matrix to Matrix Market format using SciPy.
    Ensures 'coordinate real symmetric' format.
    """
    print(f"Exporting MTX to {filename}...")
    
    # 1. Convert to Sparse Coordinate format
    # This ensures mmwrite uses "coordinate" format instead of "array"
    sparse_mat = sparse.coo_matrix(matrix)
    
    # 2. Write to file
    # We specify symmetry='symmetric' so SciPy checks and writes only the 
    # necessary entries (usually lower triangular) for efficient storage.
    mmwrite(filename, sparse_mat, symmetry='symmetric')
    
    print("MTX export complete.")
# ==========================================
# 3. RMSD CALCULATION STRATEGIES
# ==========================================

def compute_pairwise_rmsd(structures, include_hydrogens=True):
    """
    Strategy A: Pairwise Alignment (Slow, Robust).
    Used for highly symmetric molecules to avoid combinatorial explosion.
    """
    work_structures = prepare_rmsd_structures(structures, include_hydrogens)
    n_conformers = len(work_structures)
    rmsd_matrix = np.zeros((n_conformers, n_conformers))
    
    for i in tqdm(range(n_conformers), desc="Pairwise RMSD"):
        for j in range(i + 1, n_conformers):
            # AlignMol aligns j onto i. Returns best RMSD.
            val = AllChem.AlignMol(work_structures[j]['mol'], work_structures[i]['mol'])
            rmsd_matrix[i, j] = val
            rmsd_matrix[j, i] = val
            
    return rmsd_matrix

def compute_symmetry_aware_rmsd(structures, include_hydrogens=True):
    # 1. Prepare
    work_structures = prepare_rmsd_structures(structures, include_hydrogens)
    n = len(work_structures)
    
    # 2. Build Multi-Conformer Mol
    mol_multi = build_multi_conformer_molecule(work_structures)
    
    print("Computing BestRMS Matrix (Symmetry Aware)...")
    
    # 3. Pre-calculate Maps
    # We find all automorphism permutations.
    # uniquify=False is MANDATORY to get symmetry mappings.
    matches = mol_multi.GetSubstructMatches(mol_multi, uniquify=False)
    
    # Convert to format for GetAllConformerBestRMS: List of Lists of (Probe, Ref)
    custom_maps = []
    for match in matches:
        atom_map = []
        for query_idx, match_idx in enumerate(match):
            atom_map.append((query_idx, match_idx))
        custom_maps.append(atom_map)
        
    print(f"Using {len(custom_maps)} symmetry mappings per pair.")
        
    # 4. Run Calculation (C++ Optimized)
    rms_list = AllChem.GetAllConformerBestRMS(mol_multi, map=custom_maps)
    mat = np.zeros(n**2).reshape(n,n)
    k = 0
    for i in range(n):
        for j in range(i):
            mat[i,j] = rms_list[k]
            mat[j,i] = rms_list[k]
            k+=1
    return mat

def compute_rmsd_matrix(structures, symmetry_threshold=1000, include_hydrogens=True):
    """
    Main Coordinator: Selects the best RMSD strategy.
    """
    if not structures:
        return np.array([])
    
    # Decide complexity based on the graph we will actually use
    ref_mol = structures[0]['mol']
    if not include_hydrogens:
        ref_mol = Chem.RemoveHs(ref_mol)
        
    complexity = get_symmetry_complexity(ref_mol)
    
    method_name = "Heavy-Atom" if not include_hydrogens else "All-Atom"
    print(f"RMSD Mode: {method_name} | Symmetry Complexity: {complexity}")

    if complexity > symmetry_threshold:
        print(f"-> High symmetry detected. Using Iterative Pairwise alignment.")
        return compute_pairwise_rmsd(structures, include_hydrogens)
    
    try:
        print(f"-> Standard symmetry. Using Batch Matrix alignment.")
        return compute_symmetry_aware_rmsd(structures, include_hydrogens)
    except ValueError as e:
        print(f"-> Batch alignment failed ({e}). Fallback to Pairwise.")
        return compute_pairwise_rmsd(structures, include_hydrogens)

# ==========================================
# 4. FILTERING / CLUSTERING
# ==========================================

def filter_conformers(structures, rmsd_matrix, threshold=0.5):
    """
    Greedy clustering to remove duplicates.
    """
    n_conformers = len(structures)
    is_assigned = np.zeros(n_conformers, dtype=bool)
    clusters = {}
    unique_indices = []
    
    print(f"Filtering conformers (Threshold: {threshold} Å)...")
    
    for i in range(n_conformers):
        if is_assigned[i]:
            continue
            
        # Found a new cluster leader
        unique_indices.append(i)
        clusters[i] = [i]
        is_assigned[i] = True
        
        # Absorb neighbors
        for j in range(i + 1, n_conformers):
            if not is_assigned[j]:
                if rmsd_matrix[i, j] < threshold:
                    clusters[i].append(j)
                    is_assigned[j] = True
    
    unique_structures = [structures[i] for i in unique_indices]
    return unique_structures, clusters

def print_filtering_report(clusters):
    print(f"\n--- Filtering Report ---")
    print(f"Total Unique Clusters: {len(clusters)}")
    for leader, members in clusters.items():
        duplicates = [m for m in members if m != leader]
        if duplicates:
            print(f"  Conformer {leader:3d} kept. Merged: {duplicates}")
        else:
            print(f"  Conformer {leader:3d} kept. (Unique)")
    print("-" * 25)

# ==========================================
# 5. MAIN EXECUTION
# ==========================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    parser.add_argument('--threshold', type=float, default=0.25, help="RMSD threshold for filtering (default: 0.25 Å)")
    parser.add_argument('--include-hydrogens', action='store_true', help="Include hydrogens in RMSD calculation (default: False)")  
    args=parser.parse_args()

    # Settings
    XYZ_FILE = args.file  
    RMSD_THRESHOLD = args.threshold         # Angstroms
    INCLUDE_H = args.include_hydrogens      # Set False for heavy-atom only RMSD
    
    # 1. Load
    structures = load_xyz_ensemble(XYZ_FILE, charge=0)
    
    if structures:
        # 2. Compute Matrix
        matrix = compute_rmsd_matrix(
            structures, 
            symmetry_threshold=1000, 
            include_hydrogens=INCLUDE_H
        )
        
        # 3. Filter
        unique_confs, cluster_data = filter_conformers(
            structures, 
            matrix, 
            threshold=RMSD_THRESHOLD
        )
        
        # 4. Report
        print_filtering_report(cluster_data)
        
        print(f"Done. Kept {len(unique_confs)} unique structures.")
        
        # Save unique structures back to XYZ
        with open("unique_conformers.xyz", "w") as f:
            for s in unique_confs:
                f.write(Chem.MolToXYZBlock(s['mol']))
        
        export_clusters_to_xyz(structures, cluster_data, output_dir="clusters")
        export_matrix_to_mtx(matrix, "rmsd_matrix.mtx")