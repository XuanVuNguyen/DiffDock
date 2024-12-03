import os
import glob
from typing import Optional

import pandas as pd
import numpy as np

from rdkit import Chem

import argparse

def get_com_sdf(sdf_file):
    coords = []
    with Chem.SDMolSupplier(sdf_file) as supp:
        for mol in supp:
            if mol is None:
                return None
            for i, atom in enumerate(mol.GetAtoms()):
                positions = mol.GetConformer().GetAtomPosition(i)
                # print(atom.GetSymbol(), positions.x, positions.y, positions.z)
                coords.append([positions.x, positions.y, positions.z])
        
            break
    coords = np.array(coords)
    return np.mean(coords, axis=0)

def get_com_mol2(mol2_file):
    with open(mol2_file) as f:
        lines = f.readlines()
        record = False
        coords = []
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                record = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                record = False
                break
            if record:
                info = line.strip().split()
                x, y, z = info[2:5]
                coords.append([float(x), float(y), float(z)])
    coords = np.array(coords)
    return np.mean(coords, axis=0)

def posebuster_validate(pose, protein, log_file="_pose_buster.log"):
    os.system(f"bust {pose} -p {protein} > {log_file}")
    with open(log_file, "r") as f:
        for line in f:
            if line.split()[2] == "passes":
                return True
            return False

def pocket_validate(pose_file:str, pocket_reference: str, distance_threshold: float=15):
    if pocket_reference.endswith(".mol2"):
        ref_coord = get_com_mol2(pocket_reference)
    elif pocket_reference.endswith(".sdf"):
        ref_coord = get_com_sdf(pocket_reference)
    lig_coord = get_com_sdf(pose_file)
    
    if lig_coord is None:
        return False
    
    dist = np.linalg.norm(ref_coord - lig_coord)
    
    return dist < distance_threshold

def pose_validate(pose_file:str, pocket_reference: str, protein_file:str, distance_threshold: float=15):
    return pocket_validate(pose_file, pocket_reference, distance_threshold) and posebuster_validate(pose_file, protein_file)
    

def get_pose_by_rank(ligand_dir, rank):
    ligand_files = glob.glob(os.path.join(ligand_dir, "rank*_confidence*.sdf"))
    ranks = [int(os.path.basename(ligand_file).split("_")[0].replace("rank", "")) for ligand_file in ligand_files]
    rank_dict = {r: ligand_file for r, ligand_file in zip(ranks, ligand_files)}
    return rank_dict[rank]

def energy_minimize(posedir: str, pocketfile: str, flexdist: Optional[float]):
    pose_dirs = glob.glob(os.path.join(posedir, "*"))
    pose_dirs = [pose_d for pose_d in pose_dirs if os.path.isdir(pose_d)]
    pose_dirs = sorted(pose_dirs, key=lambda x: int(os.path.basename(x).split("_")[1]))
    
    diffdock_input = pd.read_csv(os.path.join(posedir, "diffdock_input.csv"))
    diffdock_input = diffdock_input.set_index("complex_name")
    
    if os.path.exists(os.path.join(posedir, "bad_ligand.txt")):
        os.remove(os.path.join(posedir, "bad_ligand.txt"))
        
    for pose_d in pose_dirs:
        target_file = diffdock_input.loc[os.path.basename(pose_d), "protein_path"]
        
        pocket_correct = False
        for rank in range(1, 11):
            ligand_file = get_pose_by_rank(pose_d, rank)
            if pose_validate(ligand_file, pocketfile, target_file):
                pocket_correct = True
                break
        
        if not pocket_correct:
            with open(os.path.join(posedir, "bad_ligand.txt"), "a") as f:
                f.write(f"{pose_d}\n")
            continue
        
        
        # ligand_file = os.path.join(pose_d, "rank1.sdf")
        old_minimized_file = glob.glob(os.path.join(pose_d, "minimized*.sdf"))
        for old_file in old_minimized_file:
            os.remove(old_file)
        
        output_file = os.path.join(pose_d, f"minimized_rank{rank}.sdf")
        
        if flexdist > 0:
            os.system(f"./gnina -r {target_file} -l {ligand_file} --minimize --flexdist_ligand {ligand_file} --flexdist {flexdist} -o {output_file}")
        else:
            os.system(f"./gnina -r {target_file} -l {ligand_file} --minimize -o {output_file}")


    

# target_file = "mchr1/mchr1_data/mchr1_esmflow/model_2.pdb"
# result_ligand_dir = "results_October/rosettafold_active"
# pocket_reference = "kalasanty_results/rosettafold_active/pocket0.mol2"
# flex_dist = 3.5

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--pose", type=str, required=True, default="results_October/rosettafold_active")
    parser.add_argument("-p", "--pocket", type=str, required=True, default="kalasanty_results/rosettafold_active/pocket0.mol2")
    parser.add_argument("-f", "--flexdist", type=float, required=False, default=3.5)
    args = parser.parse_args()
    energy_minimize(args.pose, args.pocket, args.flexdist)