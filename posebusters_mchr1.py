import os
import glob as glob

from tqdm import tqdm
import pandas as pd

def pose_buster(pose_path, protein_path, out_file):
    os.system(f"bust {pose_path} -p {protein_path} >> {out_file}")
    
def _get_rank(ligand_path):
    basename = os.path.basename(ligand_path)
    rank = basename.lsplit("_", 1)[0]
    rank = rank.replace("rank", "")
    return int(rank)

pose_protein_path_pairs = pd.read_csv("protein_ligands_pair.csv")
for row in tqdm(pose_protein_path_pairs.itertuples(), total=len(pose_protein_path_pairs)):
    protein_path = row.protein_path
    pose_dir = row.ligand_pose_directory
    
    log_file = os.path.join(pose_dir, "posebusters.log")
    if os.path.exists(log_file):
        os.remove(log_file)
    os.system(f"touch {log_file}")
    
    ligand_dirs = glob.glob(os.path.join(pose_dir, "*"))
    ligand_dirs = [p for p in ligand_dirs if os.path.isdir(p)]
    ligand_dirs = sorted(ligand_dirs, key=lambda x: int(x.rsplit("_", 1)[-1]))
    
    for lig_dir in ligand_dirs:
        rank1_pose = os.path.join(lig_dir, "rank1.sdf")
        pose_buster(rank1_pose, protein_path, log_file)
    