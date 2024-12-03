import os
from glob import glob  
import pandas as pd

def TMscore(protein, template, tmp_file: str = "template_modeling/tmscore_tmp.txt"):
    os.system(f"template_modeling/TMscore -seq {protein} {template} > {tmp_file}")
    with open(tmp_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("TM-score"):
                return float(line.split()[2])
    return None

def get_pdb_files(dir_path: str, index_id = 1):
    paths =  glob(f"{dir_path}/*.pdb")
    return sorted(paths, key=lambda x: int(x.split("/")[-1].split(".")[0].split("_")[index_id]))

alphaflow_paths = get_pdb_files("mchr1/mchr1_data/mchr1_alphaflow")
esmflow_paths = get_pdb_files("mchr1/mchr1_data/mchr1_esmflow")
af2_inactive_paths = get_pdb_files("mchr1/mchr1_data/mchr1_af2_inactive")
af2_active_paths = get_pdb_files("mchr1/mchr1_data/mchr1_af2_active", 2)
rosettafold_paths = [
    "mchr1/mchr1_data/all_target_structures/Q99705_model.active.crderr.pdb",
    "mchr1/mchr1_data/all_target_structures/Q99705_model.inactive.crderr.pdb"
]

protein_paths = alphaflow_paths + esmflow_paths + af2_inactive_paths + af2_active_paths + rosettafold_paths


templates = {
    "inactive": "template_modeling/7ul5_A_inactive.pdb",
    'active': "template_modeling/7wic_R_active.pdb"
}
output = "template_modeling/tmscore_all_models.csv"

if __name__ == "__main__":
    with open(output, "w") as f:
        f.write("protein,inactive_tmscore,active_tmscore\n")
        
    for path in protein_paths:
        inactive_tmscore = TMscore(path, templates['inactive'])
        active_tmscore = TMscore(path, templates['active'])
        with open(output, "a") as f:
            f.write(f"{path},{inactive_tmscore},{active_tmscore}\n")