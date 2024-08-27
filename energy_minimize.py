import os
import glob

target_file = "mchr1/mchr1_data/mchr1_active/mchr1_active_2_0.79.pdb"
result_ligand_dir = "results/alphafold_gpcr_active_2"

for ligand_dir in glob.glob(os.path.join(result_ligand_dir, "*")):
    ligand_file = os.path.join(ligand_dir, "rank1.sdf")
    output_file = os.path.join(ligand_dir, "minimized.sdf")
    os.system(f"./gnina -r {target_file} -l {ligand_file} --minimize -o {output_file}")