import os
from glob import glob

import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level = logging.INFO)

model_dir = "mchr1/mchr1_data/mchr1_af2_active"
ligand_csv = "mchr1/ligand_data/242310_sampled_dataset.csv"
output_dir_format = "results_October/af2_active_{}/"

def get_pdb_files(dir_path: str, index_id = 1):
    paths =  glob(f"{dir_path}/*.pdb")
    return sorted(paths, key=lambda x: int(x.split("/")[-1].split(".")[0].split("_")[index_id]))

if __name__ == "__main__":
    protein_paths = get_pdb_files(model_dir, 2)
    logger.info(f"The following pdb files are found for docking: ")
    for p in protein_paths:
        logger.info(p)
    
    tmp_diffdock_input = "mchr1/diffdock_input.csv"
    if os.path.exists(tmp_diffdock_input):
        os.remove(tmp_diffdock_input)
    
    for i, protein_path in enumerate(protein_paths):
        out_dir = output_dir_format.format(i)
        logger.info(f"Running DiffDock for: {protein_path}")
        logger.info("Preparing input file...")
        os.system(f"python mchr1/diffdock_input_prepare.py -r {protein_path} -l {ligand_csv} -o {tmp_diffdock_input}")
        logger.info("Running DiffDock inference...")
        os.system(f"python -m inference --config default_inference_args.yaml --protein_ligand_csv {tmp_diffdock_input} --out_dir {out_dir}")
        os.system(f"cp -r {tmp_diffdock_input} {out_dir}")
        