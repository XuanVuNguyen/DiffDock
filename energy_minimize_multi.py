import os 
import glob
import re
from tqdm import tqdm

import pandas as pd

from energy_minimize import energy_minimize

posedirs = glob.glob("results_October/af2_active_*")
posedirs = [posedir for posedir in posedirs if os.path.isdir(posedir)]
posedirs = sorted(posedirs, key=lambda x: int(re.search(r"af2_active_(\d+)", x).group(1)))
print("Number of found pose dirs: ", len(posedirs))

pocket_csv = "kalasanty_pockets.csv"
flexdist = 5

if __name__ == "__main__":
    pocket_df = pd.read_csv(pocket_csv)
    pocket_df = pocket_df.set_index("protein_structure")

    for i, posedir in enumerate(posedirs):
        diffdock_input_csv = os.path.join(posedir, "diffdock_input.csv")
        diffdock_input = pd.read_csv(diffdock_input_csv)
        protein_structure = diffdock_input["protein_path"].tolist()[0]
        pocket_reference = pocket_df.loc[protein_structure, "pocket"]
        if pocket_reference is not None:
            energy_minimize(posedir, pocket_reference, flexdist)