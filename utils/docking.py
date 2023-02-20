import os
import shutil
import pandas as pd

from tqdm import tqdm
from process import get_smiles_selfies

enz_name = "C7RAM1"
sub_name = "2OATA"
sub_smiles_seq = "C(=O)C(=O)O"

def get_complex(
    database:str,
    output:str,
    temp_dir: str = None,
):

    if not os.path.exists(output):
        os.makedirs(output)
    else:
        shutil.rmtree(output)

    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    else:
        shutil.rmtree(temp_dir)

    data = pd.read_csv(database, sep=";").reset_index()
    for index, row in tqdm(data.iterrows(), total=data.shape[0]):
        print(row.enzyme, row.substrat)

        enz_path = os.path.join(
            os.path.abspath("data/raw/enzymes"),
            f"{row.enzyme}",
            f"{row.enzyme}_relaxed.pdb",
        )
        out_sub_path = os.path.join(temp_dir, f"{row.substrat}.sdf")

        # convert smiles to sdf
        if not os.path.exists(out_sub_path):
            open(f"{row.substrat}.smi", "w").writelines(f"{row.smiles}\n")

            os.system(f"""
                obabel \
                {row.substrat}.smi \
                -O {out_sub_path} \
                --gen3d
            """)

        # docking with gnina
        os.system(f"""
            gnina/build/bin/gnina \
            -r {enz_path} \
            -l {out_sub_path} \
            --autobox_ligand {row.enzyme}.pdb \
            -o docked_{row.enzyme.lower()}_{row.substrat.replace(",", "_").lower()}.pdb \
            --exhaustiveness 64
        """)        
        break

    shutil.rmtree(temp_dir)


get_complex(
    database="data/raw/ta_dataset_v3.csv",
    output="data/raw/complex",
    temp_dir="data/raw/temp",
)