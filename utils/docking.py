import os
import shutil
import pandas as pd

from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdmolfiles

def combine_sdf_pdb(sdf_file, pdb_file):
    molecule = Chem.MolFromMolFile(sdf_file)
    conf = molecule.GetConformer()
    positions = [conf.GetAtomPosition(i) for i in range(molecule.GetNumAtoms())]
    print(positions)


def get_complex(
    database:str,
    output:str,
):

    # at the moment, this code only takes into account the presence
    # of the enzyme and the amine donor, which is not the case in a
    # reaction involving the transaminase. i need to understand the
    # reaction well to know how to integrate the amine acceptor.
    
    if not os.path.exists(output):
        os.makedirs(output)
    else:
        shutil.rmtree(output)

    data = pd.read_csv(database, sep=";").reset_index()
    for index, row in tqdm(data.iterrows(), total=data.shape[0]):
        print(row.enzyme, row.substrat)

        enz_path = os.path.join(
            os.path.abspath("data/raw/enzymes"),
            f"{row.enzyme}",
            f"{row.enzyme}_relaxed.pdb",
        )

        # convert smiles to sdf
        open(f"{row.substrat}.smi", "w").writelines(f"{row.smiles}\n")

        os.system(f"""
            obabel \
            {row.substrat}.smi \
            -O {f"{row.substrat}.sdf"} \
            --gen3d
        """)
        os.remove(f"{row.substrat}.smi")

        # docking with gnina
        os.system(f"""
            gnina/build/bin/gnina \
            -r {enz_path} \
            -l {row.substrat}.sdf \
            --autobox_ligand {enz_path} \
            -o {output}/docked_{row.enzyme.lower()}_{row.substrat.replace(",", "_").lower()}.pdb \
            --exhaustiveness 64
        """)        
        os.remove(f"{row.substrat}.sdf")



if __name__ == "__main__":
    # get_complex(
    #     database="data/raw/ta_dataset_v3.csv",
    #     output="data/raw/complex",
    # )

    combine_sdf_pdb(
        sdf_file="examples/sub.sdf",
        pdb_file="examples/enz.pdb",
    )